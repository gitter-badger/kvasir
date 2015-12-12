#!/usr/bin/env python
# by Kevin Bonham, PhD (2015)
# for Dutton Lab, Harvard Center for Systems Biology, Cambridge MA
# Unless otherwise indicated, licensed under GNU Public License (GPLv3)

'''
Must have Mongod running, in terminal: `mongod --dbpath path/to/db`
'''

import os, re
from subprocess import Popen, PIPE
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from bson.objectid import ObjectId
from pymongo.cursor import Cursor
from KvasirHGT import core_hgt_groups
import KvDataStructures as kv
from DataImport import get_dna_seq
from KvasirPhylo import get_list_from_tree

def make_blast_db(source, name=None, remove_source=True):
    """
    Produces BLAST database from `source`
    Optional - provide name (defaults to `source`)
    Set remove_source=False to keep fasta file (if created)
    
    Source types:
    - fasta file (use path, must end with `.fna`)
    - Mongo collection (use name of collection)
    - list of dicts containing at least keys `species`, `_id`, `dna_seq`
    - Mongo cursor eg. `collection.find({'key':value})`
    """
    # If there's no directory for blast db's, create one
    if not os.path.isdir('blast_databases/'):
        os.makedirs('blast_databases')
    
    output_fasta = None
    
    if os.path.isfile(source):
        # Input is fasta file?
        if source.endswith('.fna'):
            output_fasta = source
            if not name:
                name = os.path.basename(source)[:-4]
            remove_source = False
        else:
            print "Not a valid file type, use .fna"

    else:
        output_fasta = '{0}_all.fasta'.format(kv.db.name)     
        genes = None
        with open(output_fasta, 'w+') as output_handle:
            if source in kv.get_collections():
                genes = kv.get_collection(source).find()
                if not name:
                    name = source
            elif type(source) == list:
                genes = source
            elif type(source) == Cursor:
                genes = source
        
            for gene in genes:
                output_handle.write('>{0}|{1}\n{2}\n'.format(
                    gene['species'],
                    gene['_id'],
                    gene['dna_seq'],
                    )
                )

    while not name:
        name = str(raw_input("enter name for BLAST database: "))

    # calls makeblastdb from shell
    print "making a database!"
    Popen(
        ['makeblastdb',
        '-in', output_fasta,
        '-dbtype', 'nucl',
        '-out', 'blast_databases/{0}'.format(name),
        '-title', name,
        ]
    ).wait() # waits for this operation to terminate before moving on

    if remove_source:
        os.remove(output_fasta)

def make_prot_blast_db(source, name=None, remove_source=True):
    """
    Produces BLAST database from `source`
    Optional - provide name (defaults to `source`)
    Set remove_source=False to keep fasta file (if created)
    
    Source types:
    - fasta file (use path, must end with `.faa`)
    - Mongo collection (use name of collection)
    - list of dicts containing at least keys `species`, `_id`, `aa_seq`
    - Mongo cursor eg. `collection.find({'key':value})`
    """
    # If there's no directory for blast db's, create one
    if not os.path.isdir('blast_databases/'):
        os.makedirs('blast_databases')
    
    output_fasta = None
    
    if os.path.isfile(source):
        # Input is fasta file?
        if source.endswith('.faa'):
            output_fasta = source
            if not name:
                name = os.path.basename(source)[:-4]
            remove_source = False
        else:
            print "Not a valid file type, use .faa"

    else:
        output_fasta = '{0}_all.fasta'.format(kv.db.name)     
        genes = None
        with open(output_fasta, 'w+') as output_handle:
            if source in kv.get_collections():
                genes = kv.get_collection(source).find()
                if not name:
                    name = source
            elif type(source) == list:
                genes = source
            elif type(source) == Cursor:
                genes = source
        
            for gene in genes:
                output_handle.write('>{0}|{1}\n{2}\n'.format(
                    gene['species'],
                    gene['_id'],
                    gene['aa_seq'],
                    )
                )

    while not name:
        name = str(raw_input("enter name for BLAST database: "))

    # calls makeblastdb from shell
    print "making a database!"
    Popen(
        ['makeblastdb',
        '-in', output_fasta,
        '-dbtype', 'prot',
        '-out', 'blast_databases/{0}'.format(name),
        '-title', name,
        ]
    ).wait() # waits for this operation to terminate before moving on

    if remove_source:
        os.remove(output_fasta)

def blast_vs_fasta(query, subject):
    """
    Blast `query` against `subject`. Both must be paths to fasta file
    Returns list of lists, each `[sseqid, qseqid, pident, qlen, length]`   
    """
    out = Popen(
        ['blastn',
        '-query', query,
        '-subject', subject,
        '-outfmt', '10 sseqid qseqid pident qlen length',
        ], stdout=PIPE
    ).communicate()[0]
    if out:
        result = re.split(r'[,\n]', out.rstrip())[0]
        return [result[i:i+4] for i in range(len(result))[0::4]]

def blast_vs_db(query, db):
    """
    Blast `subject` (fasta file) against `db` (blast db). 
    Returns list of lists, each `[qseqid, sseqid, pident, length]`   
    """
    out = Popen(
        ['blastn',
        '-query', query,
        '-db', db,
        '-outfmt', '10 qseqid sseqid pident length',
        ], stdout=PIPE
    ).communicate()[0]
    if out:
        result = re.split(r'[,\n]', out.rstrip())
        return [result[i:i+4] for i in range(len(result))[0::4]]

def blastp_vs_db(query, db):
    """
    Blast `subject` (fasta file) against `db` (blast db). 
    Returns list of lists, each `[qseqid, sseqid, pident, length]`
    NOTE: Will not work if qseqid or sseqid has a "," in it   
    """
    out = Popen(
        ['blastp',
        '-query', query,
        '-db', db,
        '-outfmt', '10 qseqid sseqid pident length qlen',
        ], stdout=PIPE
    ).communicate()[0]
    if out:
        result = re.split(r'[,\n]', out.rstrip())
        return [result[i:i+5] for i in range(len(result))[0::5]]

def core_hgt_blast(perc_identity='99'):
    """
    Blasts all core genomes against core db
    - Set `perc_identity` if desired (default = 99)
    """
    if not os.path.isdir('blast_results/core/'):
        os.makedirs('blast_results/core/')
    for species in kv.get_collection('core').distinct('species'):
        query_fasta = 'blast_results/core/{}_tmp.fna'.format(species)
        
        with open(query_fasta, 'w+') as query_handle:
            for query in kv.get_collection('core').find({'species':species}):
                if query['type'] == 'gene':
                    query_handle.write('>{0}|{1}\n{2}\n'.format(
                        query['species'],
                        query['_id'],
                        query['dna_seq']
                        )
                    )
        
        print 'Blasting {0}'.format(species)
        out = Popen(
            ['blastn',
            '-query', query_fasta,
            '-db', 'blast_databases/core',
            '-outfmt', '5',
            '-out', 'blast_results/core/{}_{}_blast.xml'.format(species, perc_identity),
            '-perc_identity', perc_identity
            ],
            stdout=PIPE
        ).communicate()[0]

        os.remove(query_fasta)

def blast_to_db(db='core', perc_identity='99'):
    blast_dir = 'blast_results/{}/'.format(db)
    for f in os.listdir(blast_dir):
        if f.endswith('{}_blast.xml'.format(perc_identity)):
            file_handle = 'blast_results/{}/{}'.format(db,f)
            with open(file_handle, 'r') as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                hits_dict = {}
                for blast_record in blast_records:
                    query_parse = re.search(r'(\w+)\|(\w+)', blast_record.query)
                    query_genus_parse = re.match(r'([A-Za-z]+)_', blast_record.query)
                    query_genus = query_genus_parse.group(1)
                    query_name = query_parse.group(1)
                    query_id = query_parse.group(2)

                    hits_dict[query_id] = []

                    for alignment in blast_record.alignments:
                        hit_parse = re.search(r'(\w+)\|(\w+)', alignment.hit_def)
                        hit_genus_parse = re.match(r'([A-Za-z]+)_', alignment.hit_def)
                        hit_genus = hit_genus_parse.group(1)

                        hit_name = hit_parse.group(1)
                        hit_id = hit_parse.group(2)
                        if query_name == hit_name:
                            pass
                        elif query_genus == hit_genus:
                            print "Oops! {} and {} are the same genus, skipping...".format(query_name, hit_name)
                            pass
                        elif kv.get_mongo_record(hit_name, hit_id)['type'] == '16S':
                            print 'Skipping 16S hit'
                        else:
                            print '=======\nhit for {0} detected:\nspecies: {1}\n======='.format(query_name, hit_name)
                            hits_dict[query_id].append((hit_name, hit_id))
                    
                print 'Updataing mongoDB with hits'
                hits_collection = kv.get_collection('hits')
                hits_collection.update_one(
                    {'species':query_name},
                    {'$set':{'{}_hits_{}'.format(db, perc_identity):{x:hits_dict[x] for x in hits_dict if hits_dict[x]}}},
                    upsert=True
                    ) 

def hits_reset():
    kv.remove_collection('hits')

def other_blast():
    groups_list = core_hgt_groups()
    groups_list.sort(key=len, reverse=True)
    
    for i in range(len(groups_list)):
        group_hits = []
        kv.make_id_list_fasta(groups_list[i], 'core')
        results = blast_vs_db('tmp.fna', 'blast_databases/other')

        hits_collection = kv.get_collection('hits')
        if results:
            for j in range(len(results)):
                group_hits.append(results[j])
        hits_collection.insert_one({'group':(i+1), 'group_hits':group_hits})

def tree_from_gb(some_gb):
    with open(some_gb, 'r') as input_handle, open('species_list.txt', 'w+') as output_handle:
        gbk = SeqIO.parse(input_handle, 'gb')
        s_list = []
        for record in gbk:
            new_name = record.annotations['source'].replace(' ', '_').replace(',', '')
            if new_name in s_list:
                pass
            else:
                s_list.append(new_name)
                output_handle.write('{}\n'.format(new_name))
                

def parse_faa(faa):
    '''
    Do not use.
    '''
    with open(faa, 'r') as input_handle:
        in_faa = SeqIO.parse(input_handle, 'fasta')
        for record in in_faa:
            d = record.description
            a = d.split('|')
            species, contig = a
            kv.get_collection('all_cheese').insert_one({
                    'species':species,
                    'contig':contig,
                    'seq':str(record.seq)
                    }
                )

def scratch_import(gbk):
    '''
    Used once.
    '''
    with open(gbk, 'r') as open_file:
        collection = kv.get_collection('all_cheese')
        for record in SeqIO.parse(open_file, 'gb'):
            current_contig = record.name
            current_species = record.annotations['source']
            for feature in record.features:
                if feature.type == 'CDS':
                    parsed_location = kv.get_gene_location(feature.location)
                    gene_record = {
                        'species':current_species.replace(' ', '_').replace(',',''),
                        'location':{
                            'contig':current_contig,
                            'start':parsed_location[0],
                            'end':parsed_location[1],
                            'strand':parsed_location[2],
                            'index':None
                        },
                        'annotation':feature.qualifiers['product'][0],
                        'dna_seq':get_dna_seq(feature, record),
                        'aa_seq':feature.qualifiers['translation'][0],
                        'type':'gene'
                        }
                    collection.insert_one(gene_record)

def parse_blast(blast_results):
    
    for result in blast_results:
        query = result[0]
        sspecies, scontig = result[1].split('|')
        pident = float(result[2])
        length = float(result[3])
        qlen = float(result[4])

        kv.get_collection("all_cheese").insert_one({
            'type':'blast_result',
            'query':query,
            'query_length':qlen,
            'subject_species':sspecies,
            'subject_contig':scontig,
            'perc_identity':pident,
            'match_length':length
            })
        # if length > 0.8 * qlen and pident > 50:
        #     print "===\n{}, {}\n{}, {}\n{}, {}\n===\n".format(
        #         query,
        #         qlen,
        #         sspecies,
        #         scontig,
        #         pident,
        #         length
        #         )

def make_draw_template():
    
    with open('draw_template.txt', 'w+') as o:
        l = get_list_from_tree('/Users/KBLaptop/computation/tmp/blast_tests/ML_tree.newick')
        collection = kv.get_collection('all_cheese')
        queries = ['Transcriptional_regulator_AraC_family_CDS_translation','iron-chelator_utilization_protein_CDS_translation','Putative_transport_protein/putative_regulator_CDS_translation','Vitamin_B12_ABC_transporter_B12-binding_component_BtuF_CDS_translation','Vitamin_B12_ABC_transporter_B12-binding_component_BtuF_CDS_translation','ABC-type_Fe3+-siderophore_transport_system_permease_component_CDS_translation','ABC-type_Fe3+-siderophore_transport_system_permease_2_component_CDS_translation','ABC-type_Fe3+-siderophore_transport_system_ATPase_component_CDS_translation','Transport_ATP-binding_protein_CydCD_CDS_translation','ABC_transporter_ATP-binding_protein_CDS_translation','Transcriptional_regulator_AraC_family_CDS_translation',]

        for species in l:
            print species
            o.write('{}\n'.format(species))
            
            for query in queries:
                results = collection.find({
                    'type':'blast_result',
                    'subject_species':species,
                    'query':query
                })
                max_pident = 0
                for result in results:
                    match_length = float(result['match_length'])
                    query_length = float(result['query_length'])
                    pident = float(result['perc_identity'])
                    if match_length > 0.8 * query_length:
                        if result['perc_identity'] > max_pident:
                            max_pident = result['perc_identity']
                
                print max_pident
                o.write('{}\n'.format(max_pident))

if __name__ == '__main__':
    # import sys
    kv.mongo_init('pacbio2')
    os.chdir('/Users/KBLaptop/computation/tmp/blast_tests/')
    # kv.mongo_init(sys.argv[1])
    # os.chdir('output/{}/'.format(sys.argv[1]))
    # make_blast_db('core')
    make_prot_blast_db('other')
    # hits_reset()
    # hgt_blast(perc_identity='90')
    # hgt_blast(perc_identity='95')
    # hgt_blast(perc_identity='99')
    # blast_to_db(perc_identity='90')
    # blast_to_db(perc_identity='95')
    # blast_to_db(perc_identity='99')
    
    # other_blast()
    # scratch_import("IMG_cheese-milk-dairy_30894.gbk")
    # blast_results = blastp_vs_db("Arthro RUSTI subregion 1 translation.fasta", "other")
    # parse_blast(blast_results)

    # tree_from_gb("IMG_cheese-milk-dairy_30894.gbk")
    # make_draw_template()
