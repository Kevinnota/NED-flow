#!/usr/bin/env python

import argparse
from tqdm import tqdm
import re
import os

parser = argparse.ArgumentParser(prog='ned-check.py', description='Checks if all mappings are done and checks for missing files')
parser.add_argument('--path', "-p", help='path to the direcotry NED-flow was run in', default='.')
parser.add_argument('--path_to_references', "-pr", help='path to the references', default='/mnt/sequencedb/bowtie2_RefSeq/800_ncbi_genomes/')
parser.add_argument('--bam', "-b", help='path to bam direcotry' ,default = 'bams')	
parser.add_argument('--version',  help='print version',action='store_true', default=False)
args=parser.parse_args()

version="NCBI Genome indexer: version - v1.3 23 Jul-2025"

def read_run_logs():
	path = args.path.strip('/')
	file = open(f'{args.path}/run_logs/database_version.txt')
	db2version = {}
	db_set_list = set()

	print(f'\nReads mapped to to fillowing databases with versions\n')
	for row in file:
		if '#' in row:
			continue
		db = re.sub(' .*', '', row.strip())
		db_set_list.add(db)
		version = re.sub('database_version=|;', '', row.strip().split(' ')[1])
		print(f'\t{db}, version {version}')
		db2version[db]=version
	file.close()
	return db2version, db_set_list
	
def read_database_summary():
	db_summary_dict={}
	for db in db2version.keys():
		version_formatted = "{:04}".format(int(db2version[db]))
		db_summary = open(f"{args.path_to_references}/{db}/database_version_summary/summary_version_{version_formatted}.tsv")
		for row in db_summary:
			if 'assembly_accession' in row:
				continue
			
			assembly_accession, taxa_name, phylum_name, order_name, family_name, genus_name, taxid, Assembly_level, Genome_representation = row.strip().split('\t')
			db_summary_dict[assembly_accession] = [db, taxa_name]
		db_summary.close()

	return db_summary_dict

def check_if_assembly_is_missing():
	path = args.path.strip('/')
	assemblies = os.listdir(f'{path}/{args.bam}')
	set_assemblies=set()
	for assembly in assemblies:
		assembly = assembly.strip('_mapped.bam')
		try :
			assembly in db_summary_dict
			set_assemblies.add(assembly)
		except:
			print(f'{assembly} not in db')
	
	#print(set(db_summary_dict.keys()))
	#print(set_assemblies)
	#quit()
	missing = set(db_summary_dict.keys()) - set_assemblies

	print(f'\nThe following assemlies are in the database but not mapped to:')
	for db in db_set_list:
		print(f'\n{db}:')
		i=0
		for assembly in missing:
			DB, taxa_name = db_summary_dict[assembly]
			if DB == db:
				print(f'\t{assembly}\t{db}\t{taxa_name}')
				i+=1
		if i == 0:
			print(f'\tthere were no missing assemblies')
	print(f'\n')

if __name__ == '__main__':
    if args.version == True :
        print(version)
        quit()

    db2version, db_set_list = read_run_logs()
    db_summary_dict = read_database_summary()
    check_if_assembly_is_missing()



