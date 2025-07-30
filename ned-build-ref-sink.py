#!/usr/bin/env python

import argparse
import os
import re
import subprocess
import glob

parser = argparse.ArgumentParser(prog='Sink index', description='')
parser.add_argument('--path_to_ref', "-p", help='path to the reference directory [input is directory]', default='.')
parser.add_argument('--version',  help='print version',action='store_true', default=False)

args=parser.parse_args()

version="Bacterial sink: version - v1.0 25 March-2025"
version="sink index: version - v1.0 18 April-2025"


path_split = os.path.abspath(args.path_to_ref).split('/')
if 'fungi' in path_split:
	sink_group = 'fungi'
elif 'bacteria' in path_split:
	sink_group = 'bacteria'
elif 'archaea' in path_split:
	sink_group = 'archaea'
else:
	print(f'sink olny supported for fungi, bacteria, and archaea | this has to be present in the path')
	quit()


def check_database_version():
	db_version_summaries = os.listdir(f'{args.path_to_ref}/database_version_summary/')
	database_version = re.sub('summary_version_|.tsv','' ,db_version_summaries[-1])
	try:
		last_version = open('sink_version', 'r')
		if last_version == database_version:
			print('sink is up to data')
			quit()
	except:
		pass

def read_database_summary():
	db_version_summaries = os.listdir(f'{args.path_to_ref}/database_version_summary/')
	database_version = re.sub('summary_version_|.tsv','' ,db_version_summaries[-1])
	print(f"database version: {database_version}")
	file = open(f'{args.path_to_ref}/database_version_summary/{db_version_summaries[-1]}')
	
	phylum_dict = {}

	i = 0
	for row in file:
		if 'assembly_accession' in row:
			continue
		assembly_accession, taxa_name, phylum_name, order_name, family_name, genus_name, taxid, Assembly_level, Genome_representation =row.strip().split('\t')
		if phylum_name not in phylum_dict:
			phylum_dict[phylum_name] = []
		phylum_dict[phylum_name].append(glob.glob(f"{assembly_accession}/*.gz")[0])
	
	file.close()
	return phylum_dict, database_version
		

def index_bacteria_sinks(phylum_dict):
	if not os.path.isdir(f'{args.path_to_ref}/sink_phylum_indexes'):
		os.makedirs(f'{args.path_to_ref}/sink_phylum_indexes')

	# cat'ing fasta files	
	for phylum in phylum_dict:
		#print(phylum_dict[phylum])
		if phylum == 'NA':
			continue
		if 'Candidatus' in phylum:
			continue

		print(f"writing {phylum}")

		if not os.path.isdir(f"{args.path_to_ref}/sink_phylum_indexes/{phylum}_sink"):
			os.makedirs(f"{args.path_to_ref}/sink_phylum_indexes/{phylum}_sink")
		subprocess.run(['cat'] + phylum_dict[phylum], stdout=open(f"{args.path_to_ref}/sink_phylum_indexes/{phylum}_sink/{phylum}_{sink_group}.fna.gz", 'w'))

if __name__ == '__main__':
    if args.version == True :
        print(version)
        quit()

    check_database_version()
    phylum_dict, database_version = read_database_summary()
    index_bacteria_sinks(phylum_dict)
    file = open('sink_version', 'w')
    file.write(database_version)
    file.close

