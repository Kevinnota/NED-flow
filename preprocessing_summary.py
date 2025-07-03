#!/usr/bin/env python

import argparse
import re
import json

parser = argparse.ArgumentParser(prog='Taxon classifier', description='', usage='%(prog)s [options]')
parser.add_argument('--fastp_json', "-fj", nargs='+', help='input fastp json files')
parser.add_argument('--dust_summary', "-ds", nargs='+', help='output from dustmasker filter')

parser.add_argument('--version',  help='print version',action='store_true', default=False)
args=parser.parse_args()

#version="taxon classifier - v1.0 18 jun 2025"
### New tool to calculate library filtering stats from fastp output
version="taxon classifier - v1.0 18 jun 2025"

def read_dust_filter_file():
	lib2dust_removed_dict={}
	for dust_summary in args.dust_summary:
		file = open(dust_summary, 'r')
		try:
			for row in file:
				if "library =" in row:
					library = re.sub('library = ', '', row.strip())
					
				if 'number sequence removed' in row:
					n_seq_removed = re.sub('number sequence removed|, all.*', '', row.strip())
			#print(f"{library}, {n_seq_removed}")
			lib2dust_removed_dict[library] = n_seq_removed
			pass
		except:
			pass

	return lib2dust_removed_dict
	
	
	print(lib2dust_removed_dict)
def read_fastp_json_fils():
	header = ['library',
					'N_reads_before_filter',
		 			'N_reads_after_filter',
		 			'duplication_rate',
		 			'low_quality_reads',
		 			'too_many_N_reads',
		 			'low_complexity_reads', 
					'too_short_reads',
					'low_dust_score'] 
	print('\t'.join(header))
	for log in args.fastp_json:

		libary_name = re.sub('.*/|s_._|_S.*', '', re.sub('_report.json', '', log))
		
		try:
			dust_removed=lib2dust_removed_dict[libary_name]
		except:
			dust_removed="NA"
			pass

		log_file = json.load(open(log, 'r'))
		value_list=[libary_name,
					log_file['summary']['before_filtering']['total_reads'],
		 			log_file['summary']['after_filtering']['total_reads'],
		 			log_file['duplication']['rate'],
		 			log_file['filtering_result']['low_quality_reads'],
		 			log_file['filtering_result']['too_many_N_reads'],
		 			log_file['filtering_result']['low_complexity_reads'], 
					log_file['filtering_result']['too_short_reads'],
					dust_removed]

		print("\t".join(str(v) for v in value_list))
		

if __name__ == '__main__':
    if args.version == True :
        print(version)
        quit()
    lib2dust_removed_dict = read_dust_filter_file()
    read_fastp_json_fils()
