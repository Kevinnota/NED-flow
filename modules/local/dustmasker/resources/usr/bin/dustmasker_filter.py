#!/usr/bin/env python

import argparse
import re
import gzip
import io

parser = argparse.ArgumentParser(prog='dustmasker filter', description='')
parser.add_argument('--fastq', "-fq", help='input fastq')
parser.add_argument('--dust_threshold', "-dt", help='number of low complexity bases to be removed', default=10, type=int)
parser.add_argument('--dustmasker_acclist', "-d", help='acclist obtained from dustmasker')

parser.add_argument('--version',  help='print version',action='store_true', default=False)
args=parser.parse_args()

version="taxon classifier - v2.1 27 Oct 2024"

def read_dust_text():
    file = open(args.dustmasker_acclist)
    
    read_dict = {}

    for row in file:
        read, start, end = row.strip().split()
        read = re.sub(">", "", read)
        if read not in read_dict:
            read_dict[read]=0
        read_dict[read] += int(end)-int(start)
    
    for read, value in list(read_dict.items()):
        if value <= args.dust_threshold:
            del read_dict[read]

    read_to_filter = set(read_dict.keys())
    return read_to_filter
        
def rewrite_fastq():
    fastq = io.TextIOWrapper(gzip.open(args.fastq, 'rb'))
    fastq_out = gzip.open(re.sub('fastq.gz', 'dust.fastq.gz', args.fastq), 'wt')
    i=0

    while True:
        header = fastq.readline()   # Read the header line (first line of the group)
        library = re.sub('.*:', '', header)
        header_no_RG = re.sub("@", '' ,header.split('\t')[0].strip())
        if not header:  # End of file check
            break
        header = re.sub('\t', ':', header)
        sequence = fastq.readline()  # Read the sequence line (second line)
        plus = fastq.readline()      # Read the plus separator line (third line)
        quality = fastq.readline()   # Read the quality score line (fourth line)

        
        if header_no_RG in read_to_filter:
            i+=1
            continue
        
        fastq_out.write(header)
        fastq_out.write(sequence)
        fastq_out.write(plus)
        fastq_out.write(quality)

    fastq.close()
    fastq_out.close()

    print(f'library = {library}')
    print(f"dustmasker.py version {version}")
    print(f"dustmasker.py dust_threshold {args.dust_threshold}")
    print(f"number sequence removed {i}, all found: {i==len(read_to_filter)}")
    
if __name__ == '__main__':
    if args.version == True :
        print(version)
        quit()
    read_to_filter = read_dust_text()
    rewrite_fastq()
