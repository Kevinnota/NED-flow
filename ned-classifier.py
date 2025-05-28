#!/usr/bin/env python
 
import argparse
from tqdm import tqdm
import simplejson as json
import gzip
import pysam
import re
from scipy import stats
from Bio.Seq import Seq
import os
from scipy.stats import linregress
import numpy as np
from intervaltree import IntervalTree
from collections import Counter
import statistics
from sklearn.linear_model import LinearRegression
import glob
#from concurrent.futures import ThreadPoolExecutor
#import concurrent.futures
from multiprocessing import Pool
import sys
import time
import signal

#trouble shooting
#from pympler.asizeof import asizeof

def parse_tuple(value):
    try:
        # Split the string by commas, convert to float, and return as a tuple
        return tuple(map(float, value.split(',')))
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Invalid tuple format: {value}. Expected format: 0.05,0.05")

parser = argparse.ArgumentParser(prog='Taxon classifier', description='', usage='%(prog)s [options]')
parser.add_argument('--bam', "-b", nargs='+', help='bam(s)')
parser.add_argument('--library_type', "-lt", help='library type, single or double stranded, default SS', choices=['SS', "ss", "DS", "ds"], default='SS')
parser.add_argument('--path_to_references', "-pr", help='path to the references', default='/mnt/sequencedb/bowtie2_RefSeq/800_ncbi_genomes/')
parser.add_argument('--db_summary', "-db", nargs='+', help='database summary files - manual path selection')
parser.add_argument('--db_log', "-dl", help='database summary log file written by nf - automatically selects the correct database summary file')

parser.add_argument('--json', "-j", help='json file with RG')
parser.add_argument('--damage', "-d", help='number of nucleotides to calculate damage on, default is first', type=int, default=3)
parser.add_argument('--damage_threshold', "-dt", help='Damage threshold, (p5,p3) 0.05,0.05', type=parse_tuple, default=(0.05, 0.05))
parser.add_argument('--binomial_ci', help='number of nucleotides to calculate damage on, default is first', action='store_true', default=False)

parser.add_argument('--write_taxonomy_split_bam', '-ws',  help='write bam split files for genus and family if there are more than 100 reads assigned', action='store_true', default=False)
parser.add_argument('--split_bam_out', '-wso', help='split_bam_out', default='split_bams')

parser.add_argument('--taxa', '-tz', help='run the classifier only on a taxonomic level, such as bovids', default=None)
parser.add_argument('--print_r2_table', '-r2', help='write r2 table for each assembly', action='store_true', default=False)
parser.add_argument('--r2_table_out', '-r2o', help='r2 output dir', default='r2_tables')

parser.add_argument('--threads', '-t', help='print version', default=8, type=int)
parser.add_argument('--exclude_assembly', '-e', help='print version',nargs='+', default='')

parser.add_argument('--output', '-o', help='print version', default='summary_data.tsv')
parser.add_argument('--min_distance_sink', "-mds",  help='maximum distance from sink, default 0.90', type=float, default=0.90)
parser.add_argument('--min_distance', "-md",  help='maximum distance from reference, default 0.95', type=float, default=0.95)
parser.add_argument('--max_distance', "-Md",  help='maximum distance from reference, default 1', type=float, default=1.0)
parser.add_argument('--ignore_small_contigs', "-isc",  help='this flag will ignore contigs equal or shorter than set value, default 0', type=int, default=0)
parser.add_argument('--min_read_length', "-ml",  help='min read_length', type=int, default=35)

#tmp flag
parser.add_argument('--bact_sink', "-bs",  help='Whether bact_sink needs was used to map, put False if bacteria were mapped independently (default=False).', action='store_true', default=False)

parser.add_argument('--version',  help='print version',action='store_true', default=False)

args=parser.parse_args()

#version="taxon classifier - v2 15 Aug 2024"
### standard 5% threshold 
version="taxon classifier - v2.1 27 Oct 2024"
    # New, coverage filter Posion distribution^2
    # R2 contig filter 
    # Takes unfiltered files 
    # fcs_report.txt to remove bad contigs
    # r2 filter is done on only contigs that passes filtering.

version="taxon classifier - v2.2 5 Nov 2024"
    # multi threaded bam precessing
    # add mean read length

version="taxon classifier - v2.3 14 Nov 2024"
    # r2 for genus+family assignments
    # can read is log file that is written by nf for getting the right version of the database

    # extra columns:
    #   total number of reads mapping
    #   reads removed because of coverage
    #   reads removed on bad contigs
    #   order_count + phylum_count per assembly

    # it can now run for a family, order, or genus; this is just to get a quick idea of the table, don't think this should be used other than trouble shooting
version="taxon classifier - v2.4 11 March 2025"
    #new columns in the table
    #   order and phylum count per assembly
    #   order deamination stats
    #   order r2 
    #   calculating number of deaminate fragments 
    #   defult 3 nt for deamination rate

version="taxon classifier - v2.5 26 March 2025"
    # can handel bacterial sink

version="taxon classifier - v2.6 06 April 2025"
    # can handel mapped reads from two or more bam directories
    # reads the bam files only onces and stores the relevant informaton for the specific purpose

version="taxon classifier - v2.7 20 April 2025"
    # add bacterial sink column to the table
    # updated assign_taxonomy() function. Lot of unessesary code lines.
    # as soon as a read hit the fcs_reported area in one assembly, the whole read is excluded from taxonomic assignment <- not sure if this is good idea (swithed it of)
    #   rational is that human DNA is common contaminat, so human reads are being removed.

version="taxon classifier - v1 28 May 2025"
    # first version on GitHub

###################
# working on now
#
#   coverage filter between libraries - block overlap filtering - reasoning, unlikely a reads will overlap the same region in multiple samples 
#   its different from coverage filer per libaray - but no sure how its different
#   the problem is that a reads can not overlap any other reads but still be in a similar block as other libraries
#   
#   -> add a column about the total number of reads for a sample?
#   -> maybe the total nmber of reads mapped for the sample?
#   -> and if it's a bacterial sink? Or the DB name?
#
#
#

###################
### To add ####
### - double stranded deamination counts. So look at C->T's and G->A's  

NUCL=set(["A", "T", "C", "G"])

if args.print_r2_table == True:
    if not os.path.exists(args.r2_table_out):
        os.makedirs(args.r2_table_out)

complement = str.maketrans("ACGTacgt", "TGCAtgca")

def make_list_of_bam_files(bams, fcs):
    input_is_dir = False
    bams_list = []
    bam_file_dict={}

    for dir in bams:
        if os.path.isdir(dir):
            input_is_dir=True
            for file in os.listdir(dir):
                if 'bam' in file:
                    bams_list.append(dir+'/'+file)
                 
    if input_is_dir == True:
        bams = bams_list
    
    list_bams_to_preccess=[]
    if args.taxa == None:
        list_bams_to_preccess = bams
    else :    
        for bam in bams:
            assembly=re.sub('.*/|_[a-zA-Z].*', '', bam)
            if args.bact_sink == True:
                if 'concat' in bam:
                    #print(bam)
                    list_bams_to_preccess.append(bam)

            if assembly not in accession_id2taxa:
                continue
            if args.taxa in accession_id2taxa[assembly]:
                list_bams_to_preccess.append(bam)
    
    list_bams_to_preccess_merge=[]
    for bam in list_bams_to_preccess:
        assembly=re.sub('.*/|_[a-zA-Z].*', '', bam)
        if assembly not in bam_file_dict:
            bam_file_dict[assembly]=[]
        bam_file_dict[assembly].append(bam)

    for assembly in bam_file_dict:
        list_bams_to_preccess_merge.append(bam_file_dict[assembly])

    if fcs == True:
        list_of_assemblies=set()
        for bam in list_bams_to_preccess:
            list_of_assemblies.add(re.sub('.*/|_[a-zA-Z].*', '', bam))

        return list_of_assemblies

    else:
        return list_bams_to_preccess_merge

def read_bam_file(bams, for_what):
    bam_dict={}
    total_readcount_dict={}
    HEADER_DONE=False

    if re.search(r'concat|sink', bams[0]):
        type='bact_sink'
    else:
        type='ref_genome'

    for bam in bams:

        assembly=re.sub('.*/|_[a-zA-Z].*', '', bam)
        save = pysam.set_verbosity(0)
        samfile = pysam.AlignmentFile(bam, "rb", check_sq=False)
        pysam.set_verbosity(save)
        if assembly not in total_readcount_dict:
            total_readcount_dict[assembly]={}
        if HEADER_DONE == False and (for_what == 'read_filtering' or for_what == 'r2_stats'):
            header = samfile.header.to_dict()
            HEADER_DONE=True

            contig2lenght={}
            contig2count={}
            genome_length_unfiltered=0

            short_contigs=set()
            for contig in header['SQ']:
                if contig['LN'] <= args.ignore_small_contigs:
                    short_contigs.add(contig['SN'])
                    continue
                contig2lenght[contig['SN']]=contig['LN']
                contig2count[contig['SN']]=0
                genome_length_unfiltered+=contig['LN']

        fcs_removed_reads=set()
        for read in samfile:
            if read.query_length<args.min_read_length:
                continue
            #if read.mapping_quality <20:
            #    continue
            if for_what == 'read_filtering':
                if read.reference_name in short_contigs:
                    continue
            if "I" in read.cigarstring or "D" in read.cigarstring:
                continue
            
            # this removes Ns that where in the reference 
            # its a harsh filter maybe since this would kick out reads for one taxa but keep it for another
            if [ Nuc for Nuc in re.sub('[0-9]{1,}', '', read.get_tag('MD')) if Nuc not in ['A', 'T', 'G', 'C']]:
                continue
            
            # this removes Ns in the read
            if any(base not in NUCL for base in read.query_sequence.upper()):
                continue
            library = re.sub('.*:', '', read.query_name)
            if library not in total_readcount_dict[assembly]:
                total_readcount_dict[assembly][library] = 0
            total_readcount_dict[assembly][library]+=1

            distance=1-(read.get_tag("NM")/len(read.query))

            if type !='bact_sink':
                if not args.min_distance <= distance <= args.max_distance :
                    continue
            elif type =='bact_sink':
                if not args.min_distance_sink <= distance <= args.max_distance :
                    continue
            if assembly in fcs_info_dict:
                if read.query_name in fcs_info_dict[assembly]:
                    continue

            if assembly+read.reference_name in fcs_interval_tree_dict:
                overlapping_ranges=fcs_interval_tree_dict[assembly+read.reference_name].overlap(read.reference_start, read.reference_end)
                if len(overlapping_ranges) != 0 : 
                    fcs_removed_reads.add(read.query_name)
                    continue

            if 'concat' in bam:
                if read.reference_name in fcs_interval_tree_dict:
                    overlapping_ranges=fcs_interval_tree_dict[read.reference_name].overlap(read.reference_start, read.reference_end)
                    if len(overlapping_ranges) != 0 : 
                        fcs_removed_reads.add(read.query_name)
                        continue

            if for_what == 'read_filtering':
                bam_dict[read.query_name] = (read.query_sequence, read.reference_name, read.reference_start, read.reference_end)
            
            elif for_what == 'deamination':
                seq = read.query_sequence
                ref = read.get_reference_sequence()
                if read.is_reverse:
                    ref=Seq(read.get_reference_sequence()).reverse_complement()
                    seq = read.get_forward_sequence()

                bam_dict[read.query_name] = (seq, ref)

            elif for_what == 'r2_stats':
                bam_dict[read.query_name] = (read.query_sequence, read.reference_name)

        samfile.close()
    #print(f"{assembly} nreads={len(total_readcount_dict)}")
    if for_what == 'read_filtering':
        return bam_dict, total_readcount_dict, genome_length_unfiltered, contig2count, contig2lenght, header, type, fcs_removed_reads
    if for_what == 'deamination':
        return bam_dict, type
    if for_what == 'r2_stats':
        return bam_dict, header, type

def reverse_complement(sequence):
    """
    Returns the reverse complement of a DNA sequence.
    """
    return sequence.translate(complement)[::-1]

def read_fcs_reports():
    list_of_assemblies = make_list_of_bam_files(args.bam, True)
    progress=tqdm(total=len(list_of_assemblies), desc='total number of fcs_report read')
    fcs_info_dict={}
    fcs_interval_tree_dict = {}
    i=0
    for assembly in list_of_assemblies:
        #assembly=re.sub('.*/|_[a-zA-Z].*', '', bam)
        file_path = os.path.join(args.path_to_references, "*", assembly, "*fcs_report.txt")
        file = glob.glob(file_path)
        if len(file) == 0:
            progress.update()
            continue
        fcs_report = open(file[0], 'r')

        range_dict={}
        if assembly not in fcs_info_dict:
            fcs_info_dict[assembly] = set()
        for row in fcs_report:
            if "#" in row:
                continue
            seq_id, start_pos, end_pos, seq_len, action, contam_type, coverage, contam_details = row.strip().split('\t')
            seq_id = assembly+seq_id
            if start_pos == "1" and end_pos == seq_len:
                fcs_info_dict[assembly].add(seq_id)
            else:
                if seq_id not in range_dict:
                    range_dict[seq_id]=[]
                range_dict[seq_id].append((start_pos,end_pos))                

        for seq_id in range_dict:
            tree = IntervalTree()
            for start_pos, end_pos in range_dict[seq_id]:
                
                tree.addi(int(start_pos), int(end_pos))
            fcs_interval_tree_dict[seq_id] = tree
       
        progress.update()    
    progress.close()

    if args.bact_sink==True:
        list_of_assemblies = os.listdir(f"{args.path_to_references}/bacteria/")

        progress=tqdm(total=len(list_of_assemblies), desc='total number of fcs_report read back sink')

        for assembly in list_of_assemblies:
            assembly=re.sub('.*/|_[a-zA-Z].*', '', assembly)
            file_path = os.path.join(args.path_to_references, "*", assembly, "*fcs_report.txt")
            file = glob.glob(file_path)
            if len(file) == 0:
                progress.update()
                continue
            fcs_report = open(file[0], 'r')

            range_dict={}
            if assembly not in fcs_info_dict:
                fcs_info_dict[assembly] = set()
            for row in fcs_report:
                if "#" in row:
                    continue
                seq_id, start_pos, end_pos, seq_len, action, contam_type, coverage, contam_details = row.strip().split('\t')
                if start_pos == "1" and end_pos == seq_len:
                    fcs_info_dict[assembly].add(seq_id)
                else:
                    if seq_id not in range_dict:
                        range_dict[seq_id]=[]
                    range_dict[seq_id].append((start_pos,end_pos))                

            for seq_id in range_dict:
                tree = IntervalTree()
                for start_pos, end_pos in range_dict[seq_id]:
                    
                    tree.addi(int(start_pos), int(end_pos))
                fcs_interval_tree_dict[seq_id] = tree
           
            progress.update()    
        progress.close()


    return fcs_info_dict, fcs_interval_tree_dict

def read_database_summary():
    accession_id2taxa={}
    try:
        if args.db_summary != None:
            database_files = args.db_summary
        elif args.db_log != None:
            print(f"database versions for ")
            database_files=[]
            
            db_files = open(args.db_log)
            for row in db_files:
                db, version, date = row.strip().split(' ')
                version = int(re.sub(";", "", version.split('=')[1]))
                print(f"{db}, V{version:04}")
                database_files.append(f"{args.path_to_references}/{db}/database_version_summary/summary_version_{version:04}.tsv")
            db_files.close()

            if args.bact_sink == True:
                version=2
                database_files.append(f"{args.path_to_references}/bacteria/database_version_summary/summary_version_{version:04}.tsv")
        pass
    
    except:
        print(f"database file not found, either specify -db or -dl")
        quit()

    #print(database_files)
    for file in database_files:
        summary_file = open(file)
        db=re.sub('.*/', '', re.sub('/database_version_summary.*','', file))
        for row in summary_file:
            assembly_accession, taxa_name, phylum_name, order_name, family_name, genus_name, taxid, Assembly_level, Genome_representation = row.strip().split('\t')
            accession_id2taxa[assembly_accession] = [taxid, taxa_name, phylum_name, order_name, family_name, genus_name, db]
            
            if 'bacteria' in file and args.bact_sink == True :
                if phylum_name in accession_id2taxa:
                    continue
                accession_id2taxa[phylum_name] = ['bact_sink', phylum_name, phylum_name, phylum_name, phylum_name, phylum_name, db]
               
    return accession_id2taxa

def read_RG_json():
    with gzip.open(args.json, 'r') as json_file:
        read_dict = json.load(json_file)
    return read_dict

def process_bam_files(bam):
    assembly_not_reported=None

    assembly=re.sub('.*/|_[a-zA-Z].*', '', bam[0])
    

    contigs_to_ignore=set()
    read2taxon_dict={}
    read_dict={}
    total_95_readcount_dict={}


    total_95_readcount_dict[assembly]={}
    read_count_removed_contigs={}
    read_count_removed_contigs[assembly]={}

    read_count_high_coverage={}
    read_count_high_coverage[assembly]={}

    read_ranges={}

    if assembly not in accession_id2taxa:
        total_readcount_dict={}
        genome_length=0
        genome_length_unfiltered = 0
        assembly_not_reported = assembly
        fcs_removed_reads=set()
        return read2taxon_dict, read_dict, total_95_readcount_dict, contigs_to_ignore, assembly, total_readcount_dict, read_count_removed_contigs, read_count_high_coverage, genome_length, genome_length_unfiltered, assembly_not_reported, fcs_removed_reads   
    
    taxid, taxa_name, phylum_name, order_name, family_name, genus_name, db = accession_id2taxa[assembly]

    bam_dict, total_readcount_dict, genome_length_unfiltered, contig2count, contig2lenght, header, type, fcs_removed_reads = read_bam_file(bam, 'read_filtering')

    total_bases_sequenced=0
    mean_readlength=[]
    read_count={}

    read_to_exclude_coverage=set()
    list_contigs_removed_elbow=set()
    list_contigs_removed_ratio=set()
    
    #total_readcount_dict[assembly]={}
    if assembly not in accession_id2taxa:
        genome_length_unfiltered = 0 
        genome_length=0
        assembly_not_reported = assembly
        return read2taxon_dict, read_dict, total_95_readcount_dict, contigs_to_ignore, assembly, total_readcount_dict, read_count_removed_contigs, read_count_high_coverage, genome_length, genome_length_unfiltered, assembly_not_reported, fcs_removed_reads
 
    if ( type != 'bact_sink' ) and ( db != 'bacteria' ) :

        total_bases_sequenced_dict={}
        for read in bam_dict:

            query_sequence, reference_name, reference_start, reference_end = bam_dict[read]

            library = re.sub('.*:', '', read)

            if library not in read_count:
                read_count[library]=0
            read_count[library]+=1
            #contig2count[reference_name]+=1
            
            if library not in total_bases_sequenced_dict:
                total_bases_sequenced_dict[library] = 0
            total_bases_sequenced_dict[library]+=len(query_sequence)

            mean_readlength.append(len(query_sequence))
            if reference_name not in read_ranges:
                read_ranges[reference_name]={}
            # add library after read.reference_name
            if library not in read_ranges[reference_name]:
                read_ranges[reference_name][library]=[]
            read_ranges[reference_name][library].append((reference_start, reference_end))

        #print('genome_length=', genome_length, 'mean_read_length=', statistics.mean(mean_readlength),'total_bases_sequenced=', total_bases_sequenced)
        if genome_length_unfiltered == 0:
            #print(assembly, 'genome length == 0')
            genome_length=0
            assembly_not_reported = assembly
            return read2taxon_dict, read_dict, total_95_readcount_dict, contigs_to_ignore, assembly, total_readcount_dict, read_count_removed_contigs, read_count_high_coverage, genome_length,genome_length_unfiltered, assembly_not_reported, fcs_removed_reads
       
        
        ########################
        #  Coverage filtering  #
        ########################

        np.random.seed(42)
        thresholds_dict={}

        for library in read_count:
            if read_count[library] != 0:
                coverage=total_bases_sequenced_dict[library]/genome_length_unfiltered
                
                poisson_distribution = Counter(np.random.poisson(lam=(coverage+coverage), size=read_count[library]))
                sorted_keys = sorted(poisson_distribution.keys())
                cumulative_counts = np.cumsum([poisson_distribution[key] for key in sorted_keys])
                normalized_cumsum=cumulative_counts/read_count[library]
                index = np.searchsorted(normalized_cumsum, 0.99)
                max_overlap = sorted_keys[index] ## calculates the coverage threshold
                thresholds_dict.update({library:max_overlap})
            else:

                max_overlap = 0
                thresholds_dict.update({library:max_overlap})
            #print(f'asslebly {assembly}; library:max_overlap; {library}:{max_overlap}')
        
        interval_trees = {}

        for reference_name in read_ranges:
            for library, ranges in read_ranges[reference_name].items():
            
                sorted_intervals = sorted(ranges, key=lambda x: (x[0], x[1]))
                interval_with_overlap=set()
            
                for i in range(len(sorted_intervals)):
                    interval_with_overlap.add(sorted_intervals[i]) # its adding all ranges, its actually faster than the speed up I did. 
      
                tree = IntervalTree()
                for min_val, max_val in interval_with_overlap:
                    tree.addi(min_val, max_val)
                if reference_name not in interval_trees:
                    interval_trees[reference_name] = {}
                interval_trees[reference_name][library] = tree

        for reference_name in contig2count: # making sure all contigs are 0
            contig2count[reference_name] = 0

        error_given_before=False
        for read in bam_dict:
            query_sequence, reference_name, reference_start, reference_end = bam_dict[read]

            library = re.sub('.*:', '', read)
            # looking for overlapping reads
            #if reference_name in interval_trees:
            overlapping_ranges = interval_trees[reference_name][library].overlap(reference_start-1, reference_end+1)
            count = len(overlapping_ranges)-1# Subtract 1 to account for the read itself

            if count <= thresholds_dict[library]:
                try:
                    contig2count[reference_name]+=1
                except:
                    if error_given_before == False:
                        error_given_before=True
                        print(f"ERROR - assembly:{assembly} has multiple ref genome versions")
                        print(f"This probably means one of the two or more runs was done to a different version - so UPDATE")
                        #Pool.terminate()
            else:
                read_to_exclude_coverage.add(read)
        
        ### r2 filter should be after covering filter, if some regions have tons of data, could change r2 even though all of it is from a few reagions
        
        if sum(contig2count.values()) == 0:
            genome_length=0
            assembly_not_reported = assembly
            return read2taxon_dict, read_dict, total_95_readcount_dict, contigs_to_ignore, assembly, total_readcount_dict, read_count_removed_contigs, read_count_high_coverage, genome_length, genome_length_unfiltered, assembly_not_reported, fcs_removed_reads
        
        ########################
        #  Contig filtering    #
        ########################

        values = list(contig2lenght.values())
        int_values = [int(value) for value in values]
        x = np.array(int_values).reshape((-1, 1))
        sorted_indices = np.argsort(x.flatten())

        values = list(contig2count.values())
        int_values = [int(value) for value in values]
        y = np.array(int_values)    

        contig_sorted = np.array(list(contig2lenght.keys()))[sorted_indices]
        x_sorted = x[sorted_indices]
        y_sorted = y[sorted_indices]
        

        i=0
       
        model = LinearRegression(fit_intercept=False).fit(x, y)
        y_pred = model.predict(x_sorted)
        
        if sum(int_values) != 0:
            ratio = [y_sorted[i] / y_pred[i] for i in range(len(y_pred))]
        
        ratio_cumsum = np.cumsum(ratio, axis=0)
        prop_ratio = ratio_cumsum/sum(ratio)
        
        elbow_length_filter = x_sorted[find_elbow(x_sorted, prop_ratio)][0]
        contig2ratio={}
        r=0
        elbow=True
        
        if elbow_length_filter > 2000000:  # threshold set only if the elbow point is smaller than 2M
            elbow=False
        
        #print(f'asslebly {assembly}; elbow_length_filter= {elbow_length_filter}')
        elbow_filtered_ratios = []
        elbow_contig_filtered = []
        for contig in contig_sorted:
            if elbow==True and contig2lenght[contig] <= elbow_length_filter :
                r+=1
                continue    
            contig2ratio[contig]=ratio[r]
            elbow_filtered_ratios.append(ratio[r])
            elbow_contig_filtered.append(contig2lenght[contig])
            r+=1

        elbow_filtered_ratios = np.array(elbow_filtered_ratios)
        elbow_contig_filtered = np.array(elbow_contig_filtered)
        
        std = np.std(elbow_filtered_ratios)
        if sum(elbow_contig_filtered)==0:
            genome_length=0
            assembly_not_reported = assembly
            return read2taxon_dict, read_dict, total_95_readcount_dict, contigs_to_ignore, assembly, total_readcount_dict, read_count_removed_contigs, read_count_high_coverage, genome_length, genome_length_unfiltered, assembly_not_reported, fcs_removed_reads
        
        ratio_threshold = np.average(elbow_filtered_ratios, weights=elbow_contig_filtered)+(3*std)
            
        taxa_name = re.sub(r"\'|\.", "", re.sub(" ", "_", taxa_name))
        
        for X, Y, contig, RATIO, in zip(x_sorted, y_sorted, contig_sorted, ratio):
            if elbow==True and (X[0] <= elbow_length_filter) :
                list_contigs_removed_elbow.add(contig)

            elif RATIO >= ratio_threshold :
                list_contigs_removed_ratio.add(contig)

    genome_length=0
    i=0
    short_contigs=set()

    for contig in header['SQ']:
        try:
            if elbow==True and contig['LN'] <= elbow_length_filter :
                continue
            if contig['SN'] in list_contigs_removed_ratio or contig['SN'] in list_contigs_removed_elbow:
                continue
            genome_length+=contig['LN']
        except:
            pass

    for read in bam_dict:
        query_sequence, reference_name, reference_start, reference_end = bam_dict[read]

        library = re.sub('.*:', '', read)

        if read in read_to_exclude_coverage:
            if library not in read_count_high_coverage[assembly]:
                read_count_high_coverage[assembly][library]=0
            read_count_high_coverage[assembly][library]+=1
            continue

        if reference_name in list_contigs_removed_elbow:
            #print(read.reference_name)
            if library not in read_count_removed_contigs[assembly]:
                read_count_removed_contigs[assembly][library]=0
            read_count_removed_contigs[assembly][library]+=1
            continue

        if reference_name in list_contigs_removed_ratio:
            if library not in read_count_removed_contigs[assembly]:
                read_count_removed_contigs[assembly][library]=0
            read_count_removed_contigs[assembly][library]+=1
            continue

        read_dict[read]=library
        if library not in total_95_readcount_dict[assembly]:
            total_95_readcount_dict[assembly][library]=0
        total_95_readcount_dict[assembly][library]+=1
        if read not in read2taxon_dict:
            read2taxon_dict[read]=[]
        read2taxon_dict[read].append(assembly)

    #samfile.close()

    contigs_to_ignore=list_contigs_removed_ratio|list_contigs_removed_elbow
    
    try:
        os.makedirs(f"bams_gzset")
    except:
        pass
    
    file = gzip.open(f"bams_gzset/{assembly}.set.gz", 'wt')
    for item in list_contigs_removed_ratio:
        file.write(f"{item}\tR\n")
    for item in list_contigs_removed_elbow:
        file.write(f"{item}\tE\n")
        
    file.close()
    
    if len(read_dict) < 100:
        assembly_not_reported = assembly
    
    contigs_to_ignore=set()

    return read2taxon_dict, read_dict, total_95_readcount_dict, contigs_to_ignore, assembly, total_readcount_dict, read_count_removed_contigs, read_count_high_coverage, genome_length, genome_length_unfiltered, assembly_not_reported, fcs_removed_reads


def filter_bam_files():
    
    list_bams_to_preccess = make_list_of_bam_files(args.bam, False)

    
    with Pool(processes=args.threads) as pool:
        results = list(tqdm(pool.imap_unordered(process_bam_files, list_bams_to_preccess, chunksize=1), 
                        total=len(list_bams_to_preccess), 
                        desc="Bam files processed"))
        
    # Aggregate results from all processes
    read2taxon_dict, read_dict, total_95_readcount_dict, contigs_to_ignore, total_readcount_dict, read_count_removed_contigs, read_count_high_coverage, assembly2genome_lenght, assembly2genome_lenght_unflitered, assembly_not_reported, fcs_removed_reads = {}, {}, {}, set(), {}, {}, {}, {}, {}, set(), set()
    assembly_list = []
    progress=tqdm(total=len(list_bams_to_preccess), desc='Combining results')

    for result in results:
        r2taxon, rdict, trcount, cignore, assembly, ctotal_readcount_dict, cread_count_removed_contigs, cread_count_high_coverage, genome_length, genome_length_unfiltered, cassembly_not_reported, cfcs_removed_reads = result
        for read in r2taxon:
            if read not in read2taxon_dict:
                read2taxon_dict[read] = []
            read2taxon_dict[read].append(assembly)
        
        assembly_list.append(assembly)
        read_dict.update(rdict)
        total_95_readcount_dict.update(trcount)
        contigs_to_ignore.update(cignore)
        total_readcount_dict.update(ctotal_readcount_dict)
        read_count_removed_contigs.update(cread_count_removed_contigs)
        read_count_high_coverage.update(cread_count_high_coverage)
        assembly2genome_lenght.update({assembly:genome_length})
        assembly2genome_lenght_unflitered.update({assembly:genome_length_unfiltered})
        assembly_not_reported.add(cassembly_not_reported)
        fcs_removed_reads.update(cfcs_removed_reads)
        
        progress.update()
    progress.close()
    
    return read2taxon_dict, read_dict, total_95_readcount_dict, contigs_to_ignore, total_readcount_dict, read_count_removed_contigs, read_count_high_coverage, assembly2genome_lenght, assembly2genome_lenght_unflitered, assembly_not_reported, assembly_list, fcs_removed_reads

def assign_taxonomy():
    taxon_dict={'genus':{},
                'family':{},
                'order':{},
                'phylum':{},
                'back_sink':{}}
    i=0
    assembly_mapped_dict={}
    library_list=set()
    progress=tqdm(total=len(read2taxon_dict), desc='Reads taxonomy assigned to')
    print_error=True
    for read in read2taxon_dict:
        #if read in fcs_removed_reads:
        #    continue

        phylum_list=set()
        order_list=set()
        family_list=set()
        genus_list=set()
        bact_sink_list=set()

        library = read_dict[read]
        library_list.add(library)
        if library not in taxon_dict['genus']:
            taxon_dict['genus'][library]={}
        if library not in taxon_dict['family']:
            taxon_dict['family'][library]={}
        if library not in taxon_dict['order']:
            taxon_dict['order'][library]={}
        if library not in taxon_dict['phylum']:
            taxon_dict['phylum'][library]={}
        if library not in taxon_dict['back_sink']:
            taxon_dict['back_sink'][library]={}
        
        assembly_list = []
        phylum_count_list, order_count_list = [], []


        for assembly in read2taxon_dict[read]:
            if assembly not in accession_id2taxa:
                if print_error == True:
                    print_error = False
                continue
            
            taxid, taxa_name, phylum_name, order_name, family_name, genus_name, db = accession_id2taxa[assembly]
            
            if taxid=='bact_sink':
                bact_sink_list.add(taxa_name)
            else:
                bact_sink_list.add(phylum_name)

            phylum_list.add(phylum_name)
            phylum_count_list.append(phylum_name)
            order_list.add(order_name)
            order_count_list.append(order_name)
            
            family_list.add(family_name)
            genus_list.add(genus_name)
            assembly_list.append(assembly)

        for assembly in assembly_list:
            if assembly not in assembly_mapped_dict.keys():
                assembly_mapped_dict[assembly]={}
            if library not in assembly_mapped_dict[assembly]:
                assembly_mapped_dict[assembly][library] = { 'genus':set(),
                                                            'family':set(),
                                                            'order':set(),
                                                            'phylum':set(),
                                                            'bact_sink':set()}

        if len(genus_list) == 1:
            genus_list = next(iter(genus_list)) # since the genus_list is a set() list, next(inter) is used to get the actual value
            if genus_list not in taxon_dict['genus'][library]:
                taxon_dict['genus'][library][genus_list] = 0
            taxon_dict['genus'][library][genus_list]+=1
            
            for assembly in assembly_list:
            #   if assembly not in assembly_mapped_dict.keys():
            #        assembly_mapped_dict[assembly]={}
            #    if library not in assembly_mapped_dict[assembly]:
            #        assembly_mapped_dict[assembly][library] = { 'genus':set(),
            #                                                    'family':set(),
            #                                                    'order':set(),
            #                                                    'phylum':set()}
                assembly_mapped_dict[assembly][library]['genus'].add(read)
    
        elif len(family_list) == 1:
            family_list = next(iter(family_list))
            if family_list not in taxon_dict['family'][library]:
                taxon_dict['family'][library][family_list] = 0
            taxon_dict['family'][library][family_list]+=1
                    
            for assembly in assembly_list:
                #if assembly not in assembly_mapped_dict.keys():
                #    assembly_mapped_dict[assembly]={}
                #if library not in assembly_mapped_dict[assembly]:
                #    assembly_mapped_dict[assembly][library] = {'genus':set(),
                #                                               'family':set(),
                #                                                'order':set(),
                #                                                'phylum':set()}
                assembly_mapped_dict[assembly][library]['family'].add(read)

                i+=1
    
        elif len(order_list) == 1:
            order_list = next(iter(order_list))
            if order_list not in taxon_dict['order'][library]:
                taxon_dict['order'][library][order_list] = 0
            taxon_dict['order'][library][order_list]+=1
            
            for assembly in assembly_list:
                #if assembly not in assembly_mapped_dict.keys():
                #    assembly_mapped_dict[assembly]={}
                #if library not in assembly_mapped_dict[assembly]:
                #    assembly_mapped_dict[assembly][library] = {'genus':set(),
                #                                               'family':set(),
                #                                                'order':set(),
                #                                                'phylum':set()}
                assembly_mapped_dict[assembly][library]['order'].add(read)

        #if( len(order_list) == 2 ):
            #order_count = Counter(order_count_list)
            #phylum_count = Counter(phylum_count_list)
            #print(read)
            #print(f"order_count {order_count}")
            #print(f"phylum_count {phylum_count}")
            #print(f'order assemies {read2taxon_dict[read]}')
            #print('\n')


        elif len(phylum_list) == 1:
            phylum_list = next(iter(phylum_list))
            if phylum_list not in taxon_dict['phylum'][library]:
                taxon_dict['phylum'][library][phylum_list] = 0
            taxon_dict['phylum'][library][phylum_list]+=1
            
            for assembly in assembly_list:
                #if assembly not in assembly_mapped_dict.keys():
                #    assembly_mapped_dict[assembly]={}
                #if library not in assembly_mapped_dict[assembly]:
                #    assembly_mapped_dict[assembly][library] = {'genus':set(),
                #                                               'family':set(),
                #                                                'order':set(),
                #                                                'phylum':set()}
                assembly_mapped_dict[assembly][library]['phylum'].add(read)            

        # reads only the total number of bacteria that mapped to the reference genome
        elif len(bact_sink_list) != 0 :
            if 'bacteria' not in taxon_dict['back_sink'][library]:
                taxon_dict['back_sink'][library]['bacteria'] = 0
            taxon_dict['back_sink'][library]['bacteria']+=1

            for assembly in assembly_list:
                assembly_mapped_dict[assembly][library]['bact_sink'].add(read)
            
            #print(f'phylum_names = {bact_sink_list}')
            #print(f'taxon_dict = {taxon_dict}')
            #quit()
        
        #phylum_count = Counter(phylum_count_list)
        #print(f"phylum_count {phylum_count}")
        

        progress.update()
    #print(f"assembly {assembly}; assembly_mapped_dict {assembly_mapped_dict}")

    return taxon_dict, assembly_mapped_dict, library_list #read2taxa

def get_deamination_stat(bam):
    deamination_count={}
    
    #save = pysam.set_verbosity(0)
    #samfile = pysam.AlignmentFile(bam, "rb", check_sq=False)
    #pysam.set_verbosity(save)
    assembly=re.sub('.*/|_[a-zA-Z].*', '', bam[0])
    if assembly in args.exclude_assembly:
        return 
    i=0
    read2deam={}

    genus_read_length={assembly:{}}
    family_read_length={assembly:{}}
    order_read_length={assembly:{}}
    phylum_read_length={assembly:{}}

    bam_dict, type = read_bam_file(bam, 'deamination')
    
    if len(bam_dict) == 0:
        return
    if  type == 'bact_sink' :
        return 

    for read in bam_dict:
        
        seq, ref = bam_dict[read]
            
        c_in_5, deam_5, c_in_3, deam_3 = [False] * 4
        # count c's

        if args.library_type.upper() == 'SS':

            if "C" in ref[0:args.damage].upper():
                c_in_5 =  True 
                for ref_nt, seq_nt, in zip(ref[0:args.damage:].upper(), seq[0:args.damage:].upper()):
                    if ref_nt == "C" and seq_nt == 'T':
                        deam_5=True

            if "C" in ref[-args.damage:].upper():
                c_in_3 =  True 
                for ref_nt, seq_nt, in zip(ref[-args.damage:].upper(), seq[-args.damage:].upper()):
                    if ref_nt == "C" and seq_nt == 'T':
                        deam_3=True
        
        if args.library_type.upper() == 'DS' :
            if "G" in ref[0:args.damage].upper():
                c_in_5 =  True 
                for ref_nt, seq_nt, in zip(ref[0:args.damage:].upper(), seq[0:args.damage:].upper()):
                    if ref_nt == "C" and seq_nt == 'T':
                        deam_5=True

            if "C" in ref[-args.damage:].upper():
                c_in_3 =  True 
                for ref_nt, seq_nt, in zip(ref[-args.damage:].upper(), seq[-args.damage:].upper()):
                    if ref_nt == "G" and seq_nt == 'A':
                        deam_3=True            

        read2deam[read] = [c_in_5, deam_5, c_in_3, deam_3]
        
        ## get median read lengths

        library = re.sub('.*:', '', read)
        if assembly not in assembly_mapped_dict:
            continue
        if library not in assembly_mapped_dict[assembly]:
            continue
        
        if read in assembly_mapped_dict[assembly][library]['genus']:
            if library not in genus_read_length[assembly]:
                genus_read_length[assembly][library]=[]
            genus_read_length[assembly][library].append(len(read))

        if read in assembly_mapped_dict[assembly][library]['family']:
            if library not in family_read_length[assembly]:
                family_read_length[assembly][library]=[]
            family_read_length[assembly][library].append(len(read))

        if read in assembly_mapped_dict[assembly][library]['order']:
            if library not in order_read_length[assembly]:
                order_read_length[assembly][library]=[]
            order_read_length[assembly][library].append(len(read))


    for assembly in genus_read_length:
        for library in genus_read_length[assembly]:
            genus_read_length[assembly][library]=str(f"{np.median(genus_read_length[assembly][library]):.0f}")
  
    for assembly in family_read_length:
        for library in family_read_length[assembly]:
            family_read_length[assembly][library]=str(f"{np.median(family_read_length[assembly][library]):.0f}")
    
    for assembly in order_read_length:
        for library in order_read_length[assembly]:
            order_read_length[assembly][library]=str(f"{np.median(order_read_length[assembly][library]):.0f}")
    

    #progress.update()
    #samfile.close()
    library_l=set()
    #calculate deam stats
    deamination_count[assembly]={}
    for read_name in read2deam:
        library = re.sub('.*:', '', read_name)
        if assembly not in assembly_mapped_dict:
            continue
        if library not in assembly_mapped_dict[assembly]:
            continue
        
        if read_name in assembly_mapped_dict[assembly][library]['genus']:
            assignment='genus'

        elif read_name in assembly_mapped_dict[assembly][library]['family']:
            assignment='family'
        
        elif read_name in assembly_mapped_dict[assembly][library]['order']:
            assignment='order'
        
        
        else:
            continue     
            
        if library not in deamination_count[assembly]:
            deamination_count[assembly][library]={'genus':{}, 'family':{}, 'order':{}}
            for taxa_level in ['genus', 'family', 'order']:
                deamination_count[assembly][library][taxa_level]={'sum_5p':0, 
                                                                'sum_3p':0,
                                                                'deam_5p':0,
                                                                'deam_3p':0,
                                                                'cond_3_sum_5p':0,
                                                                'cond_5_sum_3p':0,
                                                                'con_deam':0,
                                                                'deam_count':0}
                    
        if read2deam[read_name][0]==True:
            deamination_count[assembly][library][assignment]['sum_5p']+=1
        if read2deam[read_name][2]==True:
            deamination_count[assembly][library][assignment]['sum_3p']+=1
        if read2deam[read_name][1]==True:
            deamination_count[assembly][library][assignment]['deam_5p']+=1
        if read2deam[read_name][3]==True:
            deamination_count[assembly][library][assignment]['deam_3p']+=1

        if read2deam[read_name][1]==True | read2deam[read_name][3]==True:
            deamination_count[assembly][library][assignment]['deam_count']+=1

            # 5p deaminated, and C in p3            
        if read2deam[read_name][1]==True and read2deam[read_name][2] == True:
            deamination_count[assembly][library][assignment]['cond_5_sum_3p']+=1
            # 3p deaminated, and C in p5            
        if read2deam[read_name][3]==True and read2deam[read_name][0] == True:
            deamination_count[assembly][library][assignment]['cond_3_sum_5p']+=1
            # 5p deaminated, p3 deaminated           
        if read2deam[read_name][1]==True and read2deam[read_name][3] == True:
            deamination_count[assembly][library][assignment]['con_deam']+=1     
    
    return deamination_count, genus_read_length, family_read_length, order_read_length

def deamination_stat_multi_thread():
    
    #list_bams_to_preccess=[]
    #if args.taxa == None:
    #    for bam in args.bam:
    #        assembly=re.sub('.*/|_[a-zA-Z].*', '', bam)
    #        if assembly in assembly_not_reported:
    #            continue            
    #        list_bams_to_preccess.append(bam)
        
    #else :
    #    for bam in args.bam:
    #        assembly=re.sub('.*/|_[a-zA-Z].*', '', bam)
    #        if assembly not in accession_id2taxa:
    #            continue
    #        if assembly in assembly_not_reported:
    #            continue
    #        if args.taxa in accession_id2taxa[assembly]:
    #            list_bams_to_preccess.append(bam)

    list_bams_to_preccess = make_list_of_bam_files(args.bam, False)

    with Pool(processes=args.threads) as pool:
        results = list(tqdm(pool.imap_unordered(get_deamination_stat, list_bams_to_preccess), total=len(list_bams_to_preccess), desc="Bam files processed deamination"))

    # Aggregate results from all processes
    c_deamination_count, c_genus_read_length, c_family_read_length, c_order_read_length = {}, {}, {}, {}

    for result in results:  
        if result == None:
            continue
        deamination_count, genus_read_length, family_read_length, order_read_length = result
        c_deamination_count.update(deamination_count)
        c_genus_read_length.update(genus_read_length)
        c_family_read_length.update(family_read_length)
        c_order_read_length.update(order_read_length)

    return c_deamination_count, c_genus_read_length, c_family_read_length, c_order_read_length

def binomial_ci(x, n, alpha=0.05):
    #Clopper-Pearson interval for binomial distribution
    #x is number of successes, n is number of trials, alpha the percent
    if n==0:
        return ["N/A","N/A"]
    lower = stats.beta.interval(1-alpha, x,n-x+1)[0] if x != 0 else 0
    upper = stats.beta.interval(1-alpha, x+1,n-x)[1] if x != n else 1
    conf_int = [f"{round(lower,3):.3f}", f"{round(upper,3):.3f}"]
    return conf_int

def count_total_number_taxon():

    taxa_count_dict={}
    for assembly in total_95_readcount_dict:
        if assembly not in accession_id2taxa:
            continue
        taxid, taxa_name, phylum_name, order_name, family_name, genus_name, db = accession_id2taxa[assembly]
        for library in total_95_readcount_dict[assembly]:
            #print(library)
            if library not in taxa_count_dict:
                taxa_count_dict[library]={}
            
            if genus_name not in taxa_count_dict[library].keys():
                taxa_count_dict[library][genus_name]=set()

            if family_name not in taxa_count_dict[library].keys():
                taxa_count_dict[library][family_name]=set()
            
            if order_name not in taxa_count_dict[library].keys():
                taxa_count_dict[library][order_name]=set()

            if phylum_name not in taxa_count_dict[library].keys():
                taxa_count_dict[library][phylum_name]=set()


            if assembly not in assembly_mapped_dict:
                continue
            if library not in assembly_mapped_dict[assembly]:
                continue

            taxa_count_dict[library][genus_name].update(assembly_mapped_dict[assembly][library]['genus'])
            taxa_count_dict[library][family_name].update(assembly_mapped_dict[assembly][library]['family'])
            taxa_count_dict[library][order_name].update(assembly_mapped_dict[assembly][library]['order'])
            taxa_count_dict[library][phylum_name].update(assembly_mapped_dict[assembly][library]['phylum'])

    summary_taxa_count_dict={}
    for library in taxa_count_dict:
        if library not in summary_taxa_count_dict:
            summary_taxa_count_dict[library]={}
        for taxa in taxa_count_dict[library]:
            summary_taxa_count_dict[library][taxa]=len(taxa_count_dict[library][taxa])
    
    return summary_taxa_count_dict

def get_R2_values(bam):
    R2_values={}
    #save = pysam.set_verbosity(0)
    #samfile = pysam.AlignmentFile(bam, "rb", check_sq=False)
    #pysam.set_verbosity(save)
    assembly=re.sub('.*/|_[a-zA-Z].*', '', bam[0])

    try:
        taxid, taxa_name, phylum_name, order_name, family_name, genus_name, db = accession_id2taxa[assembly]
    except:
        return 

    bam_dict, header, type = read_bam_file(bam, 'r2_stats')
    if len(bam_dict) == 0:
        return
    if ( type == 'bact_sink' ) or ( db == 'bacteria' ) :
        return 

    #header = samfile.header.to_dict()
    #print(assembly)
    gzset=gzip.open(f"bams_gzset/{assembly}.set.gz", 'rt')
    contigs_to_ignore = set(contg.strip().split('\t')[0] for contg in gzset)
    gzset.close()

    filter_dict={}
    gzset=gzip.open(f"bams_gzset/{assembly}.set.gz", 'rt')
    for row in gzset:
        contig, filter= row.strip().split('\t')
        filter_dict[contig]=filter
    gzset.close()
    contig_to_lenght={}
    contig_to_lenght_unfiltered = {}

    for contig in header['SQ']:
        if contig['SN'] in contigs_to_ignore:
            contig_to_lenght_unfiltered[contig['SN']]=contig['LN'] 
            continue
        contig_to_lenght[contig['SN']]=contig['LN'] 
        contig_to_lenght_unfiltered[contig['SN']]=contig['LN']

    contig_to_count={}
    library_list=set()

    for read in bam_dict:
        #if read.query_length<args.min_read_length:
        #    continue
        #if read.mapping_quality <20:
        #    continue
        query_sequence, reference_name = bam_dict[read]
        library = re.sub('.*:', '', read)
        library_list.add(library)
        if assembly not in assembly_mapped_dict:
            continue
        if library not in assembly_mapped_dict[assembly]:
            continue
        
        
        if library not in contig_to_count:
            contig_to_count[library]={}
        if 'all' not in contig_to_count[library]:
            contig_to_count[library]['all']={}
        if reference_name not in contig_to_count[library]['all']:
            contig_to_count[library]['all'][reference_name]=0

        contig_to_count[library]['all'][reference_name]+=1

        if reference_name in contigs_to_ignore:
            continue

        if read in assembly_mapped_dict[assembly][library]['genus']:
            assignment='genus'

        elif read in assembly_mapped_dict[assembly][library]['family']:
            assignment='family'

        elif read in assembly_mapped_dict[assembly][library]['order']:
            assignment='order'

        else:
            continue     

        if assignment not in contig_to_count[library]:
            contig_to_count[library][assignment]={}

        if reference_name not in contig_to_count[library][assignment]:
            contig_to_count[library][assignment][reference_name]=0
        contig_to_count[library][assignment][reference_name]+=1
    library_list_100=set()
    
    
    for library in library_list:
     
        contig_count_dict={}
        for assignment in ['all', 'genus', 'family', 'order']:
            contig_length=[]
            contig_count=[]
            
            if library not in contig_to_count:
                continue
            if assignment not in contig_to_count[library]:
                continue

            if sum(list(contig_to_count[library][assignment].values()))<100:
                continue

            library_list_100.add(library+':'+assignment)
            #print(f"library: {library}; assignment: {assignment};  numer of reads : {sum(list(contig_to_count[library][assignment].values()))}")
            if assignment == 'all':
                continue

            for contig in contig_to_lenght:
                contig_length.append(contig_to_lenght[contig])
                
                if contig not in contig_to_count[library][assignment]:
                    contig_count.append(0)
                else:
                    contig_count.append(contig_to_count[library][assignment][contig])
            
            x = np.array(contig_length).reshape((-1, 1))
            y = np.array(contig_count)    
            model = LinearRegression(fit_intercept=False).fit(x, y)
            r_squared = model.score(x, y)

            #print(f"assembly = {assembly}; library = {library}; r_squared = {r_squared}")
            try:
                if sum(contig_count) < 100:
                    raise ValueError
                x = np.array(contig_length).reshape((-1, 1))
                y = np.array(contig_count)    
                model = LinearRegression(fit_intercept=False).fit(x, y)
                r_squared = model.score(x, y)
                #slope, intercept, r_value, p_value, std_err = linregress(contig_length, contig_count)
                #r_squared = r_value ** 2
                #if p_value >= 0.05:
                #    r_squared=0
            
            except:
                r_squared=-1
                #p_value=1

            if assembly not in R2_values:
                R2_values[assembly]={}
            if library not in R2_values[assembly]:
                R2_values[assembly][library]={}
            if  assignment not in R2_values[assembly][library]:
                R2_values[assembly][library][assignment]=r_squared

            contig_count_dict[assignment]=contig_count
        
        try:
            contig_count = [a + b for a, b in zip(contig_count_dict['genus'], contig_count_dict['family'])]
            if sum(contig_count) < 100:
                continue
            x = np.array(contig_length).reshape((-1, 1))
            y = np.array(contig_count) 
            model = LinearRegression(fit_intercept=False).fit(x, y)
            r_squared = model.score(x, y)
            #slope, intercept, r_value, p_value, std_err = linregress(contig_length, contig_count)
            #r_squared = r_value ** 2
            #if p_value >= 0.05:
            #    r_squared=0
        except:
            r_squared=-1
            #p_value=-1

        if assembly not in R2_values:
            R2_values[assembly] = {}
        if library not in R2_values[assembly]:
            R2_values[assembly][library]={}
        R2_values[assembly][library]['genus+family']=r_squared

    try:
        taxid, taxa_name, phylum_name, order_name, family_name, genus_name, db = accession_id2taxa[assembly]
    except:
        return

    if args.print_r2_table == True:
        taxa_name = re.sub(" ", "_", taxa_name)
        file = gzip.open(f"{args.r2_table_out}/{taxa_name}_{assembly}_r2.tsv.gz", 'wt')

        library_list_100 = sorted(library_list_100)
        
        updated_list = [re.sub(":", "_", item) for item in library_list_100]

        file.write('\t'.join(['contig', 'contig_length', 'filtered', 'filter_class']+list(updated_list))+'\n')

        
        for contig in contig_to_lenght_unfiltered:
            line=[contig, contig_to_lenght_unfiltered[contig]]
            if contig in contigs_to_ignore:
                line.append("TRUE")
                line.append(filter_dict[contig])
            else:
                line.append("FALSE")
                line.append("FALSE")

            for item in library_list_100:
                library, assignment = item.split(":", 1)
                try :
                    line.append(contig_to_count[library][assignment][contig])
                except :
                    line.append(0)
            line = [str(item) for item in line]

            file.write('\t'.join(line)+'\n')
        file.close()

    return R2_values

def R2_stat_multi_thread():
    list_bams_to_preccess = make_list_of_bam_files(args.bam, False)

    with Pool(processes=args.threads) as pool:
        results = list(tqdm(pool.imap_unordered(get_R2_values, list_bams_to_preccess), total=len(list_bams_to_preccess), desc="Bam files processed r2"))

    # Aggregate results from all processes
    c_R2_values = {}
    
    for result in results:
        try :
            c_R2_values.update(result)
        except:
            #print("result: {result}")
            pass
            
    return c_R2_values

def write_taxonomy_split_bam(bams):
    files_made=False
    for bam in bams :
        save = pysam.set_verbosity(0)
        samfile = pysam.AlignmentFile(bam, "rb")
        pysam.set_verbosity(save)
        assembly=re.sub('.*/|_[a-zA-Z].*', '', bam)
        if assembly in assembly_not_reported:
            return    
        if assembly not in assembly_mapped_dict:
            return
        
        total_95_readcount_dict[assembly]={}
        try:
            taxid, taxa_name, phylum_name, order_name, family_name, genus_name, db = accession_id2taxa[assembly]
        except:
            return
        taxa_name = re.sub(r"\'|\.", "",re.sub(' ', '_', taxa_name))

        if files_made == False:
            files_made = True
            bam_dict_genus, bam_dict_family, bam_dict_order = {}, {}, {}
            
            for library in library_list:
                if library not in assembly_mapped_dict[assembly]:
                    continue
                if 'genus' in assembly_mapped_dict[assembly][library]:
                    if len(assembly_mapped_dict[assembly][library]['genus'])>=100:
                        if not os.path.exists(f"{args.split_bam_out}/{taxa_name}_{assembly}"):
                            os.makedirs(f"{args.split_bam_out}/{taxa_name}_{assembly}")
                        bam_dict_genus[library]=pysam.AlignmentFile(f"{args.split_bam_out}/{taxa_name}_{assembly}/{assembly}_{library}_genus.bam", "wb", header=samfile.header)

                if 'family' in assembly_mapped_dict[assembly][library]:
                    if len(assembly_mapped_dict[assembly][library]['family'])>=100:
                        if not os.path.exists(f"{args.split_bam_out}/{taxa_name}_{assembly}"):
                            os.makedirs(f"{args.split_bam_out}/{taxa_name}_{assembly}")
                        bam_dict_family[library]=pysam.AlignmentFile(f"{args.split_bam_out}/{taxa_name}_{assembly}/{assembly}_{library}_family.bam", "wb", header=samfile.header)

                if 'order' in assembly_mapped_dict[assembly][library]:
                    if len(assembly_mapped_dict[assembly][library]['order'])>=100:
                        if not os.path.exists(f"{args.split_bam_out}/{taxa_name}_{assembly}"):
                            os.makedirs(f"{args.split_bam_out}/{taxa_name}_{assembly}")
                        bam_dict_order[library]=pysam.AlignmentFile(f"{args.split_bam_out}/{taxa_name}_{assembly}/{assembly}_{library}_order.bam", "wb", header=samfile.header)

        #genus_bam = pysam.AlignmentFile(f"split_bams/"+assembly+"_"+taxa_name+'_genus.bam', "wb", template=samfile)
        #family_bam = pysam.AlignmentFile("split_bams/"+assembly+"_"+taxa_name+'_family.bam', "wb", template=samfile)
        
        for read in samfile:
            if read.query_length<args.min_read_length:
                continue
            #if read.mapping_quality <20:
            #    continue
            library=re.sub('.*:', '', read.query_name)

            if assembly not in assembly_mapped_dict:
                return
            if library not in assembly_mapped_dict[assembly]:
                continue
            
            if read.query_name in assembly_mapped_dict[assembly][library]['genus'] and len(assembly_mapped_dict[assembly][library]['genus'])>=100:
                bam_dict_genus[library].write(read)

            elif read.query_name in assembly_mapped_dict[assembly][library]['family'] and len(assembly_mapped_dict[assembly][library]['family'])>=100:
                bam_dict_family[library].write(read)
        
            elif read.query_name in assembly_mapped_dict[assembly][library]['order'] and len(assembly_mapped_dict[assembly][library]['order'])>=100:
                bam_dict_order[library].write(read)

    samfile.close()
    for file in bam_dict_genus.values():
        file.close()
    for file in bam_dict_family.values():
        file.close()

def write_split_bams_multi_thread():
    if not os.path.exists('split_bams'):
        os.makedirs('split_bams')

    list_bams_to_preccess = make_list_of_bam_files(args.bam, False)

    with Pool(processes=args.threads) as pool:
        results = list(tqdm(pool.imap(write_taxonomy_split_bam, list_bams_to_preccess), total=len(list_bams_to_preccess), desc="Bam files to write"))

def find_elbow(x, y, threshold=0):
    """
    Finds the index of the elbow point in a curve defined by x and y.
    
    Parameters:
    x (array-like): x-coordinates of the points.
    y (array-like): y-coordinates of the points.
    threshold (float): Threshold for deciding if the elbow is significant.
    
    Returns:
    int or None: Index of the elbow point, or None if no significant elbow.
    """
    # Convert inputs to numpy arrays
    x = np.array(x).flatten()
    y = np.array(y)
    
    # Start and end points
    #print(x[0], y[0])
    start = np.array([x[0], y[0]])
    end = np.array([x[-1], y[-1]])
    
    # Line vector and its normalized version
    line_vec = end - start
    line_vec_norm = line_vec / np.linalg.norm(line_vec)
    
    # Distances to the line
    vec_from_start = np.column_stack((x - start[0], y - start[1]))
    proj_lengths = np.dot(vec_from_start, line_vec_norm)
    proj_points = np.outer(proj_lengths, line_vec_norm)
    dists = np.sqrt(np.sum((vec_from_start - proj_points)**2, axis=1))
    
    # Find the index of the maximum distance (elbow point)
    elbow_idx = np.argmax(dists)
    
    # Decide if the elbow is significant
    if dists[elbow_idx] < threshold * np.max(y):
        return None  # No significant elbow
    else:
        return elbow_idx

def write_table():
    file = open(args.output, 'w')
    levels = ['genus','family', 'order']

    file.write('\t'.join(['library', 'assembly', 'genome_length_unflitered', 'genome_length_filtered', 'phylum_name', 'order_name', 'family_name', 'genus_name', 'taxa_name', 
                          'mapped_read_count', 'removed_high_coverage', 'removed_excluded_contig', 'read_count_mapped_95%', 'bact_sink_count', 'genus_count_assembly', 'genus_deam_count', 'sum_genus_count', 'family_count_assembly', 'family_deam_count', 'sum_family_count','order_count_assembly', 'order_deam_count', 'sum_order_count', 'phylum_count_assembly', 'sum_phylum_count',
                         'median_read_length_genus', 'median_read_length_family', 'median_read_length_order',
                         'ancientness_genus', 'ancientness_family', 'ancientness_order',
                         'r2_genus', 'r2_family', 'r2_order', 'r2_genus+family',
                         'deam5p_frac_genus', 'deam5p_genus', 'ref_C_sum_5p_genus',
                         'deam3p_frac_genus', 'deam3p_genus', 'ref_C_sum_3p_genus',
                         '5p_con_3p_genus', '3p_con_5p_genus', 'con_deam_genus', 'ref_C_sum_5p_con_3p_genus', 'ref_C_sum_3p_con_5p_genus',
                         'deam5p_frac_family', 'deam5p_family', 'ref_C_sum_5p_family',
                         'deam3p_frac_family', 'deam3p_family','ref_C_sum_3p_family',
                         '5p_con_3p_family', '3p_con_5p_family', 'con_deam_family', 'ref_C_sum_5p_con_3p_family', 'ref_C_sum_3p_con_5p_family',
                         'deam5p_frac_order', 'deam5p_order', 'ref_C_sum_5p_order',
                         'deam3p_frac_order', 'deam3p_order','ref_C_sum_3p_order',
                         '5p_con_3p_order', '3p_con_5p_order', 'con_deam_order', 'ref_C_sum_5p_con_3p_order', 'ref_C_sum_3p_con_5p_order'
                         ])+'\n')
    
    
    for assembly in assembly_list:
        if assembly not in accession_id2taxa:
            continue
        
        if assembly in assembly_not_reported:
                continue
        
        taxid, taxa_name, phylum_name, order_name, family_name, genus_name, db = accession_id2taxa[assembly]
        if assembly not in assembly_mapped_dict:
            continue

        for library in assembly_mapped_dict[assembly]:
            if library not in assembly_mapped_dict[assembly]: 
                continue   
            
            if library in total_readcount_dict[assembly]:
                total_mapped = str(total_readcount_dict[assembly][library])
            else:
                total_mapped='0'
            if library in total_95_readcount_dict[assembly]:
                all_count_95 = str(total_95_readcount_dict[assembly][library])  
            else:
                all_count_95='0'

            if library in read_count_high_coverage[assembly]:
                reads_removed_coverage = str(read_count_high_coverage[assembly][library])
            else: 
                reads_removed_coverage = '0'

            if library in read_count_removed_contigs[assembly]:
                reads_removed_contigs = str(read_count_removed_contigs[assembly][library])
            else: 
                reads_removed_contigs = '0'

            genus_count = str(len(assembly_mapped_dict[assembly][library]['genus']))
            family_count = str(len(assembly_mapped_dict[assembly][library]['family']))
            order_count = str(len(assembly_mapped_dict[assembly][library]['order']))
            phylum_count = str(len(assembly_mapped_dict[assembly][library]['phylum']))
            bact_sink_count = str(len(assembly_mapped_dict[assembly][library]['bact_sink']))

            try:            
                genus_deam_count = str(deamination_count[assembly][library]['genus']['deam_count'])
            except:
                genus_deam_count = '0'
            
            try:
                family_deam_count = str(deamination_count[assembly][library]['family']['deam_count'])
            except:
                family_deam_count = '0'
            
            try:
                order_deam_count = str(deamination_count[assembly][library]['order']['deam_count'])
            except:
                order_deam_count = '0'
            
            sum_genus_count='0'
            if genus_name in summary_taxa_count_dict[library]:
                sum_genus_count = str(summary_taxa_count_dict[library][genus_name])
            
            sum_family_count='0'
            if family_name in summary_taxa_count_dict[library]:
                sum_family_count = str(summary_taxa_count_dict[library][family_name])
            
            sum_order_count='0'
            if order_name in summary_taxa_count_dict[library]:
                sum_order_count = str(summary_taxa_count_dict[library][order_name])
            
            sum_phylum_count='0'
            if phylum_name in summary_taxa_count_dict[library]:
                sum_phylum_count = str(summary_taxa_count_dict[library][phylum_name])
            
            if assembly in genus_read_length:
                if library in genus_read_length[assembly]:
                    median_read_count_genus=genus_read_length[assembly][library]
                else:
                    median_read_count_genus=str(0)
            else:
                median_read_count_genus=str(0)

            if assembly in family_read_length:
                if library in family_read_length[assembly]:
                    median_read_count_family=family_read_length[assembly][library]
                else:
                    median_read_count_family=str(0)
            else:
                median_read_count_family=str(0)

            if assembly in order_read_length:
                if library in order_read_length[assembly]:
                    median_read_count_order=order_read_length[assembly][library]
                else:
                    median_read_count_order=str(0)
            else:
                median_read_count_order=str(0)

            deamination_list=[]
            ancientness_list=[]
            r2_list=[]
            for level in levels:
                #print(assembly, library, level)
                try:
                    r2 = R2_values[assembly][library][level]
                    r2_list.append(str(f"{r2:.2f}"))
                except:
                    r2_list.append(str(f"{-1}"))

                try:    
                    #if library in deamination_count[assembly]:
                    sum_5p = deamination_count[assembly][library][level]['sum_5p']
                    deam_5p = deamination_count[assembly][library][level]['deam_5p']
                except:
                    sum_5p = 0
                    deam_5p = 0

                if sum_5p != 0:
                    frac_5p_deam = f"{deam_5p/sum_5p:.3f}"
                else:
                    frac_5p_deam = 0


                try:
                    #if library in deamination_count[assembly]:
                    sum_3p = deamination_count[assembly][library][level]['sum_3p']
                    deam_3p = deamination_count[assembly][library][level]['deam_3p']
                except:
                    sum_3p=0
                    deam_3p=0

                if sum_3p != 0:
                    frac_3p_deam = f"{deam_3p/sum_3p:.3f}"   
                else:
                    frac_3p_deam = 0

                try:
                #if library in deamination_count[assembly]:
                    C_sum_5p_con_3p = deamination_count[assembly][library][level]['cond_3_sum_5p']
                    C_sum_3p_con_5p = deamination_count[assembly][library][level]['cond_5_sum_3p']
                    deam_read_count = deamination_count[assembly][library][level]['con_deam']
                    if C_sum_5p_con_3p != 0:
                        con_5p_3p = f"{deam_read_count/C_sum_5p_con_3p:.3f}"
                    else:
                        con_5p_3p = 0
                    if  C_sum_3p_con_5p != 0:
                        con_3p_5p = f"{deam_read_count/C_sum_3p_con_5p:.3f}"
                    else:
                        con_3p_5p = 0
                except:
                    con_3p_5p=0
                    C_sum_3p_con_5p=0
                    con_5p_3p=0
                    C_sum_5p_con_3p=0
                    deam_read_count=0
                deamination_list.extend( list(map(str, [frac_5p_deam, deam_5p, sum_5p, frac_3p_deam, deam_3p, sum_3p, con_5p_3p, con_3p_5p, deam_read_count, C_sum_5p_con_3p, C_sum_3p_con_5p])))

                ancientness = '-'           

                try:
                    if args.binomial_ci == True:
                        test51,test31 = float(binomial_ci(deam_5p, sum_5p)[0]),float(binomial_ci(deam_3p, sum_3p)[0])
                        
                    else:
                        try:
                            test51=deam_5p/sum_5p
                        except :
                            test51 = 0
                        try:
                            test31=deam_3p/sum_3p
                        except :
                            test31 = 0
                    
                    threshold_5p = args.damage_threshold[0]
                    threshold_3p = args.damage_threshold[1]
                    
                    if test51 and test31:
                        if test51 >= threshold_5p or test31 >= threshold_3p:
                            ancientness = '+'
                        if test51 >= threshold_5p and test31 >= threshold_3p:
                            ancientness = '++'
                    #print(f"level= {level}; genus_count= {genus_count}; test51 = {test51}; test31 = {test31}, ancientness {ancientness}")
                
                except ValueError:
                    pass

                #print(f'ancientness_list; {ancientness_list}')
                ancientness_list.append(ancientness)
            try:
                if R2_values[assembly][library]['genus+family'] != -1 :
                    r2_list.append(str(f"{R2_values[assembly][library]['genus+family']:.2f}"))
                else:
                    r2_list.append(f"-1")
            except:
                r2_list.append("NA")

            #print(f'ancientness_list; {ancientness_list}')
            #print([library, assembly, str(assembly2genome_lenght[assembly]), phylum_name, order_name, family_name, genus_name, taxa_name, total_mapped,reads_removed_coverage, reads_removed_contigs, all_count_95, genus_count, sum_genus_count, family_count, sum_family_count, sum_order_count, sum_phylum_count, median_read_count_genus, median_read_count_family]+ancientness_list+r2_list+deamination_list)
            file.write('\t'.join([library, assembly, str(assembly2genome_lenght_unflitered[assembly]), str(assembly2genome_lenght[assembly]), phylum_name, order_name, family_name, genus_name, taxa_name, total_mapped,reads_removed_coverage, reads_removed_contigs, all_count_95, bact_sink_count, genus_count, genus_deam_count, sum_genus_count, family_count, family_deam_count, sum_family_count, order_count, order_deam_count, sum_order_count, phylum_count, sum_phylum_count, median_read_count_genus, median_read_count_family, median_read_count_order]+ancientness_list+r2_list+deamination_list)+'\n')

if __name__ == '__main__':
    if args.version == True :
        print(version)
        quit()

    accession_id2taxa = read_database_summary()

    fcs_info_dict, fcs_interval_tree_dict = read_fcs_reports()

    read2taxon_dict, read_dict, total_95_readcount_dict, contigs_to_ignore, total_readcount_dict, read_count_removed_contigs, read_count_high_coverage, assembly2genome_lenght, assembly2genome_lenght_unflitered, assembly_not_reported, assembly_list, fcs_removed_reads = filter_bam_files()
    
    
    taxon_dict, assembly_mapped_dict, library_list = assign_taxonomy()
    
    del read_dict, read2taxon_dict # remove ram heavy dictonaries 
    
    deamination_count, genus_read_length, family_read_length, order_read_length = deamination_stat_multi_thread()
    summary_taxa_count_dict = count_total_number_taxon()

    R2_values = R2_stat_multi_thread()
    write_table()

    if args.write_taxonomy_split_bam == True:
        write_split_bams_multi_thread()
    
