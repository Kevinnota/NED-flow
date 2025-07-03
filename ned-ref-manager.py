#!/usr/bin/env python

import argparse
from tqdm import tqdm
import os
import re
import datetime
from ete3 import NCBITaxa
import subprocess
import sys
import shutil
import glob
from colorama import Fore, Style
import fnmatch
import requests

parser = argparse.ArgumentParser(prog='NCBI Genome indexer', description='')
parser.add_argument('--database', "-db", help='the genbank databases [archaea, bacteria, fungi, invertebrate, vertebrate_mammalian, vertebrate_other, plant, protozoa, viral]')
parser.add_argument('--path_to_ref_db', "-p", help='path to the reference directory', default='.')
parser.add_argument('--assembly_list', "-al", help='list of assemblies to download')
parser.add_argument('--version',  help='print version',action='store_true', default=False)
args=parser.parse_args()

#version="NCBI Genome indexer: version - v1.0 05 Jul-2024"
#version="NCBI Genome indexer: version - v1.1 05 Aug-2024" 
            ## V1.1 
            # It will remove genomes that are no longer representative for the species.
            # It will also remove genomes that are excluded from refseq if there are any. 
version="NCBI Genome indexer: version - v1.2 24 Aug-2024" 
            # excludes hybrids

ftp_server = 'ftp.ncbi.nlm.nih.gov'


############
# wish list
#   -> downloading specific assemblies that are not reference genomes
#       flage - include GCA_.... (wolf, instead of only take dog)
############
## count how many genomes were successfully downloaded
N_downloaded_success=0

## count how many genomes were not downloaded
N_downloaded_failed=0

## how many genomes were updated
N_genomes_updated=0

def get_right_directory():
    
    if args.path_to_ref_db == '.':
        path = os.getcwd()
    else:
        path = args.path_to_ref_db.split('/')

    if path.split('/')[-1] == args.database:
        pass
    elif path.split('/')[-1] != args.database: 
        args.path_to_ref_db = f'{path}/{args.database}'

def read_assembly_list():
    if ',' in args.assembly_list:
        assembly_list = args.assembly_list.split(',')
    else:
        assembly_list=[]
        file = open(args.assembly_list)
        for row in file:
            assembly_list.append(row.strip())
    return assembly_list

def get_the_assembly_summary():
    print('Download assembly summary from NCBI - might take a few minutes')
    ftp_path = 'rsync://ftp.ncbi.nlm.nih.gov/genomes/genbank/'+args.database+'/assembly_summary.txt'
    rsync_command = ['rsync', '-a', '--quiet', ftp_path, f'{args.path_to_ref_db}/{args.database}_assembly_summary.txt']
    try:
        subprocess.run(rsync_command, check=True)
    except subprocess.CalledProcessError as e:
            print(f"Error downloading files: {e}")

def read_database_version():
    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%Y-%m-%d")
    if not os.path.isdir(f"{args.path_to_ref_db}/log_files"):
        os.makedirs(f"{args.path_to_ref_db}/log_files")
    
    if not os.path.isfile(f"{args.path_to_ref_db}/log_files/database_version"):
        file = open(f"{args.path_to_ref_db}/log_files/database_version", 'w')
        file.write('database_version=1; date='+formatted_datetime+'\n')
        file.close()
        database_version=1

    else:
        file = open(f"{args.path_to_ref_db}/log_files/database_version", 'r')
        database = file.readlines()[0]
        database_version = int(re.sub(';.*', '', re.sub('database_version=', '', database)))
        database_version+=1
        file.close()

        file = open(f"{args.path_to_ref_db}/log_files/database_version", 'w')
        file.write('database_version='+str(database_version)+'; date='+formatted_datetime+'\n')
        file.close()

    return database_version

def read_and_update_download_history():    
    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%Y%m%d")
    date_start_download=formatted_datetime

    if not os.path.isdir(f"{args.path_to_ref_db}/log_files"):
        os.makedirs(f"{args.path_to_ref_db}/log_files")

    if not os.path.isfile(f"{args.path_to_ref_db}/log_files/download_history.tsv"):
        history_file = open(f"{args.path_to_ref_db}/log_files/download_history.tsv", 'w')
        history_file.write('\t'.join(['ref_db_status', 'assembly_accession', 'taxid', 'download_date', 'assembly_level', 'seq_rel_date'])+'\n')
        history_file.close()

    if not os.path.isfile(f"{args.path_to_ref_db}/log_files/download.log"):
        log_file = open(f"{args.path_to_ref_db}/log_files/download.log", 'w')
        log_file.write('\t'.join(['database_version','run_start_date', 'N_genomes_downloaded_before', 'N_genomes_downloaded', 'N_genomes_updated', 'total_N_genomes' ,'N_genomes_failed', 'N_genomes_removed'])+'\n')
        log_file.close()

    return date_start_download

def get_common_names_gbif(scientific_name):
    url = "https://api.gbif.org/v1/species/match"
    params = {"name": scientific_name}
    response = requests.get(url, params=params).json()
    if 'usageKey' in response:
        species_key = response['usageKey']
        vernacular_url = f"https://api.gbif.org/v1/species/{species_key}/vernacularNames"
        vernaculars = requests.get(vernacular_url).json()

        for v in vernaculars['results']:
            if v.get('language') == 'eng' and 'vernacularName' in v:
                return v['vernacularName']
        return "None" 
        

def read_genbank_assembly_summary(): 
    ## counts how many genomes were already downloaded.
    N_already_downloaded_before=0   
    ## creates a list with genomes that need to be updated because a new version is available.
    required_to_update = []
    ## creates a list with genomes that are new and need to be downloaded.
    new_need_to_download = []
    #for file in args.assembly_summary:
    assembly_summary = open(f"{args.database}_assembly_summary.txt", 'r')

    total_number_bases_download = 0
    taxon_names_update=[]
    taxon_names_to_download=[]
    assemblies_to_remove=[]
    list_taxa_to_remove=[]

    representative_dict={}
    excluded_from_refseq_dict={}
    taxonomy_dict={}
    assembly2taxa_name={}
    not_ref_genome = set()
    for row in assembly_summary:
        row_list = row.strip().split('\t')
        if '##' in row:
            continue
        if '#assembly_accession':
            i=0   
            for item in row_list:
                if 'assembly_accession' in item:
                    assembly_accession_index=i

                elif 'taxid' in item:
                    taxid=i
                elif 'refseq_category' in item:
                    refseq_category=i
                elif 'assembly_level' in item:
                    assembly_level=i 
                elif 'seq_rel_date' in item:
                    seq_rel_date=i
                elif 'ftp_path' in item:
                    ftp_path=i
                elif 'genome_size' in item:
                    genome_size=i
                elif 'organism_name' in item:
                    organism_name=i
                elif 'excluded_from_refseq' in item:
                    excluded_from_refseq=i

                i+=1
        if args.assembly_list != None:
            if re.sub('\\.[0-9].*', '', row_list[assembly_accession_index]) not in assembly_list:
                continue

        #make dict with reference genomes that are representative and which are not
        representative_dict[re.sub('\\.[0-9].*', '', row_list[assembly_accession_index])] = row_list[refseq_category]
        excluded_from_refseq_dict[re.sub('\\.[0-9].*', '', row_list[assembly_accession_index])] = row_list[excluded_from_refseq]
        
        assembly_accession_no_version = re.sub("\\.[1-9]*.", '', row_list[assembly_accession_index])
        taxonomy_dict[assembly_accession_no_version] = [assembly_accession_no_version, row_list[taxid], row_list[assembly_accession_index], row_list[assembly_level], row_list[seq_rel_date], row_list[ftp_path]+'/']
        
        assembly2taxa_name[re.sub('\\.[0-9].*', '', row_list[assembly_accession_index])] = row_list[organism_name]
        #print(assembly2taxa_name)
        
        #if row_list[refseq_category] != 'representative genome':
        if args.assembly_list == None:
            if row_list[refseq_category] != 'reference genome':
                not_ref_genome.add(re.sub("\\.[1-9]*.", '', row_list[assembly_accession_index]))
                continue

        if ' x ' in row_list[organism_name]:
            #assemblies_to_remove.append(re.sub('\\.[0-9].*', '', row_list[assembly_accession_index]))
            #list_taxa_to_remove.append(row_list[organism_name])
            continue
        if row_list[excluded_from_refseq] not in ['na', 'derived from single cell']:
            continue

        if row_list[ftp_path] == "na":
            print('SKIP:', row_list[assembly_accession_index], 'does not have a ftp path')
            continue
        
        if os.path.isdir(f"{args.path_to_ref_db}/{assembly_accession_no_version}"):
            files_in_dir = os.listdir(f"{args.path_to_ref_db}/{assembly_accession_no_version}")
      

            if any(fnmatch.fnmatch(filename, row_list[assembly_accession_index] + "*fna.gz") for filename in files_in_dir):
            #if any(row_list[assembly_accession_index] in filename for filename in files_in_dir):
                N_already_downloaded_before+=1
            else:
                required_to_update.append([assembly_accession_no_version, row_list[taxid], row_list[assembly_accession_index], row_list[assembly_level], row_list[seq_rel_date], row_list[ftp_path]+'/'])
                taxon_names_update.append(f"{assembly_accession_no_version}\t{row_list[organism_name]}")
        else:
            new_need_to_download.append([assembly_accession_no_version, row_list[taxid], row_list[assembly_accession_index], row_list[assembly_level], row_list[seq_rel_date], row_list[ftp_path]+'/'])
            total_number_bases_download+=int(row_list[genome_size])
            taxon_names_to_download.append(f"{assembly_accession_no_version}\t{row_list[organism_name]}")

    ## need to check what this does 230-271
    try:
        listed_file = os.listdir(args.path_to_ref_db)
    except:
        listed_file = []
    
    for assembly_accession in listed_file:
        #print(assembly_accession)
        matching_files = glob.glob(f'{args.path_to_ref_db}/{assembly_accession}/*_assembly_report.txt')
        #print(matching_files)
        if "GCA" not in assembly_accession:
            continue
        if assembly_accession not in assembly2taxa_name:
            file = open(matching_files[0],  'r')
            #print(f'lalalal_______________')
            #print(file)
            for row in file:
                if "Organism name" in row:
                    taxon_name=re.sub(".*:  | \\(.*", "", row.strip())
                    break
        if assembly_accession not in representative_dict:
            print(assembly_accession, taxon_name, "<---")
            #assemblies_to_remove.append(assembly_accession)
            continue
        #if representative_dict[assembly_accession] != 'representative genome':
        if representative_dict[assembly_accession] != 'reference genome':
            #print(representative_dict[assembly_accession])
            #print(assembly2taxa_name[assembly_accession])
            try:
                taxon_names_to_download.remove(f"{assembly_accession}\t{assembly2taxa_name[assembly_accession]}")
            except:
                ""
            
            try:
                taxon_names_to_download.remove(re.match(r"(\w+\s\w+)", f"{assembly_accession}\t{assembly2taxa_name[assembly_accession]}").group())
            except:
                ""

            #if assembly_accession in assembly2taxa_name:
            #    taxon_names_update.append(f"{assembly_accession}\t{assembly2taxa_name[assembly_accession]}")
            #else:
            #    taxon_names_update.append(taxon_name)
            assemblies_to_remove.append(assembly_accession)

            #print(taxonomy_dict[assembly_accession])

        if excluded_from_refseq_dict[assembly_accession] not in ['na', 'derived from single cell']: 
            #print(excluded_from_refseq_dict[assembly_accession] )
            list_taxa_to_remove.append(f"{assembly_accession}\t{assembly2taxa_name[assembly_accession]}")
            assemblies_to_remove.append(assembly_accession) 

    if len(taxon_names_to_download) <= 250:
        print('')
        print(Fore.GREEN +"taxa that are new to download"+Style.RESET_ALL)
        taxon_list = [re.sub(".*\t", '', item) for item in taxon_names_update]
        n_new_downloads=0
        for item in taxon_names_to_download:
            if re.sub(".*\t", '', item) in taxon_list:
                continue
            n_new_downloads+=1
           
            print(f'{item} \t {get_common_names_gbif(item.split('\t')[1])}')
        print("")
    else:
        n_new_downloads = len(taxon_names_to_download)
        
    if len(taxon_names_update) <= 250:
        print(Fore.YELLOW+"taxa that need to be updated"+Style.RESET_ALL)
        for item in taxon_names_update:
            print(f'{item} \t {get_common_names_gbif(item.split('\t')[1])}')
        print("")

    if len(list_taxa_to_remove) <= 250:
        print(Fore.RED +"taxa that need to be removed"+Style.RESET_ALL)
        taxon_list = [re.sub(".*\t", '', item) for item in taxon_names_update]
        k=0
        n_remove=0
        for item in list_taxa_to_remove:
            if re.sub(".*\t", '', item) in taxon_list:
                k+=1
                continue
            print(f'{item} \t {get_common_names_gbif(item.split('\t')[1])}')
            n_remove+=1
        if k == len(list_taxa_to_remove):
            print(f"No assemblies to remove")
        print("")

    #print("")
    #print(assemblies_to_remove)
    #quit()
    response = input('\n'+str(n_new_downloads)+" new genomes, "+str(len(taxon_names_update))+" to be updated, and "+str(n_remove)+" assemblies should be removed\nRequires and estimate of "+ str(int((total_number_bases_download*0.34)/(1024*1024*1024))) +" gb of free disk space\n\nDo you want to continue (Yes/No): ")
    if response not in ['y', 'Y', 'yes', 'Yes'] :
        sys.exit(1) 

    return required_to_update, new_need_to_download, N_already_downloaded_before, assemblies_to_remove, taxonomy_dict

def download_new_genomes(N_downloaded_success, N_downloaded_failed):
    history_file = open(f'{args.path_to_ref_db}/log_files/download_history.tsv', 'a')
    
    ## create lines for for writing a log file
    for item in new_need_to_download:
        current_datetime = datetime.datetime.now()
        formatted_datetime = current_datetime.strftime("%Y-%m-%d")
        assembly_accession_no_version, taxid, assembly_accession, assembly_level, seq_rel_date, ftp_path = item
        history_file.write('\t'.join(['first_download', str(assembly_accession), str(taxid), str(formatted_datetime), str(assembly_level), str(seq_rel_date)])+'\n')
        
        os.makedirs(f"{args.path_to_ref_db}/{assembly_accession_no_version}")

        ftp_path = re.sub("https", "rsync", ftp_path)
        rsync_command = ['rsync', '-a', '--progress', '--exclude=*_from_genomic*', '--exclude=*_cds_*', '--exclude=*_rna_*', '--include=*_genomic.fna.gz', '--include=*assembly_report.txt', '--include=*_fcs_report.txt', '--exclude=*', ftp_path, f"{args.path_to_ref_db}/{assembly_accession_no_version}"]

        # Execute the rsync command
        try:
            subprocess.run(rsync_command, check=True)
            print("Files downloaded successfully.")
            N_downloaded_success+=1
            log=open(f'{args.path_to_ref_db}/{assembly_accession_no_version}/{assembly_accession_no_version}_download.log', 'a')
            log.write("\t".join(['first_download', assembly_accession, formatted_datetime])+'\n')
            log.close

        except subprocess.CalledProcessError as e:
            print(f"Error downloading files: {e}")
            N_downloaded_failed+=1
    return N_downloaded_success, N_downloaded_failed

def update_genomes(N_genomes_updated, N_downloaded_failed):
    history_file = open(f'{args.path_to_ref_db}/log_files/download_history.tsv', 'a')
    
    ## create lines for for writing a log file
    for item in required_to_update:
        current_datetime = datetime.datetime.now()
        formatted_datetime = current_datetime.strftime("%Y-%m-%d")
        assembly_accession_no_version, taxid, assembly_accession, assembly_level, seq_rel_date, ftp_path = item
        history_file.write('\t'.join(['updated', str(assembly_accession), str(taxid), str(formatted_datetime), str(assembly_level), str(seq_rel_date)])+'\n')
        
        for file in os.listdir(assembly_accession_no_version):
            if 'log' in file:
                continue
            else:
                os.remove(assembly_accession_no_version+'/'+file)
        
        ftp_path = re.sub("https", "rsync", ftp_path)
        rsync_command = ['rsync', '-a', '--progress', '--include=*fna.gz', '--include=*assembly_report.txt', '--exclude=*', ftp_path, assembly_accession_no_version]

        # Execute the rsync command
        try:
            subprocess.run(rsync_command, check=True)
            print("Files downloaded successfully.")
            N_genomes_updated+=1
            log=open(assembly_accession_no_version+'/'+assembly_accession_no_version+'_download.log', 'a')
            log.write("\t".join(['updated', assembly_accession, formatted_datetime])+'\n')
            log.close()
        except subprocess.CalledProcessError as e:
            print(f"Error downloading files: {e}")
            N_downloaded_failed+=1
    
    return N_genomes_updated, N_downloaded_failed
            
def write_log_file(date_start_download, database_version):

    log_file = open(f'{args.path_to_ref_db}/log_files/download.log', 'a')
    log_file.write('\t'.join(map(str,[database_version, date_start_download , N_already_downloaded_before, N_downloaded_success, N_genomes_updated, N_already_downloaded_before+N_downloaded_success, N_downloaded_failed, N_removed]))+'\n')
    log_file.close()

def write_db_version_summary(database_version):
    ncbi = NCBITaxa()

    if not os.path.isdir(f"{args.path_to_ref_db}/database_version_summary"):
        os.makedirs(f"{args.path_to_ref_db}/database_version_summary")
    
    file = open(f'{args.path_to_ref_db}/database_version_summary/summary_version_{database_version:04}.tsv', 'w')  
    file.write("\t".join(['assembly_accession', 'taxa_name', 'phylum_name', 'order_name', 'family_name', 'genus_name', 'taxid', 'Assembly_level', 'Genome_representation'])+'\n')
    listed_file = os.listdir(args.path_to_ref_db)
    
    for assembly_accession in listed_file:
        if "GCA" not in assembly_accession:
            continue

        files = os.listdir(f'{args.path_to_ref_db}/{assembly_accession}')
        try :
            assembly_report = [file for file in files if "assembly_report.txt" in file][0]
        except:
            print(assembly_accession, 'has no assembly_report.txt')
            continue
        
  
        assembly_report_file = open(f'{args.path_to_ref_db}/{assembly_accession}/{assembly_report}')
        for row in assembly_report_file:
            if "Taxid" in row:
                taxid = re.sub('# Taxid: ', '', row).strip()
            elif 'Assembly level' in row:
                Assembly_level = re.sub('# Assembly level: ', '', row).strip()
            elif 'Genome representation' in row:
                Genome_representation = re.sub('# Genome representation: ', '', row).strip()

        # assembly_accession from file name

        try:
            taxon_name = ncbi.get_taxid_translator([int(taxid)])[int(taxid)]
        except :
            ncbi.update_taxonomy_database()
            try :
                taxon_name = ncbi.get_taxid_translator([int(taxid)])[int(taxid)]
            except:
                taxid='1'
        lineage = ncbi.get_lineage(int(taxid))
    
        phylum_name='NA'
        order_name='NA'
        family_name='NA'
        genus_name='NA'

        for item in lineage :
            if ncbi.get_rank([item])[item] in ['genus', 'family', 'order', 'phylum'] :
                if ncbi.get_rank([item])[item] == 'phylum':
                    phylum_name = ncbi.get_taxid_translator([int(item)])[int(item)]
                if ncbi.get_rank([item])[item] == 'order':
                    order_name = ncbi.get_taxid_translator([int(item)])[int(item)]
                if ncbi.get_rank([item])[item] == 'family':
                    family_name = ncbi.get_taxid_translator([int(item)])[int(item)]
                if ncbi.get_rank([item])[item] == 'genus':
                    genus_name = ncbi.get_taxid_translator([int(item)])[int(item)]

        file.write("\t".join([assembly_accession, 
                              taxon_name,
                              phylum_name,
                              order_name, 
                              family_name, 
                              genus_name,
                              taxid, 
                              Assembly_level, 
                              Genome_representation])+'\n')
    file.close()

def remove_assemblies():
    history_file = open(f'{args.path_to_ref_db}/log_files/download_history.tsv', 'a')
    N_removed=0
    for assembly in assemblies_to_remove:
        current_datetime = datetime.datetime.now()
        formatted_datetime = current_datetime.strftime("%Y-%m-%d")
        assembly_accession_no_version, taxid, assembly_accession, assembly_level, seq_rel_date, ftp_path = taxonomy_dict[assembly]
        history_file.write('\t'.join(['removed', str(assembly_accession), str(taxid), str(formatted_datetime), str(assembly_level), str(seq_rel_date)])+'\n')
        if os.path.isdir(assembly):
            shutil.rmtree(assembly)
            N_removed+=1

    return N_removed

if __name__ == '__main__':
    if args.version == True :
        print(version)
        quit()

    print("start main\n")
    get_right_directory()
    if args.assembly_list != None:
        print(f'Running contume reference assembly mode')
        assembly_list = read_assembly_list()
        

    ## Download the ncbi genome list and creates a list with all genomes that need to be download and updated
    get_the_assembly_summary()
    required_to_update, new_need_to_download, N_already_downloaded_before, assemblies_to_remove, taxonomy_dict = read_genbank_assembly_summary()
    
    ## Creates all log files
    database_version = read_database_version()
    date_start_download = read_and_update_download_history()
    
    ## Download and update the genomes
    N_downloaded_success, N_downloaded_failed = download_new_genomes(N_downloaded_success, N_downloaded_failed)
    N_genomes_updated, N_downloaded_failed = update_genomes(N_genomes_updated, N_downloaded_failed)
    
    ## Remove assemblies
    N_removed = remove_assemblies()

    # Write log files
    write_db_version_summary(database_version)
    write_log_file(date_start_download, database_version)

