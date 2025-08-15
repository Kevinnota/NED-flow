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
import gzip

parser = argparse.ArgumentParser(prog='ned-ref-manager.py', description='Downloads and manages reference genomes for NED-flow')
parser.add_argument('--database', "-db", nargs='+', help='the genbank databases [all, vertebrates, archaea, bacteria, fungi, invertebrate, vertebrate_mammalian, vertebrate_other, plant, protozoa, viral]')
parser.add_argument('--path_to_ref_db', "-p", help='path to the reference directory, default is the curent working directory', default='.')
parser.add_argument('--assembly_list', "-al", help='list of assemblies to download')

parser.add_argument('--check-db', action='store_true', help='tool to check the status of the database')
parser.add_argument('--version',  help='print version',action='store_true', default=False)
args=parser.parse_args()

#version="ned-ref-manager.py: version - v1.0 05 Jul-2024"
#version="ned-ref-manager.py: version - v1.1 05 Aug-2024" 
            ## V1.1 
            # It will remove genomes that are no longer representative for the species.
            # It will also remove genomes that are excluded from refseq if there are any. 
#version="ned-ref-manager.py: version - v1.2 24 Aug-2024" 
            # excludes hybrids
version="ned-ref-manager.py: version - v1.3 23 Jul-2025"
            # fixed bug were new reference genomes of species that already excist in the database are now shown as updates.
            # for genomes that do not have a ftp link in the assembly_summary.txt, the link is created from availible information.
            # made a function to check if the database is correct - remove empty directories and gives stats

version="ned-ref-manager.py: version - v2.0 14 Aug-2025"
            # This version has an --all option, to update and download everything in one go
            # It is also resturctured and creates database structre

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
    
    # makes a database list so it can download the assembly stats.
    database_list = args.database
    if args.assembly_list != None:
        args.database = 'all'
    
    if 'all' in args.database:
        database_list = ['invertebrate', 'vertebrate_mammalian', 'vertebrate_other', 'plant']
    
    elif 'vertebrates' in args.database:
        database_list = ['vertebrate_mammalian', 'vertebrate_other']

    print('Download assembly summary from NCBI - might take a few minutes')
    try:
        os.makedirs(f"{args.path_to_ref_db}/tmp_download")
    except:
        pass

    for database in database_list:
        ftp_path = 'rsync://ftp.ncbi.nlm.nih.gov/genomes/genbank/'+database+'/assembly_summary.txt'
        
        rsync_command = ['rsync', '-a', '--no-motd', '--quiet', ftp_path, f'{args.path_to_ref_db}/tmp_download/{database}_assembly_summary.txt']
        
        for attempt in range(5):
            try:
                subprocess.run(rsync_command, check=True, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError as e:
                    print(f"Error downloading files: {e}")

    return database_list

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
                vernacularName = v['vernacularName']
                if ',' in v['vernacularName']:
                    vernacularName=re.sub(',.*', '', vernacularName)
                    return vernacularName
                elif vernacularName != "" :
                    return vernacularName
        return "" 
        

def read_genbank_assembly_summary(): 
    ## counts how many genomes were already downloaded.
    N_already_downloaded_before=0   
    ## creates a list with genomes that need to be updated because a new version is available.
    required_to_update = []
    ## creates a list with genomes that are new and need to be downloaded.
    new_need_to_download = []


    total_number_bases_download = 0
    taxon_names_update=[]
    taxon_names_to_download=[]
    assemblies_to_remove=[]
    list_taxa_to_remove=[]
    list_no_longer_ref_genome=[]

    representative_dict={}
    excluded_from_refseq_dict={}
    taxonomy_dict={}
    assembly2taxa_name={}
   
    not_ref_genome = set()

    for database in database_list: 
        assembly_summary = open(f"tmp_download/{database}_assembly_summary.txt", 'r')
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
                    elif 'asm_name' in item:
                        asm_name=i
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
            
            ## takes only assemblies that are given in a list
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

            ## this makes an ftp link if the link is not given in the assembly_summary file
            if row_list[ftp_path] == "na":
                #print('SKIP:', row_list[assembly_accession_index], 'does not have a ftp path')
                assembly_number_no_version=re.sub('GCA_|\\.[0-9].*', '', row_list[assembly_accession_index])
                
                triplets = [assembly_number_no_version[i:i+3] for i in range(0, len(assembly_number_no_version), 3)]
                assembly_name = re.sub(' ', '_', row_list[asm_name])
                ftp_path_sting=f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{triplets[0]}/{triplets[1]}/{triplets[2]}/{row_list[assembly_accession_index]}_{assembly_name}"
                
                
                row_list[ftp_path]=ftp_path_sting
                row_list[ftp_path]
                #print(row_list)
                #quit()
                #continue

            if os.path.isdir(f"{args.path_to_ref_db}/{assembly_accession_no_version}"):
                files_in_dir = os.listdir(f"{args.path_to_ref_db}/{assembly_accession_no_version}")
                assebly_name=re.sub(r'\+', '', row_list[asm_name])
                assebly_name=re.sub(' _', '_', assebly_name)
                assebly_name=re.sub('_ ', '_', assebly_name)
                assebly_name=re.sub(' ', '_', assebly_name)

                if any(filename == f"{row_list[assembly_accession_index]}_{assebly_name}_genomic.fna.gz" for filename in files_in_dir):
                #if any(fnmatch.fnmatch(filename, row_list[assembly_accession_index] + "*fna.gz") for filename in files_in_dir):
                #if any(row_list[assembly_accession_index] in filename for filename in files_in_dir):
                    N_already_downloaded_before+=1
                else:
                    #print('')
                    #print(row_list[asm_name])
                    #print(f"{row_list[assembly_accession_index]}_{assebly_name}_genomic.fna.gz")
                    #print('')
                    #print(files_in_dir)
                    #print('-----------_------------_------------')
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
            #print(f'{assembly_accession}, {matching_files}, "\n_________\n"')
            if "GCA" not in assembly_accession or len(matching_files)==0:
                continue
            if assembly_accession not in assembly2taxa_name:
                file = open(matching_files[0],  'r')
                for row in file:
                    if "Organism name" in row:
                        taxon_name=re.sub(".*:  | \\(.*", "", row.strip())
                        break

            ## this shows that genomes that are represed - they are no longer in the list of availible genomes 
            if assembly_accession not in representative_dict:
                #print(assembly_accession, taxon_name, "<---")
                list_taxa_to_remove.append(f"{assembly_accession}\t{taxon_name}")
                assemblies_to_remove.append(assembly_accession) 
                continue
            
            ## this gives a list of taxa that used to be reference genome but are no longer - for these genomes there will be an 'update'
            if representative_dict[assembly_accession] != 'reference genome':

                try:
                    taxon_names_to_download.remove(f"{assembly_accession}\t{assembly2taxa_name[assembly_accession]}")
                except:
                    ""
            
                try:
                    taxon_names_to_download.remove(re.match(r"(\w+\s\w+)", f"{assembly_accession}\t{assembly2taxa_name[assembly_accession]}").group())
                except:
                    ""

                #print(f"{assembly_accession}, {assembly2taxa_name[assembly_accession]}")
                list_no_longer_ref_genome.append(assembly2taxa_name[assembly_accession])
                assemblies_to_remove.append(assembly_accession)


            if excluded_from_refseq_dict[assembly_accession] not in ['na', 'derived from single cell']: 
                list_taxa_to_remove.append(f"{assembly_accession}\t{assembly2taxa_name[assembly_accession]}")
                assemblies_to_remove.append(assembly_accession) 
            assembly_summary.close()
    if len(taxon_names_to_download) <= 250:
        print('')
        print(Fore.GREEN +"taxa that are new to download"+Style.RESET_ALL)
        taxon_list = [re.sub(".*\t", '', item) for item in taxon_names_update]
        n_new_downloads=0
        for item in taxon_names_to_download:
            if re.sub(".*\t", '', item) in list_no_longer_ref_genome:
                continue
            
            if re.sub(".*\t", '', item) in taxon_list:
                continue
            n_new_downloads+=1
           
            common_name = get_common_names_gbif(item.split('\t')[1])
            #print(f"\n{common_name}")
            if common_name != "" and common_name != None :
                print(f'{item} ({common_name})')
            else :
                print(f'{item}')
        print("")

    else:
        n_new_downloads = len(taxon_names_to_download)

    
    if len(taxon_names_update) <= 250 and not len(taxon_names_update)==0:
        print(Fore.YELLOW+"taxa that need to be updated"+Style.RESET_ALL)
        print('\t'+Fore.YELLOW+"new version"+Style.RESET_ALL)
        for item in taxon_names_update:
            common_name = get_common_names_gbif(item.split('\t')[1])
            if common_name != "" and common_name != None :
                print(f'{item} ({common_name})')
            else :
                print(f'{item}')
        print("")

        print('\t'+Fore.YELLOW+"new genome"+Style.RESET_ALL)
        for item in taxon_names_to_download:
            if re.sub(".*\t", '', item) not in list_no_longer_ref_genome:
                continue
            common_name = get_common_names_gbif(item.split('\t')[1])
            if common_name != "" and common_name != None :
                print(f'{item} ({common_name})')
            else :
                print(f'{item}')
        print("")

    if len(list_taxa_to_remove) <= 250:
        if len(taxon_names_update) != 0:
            print(Fore.RED +"taxa that need to be removed"+Style.RESET_ALL)

        taxon_list = [re.sub(".*\t", '', item) for item in taxon_names_update]
        k=0
        n_remove=0
        for item in list_taxa_to_remove:
            if re.sub(".*\t", '', item) in taxon_list:
                k+=1
                continue
            if re.sub(".*\t", '', item) in list_no_longer_ref_genome:
                continue

            common_name = get_common_names_gbif(item.split('\t')[1])
            if common_name != "" and common_name != None :
                print(f'{item} ({common_name})')
            else :
                print(f'{item}')
            n_remove+=1
        #if k == len(list_taxa_to_remove):
        #    print(f"No assemblies to remove")
        print("")

    #print(assemblies_to_remove)
    #print(f"\nNumber of genomes to remove {len(assemblies_to_remove)}\n")

    response = input(
    f"\n{n_new_downloads} new genomes, {len(taxon_names_update)} to be updated, and {n_remove} assemblies should be removed\n"
    f"Requires an estimate of {Fore.RED}{int((total_number_bases_download * 0.34) / (1024 * 1024 * 1024))} gb{Style.RESET_ALL} of free disk space\n\n"
    "Do you want to continue (Yes/No): "
    )
    
    if response not in ['y', 'Y', 'yes', 'Yes'] :
        shutil.rmtree(f"{args.path_to_ref_db}/tmp_download")
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
        ftp_path = re.sub('/$', '', ftp_path)
        assembly_name = re.sub('.*/', '', ftp_path)

        #rsync_command = ['rsync', '-a', '--no-motd', '--progress', '--exclude=*_from_genomic*', '--exclude=*_cds_*', '--exclude=*_rna_*', '--include=*_genomic.fna.gz', '--include=*assembly_report.txt', '--include=*_fcs_report.txt', '--exclude=*', ftp_path, f"{args.path_to_ref_db}/{assembly_accession_no_version}"]
        rsync_command = ['rsync', '-a', '--no-motd', '--progress', f'{ftp_path}/{assembly_name}_genomic.fna.gz', f'{ftp_path}/{assembly_name}_fcs_report.txt', f'{ftp_path}/{assembly_name}_assembly_report.txt', f"{args.path_to_ref_db}/{assembly_accession_no_version}"]

        # Execute the rsync command
        for attempt in range(5):
            try:
                subprocess.run(rsync_command, check=True)
                print("Files downloaded successfully.")
                N_downloaded_success+=1
                log=open(f'{args.path_to_ref_db}/{assembly_accession_no_version}/{assembly_accession_no_version}_download.log', 'a')
                log.write("\t".join(['first_download', assembly_accession, formatted_datetime])+'\n')
                log.close
                break

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
        ftp_path = re.sub('/$', '', ftp_path)
        assembly_name = re.sub('.*/', '', ftp_path)
        #print(f'{ftp_path}/{assembly_name}_genomic.fna.gz')
        #print(assembly_name)
        #print(ftp_path)
        #quit()
        #rsync_command = ['rsync', '-a', '--no-motd', '--progress', '--include=*fna.gz', '--include=*assembly_report.txt', '--exclude=*', ftp_path, assembly_accession_no_version]
        rsync_command = ['rsync', '-a', '--no-motd', '--progress', f'{ftp_path}/{assembly_name}_genomic.fna.gz', f'{ftp_path}/{assembly_name}_fcs_report.txt', f'{ftp_path}/{assembly_name}_assembly_report.txt', f"{args.path_to_ref_db}/{assembly_accession_no_version}"]
        
        # Execute the rsync command
        for attempt in range(5):
            try:
                subprocess.run(rsync_command, check=True)
                print("Files downloaded successfully.")
                N_genomes_updated+=1
                log=open(assembly_accession_no_version+'/'+assembly_accession_no_version+'_download.log', 'a')
                log.write("\t".join(['updated', assembly_accession, formatted_datetime])+'\n')
                log.close()
                break
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
        try :
            assembly_accession_no_version, taxid, assembly_accession, assembly_level, seq_rel_date, ftp_path = taxonomy_dict[assembly]
        except:
            file = open(glob.glob(f'{args.path_to_ref_db}/{assembly}/*_assembly_report.txt')[0])
            for row in file:
                if 'GenBank assembly accession:' in row:
                    assembly_accession = re.sub(': .*', '',row.strip())
                if 'Taxid:' in row:
                    taxid = re.sub(': .*', '',row.strip())
                if 'Assembly level:' in row:
                    assembly_level = re.sub(': .*', '',row.strip())
                if  'Date:' in  row:
                    seq_rel_date = re.sub(': .*', '',row.strip())
            file.close()
        
        history_file.write('\t'.join(['removed', str(assembly_accession), str(taxid), str(formatted_datetime), str(assembly_level), str(seq_rel_date)])+'\n')
        

        if os.path.isdir(assembly):
            shutil.rmtree(assembly)
            N_removed+=1
    history_file.close()
    return N_removed

def check_database_consistency():
    required_suffixes = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    required_suffixes_large = ['.1.bt2l', '.2.bt2l', '.3.bt2l', '.4.bt2l', '.rev.1.bt2l', '.rev.2.bt2l']
        
    ii=0
    current_datetime = datetime.datetime.now()
    formatted_datetime = current_datetime.strftime("%Y-%m-%d")
    files_in_dir = os.listdir(f"{args.path_to_ref_db}")

    number_of_fcs_reports_downloaded = 0
    number_of_multi_version_error = 0
    number_of_old_files = 0
    number_of_asseblies_without_fna = 0
    genome_suppressed = 0
    not_indexed = 0
    index_incomplete = 0

    list_assemblies_no_fna = []
    list_assemblies_multi = []
    no_index_list = []
    index_incomplete_list = []

    progress = tqdm(total=len(files_in_dir), desc='number of assemblies checked')

    ## reads the latest downloaded assembly summary
    assembly_summary = open(f"{next(f for f in files_in_dir if f.endswith('_assembly_summary.txt'))}", 'r')
    assembly2ftp_link={}
    assembly2genome_size={}
    for row in assembly_summary:
        row_list=row.split('\t')
        if '##' in row:
            continue
            
        i=0    
        if '#' in row:
            for item in row.split('\t'):
                if item == 'ftp_path' :
                    ftp_path_index=i
                if item == '#assembly_accession' :
                    assembly_accession_index=i
                if item == 'asm_name' :
                    asm_name_index=i
                if item == 'genome_size':
                    genome_size_index = i
                i+=1
        assembly_number_no_version = re.sub('GCA_|\\.[0-9].*', '', row_list[assembly_accession_index])
        ftp_path_string = row_list[ftp_path_index]
        
        if row_list[ftp_path_index] == '':
            triplets = [assembly_number_no_version[i:i+3] for i in range(0, len(assembly_number_no_version), 3)]
            assembly_name = re.sub(' ', '_', row_list[asm_name])

            ftp_path_string=f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{triplets[0]}/{triplets[1]}/{triplets[2]}/{row_list[assembly_accession_index]}_{assembly_name}"

        assembly2ftp_link[f'GCA_{assembly_number_no_version}'] = ftp_path_string
        assembly2genome_size[f'GCA_{assembly_number_no_version}'] = row_list[genome_size_index]

    ## checks if the ftp path is still up to date.
    for item in files_in_dir:
        #print(item)
        if "GCA" not in item:
            progress.update()
            continue
        files = os.listdir(f"{args.path_to_ref_db}/{item}")
        count = sum(f.endswith('.fna.gz') for f in files)
        
        ## check if the fna file still excist
        
        try:
            fna_file = f"{next(f for f in files if f.endswith('fna.gz'))}"
        except:
            list_assemblies_no_fna.append(item)
            #shutil.rmtree(item)

            number_of_asseblies_without_fna += 1
            continue

        assembly_accession_name = re.sub('_genomic.fna.gz', '', fna_file)

        assembly_number_no_version = re.sub('\\..*|GCA_', '', fna_file)
        triplets = [assembly_number_no_version[i:i+3] for i in range(0, len(assembly_number_no_version), 3)]

        ftp_path_fna=f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{triplets[0]}/{triplets[1]}/{triplets[2]}/{assembly_accession_name}" #/{fna_file}"
        try:
            if ftp_path_fna != assembly2ftp_link[item]:

                #print(item)
                #print(ftp_path_fna)
                #print(assembly2ftp_link[item])
                #print('')
                number_of_old_files += 1    
        except:
            #print(item)
            genome_suppressed += 1
        
        set_assembly_numbers=set()
        
        if sum(f.endswith('.fna.gz') for f in files) != 1:

            for file in files:
                assembly=re.sub('_[aA-zZ].*', '', file)
                if '.' in assembly:
                    set_assembly_numbers.add(re.sub('_[aA-zZ].*', '', file))
            if len(set_assembly_numbers) == 2:
                # get latest version of the assembly
                latest_version = max(set_assembly_numbers, key=lambda x: int(x.split(".")[1]))
            
            for file in files:
                if latest_version not in file and not "_download.log" in file:
                    os.remove(f"{args.path_to_ref_db}/{item}/{file}")
            number_of_multi_version_error+=1
            list_assemblies_multi.append(item)
        #if sum(f.endswith('.fna.gz') for f in files) == 0:
        #    print(f'{item} has no fna file')
        #    print(item)

        progress.update()
        

        ## checks if fcs files were properly downloaded, and updates the assemblies for which this was not the case.
        if sum(f.endswith('fcs_report.txt') for f in files) == 0:
            file = open(f"{args.path_to_ref_db}/{item}/{next(f for f in files if f.endswith('assembly_report.txt'))}")
            for row in file :
                if 'Assembly name' in row:
                    asm_name = re.sub('.*:| ', '', row.strip())
                if 'GenBank assembly accession' in row :
                    assembly_accession = re.sub('.*:| ', '', row.strip())
            file.close()
            
            assembly_number_no_version=re.sub('GCA_', '', item)
            triplets = [assembly_number_no_version[i:i+3] for i in range(0, len(assembly_number_no_version), 3)]
            ftp_path_sting=f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/{triplets[0]}/{triplets[1]}/{triplets[2]}/{assembly_accession}_{asm_name}/{assembly_accession}_{asm_name}_fcs_report.txt"

            ftp_path = re.sub("https", "rsync", ftp_path_sting)
            rsync_command = ['rsync', '-a', '--no-motd', ftp_path, f"{args.path_to_ref_db}/{item}"]
            for attempt in range(2):
                try:
                    subprocess.run(rsync_command, check=True, stderr=subprocess.DEVNULL)
                    log=open(f"{item}/{item}_download.log", 'a')
                    log.write("\t".join(['updated fcs report', assembly_accession, formatted_datetime])+'\n')
                    log.close()

                    number_of_fcs_reports_downloaded+=1
                except :
                    pass
                    #print(f"Error downloading files")
    

        ## checks if the file sizes is resonable
        #print(assembly2genome_size[item])
        #print(item)
        
        #file = gzip.open(glob.glob(f'{item}/{item}*_genomic.fna.gz')[0], 'rt')
        #iupac_bases = set("ACGTRYKMSWBDHVN")
        #genome_size=0
        #for row in file :
        #    if '>' in row:
        #        continue
        #    if set(row.strip().upper()).issubset(iupac_bases) != True:
        #        print(row)
        #        continue
        #    genome_size+=len(row.strip())
        #print(f'counted geonme size : {genome_size}')
        #print(f'writen genome size : {assembly2genome_size[item]}')
        #if ii == 10:
        #    quit()
        #ii += 1
        #continue
        
        file_size=os.path.getsize(glob.glob(f'{item}/{item}*_genomic.fna.gz')[0])
        try:
            estimated_file_size=int(int(assembly2genome_size[item]) * 0.34)
        except:
            continue
        if file_size > estimated_file_size*1.4:
            print(estimated_file_size)
            print(estimated_file_size*1.4)
            print(file_size)
            print(f'{item} fna file is to large {estimated_file_size} bytes')

        #elif file_size < estimated_file_size*0.85:
        #    print(file_size-estimated_file_size)
        #    print(f'{item} fna file is too small {estimated_file_size} bytes')
        
    ## check if index has been made or not, or is incomplete
        #print(glob.glob(f'{item}/{item}*_genomic.fna.gz')[0])
        file_base = re.sub('.*/', '', glob.glob(f'{item}/{item}*_genomic.fna.gz')[0])
        #print(file_base)
        #file_base = f'{re.sub('.*/','', assembly2ftp_link[item])}_genomic.fna.gz'
        #print(file_base)
        #quit()
        missing = []
        for suf_bt2, suf_bt2l in zip(required_suffixes, required_suffixes_large):
             if not any(f in files for f in [file_base + suf_bt2, file_base + suf_bt2l]):
                missing.append(file_base + suf_bt2 + " or " + suf_bt2l)
        if len(missing)!=0 and not any(filename.endswith('.tmp') for filename in files):
            no_index_list.append(item)
            not_indexed+=1
        elif len(missing)!=0 and any(filename.endswith('.tmp') for filename in files):
            for file in files:
                if file.endswith('.tmp') or file.endswith('.sa') :
                    os.remove(f'{item}/{file}')

            index_incomplete_list.append(item)
            index_incomplete+=1
  

    progress.close()
    print('')
    print(f"Checked:\n")
    print(f"\tmultiple assemblies in one dir: {number_of_multi_version_error}")
    print(f"\tfcs reports downloaded: {number_of_fcs_reports_downloaded}")
    print(f"\tftp paths that were changed: {number_of_old_files}")
    print(f"\tgenome suppressed: {genome_suppressed}")   
    print(f"\tasseblies without fna: {number_of_asseblies_without_fna}") 
    print(f"\tasseblies not indexed: {not_indexed}") 
    print(f"\tasseblies with incomplete index: {index_incomplete}") 
    print('')


    print(f'Recommendations:\n')
    if number_of_multi_version_error+number_of_fcs_reports_downloaded+number_of_old_files+genome_suppressed+not_indexed+index_incomplete != 0:
        if number_of_multi_version_error != 0:
            print(f'\t> Some asseblies had multiple versions')
            print(f'\t  check if there is only one version left for these assemblies:')
            for assembly in list_assemblies_multi:
                print(f'\t{assembly}')
            print('')

        if number_of_fcs_reports_downloaded != 0:
            print(f'\t> Some fcs reports were not present, they are downloaded now so no acction required')
            print('')

        if number_of_old_files != 0:
            print(f'\t> NCBI changed to assembly file(s)')
            print(f'\t\trun ned-ref-manager.py again, it will automatically update the new file')
            print('')

        if genome_suppressed != 0:
            print(f'\t> NCBI removed these assemblies from the list - its unfortunate, but they cant be used')
            print(f'\t\tno action required')
            print('')

        if number_of_asseblies_without_fna != 0:
            print(f'\t> Some asseblies did not have a fna file')
            for assembly in list_assemblies_no_fna:
                print(f'\t\t{assembly}')
            print(f"\t  there is noting that can be done, wait a few days and check if the file is present than")
            print('')

        if not_indexed != 0 :
            print(f'\t> Some assemblies were not indexed\n')
            for assembly in no_index_list:
                print(f'\t\t{assembly}')
            print(f'\n\t\trun nextflow ned.nf --build')
            print('')

        if index_incomplete != 0 :
            print(f'\t> Some assemblies the indexing was not completed')
            for assembly in index_incomplete_list:
                print(f'\t\t{assembly}')
            print(f'\n\tthe incomplete files were removed')
            print(f'\t\trun nextflow ned.nf --build\n')
            print('')
    else:
        print(f"Everything looks good, no action is needed")
        
    print(' ')
if __name__ == '__main__':
    if args.version == True :
        print(version)
        quit()

    print("start main\n")
    if args.check_db == True:
        check_database_consistency()
        quit()

    #get_right_directory()
    if args.assembly_list != None:
        print(f'Running contume reference assembly mode\n')
        assembly_list = read_assembly_list()
        

    ## Download the ncbi genome list and creates a list with all genomes that need to be download and updated
    database_list = get_the_assembly_summary()
    required_to_update, new_need_to_download, N_already_downloaded_before, assemblies_to_remove, taxonomy_dict = read_genbank_assembly_summary()
    quit()

    ## Creates all log files
    date_start_download = read_and_update_download_history()
    
    ## Download and update the genomes
    N_downloaded_success, N_downloaded_failed = download_new_genomes(N_downloaded_success, N_downloaded_failed)
    N_genomes_updated, N_downloaded_failed = update_genomes(N_genomes_updated, N_downloaded_failed)
    
    ## Remove assemblies
    N_removed = remove_assemblies()

    # Write log files
    database_version = read_database_version()
    write_db_version_summary(database_version)
    write_log_file(date_start_download, database_version)

