# NED-flow

NED-flow is a pipeline developed for using <ins>**N**</ins>uclear reference genomes to classify reads obtained from shotgun-sequenced ancient samples <ins>**E**</ins>nvironmental <ins>**D**</ins>NA using a Next<ins>**flow**</ins> backbone. The tool contains tools for database management & building and for preprocessing reads, mapping and taxonomic classifying. The pipeline was designed specifically for nuclear reference genomes/assemblies obtained from GenBank.

In short, NED-flow is downloading all reference genomes deposited on NCBI. Ned-ref-manager.py is written in a way that it will download all reference genomes the first time – and only updates the existing reference genomes subsequently. Only new assemblies and updated assemblies are downloaded. The make-ref-sink.py is a tool that will concatenate bacteria/fungi/archaea genomes into sinks that are used for mapping and removing potential microorganism contamination in reference genomes. The ned.nf --build is a Nextflow tool that will index only the genomes that need to be indexed. It can do this both locally and on an SGE/SLURM cluster.

> [!WARNING]
> Downloading all plant, vertebrate, and invertebrate genomes available at present requires <ins>**~23T**</ins> of free disk space. Before starting to download and build the database, make sure there is enough space available. 

>[!CAUTION]
> NED-flow does at the moment not support mtDNA or chloroplast DNA. Its for now intended to only use nuclear DNA.

## Installation

NED-flow requires Nextflow version >=22.10. Besides that it requires the instalation of a list of python libraries. All required libraries can be found in the ned-py-install.txt file and can be installed with pip.

One by one.
```
pip intsall <package_name>
```
or all at the same time.

```
pip install -r ned-py-install.txt
```

## Database management
NED-flow comes with a reference database manager to make it easier to organise reference genomes/assemblies and helps with updating and indexing. The database manager has two tools. A `ned-ref-manager.py`to download & update reference database, andcreate "sinks" for co-mapping to bacteria and a Nextflow `ned.nf --build`part for indexing.

### Getting started downloading sequences
The reference genomes/assemblies are optained from the genbank ftp server `ftp.ncbi.nlm.nih.gov/genomes/genbank/`. The `assembly_summary.txt` is downloaded and queried for the availible genomes, their status, and taxonomic information. Genbank genomes are divided into taxonomic groups, and NED-flow retains these. So the databases are split into:

*Databases intended to be used with NED-flow for taxonomy*
- plant
- vertebrate_mammalian
- vertebrate_other
- invertebrate

*databases intended as sinks for NED-flow*
>- fungi
>- archaea
>- bacteria

>[!IMPORTANT]
> viral and protazoa are not tested.
>- viral
>- protozoa

>[!IMPORTANT]
> Fungi/archaea/bacteria can be used for taxonomuc clasifying, but not recomended. NED-flow is not tested or ment for classifying these groups.

#### ned-ref-manager.py

To download the reference database for plants the easiest option is to navigate to the location where the database will be downloaded and mentaind and run the following command:
```
ned-ref-manager.py -db plant
```
The tool will start with downloading the `assembly_summary.txt` for the reference db requested. It will check if a the database has been downloaded before, or if its the first download. It will calculate the number of assembilies that are availble for downloading and it gives and estimation on the free disk space that is required. It will not check if there is enough space - it will start downloading if there is enough or not. Type `Yes or Y` if you want to download or `No or N` if you do not want to. 

The database structure is as follows:
```
plant/
	-> GCA_xxxxxxxxx/
		-> [GCA_xxxxxxxxx]_assembly_report.txt
		-> [GCA_xxxxxxxxx]_fcs_report.txt
		-> [GCA_xxxxxxxxx]_download.log
	   
```
There is one directory for each assembly. Its will contain a `assembly_report.txt` that NED-flow uses for taxonomic information, and a `_fcs_report.txt` file that contains Foreign Contamination Screening information. This file will be used to mask reagions in the assembly that will not be used for taxonomic assignment.

If you want to download only a subsection of taxa, it is possible to give ned-ref-manager.py a list of assembly accession numbers. It is still required to specify which database NED-flow should query against.

```
ned-ref-manager.py -db plant -al [example_files/ref_genome_list.tsv]
```
This can also be done if specific assemblies should be added. By default NED-flow will only download the assembly that is indicated at reference genome. For example, for Canis lupis there is only domestic dog used as reference assembly, but with this funtion it is possible to add wolf to the reference database. 

It is also possible to specify the path to the where the reference database with the `-p` flag.

```
usage: ned-ref-manager.py [-h] [--database DATABASE] [--path_to_ref_db PATH_TO_REF_DB] [--assembly_list ASSEMBLY_LIST] [--version]

Downloads and manages reference genomes for NED-flow

options:
  -h, --help            show this help message and exit
  --database DATABASE, -db DATABASE
                        the genbank databases [archaea, bacteria, fungi, invertebrate, vertebrate_mammalian, vertebrate_other, plant, protozoa, viral]
  --path_to_ref_db PATH_TO_REF_DB, -p PATH_TO_REF_DB
                        path to the reference directory, default is the curent working directory
  --assembly_list ASSEMBLY_LIST, -al ASSEMBLY_LIST
                        list of assemblies to download
  --version             print version

```

### Indexing reference asseblies.

Indexing reference asseblies is a taks that will require quite some resorces the first time. This part is done in Nextflow and is part of the main `ned.nf`. For indexing the reference assemblies run:

```
ned.nf --build --refgenomes --path_reference_dbs /path/to/ref_db/
```

```
ned-ref-manager.py -db [plant/vertebrate_mammalian/vertebrate_other/invertebrate]
```

### ned-ref-manager.py
```
~/NED-flow/ned-ref-manager.py -db plant
```

### make-ref-sink.py 

### ned.nf --build

## preproccessing and mapping


### ned.nf --preproccessing

### ned.nf --mapping

## Taxonomic classifying
### ned-classifier.py
### Output table info