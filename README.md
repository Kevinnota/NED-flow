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

*Databases intended to be used with NED-flow for taxonomy classification*
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
		-> [GCA_xxxxxxxxx.x]_[asm_name]_assembly_report.txt
		-> [GCA_xxxxxxxxx.x]_[asm_name]_fcs_report.txt
		-> [GCA_xxxxxxxxx.x]_[asm_name]_genomic.fna.gz
		-> [GCA_xxxxxxxxx]_download.log
	   
```
There is one directory for each assembly. Its will contain a `assembly_report.txt` that NED-flow uses for taxonomic information, and a `_fcs_report.txt` file that contains Foreign Contamination Screening information. This file will be used to mask reagions in the assembly that will not be used for taxonomic assignment.

If you want to download only a subsection of taxa, it is possible to give ned-ref-manager.py a list of assembly accession numbers. It is still required to specify which database NED-flow should query against [plant|vertebrate_mammalian|vertebrate_other|invertebrate].

```
ned-ref-manager.py -db plant -al [example_files/ref_genome_list.tsv]
```
This can also be done if specific assemblies should be added. By default NED-flow will only download the assembly that is indicated as reference genome. For example, for Canis lupis there is only domestic dog used as reference assembly, but with this funtion it is possible to add wolf to the reference database. 

It is also possible to specify the path to the where the reference database is/should be created with the `-p` flag.

## check the database 

During the downloading or indexing process, its possible that errors have occurred. Running `ned-ref-manager.py` with the `--check-db` flag will not update the database, but will check if there are mistakes in the database. It will check for
> - if there is one fasta file in the directory. 
> - if there are FCS reports missing. If they are missing, it will try to download them
> - if the FTP path of an assembly has changed, the normal updater will fix this.
> - if genomes are no longer available on NCBI
> - if the assembly is indexed

It will print some stats and end with a list of advice for further actions to take. 

```
ned-ref-manager.py --check-db
```

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
  --check-db CHECK_DB   tool to check the status of the database
  --version             print version

```

### Indexing reference asseblies.

Indexing reference assemblies is a task that will require quite some resources, especially the first time. This part is done in Nextflow and is part of the main `ned.nf` program. It will automatically locate assemblies that are not indexed and skip those. Therefore there is no need to specify which genome needs to be indexed but can add a wildcard. For indexing the reference assemblies, run:

```
ned.nf --build --refgenomes --executor [local/slurm/sge] --path_reference_dbs '/path/to/ref_db/GCA*'
```

Indexing on a cluster (slurm) is hard since you need to request a running time and memory. Small genomes require very little time and RAM to run, while large genomes take more resources. There is a risk in requesting too many resources that the queuing goes very slowly, and computing allocations are used up rapidly. On the other hand, requiring too few resources will result in a lot of errors with incomplete indexing attempts. With --executor sge, there is no need to specify RAM or runtime, but for SLURM, there is. So for slurm, memory and time are estimated by the input size of the FNA file. This is an estimation, so it might sometimes fail. There is also a default range between 8-128 GB of RAM and 3-24 hours in time. If an index requires more, then modify this manually or increase the max by requesting more resources by default. This can be done with the --memory/--memory-max and --time/--time-max flags. If failing persists, use `ned-ref-manager --check-db` to see which genomes failed to index.


### make-ref-sink.py 

Microbial, fungal or archeal DNA can often occur as a contaminant in reference genomes. There is an option ot make 'sinks' that download genomes and concatenate these to make sinks that reads can map against. This is a competitive mapping approach in the sense that a read maps to a sink; it will not be considered by `ned-classify.py` for taxonomic classification. Making the sinks is easy and can be done with `ned-build-ref-sink.py`. A sink can only be made for fungi, bacteria, archaea. 

```
ned-build-ref-sink.py -p [path/to/ref/db]
```

Indexing the sink can be done with `ned.nf --build` using the `--sink` flag

```
ned.nf --build --build --sink --executor [local/slurm/sge] --path_reference_dbs '/path/[bacteria/fingi/archaea]/*_sink'
```

## running NED-flow

The idea of NED-flow is that all reads are mapped indipendent to each reference genomes/assemblies in the database. This will create one bam file per mapping. To reduce the number of files, all reads are concatinated into one fastq file before mapping. This is all done automatically, but NED-flow does require the reads to have the library_id, or some kind of identifier at the end of the read header. This is used by downstream classifier for corretly assigning a read to a sample. This is also automatically done, but requires an input file with fastq_name and library identifier. 


### preproccessing and mapping


#### ned.nf --preproccessing

#### ned.nf --mapping

## Taxonomic classifying
### ned-classifier.py
### Output table info