.. NED-flow documentation master file, created by
   sphinx-quickstart on Thu Jul 31 14:57:25 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NED-flow's documentation
====================================

NED-flow is a NextFlow pipeline for analysing Nuclear (ancient) environmental DNA. 

With NED-flow it is easy to:

- download, update, and manage reference database with :code:`ned-ref-manager.py` and :code:`ned.nf --build` 
- map reads to each reference assembly independently with :code:`ned.nf --mapping`
- classify reads to the lowest common ancestor (LCA) with :code:`ned-classifier.py`, utilising data-driven quality filtering
- check database and mappings with :code:`ned-ref-manager.py --check-db` and :code:`ned-check.py`


Quickstart
------------------------------------
In this quickstart page, the basic workflow of NED-flow is explained with short descriptions of what and why it needs to be done. Every command is linked to its own page with more details and optional flags. Furthermore, there is a page with best practices and basic output visualisation in R. 


Install ned-flow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first step is downloading the NED-flow GitHub repository. All of the NED-flow suite is coded in :code:`Nexflow` and :code:`Python`. Only Nextflow needs to be installed, and a number of Python libraries. The two lines of code that need to be run are listed below. 

| Downloading/Cloning NED-flow GitHub repository:
| :code:`git clone https://github.com/kevinnota/NED-flow.git`
| :code:`cd NED-flow`

| Install Python dependencies:
| :code:`pip install -r ned-py-install.txt`


Reference database setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NED-flow uses reference genomes deposited on NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/). For NED-flow to operate, the database needs to be structured in a certain way. `ned-ref-manager.py` does this in an automated way. There are a lot of genomes available, and by default, NED-flow will download all of them. But to test if NED-flow works and get it up and running, it's recommended to start with a smaller subset of genomes. The full NED-flow database will require a lot of disk space (~22TB). It's recommended for people who work on a cluster to download and maintain it in a place that is accessible to all users. For this small download, it's totally alright to download it in the NED-flow directory. 

| If you are not in the NED-flow dir yet run:
| :code:`cd NED-flow`

| Download small subset of references:
| :code:`ned-ref-manager.py --assembly_list test_files/assmebly_list.txt`

| Ones the reference asseblies are downloaded its time to index them. This is done in Nextflow with the following command. 
| Index reference database:
| :code:`ned.nf --build --refgenomes --executor [local/slurm/sge] --path_reference_dbs '/path/to/ref_db/GCA*'`

.. admonition:: Recommendation
   
   :code:`ned.nd` is continuously adding jobs and might take a long time to run. We advise using :code:`tmux` for running :code:`ned.nd` so that it is possible to log out from the login node and have the processes continue in the background. :code:`tmux`. Check out the :code:`tmux` git page for more information https://github.com/tmux/tmux/wiki.


Preprocessing and mapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To make sure that NED-flow is doing what its suppose to do this tutorial is using two cave sediment papers from SLON et al., (2017). 

| start with making a directory for run and to store the fastq files:
| :code:`mkdir my_first_NED my_first_NED/slon_fastqs`
| :code:`cd my_first_NED`

| To download two samples run
| :code:`wget -O slon_fastqs/ERR1883475.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/005/ERR1883475/ERR1883475.fastq.gz`
| :code:`wget -O slon_fastqs/ERR1883480.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/000/ERR1883480/ERR1883480.fastq.gz`

To start :code:`NED-flow` some minimal preproccessing is done for each fastq file in the batch. This includes reads length filtering, removal of low quality reads, and reads with low complexity (dust). This is mainly done to speed up the mapping. It is not possible to add new sample to the batch, :code:`NED-flow` after the first mapping.   

| :code:`nextflow ~/NED-flow/ned.nf --preprocessing --input_fastq_tsv ~/NED-flow/example_files/slon_sample_list.tsv`
| :code:`ls -la filterd_fastq` 
.. note::
   total 485224
   drwxrwsr-x 2 kevin_nota genetics_g       100 Aug 15 14:10 ./
   drwxrwsr-x 9 kevin_nota genetics_g      4096 Aug 15 14:24 ../
   -rw-rw-r-- 1 kevin_nota genetics_g  90644805 Aug 15 14:00 ERR1883475_in_clean.dust.fastq.gz
   -rw-rw-r-- 1 kevin_nota genetics_g 406210349 Aug 15 14:10 ERR1883480_in_clean.dust.fastq.gz

Detailed documentation
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   instalation
   reference_db
   indexing_db
   preprocessing
   ned_mapper
   classify

