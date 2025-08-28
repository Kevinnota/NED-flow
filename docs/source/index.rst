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

NED-flow uses reference genomes deposited on NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/). For NED-flow to operate, the database needs to be structured in a certain way. `ned-ref-manager.py` does this in an automated way. There are a lot of genomes available, and by default, NED-flow will download all of them. But to test if NED-flow works and get it up and running, it's recommended to start with a smaller subset of genomes. The full database will require a lot of disk space (~25TB). It's recommended for people who work on a cluster to download and maintain it in a place that is accessible to all users. For this small download (19 taxa), it's totally alright to download it in the :code:`NED-flow` directory. 

| If you are not in the NED-flow dir yet run:
| :code:`cd NED-flow` 

| Make a directory for the database
| :code:`mkdir ned_ref_db ; cd ned_ref_db`

| Download small subset of references:
| :code:`ned-ref-manager.py --assembly_list ~/NED-flow/example_files/costume_refDB_list.txt`

Ones the reference asseblies are downloaded its time to index them. This is done in Nextflow with the following command. 
| Index reference database:
| :code:`ned.nf --build --refgenomes --executor [local/slurm/sge] --path_reference_dbs '/path/to/ref_db/GCA*'`

.. admonition:: Recommendation
   
   :code:`ned.nf` is continuously adding jobs on a cluster which might take a long time to run. We advise using :code:`tmux` for running :code:`ned.nd` so that it is possible to logout from the login node and have the processes continue in the background. :code:`tmux`. Check out the :code:`tmux` git page for more information https://github.com/tmux/tmux/wiki.

Preprocessing and mapping
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To make sure that NED-flow is doing what its suppose to do this tutorial is using two cave sediment papers from `Slon et al., (2017)`_. 

.. _Slon et al., (2017): https://www.science.org/doi/10.1126/science.aam9695. 

| start with making a directory for run and to store the fastq files:
| :code:`mkdir my_first_NED my_first_NED/slon_fastqs`
| :code:`cd my_first_NED`

| To download two samples run:
| :code:`wget -O slon_fastqs/ERR1883475.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/005/ERR1883475/ERR1883475.fastq.gz`
| :code:`wget -O slon_fastqs/ERR1883480.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/000/ERR1883480/ERR1883480.fastq.gz`

To start :code:`NED-flow` some minimal preproccessing is done for each fastq file in the batch. This includes reads length filtering, removal of low quality reads, and reads with low complexity (dust). This is mainly done to speed up the mapping. It is not possible to add new sample to the batch, :code:`NED-flow` after the first mapping.   

| :code:`NED-flow` NED-flow expects one fastq or unmapped bam per library. :code:`NED-flow` takes a :code:`.tsv` file with the lib_id and and file_path header. The path can be absolute or relative. An example input.tsv can be found in :code:`~/NED-flow/example_files`

| Run the NED-flow preprocessing:
| :code:`nextflow ~/NED-flow/ned.nf --preprocessing --input_fastq_tsv ~/NED-flow/example_files/slon_sample_list.tsv`

| The preprocessing results in filtered_fastq directory with one file for each sample in the input sample_list.tsv:
| :code:`ls -la filterd_fastq` 

.. code-block:: none

   total 485224
   drwxrwsr-x 2 kevin_nota genetics_g       100 Aug 15 14:10 ./
   drwxrwsr-x 9 kevin_nota genetics_g      4096 Aug 15 14:24 ../
   -rw-rw-r-- 1 kevin_nota genetics_g  90644805 Aug 15 14:00 ERR1883475_in_clean.dust.fastq.gz
   -rw-rw-r-- 1 kevin_nota genetics_g 406210349 Aug 15 14:10 ERR1883480_in_clean.dust.fastq.gz

..  warning :: Do not add new samples after mapping
   
   NED-flow is set up to be easily updated. To make it easy, it comes at the cost that no new samples can be added to the run after the first mapping has started. This is because NED-flow is concatenating all reads before mapping and then checking with the database if an assembly has been mapped to or not. NED-flow does not check which samples/reads have been mapped. For new samples, make a new directory. It is possible to classify reads over multiple run directories.

After the preprocessing, the samples are ready for mapping. It is possible to do the mapping locally, on a Slurm or an SGE cluster. The tutorial is using a Slurm cluster.

| To start the mapping:
| :code:`~/NED-flow/ned.nf --mapping --all --executor slurm --maxForks_cluster 20 --path_reference_dbs /NED-flow/ned_ref_db/`

Taxonomic assignment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The NED taxonomic classifier reads in all the :code:`.bam` files created during the mapping. It does this multi-treaded, the more treads, the faster the classifier runs. The classifier will do the following data-driven filtering in this order: (I) remove reads with a higher than expected coverage based on the genome length, number of reads and read length, (II) remove contigs that are too short and have a disproportional amount of reads, (III) remove contigs that have too many reads mapping. After this filtering, it will assign reads to genus, family, order, and phylum per assembly. 

| To classify mapped reads
| :code:`NED-flow/ned-classifier.py -b bams/ -t 1 -o slon_taxa_out.tsv -pr .`

..  warning :: Interpreting the output table

   The table that is created is a raw summary per assembly. It does not give a summary of reads for any given sample to a given taxonomic rank. With NED, it is possible to check the assembly on which the assignments are based and make a call before creating a summary using different strategies. Check the best practices page for tips on how to go about it.



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

