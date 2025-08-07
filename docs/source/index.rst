.. NED-flow documentation master file, created by
   sphinx-quickstart on Thu Jul 31 14:57:25 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NED-flow's documentation
====================================

NED-flow is a NextFlow pipeline for analysing Nuclear (ancient) environmental DNA. 

With NED-flow it is easy to:

- download, update, and manage reference database with `ned-ref-manager.py` and `ned.nf --build` 
- map reads to each reference assembly independently with `ned.nf --mapping`
- classify reads to the lowest common ancestor (LCA) with `ned-classifier.py`, utilising data-driven quality filtering
- check database and mappings with `ned-ref-manager.py --check-db` and `ned-check.py`


Quickstart
------------------------------------
In this quickstart page, the basic workflow of NED-flow is explained with short descriptions of what and why it needs to be done. Every command is linked to its own page with more details and optional flags. Furthermore, there is a page with best practices and basic output visualisation in R. 


Install ned-flow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first step is downloading the NED-flow GitHub repository. All of the NED-flow suite is coded in `Nexflow` and `Python`. So Nextflow needs to be installed and a number of Python modules. The two lines of code that need to be run are listed below. 

Downloading/Cloning NED-flow GitHub repository:

.. code-block:: bash

   git clone https://github.com/kevinnota/NED-flow.git
   cd NED-flow

Install Python dependencies:

.. code-block:: bash

   pip install -r ned-py-install.txt


Reference database setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

NED-flow uses reference genomes deposited on NCBI (https://www.ncbi.nlm.nih.gov/datasets/genome/). For NED-flow to operate, the database needs to be structured in a certain way. `ned-ref-manager.py` does this in an automated way. There are a lot of genomes available, and by default, NED-flow will download all of them. But to test if nedflow works and get it up and running, it's recommended to start with a smaller subset of genomes.

NED-flow database will requite a lot of disk space - so its recommended for people that work on a cluster to downloaded it in a place that is accessable to all users. For this small download its totally alright to download it in the NED-flow directory. 
`cd NED-flow`


Download small subset of references:

.. code-block:: bash 

   ned-ref-manager.py -db [plant|vertebrate_mammalian|vertebrate_other|invertebrate]

Index reference database:

.. code-block:: bash

   ned.nf --build --refgenomes --executor [local/slurm/sge] --path_reference_dbs '/path/to/ref_db/GCA*'

Download + index sink (optional):

.. code-block:: bash

   ned-ref-manager.py -db [bacteria|fungi|archaea]
   ned-build-ref-sink.py --path_reference_dbs [path/to/ref/db]
   ned.nf --build --sink --executor [local/slurm/sge] --path_reference_dbs '/path/[bacteria/fingi/archaea]/*_sink'

Mapping and taxonomic classification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Preprocessing reads:

.. code-block:: bash

   ned.nf --preproccessing --path_reference_dbs '/path/to/ref_db/GCA*' --input_fastq_tsv [path/to/fastq/file]

Mapping reads:

.. code-block:: bash

   ned.nf --mapping --maxForks_cluster 100 --all [--sinks] --fastq_files  

Classify reads:

.. code-block:: bash

   ned-classifier.py -b /path/to/bams -t [n_threads] --library_type [ss,SS,ds,DS] --path_to_references [path/to/reference/database] -o [output.tsv]


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

