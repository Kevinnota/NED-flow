.. _reference_db-page:

Reference Database Manager
==========================

NED-flow comes with a reference database manager to make it easier to organise reference genomes/assemblies and to help with updating and indexing. The database manager includes two tools:

- ``ned-ref-manager.py`` — used to download and update the reference database, and to create "sinks" for co-mapping to bacteria.
- Nextflow module ``ned.nf --build`` — used for indexing the database.

Downloading reference genomes
======================================

The reference genomes/assemblies are obtained from the GenBank FTP server:
``ftp.ncbi.nlm.nih.gov/genomes/genbank/``


The file ``assembly_summary.txt`` is downloaded and used to query the list of available genomes, their status, and taxonomic information. GenBank divides genomes by taxonomic groups, and NED-flow preserves this structure. The databases are split into:

**Databases intended for taxonomy classification in NED-flow:**

- ``plant``
- ``vertebrate_mammalian``
- ``vertebrate_other``
- ``invertebrate``

**Databases intended as sinks in NED-flow (for compatative-mapping):**

- ``fungi``
- ``archaea``
- ``bacteria``

.. important::

   ``viral`` and ``protozoa`` have not been tested with NED-flow for either sinks or taxonomic clasification.

.. important::

   While ``fungi``, ``archaea``, and ``bacteria`` can technically be used for taxonomic classification, this is **not recommended**. NED-flow is not tested or intended for classifying these taxonomic groups - there other dedicated tools for this.

.. Warning::
   Downloading all plant, vertebrate, and invertebrate genomes available at present requires ~23T of free disk space. Before starting to download and build the database, make sure there is enough space available. 
   for checking if NED-flow works first make a database ``--assembly_list`` flag. 

Refernece database manager
-----------------------

To download the reference database, navigate to the directory where the database should be stored, or specify ``--path_to_ref_db`` flag and run:

.. code-block:: bash

   ned-ref-manager.py -db [plant|vertebrate_mammalian|vertebrate_other|invertebrate]

This will:

- Download the ``assembly_summary.txt`` for the selected database.
- Check if a previous database exists.
- Estimate the number of available assemblies and approximate disk space required.
- Prompt the user to confirm with ``Yes/Y`` or cancel with ``No/N``.

**Database directory structure:**

.. code-block:: none

Each assembly has its own folder. The ``assembly_report.txt`` file contains taxonomic information, and the ``fcs_report.txt`` file includes Foreign Contamination Screening data, which NED-flow uses to mask regions that should be excluded from taxonomic assignment. And the ``genomic.fna.gz`` which contains their actual reference sequences.

   plant/
       └── GCA_xxxxxxxxx/
           ├── [GCA_xxxxxxxxx.x]_[asm_name]_assembly_report.txt
           ├── [GCA_xxxxxxxxx.x]_[asm_name]_fcs_report.txt
           ├── [GCA_xxxxxxxxx.x]_[asm_name]_genomic.fna.gz
           └── [GCA_xxxxxxxxx]_download.log


**Download specific assemblies:**

Is is also possible to download only a subset of assemblies using an accession list:

.. code-block:: bash

   ned-ref-manager.py -db plant --assembly_list example_files/ref_genome_list.tsv

This is useful if you want to supplement the default set of reference genomes (e.g., adding a wolf genome to complement domestic dog). Or for testing for NED-flow is working. 

**Specify a custom database location:**

Use the ``-p`` flag to define where the reference database should be stored:

.. code-block:: bash

   ned-ref-manager.py -db plant -p /path/to/database

Database Checking
-----------------

Errors can occur during download or later during the indexing. Use the ``--check-db`` flag to check the status of the database. This will validate the reference database without making updates. It will remove incomplete indexes and try to download missing files. The tool will print in the terminal recommendations for what actions should be taken.

``--check-db`` checks for:

- Presence of exactly one FASTA file per assembly directory.
- Missing FCS reports (and attempts to download them if missing).
- Changed FTP paths (automatically corrected).
- Removed or deprecated genomes.
- Indexing status.

.. code-block:: bash

   ned-ref-manager.py --check-db --path_to_ref_db [path]

Command-Line Usage
------------------

Here are all the options for ``ned-ref_manager.py``

.. code-block:: none

   usage: ned-ref-manager.py [-h] [--database DATABASE]
                             [--path_to_ref_db PATH_TO_REF_DB]
                             [--assembly_list ASSEMBLY_LIST]
                             [--check-db CHECK_DB]
                             [--version]

   Downloads and manages reference genomes for NED-flow

   options:
     -h, --help                     Show this help message and exit
     --database DATABASE, -db       GenBank database [archaea, bacteria, fungi,
                                    invertebrate, vertebrate_mammalian,
                                    vertebrate_other, plant, protozoa, viral]
     --path_to_ref_db PATH_TO_REF_DB, -p
                                    Path to the reference directory (default: current directory)
     --assembly_list ASSEMBLY_LIST, -al
                                    List of assemblies to download (TSV format)
     --check-db CHECK_DB            Check the integrity of the reference database
     --version                      Print version information

