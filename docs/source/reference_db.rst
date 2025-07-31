Reference Database Manager
==========================

NED-flow comes with a reference database manager to make it easier to organise reference genomes/assemblies and to help with updating and indexing. The database manager includes two tools:

- ``ned-ref-manager.py`` — used to download and update the reference database, and to create "sinks" for co-mapping to bacteria.
- Nextflow module ``ned.nf --build`` — used for indexing the database.

Getting Started: Downloading Sequences
--------------------------------------

The reference genomes/assemblies are obtained from the GenBank FTP server:
ftp.ncbi.nlm.nih.gov/genomes/genbank/


The file ``assembly_summary.txt`` is downloaded and queried to list available genomes, their status, and taxonomic information. GenBank divides genomes by taxonomic groups, and NED-flow preserves this structure. The databases are split into:

**Databases intended for taxonomy classification in NED-flow:**

- ``plant``
- ``vertebrate_mammalian``
- ``vertebrate_other``
- ``invertebrate``

**Databases intended as sinks in NED-flow (for co-mapping):**

- ``fungi``
- ``archaea``
- ``bacteria``

.. important::

   ``viral`` and ``protozoa`` have not been tested with NED-flow.

.. important::

   While ``fungi``, ``archaea``, and ``bacteria`` can technically be used for taxonomic classification, this is **not recommended**. NED-flow is not tested or intended for classifying these groups.

``ned-ref-manager.py``
-----------------------

To download the reference database for plants, navigate to the directory where the database should be stored and run:

.. code-block:: bash

   ned-ref-manager.py -db plant

This will:

- Download the ``assembly_summary.txt`` for the selected database.
- Check if a previous database exists.
- Estimate the number of available assemblies and approximate disk space required.
- Prompt the user to confirm with ``Yes/Y`` or cancel with ``No/N``.

**Database directory structure:**

.. code-block:: none

   plant/
       └── GCA_xxxxxxxxx/
           ├── [GCA_xxxxxxxxx.x]_[asm_name]_assembly_report.txt
           ├── [GCA_xxxxxxxxx.x]_[asm_name]_fcs_report.txt
           ├── [GCA_xxxxxxxxx.x]_[asm_name]_genomic.fna.gz
           └── [GCA_xxxxxxxxx]_download.log

Each assembly has its own folder. The ``assembly_report.txt`` file contains taxonomic information, and the ``fcs_report.txt`` file includes Foreign Contamination Screening data, which NED-flow uses to mask regions that should be excluded from taxonomic assignment.

**Download specific assemblies:**

You can download only a subset of assemblies using an accession list:

.. code-block:: bash

   ned-ref-manager.py -db plant -al example_files/ref_genome_list.tsv

This is useful if you want to supplement the default set of reference genomes (e.g., adding wolf genomes to complement domestic dog).

**Specify a custom database location:**

Use the ``-p`` flag to define where the reference database should be stored:

.. code-block:: bash

   ned-ref-manager.py -db plant -p /path/to/database

Database Checking
-----------------

If errors occurred during download or indexing, use the ``--check-db`` option. This will validate the reference database without making updates.

It checks for:

- Presence of exactly one FASTA file per assembly directory.
- Missing FCS reports (and attempts to download them if missing).
- Changed FTP paths (automatically corrected).
- Removed or deprecated genomes.
- Indexing status.

.. code-block:: bash

   ned-ref-manager.py --check-db

It will print statistics and suggest next steps to fix problems.

Command-Line Usage
------------------

.. code-block:: none

   usage: ned-ref-manager.py [-h] [--database DATABASE]
                             [--path_to_ref_db PATH_TO_REF_DB]
                             [--assembly_list ASSEMBLY_LIST]
                             [--check-db CHECK_DB]
                             [--version]

   Downloads and manages reference genomes for NED-flow

   options:
     -h, --help                      Show this help message and exit
     --database DATABASE, -db       GenBank database [archaea, bacteria, fungi,
                                    invertebrate, vertebrate_mammalian,
                                    vertebrate_other, plant, protozoa, viral]
     --path_to_ref_db PATH_TO_REF_DB, -p
                                    Path to the reference directory (default: current directory)
     --assembly_list ASSEMBLY_LIST, -al
                                    List of assemblies to download (TSV format)
     --check-db CHECK_DB            Check integrity of the reference database
     --version                      Print version information

