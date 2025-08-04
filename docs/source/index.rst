.. NED-flow documentation master file, created by
   sphinx-quickstart on Thu Jul 31 14:57:25 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NED-flow's documentation
====================================

NED-flow is a NextFlow pipeline for analysing Nuclear (ancient) environmental DNA. 

Whith NED-flow it is easy to:

- download, update, and manage reference database with ``ned-ref-manager.py`` and ``ned.nf --build`` 
- map reads to each reference assembly independently with ``ned.nf --mapping``
- classify reads to the lowest common ancestor (LCA) with ``ned-classifier.py``, utilising data-driven quality filtering
- check database and mappings with ``ned-ref-manager.py --check-db`` and ``ned-check.py``


Quickstart
------------------------------------
Install ned-flow

With these two commands:
Clone git repository

.. code-block:: bash
   git clone https://github.com/kevinnota/NED-flow.git
   cd NED-flow

Install Python dependencies:

.. code-block:: bash
   pip install -r ned-py-install.txt

Download reference database:

.. code-block:: bash
   ned-ref-manager.py -db [plant|vertebrate_mammalian|vertebrate_other|invertebrate]

Index reference database:

.. code-block:: bash
   ned.nf --build --refgenomes --executor [local/slurm/sge] --path_reference_dbs '/path/to/ref_db/GCA*'

Download + index sink (obtional):

.. code-block:: bash
   ned-ref-manager.py -db [bacteria|fungi|archaea]
   ned-build-ref-sink.py --path_reference_dbs [path/to/ref/db]
   ned.nf --build --sink --executor [local/slurm/sge] --path_reference_dbs '/path/[bacteria/fingi/archaea]/*_sink'

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
------------------------------------

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   instalation
   reference_db
   indexing_db
   preprocessing
   ned_mapper
   classify

