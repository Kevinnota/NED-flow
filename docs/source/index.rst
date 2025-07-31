.. NED-flow documentation master file, created by
   sphinx-quickstart on Thu Jul 31 14:57:25 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NED-flow's documentation
====================================

NED-flow is a pipeline developed for using **N**uclear reference genomes to classify reads obtained from shotgun-sequenced ancient samples **E**nvironmental **D**NA using a Next**flow** backbone. The tool contains tools for database management & building and for preprocessing reads, mapping and taxonomic classifying. The pipeline was designed specifically for nuclear reference genomes/assemblies obtained from GenBank.

In short, NED-flow is downloading all reference genomes deposited on NCBI. ``Ned-ref-manager.py`` is written in a way that it will download all reference genomes the first time – and only updates the existing reference genomes and downloads newly released assemblies subsequently. The ``make-ref-sink.py`` is a tool that will concatenate ``bacteria/fungi/archaea`` genomes into sinks that are used for mapping and removing potential microorganism contamination in reference genomes. The ``ned.nf`` ``--build`` is a Nextflow tool that will index only the genomes that need to be indexed. It can do this both locally and on an ``SGE/SLURM`` cluster.


.. toctree::
   :maxdepth: 2
   :caption: Contents:


Reference Database Manager
==================

   reference_db


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
