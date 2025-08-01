.. NED-flow documentation master file, created by
   sphinx-quickstart on Thu Jul 31 14:57:25 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

NED-flow's documentation
====================================

.. |D| replace:: **D**
.. |E| replace:: **E**
.. |N| replace:: **N**
.. |flow| replace:: **flow**


NED-flow is a pipeline developed for targeting |N|uclear reference genomes to taxonomically classify reads obtained from shotgun-sequenced ancient |E|nvironmental |D|NA using a Next |flow|backbone. The pipeline contains tools for database management & building and for preprocessing reads, mapping and taxonomic classifying. The pipeline was designed specifically for nuclear reference genomes/assemblies obtained from GenBank. The general mode for classifying is mapping reads to all reference genomes/assemblies indipendently rather than to large RefSeq/NT concatinated reference indexes.

In a nuttshell, NED-flow is downloading all reference genomes deposited on NCBI using the ``Ned-ref-manager.py`` which is written in a way that it will download all reference genomes the first time – and only updates the existing reference genomes and downloads newly released assemblies subsequently. The ``make-ref-sink.py`` is a tool that will concatenate ``bacteria/fungi/archaea`` genomes into sinks that are used for compatative mapping to remove potential microorganism contamination in reference genomes. The ``ned.nf`` ``--build`` is a Nextflow tool that will index only the genomes that need to be indexed. It can do this both locally and on an ``SGE/SLURM`` cluster. ``ned.nf`` ``--preprocessing`` will do standared reads lenght, quality and complexity filtering. The the reads are then mapped with ``ned.nf`` ``--mapping``, which is mapping reads with bowtie2 indipendently to all reference assemblies, or only to a specified subpart of the database such as ``plant`` or ``vertabrates``. The last part of the pipeline is ``ned-classifier.py`` which is a multi-treaded python code for assiging taxonomy. 
 

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   index
   instalation
   reference_db
   indexing_db
   ned_mapper

