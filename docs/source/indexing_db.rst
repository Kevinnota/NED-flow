.. _indexing_db-page:

Indexing reference genomes
=========================

Indexing reference assemblies is a task that will require quite some resources, especially the first time. This part is done in Nextflow and is part of the main `ned.nf` program. It will automatically locate assemblies that are not indexed and skip those. Therefore there is no need to specify which genome needs to be indexed but can add a wildcard. For indexing the reference assemblies, run:

.. code-block:: bash
	ned.nf --build --refgenomes --executor [local/slurm/sge] --path_reference_dbs '/path/to/ref_db/GCA*'


Indexing on a cluster (slurm) is hard since you need to request a running time and memory. Small genomes require very little time and RAM to run, while large genomes take more resources. There is a risk in requesting too many resources that the queuing goes very slowly, and computing allocations are used up rapidly. On the other hand, requiring too few resources will result in a lot of errors with incomplete indexing attempts. With --executor sge, there is no need to specify RAM or runtime, but for SLURM, there is. So for slurm, memory and time are estimated by the input size of the FNA file. This is an estimation, so it might sometimes fail. There is also a default range between 8-128 GB of RAM and 3-24 hours in time. If an index requires more, then modify this manually or increase the max by requesting more resources by default. This can be done with the --memory/--memory-max and --time/--time-max flags. If failing persists, use `ned-ref-manager --check-db` to see which genomes failed to index.

.. important::
	When indexing fails with ``--executor sge``. Try increasing the ``--threads`` flag. 

Make a sink 
---------------------
Microbial, fungal or archeal DNA can often occur as a contaminant in reference genomes. There is an option ot make 'sinks' that download genomes and concatenate these to make sinks that reads can map against. This is a competitive mapping approach in the sense that a read maps to a sink; it will not be considered by `ned-classify.py` for taxonomic classification. Making the sinks is easy and can be done with `ned-build-ref-sink.py`. A sink can only be made for fungi, bacteria, archaea. 

.. code-block:: bash
	ned-build-ref-sink.py --path_reference_dbs [path/to/ref/db]

Indexing the sink can be done with `ned.nf --build` using the `--sink` flag

.. code-block:: bash
	ned.nf --build --sink --executor [local/slurm/sge] --path_reference_dbs '/path/[bacteria/fingi/archaea]/*_sink'


Command-Line Usage
------------------
Here are all the options for ``ned.nf --build``

.. code-block:: none

   usage: nextflow ned.nf --build [options]

   Downloads and manages reference genomes for NED-flow

   build options:
     --executor                    	local/sge/slurm
                                    This option specifies if Nextflow will execute the computation localy or on cluster sge or slurm.
     --path_reference_dbs           Path to the reference assemblies, use as 'path/GCA_*'					
     --refgenomes                   Use when the genomes that need to be indexed are reference genomes
     --sink                         Use when the genomes that need to be indexed is a sink
	 
	slurm options (optional):
     --memory                       minimal memmory to allow for a job (default: 8 Gb )
     --memory_max                   maximal memmory to allow for a job (default: 128 Gb )
     --time                         minimal time to allow for a job (default: 3h )
     --time_max                     maximal time to allow for a job (default: 24h )

    standard options (optional):
     --threads                      The number of cores (default = 8)
     --maxForks                     Number of parallel jobs than can be run locally (default: 5)
                                    Keep in mind that every job will take --threads. so by default it will use 5 x 8 = 40 cores 
     --maxForks_cluster             Number of parallel jobs run on a cluster (default: 10).
                                    Keep in mind that every job will take --threads. so by default it will use 10 x 8 = 80 cores 
