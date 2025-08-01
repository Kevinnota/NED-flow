.. _ned_mapper-page:

NED mapping
==========

The NED-flow concept is based on mapping all reads independently to each reference genome/assembly in the database. This will create one BAM file per mapping. To reduce the number of files, all reads are concatenated into one fastq file before mapping. This is all done automatically, but NED-flow does require the reads to have the library_id, or some kind of identifier at the end of the read header. The downstream classifier uses this for correctly assigning reads to samples. This is also automatically done, but requires an input file with fastq_name and library_identifier. 

preproccessing
-----------------

NED-flow expects the reads to be pre-merged and adapter-trimmed. The preprocessing is using ``fastp`` so adapters will be looked for. But because library preps can have different sequencing adapters etc., it is safer to merge and trim reads tailored to the right library prep methods.

fastq input file:

.. code-block:: bash
	ned.nf --preproccessing --path_reference_dbs '/path/to/ref_db/GCA*' --input_fastq_tsv [path/to/fastq/file]

``input_fastq_tsv`` is a tab-separated input file which requires a header with *lib_name* and a *file_path*. The file_path is either a relative path from the directory where Nextflow is run or an absolute path. In short, the path to where the input fastq files are located. The lib_name column should contain the sample identifiers. The safest names contain no special characters. Underscores (_) and hyphens (-) are allowed, colons (:) and at (@) are not. The lib_name is independent from the file name - it's simply the sample identifier.

check ``input_fastq.tsv`` in the ``example_file`` directory


bam input (for MPI-eva internal usage):

.. code-block:: bash
	ned.nf --preproccessing --path_reference_dbs '/path/to/ref_db/GCA*' --input_fastq_tsv [path/to/tsv/file]

.. code-block:: bash
	ned.nf --preproccessing --path_reference_dbs '/path/to/ref_db/GCA*' --bams_tsv [path/to/tsv/file]

the tsv file for internal usage had to contain ``lib_id``, ``lane`` and ``run_id`` in the header. The lib_id can not be changed but has to correspont to the RG. So the indexed libary id from coreDB. The libary can be either Lib.X.XXXX or Lib_X_XXXX. 
 
.. note::
	the bam files require a RG flag - this is already in there for bams processed for th genetics department.

mapping
-----------------

.. code-block:: bash
	ned.nf --mapping --maxForks_cluster 100 --all --sinks --fastq_files  





