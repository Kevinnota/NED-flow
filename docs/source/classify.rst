.. _ned_mapper-page:

Taxonomic clasification
==========

.. code-block:: bash
	ned-classifier.py -b [path/to/mapped/bams] --path_to_references [path/to/referenge/database]


Command-Line Usage
------------------

.. code-block:: none

    usage: Taxon classifier [options]

    options:
      -h, --help            show this help message and exit
      --bam BAM [BAM ...], -b BAM [BAM ...]
                            bam(s)
      --library_type {SS,ss,DS,ds}, -lt {SS,ss,DS,ds}
                            library type, single or double stranded, default SS
      --path_to_references PATH_TO_REFERENCES, -pr PATH_TO_REFERENCES
                            path to the references
      --db_summary DB_SUMMARY [DB_SUMMARY ...], -db DB_SUMMARY [DB_SUMMARY ...]
                            database summary files - manual path selection
      --db_log DB_LOG, -dl DB_LOG
                            database summary log file written by nf - automatically selects the correct database summary file
      --json JSON, -j JSON  json file with RG
      --damage DAMAGE, -d DAMAGE
                            number of nucleotides to calculate damage on, default is first
      --damage_threshold DAMAGE_THRESHOLD, -dt DAMAGE_THRESHOLD
                            Damage threshold, (p5,p3) 0.05,0.05
      --binomial_ci         number of nucleotides to calculate damage on, default is first
      --write_taxonomy_split_bam, -ws
                            write bam split files for genus and family if there are more than 100 reads assigned
      --split_bam_out SPLIT_BAM_OUT, -wso SPLIT_BAM_OUT
                            split_bam_out
      --taxa TAXA, -tz TAXA
                            run the classifier only on a taxonomic level, such as bovids
      --print_r2_table, -r2
                            write r2 table for each assembly
      --r2_table_out R2_TABLE_OUT, -r2o R2_TABLE_OUT
                            r2 output dir
      --threads THREADS, -t THREADS
                            print version
      --exclude_assembly EXCLUDE_ASSEMBLY [EXCLUDE_ASSEMBLY ...], -e EXCLUDE_ASSEMBLY [EXCLUDE_ASSEMBLY ...]
                            print version
      --output OUTPUT, -o OUTPUT
                            print version
      --min_distance_sink MIN_DISTANCE_SINK, -mds MIN_DISTANCE_SINK
                            maximum distance from sink, default 0.90
      --min_distance MIN_DISTANCE, -md MIN_DISTANCE
                            maximum distance from reference, default 0.95
      --max_distance MAX_DISTANCE, -Md MAX_DISTANCE
                            maximum distance from reference, default 1
      --ignore_small_contigs IGNORE_SMALL_CONTIGS, -isc IGNORE_SMALL_CONTIGS
                            this flag will ignore contigs equal or shorter than set value, default 0
      --min_read_length MIN_READ_LENGTH, -ml MIN_READ_LENGTH
                            min read_length
      --bact_sink, -bs      Whether bact_sink needs was used to map, put False if bacteria were mapped independently (default=False).
      --version             print version
