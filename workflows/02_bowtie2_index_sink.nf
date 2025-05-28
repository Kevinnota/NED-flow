#!/usr/bin/env nextflow

include { BOWTIE2_INDEXER_SINK }        from '../modules/local/bowtie2_indexer_bact_sink'

workflow bowtie2_indexer_sink {

    take: 
        input_dirs
    
    main:
        BOWTIE2_INDEXER_SINK(input_dirs)

    //emit:
    //    index_files  = BOWTIE2_INDEXER.out.index_files

}