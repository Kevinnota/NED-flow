#!/usr/bin/env nextflow

include { BOWTIE2_INDEXER }        from '../modules/local/bowtie2_indexer'

workflow bowtie2_indexer {

    take: 
        input_dirs
    
    main:
        BOWTIE2_INDEXER(input_dirs)

    //emit:
    //    index_files  = BOWTIE2_INDEXER.out.index_files

}