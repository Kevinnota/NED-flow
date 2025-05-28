#!/usr/bin/env nextflow

include { MERGER_ALL }        from '../modules/local/fastq_merger_all'

workflow mergerA {

    take: 
        input_dirs
        //sample_tracker

    main:
        MERGER_ALL(input_dirs)
    
    emit:
        merged_fastq  = MERGER_ALL.out.merged_fastq
        sample_list   = MERGER_ALL.out.sample_list
}