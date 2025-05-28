#!/usr/bin/env nextflow

include { MERGER_SUB }        from '../modules/local/fastq_merger_sub'

workflow mergerB {

    take: 
        input_dirs

    main:
        MERGER_SUB(input_dirs)
    
    emit:
        merged_fastq  = MERGER_SUB.out.merged_fastq
        sample_list  = MERGER_SUB.out.sample_list
        
}