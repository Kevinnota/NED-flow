#!/usr/bin/env nextflow

include { DUSTMASKER }        from '../modules/local/dustmasker'
include { DUSTFILTER }        from '../modules/local/dustmasker'


workflow dustmasker {

    take: 
        input_fastq
    
    main:
        DUSTMASKER(input_fastq)
        DUSTFILTER(DUSTMASKER.out.fastq_in, DUSTMASKER.out.dust)

    emit:
        dusted_finished = DUSTFILTER.out.dust_filtered
        
}