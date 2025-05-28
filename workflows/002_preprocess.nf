#!/usr/bin/env nextflow

include { PREPROCESSOR }        from '../modules/local/preprocessor'

workflow preprocessing {

    take: 
        input_fastq
    
    main:
        input = input_fastq.map{
            [it, it.baseName]
        }
        
        PREPROCESSOR(input)

        

    emit:
        fastp_fastq = PREPROCESSOR.out.fastp_fastq
        
}