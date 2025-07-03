#!/usr/bin/env nextflow

include { DUSTMASKER }        from '../modules/local/dustmasker'
include { DUSTFILTER }        from '../modules/local/dustmasker'


workflow dustmasker {

    take: 
        input_fastq
          
    main:
        DUSTMASKER(input_fastq)
        DUSTFILTER(DUSTMASKER.out.fastq_in, DUSTMASKER.out.dust, input_fastq.map{it.baseName})

    emit:
        dusted_finished = DUSTFILTER.out.dust_filtered
        dustmasker_logs = DUSTFILTER.out.dustmasker_logs
}