#!/usr/bin/env nextflow

include { ID2FASTQ }        from '../modules/local/id2fastq_headers'

workflow add_id2fastq {

    take: 
        input
    
    main:
        ID2FASTQ(input)
        //input = input_fastq.map{
        //    [it, it.baseName, lib2id]
        //}
        
        
        

    emit:
        id2fastq = ID2FASTQ.out.id2fastq
        
}