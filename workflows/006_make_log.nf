#!/usr/bin/env nextflow

include { MAKE_LOG }        from '../modules/local/mapping_log'

workflow make_log {

    take: 
        sample_list
        assembly_ids

    main:
        MAKE_LOG(sample_list, assembly_ids)
    
    //emit:
    //    mapped_bam  = BOWTIE2_MAPPER.out.mapped_bam
    //    assbly_id = BOWTIE2_MAPPER.out.assbly_id
        
}