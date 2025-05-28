#!/usr/bin/env nextflow

include { BOWTIE2_MAPPER }        from '../modules/local/bowtie2_mapper_new_samples/'

workflow bowtie2_mapperB {

    take: 
        input_dirs

    main:
        BOWTIE2_MAPPER(input_dirs)
    
    emit:
        mapped_bam  = BOWTIE2_MAPPER.out.mapped_bam
        assembly_ids = BOWTIE2_MAPPER.out.assembly_ids
        
}