#!/usr/bin/env nextflow

include { BOWTIE2_MAPPER }        from '../modules/local/bowtie2_mapper'

workflow bowtie2_mapperA {

    take: 
        input_dirs
        //blah

    main:
        //input = input_dirs.map{
        //    [it, it.baseName]
        //}
        BOWTIE2_MAPPER(input_dirs)
    
    emit:
        mapped_bam  = BOWTIE2_MAPPER.out.mapped_bam
        assembly_ids = BOWTIE2_MAPPER.out.assembly_ids
        
}