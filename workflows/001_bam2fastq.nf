#!/usr/bin/env nextflow

include { BAM2FASTQ }        from '../modules/local/bam2fastq'

workflow bam2fastq {

    take: 
        bam_RG

    main:
        input = bam_RG.map{
            [it, it.baseName.getAt(0)]
        }
        BAM2FASTQ(input)

    emit:
        versions = BAM2FASTQ.out.versions
        fastq = BAM2FASTQ.out.fastq

}

