process DUSTMASKER {
    maxForks params.maxForks_PRE

    container "${ workflow.containerEngine == 'singularity'
        'https://depot.galaxyproject.org/singularity/blast:2.16.0--hc155240_3'}"

    input:
        path(input_fastq)
    
    output:
    path("*dust.txt")                                     , emit: dust
    path(input_fastq)                                     , emit: fastq_in

    script:
    """
    output_file="\$(echo ${input_fastq} | sed 's/_clean.fastq.gz/.dust.txt/')"
    dustmasker -level 1 -in <(zcat ${input_fastq} | awk 'NR % 4 == 1 {print ">" substr(\$0, 2)} NR % 4 == 2 {print}') -out \${output_file} -outfmt acclist
    
    """  

}

process DUSTFILTER {
    maxForks params.maxForks_PRE

    input:
       path(input_fastq)
       path(dust_acclist)

    publishDir "${params.out_dir}/filterd_fastq", mode: 'move', overwrite: true, pattern: '*dust.fastq.gz'
    
    output:
        path("*dust.fastq.gz")                           , emit: dust_filtered
    
    script:
    """
    
    dustmasker_filter.py -fq ${input_fastq} -d ${dust_acclist} --dust_threshold ${params.max_low_complexity_bases}
    
    """  

}
