process DUSTMASKER {
    maxForks params.maxForks

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
    maxForks params.maxForks

    input:
       path(input_fastq)
       path(dust_acclist)
       val(name)

    publishDir "${params.out_dir}/filterd_fastq", mode: 'move', overwrite: true, pattern: '*dust.fastq.gz'
    publishDir "${params.out_dir}/dust_stats", mode: 'move', overwrite: true, pattern: '*dustmasker_logs.tsv'

    output:
        path("*dust.fastq.gz")                           , emit: dust_filtered
        path("*dustmasker_logs.tsv")                     , emit: dustmasker_logs
    script:
    """
    file_out="\$(echo ${name} | sed 's/.fastq/dustmasker_logs.tsv/g' | sed 's/s_._//g')"

    dustmasker_filter.py -fq ${input_fastq} -d ${dust_acclist} --dust_threshold ${params.max_low_complexity_bases} > \${file_out}
    
    """  

}
