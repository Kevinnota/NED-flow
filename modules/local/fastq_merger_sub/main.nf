process MERGER_SUB {
    maxForks params.maxForks

    //container "${ workflow.containerEngine == 'singularity'
    //    'https://depot.galaxyproject.org/singularity/blast:2.16.0--hc155240_3'}"

    input:
        path(input_fastq)
    
    output:
    path("*merged.fastq.gz")                                     , emit: merged_fastq
    path("*sample_list.txt")                                      , emit: sample_list
    // This is necessary when new index is added since all samples need to be mapped.
    // cat $(ls | grep .dust.fastq.gz) > all_merged.fastq.gz
    // 
    
    script:
    """

    cat \$(ls ${input_fastq} | grep .dust.fastq.gz | grep -v -f "${params.out_dir}/run_logs/sample_log.txt" | sed 's|^|${input_fastq}/|') > sub_merged.fastq.gz

    ls ${input_fastq} | grep .dust.fastq.gz | grep -v -f "${params.out_dir}/run_logs/sample_log.txt" > sample_list.txt

    echo done
    """  

}
