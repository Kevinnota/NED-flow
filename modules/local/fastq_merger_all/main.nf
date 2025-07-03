process MERGER_ALL {
    maxForks params.maxForks

    //container "${ workflow.containerEngine == 'singularity'
    //    'https://depot.galaxyproject.org/singularity/blast:2.16.0--hc155240_3'}"

    input:
        path(input_fastq)
        //path(sample_run_tracker)
    
    output:
    path("*merged.fastq.gz")                                     , emit: merged_fastq
    path('list_of_samples.tsv')                                  , emit: sample_list
    //
    // This is necessary when new index is added since all samples need to be mapped
    // cat $(ls | grep .dust.fastq.gz) > all_merged.fastq.gz
    // 
    
    script:
    """

    cat \$(ls ${input_fastq} | grep .dust.fastq.gz | sed 's|^|${input_fastq}/|') > merged.fastq.gz
    ls ${input_fastq} | grep .dust.fastq.gz > list_of_samples.tsv
    """  

}
