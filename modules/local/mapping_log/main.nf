process MAKE_LOG {
    maxForks params.maxForks

    //container "${ workflow.containerEngine == 'singularity'
    //    'https://depot.galaxyproject.org/singularity/fastp:0.24.0--heae3180_1'}"

    input:
        path(input_list)
        val(assembly_ids)


    publishDir "${params.out_dir}/run_logs", mode: 'move', overwrite: true, pattern: 'sample_log.txt'
    publishDir "${params.out_dir}/run_logs", mode: 'move', overwrite: true, pattern: 'database_version.txt'

    output:
    path("sample_log.txt")                                     , emit: sample_log
    path('database_version.txt')                               , emit: database_version    

    script:
    """
    if [ -e ${params.out_dir}/run_logs/sample_log.txt ] ; then
        cat ${params.out_dir}/run_logs/sample_log.txt ${input_list} > sample_log.txt
    else 
        cat ${input_list} > sample_log.txt
    fi

    if [ "${params.invertebrate}" = "true" ] || [ "${params.all}" = "true" ]; then
        echo invertebrate \$(cat ${params.path_reference_dbs}/invertebrate/log_files/database_version) >> database_version.txt
    fi

    if [ "${params.plant}" = "true" ] || [ "${params.all}" = "true" ]; then
        echo plant \$(cat ${params.path_reference_dbs}/plant/log_files/database_version) >> database_version.txt
    fi
    
    if [ "${params.vertebrate_mammalian}" = "true" ] || [ "${params.all}" = "true" ]; then
        echo vertebrate_mammalian \$(cat ${params.path_reference_dbs}/vertebrate_mammalian/log_files/database_version) >> database_version.txt
    fi
    
    if [ "${params.vertebrate_other}" = "true" ] || [ "${params.all}" = "true" ]; then
        echo vertebrate_other \$(cat ${params.path_reference_dbs}/vertebrate_other/log_files/database_version) >> database_version.txt
    fi  
    
    sort sample_log.txt | uniq > temp.txt && mv temp.txt sample_log.txt
    
    """  

}
