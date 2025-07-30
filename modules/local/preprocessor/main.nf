process PREPROCESSOR {
    maxForks params.maxForks

    container "${ workflow.containerEngine == 'singularity'
        'https://depot.galaxyproject.org/singularity/fastp:0.24.0--heae3180_1'}"

    input:
        tuple path(input_fastq), val(name)

    output:
        path("*_clean.fastq.gz")                           , emit: fastp_fastq
        path("*_report.json")                              , emit: report_json
        path("*log")                                       , emit: fastp_log


    publishDir "${params.out_dir}/fastp_report", mode: 'move', overwrite: true, pattern: '*_report.json'
    publishDir "${params.out_dir}/fastp_log", mode: 'move', overwrite: true, pattern: '*log'

    script:
    """

    output_file="\$(echo ${input_fastq} | sed 's/\\.\\(fastq\\|fq\\).gz/_clean.fastq.gz/')"
    json_out="\$(echo ${name} | sed 's/\\.\\(fastq\\|fq\\)/_report.json/g')"
    log_out="\$(echo ${name} | sed 's/\\.\\(fastq\\|fq\\)/_fastp.log/g')"

    fastp -i ${input_fastq} -o \${output_file} --dedup --low_complexity_filter --complexity_threshold 50 --overrepresentation_analysis --length_required ${params.min_read_length} -w 4 -h report.html -j \${json_out} 2> \${log_out} 
    
    """  

}
