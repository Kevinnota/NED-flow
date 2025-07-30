process ID2FASTQ {
    maxForks params.maxForks

    input:
        tuple val(lib_name), path(input_fastq)
    
    output:
        path("${lib_name}_in.fastq.gz")                                     , emit: id2fastq

    script:
    """
    out_name="\$(echo ${lib_name}_in.fastq.gz)"
    echo \${out_name}
    echo ${lib_name}
    zcat ${input_fastq} | awk -v lib=":${lib_name}" '{ sub(/[[:space:]]+\$/, "") } NR % 4 == 1 { \$0 = \$0 lib } { print }' | gzip - > \${out_name} 


    """  

}

