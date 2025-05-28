process BOWTIE2_INDEXER {

    maxForks params.maxForks
    
    //executor 'local'
    cpus "${params.threads}"
    executor 'sge'
    penv 'smp'
    errorStrategy { 'ignore' }
    //memory '8 GB'
    //clusterOptions '-l class=*'
    

    container "${ workflow.containerEngine == 'singularity'
        'https://depot.galaxyproject.org/singularity/bowtie2:2.5.3--py38he00c5e5_0' }"

    input:
        path(input_dirs)

    publishDir "${input_dirs}", mode: 'move', overwrite: true


    when: 
        file("${input_dirs}").listFiles().findAll { it.isFile() && (it.name.endsWith('.bt2') || it.name.endsWith('.bt2l')) }.isEmpty()  
    
    script:

    """
    fasta_name="${input_dirs}/*.fna.gz"
    assembly_report="${input_dirs}/*assembly_report.txt"
    total_number_basis="\$(cat \${assembly_report} | grep -v "#" | cut -f 9 | awk '{sum += \$0} END {print sum}')"
    echo \$total_number_basis
    if [ \$total_number_basis -lt 4000000000 ] ; then

        echo \$fasta_name
        echo ${input_dirs}
        echo "make small index"
        bowtie2-build --threads ${params.threads} ${input_dirs}/*.fna.gz \$fasta_name
    
    else 
        echo "make large index"
        bowtie2-build-l --threads ${params.threads} ${input_dirs}/*.fna.gz \$fasta_name

    fi
    """  

}

//#GCA_004027775/GCA_004027775.1_BraVar_v1_BIUU_assembly_report.txt | grep -v "#" | cut -f 9 | awk '{sum += $0} END {print sum}'

