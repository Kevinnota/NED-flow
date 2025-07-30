process BOWTIE2_INDEXER {

    input:
        path(input_dirs)
    

    maxForks params.maxForks_cluster
    
    cpus params.threads
    executor params.executor
    penv { params.executor == 'sge' ? 'smp' : null }

    memory {
        def fnaFile = file(input_dirs).listFiles()?.find { it.name.endsWith('.fna.gz') }
        if (!fnaFile) {
            throw new RuntimeException("No .fna.gz file found in ${input_dirs}")
        }

        def sizeGB = fnaFile.size() / 1e9
        def maxMemGB = params.memory_max?.replaceAll(/[^\d]/, '')?.toInteger() ?: 128
        def estimatedMem = Math.max(8, Math.min((double)(sizeGB * 12), (double)maxMemGB)).round(0)
        
        return params.executor == 'slurm'
            ? (params.memory ?: "${estimatedMem} GB")
            : null
    }

    time {
    def fnaFile = file(input_dirs).listFiles()?.find { it.name.endsWith('.fna.gz') }
    if (!fnaFile) {
        throw new RuntimeException("No .fna.gz file found in ${input_dirs}")
    }

    def sizeGB = fnaFile.size() / 1e9
    def maxTimeH = params.time_max?.replaceAll(/[^\d]/, '')?.toInteger() ?: 24
    def estimatedTime = Math.max(3, Math.min((double)(sizeGB * 0.5), (double)maxTimeH)).round(0)

    return params.executor == 'slurm'
        ? (params.time ?: "${estimatedTime}h")
        : null
    }

    errorStrategy { 'ignore' }

    //clusterOptions '-l class=*'
    

    container "${ workflow.containerEngine == 'singularity'
        'https://depot.galaxyproject.org/singularity/bowtie2:2.5.3--py38he00c5e5_0' }"


    publishDir "${input_dirs}", mode: 'move', overwrite: true


    when: 
        //file("${input_dirs}").listFiles().findAll { it.isFile() && (it.name.endsWith('.bt2') || it.name.endsWith('.bt2l')) }.isEmpty()
        def files = file("${input_dirs}").listFiles().findAll { it.isFile() }

        def hasNoBt2 = files.findAll { it.name.endsWith('.bt2') || it.name.endsWith('.bt2l') }.isEmpty()
        def hasFnaGz = files.any { it.name.endsWith('.fna.gz') }

        return hasNoBt2 && hasFnaGz
    
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

