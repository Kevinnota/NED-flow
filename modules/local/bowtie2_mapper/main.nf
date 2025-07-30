process BOWTIE2_MAPPER {

    input:
        tuple path(input_fastqs), path(input_dirs), val(assembly_number)

    maxForks params.maxForks_cluster
    
    cpus params.threads
    executor params.executor
    penv { params.executor == 'sge' ? 'smp' : null }

    memory {
    def fnaFile = ["plant", "invertebrate", "vertebrate_mammalian", "vertebrate_other", "bacteria", "archaea", "bacteria/bact_sink"]
        .collect { category -> file("${params.path_reference_dbs}/${category}/${input_dirs}") }
        .findResult { dir ->
            dir.exists() ? dir.listFiles()?.find { it.name.endsWith('.fna.gz') } : null
        }

    if (!fnaFile) {
        throw new RuntimeException("No .fna.gz file found in any category under ${input_dirs}")
    }

    def sizeGB = fnaFile.size() / 1e9
    def maxMemGB = params.memory_max?.replaceAll(/[^\d]/, '')?.toInteger() ?: 128
    def estimatedMem = Math.max(8, Math.min((double)(sizeGB * 5), (double)maxMemGB)).round(0)

    return params.executor == 'slurm'
        ? (params.memory ?: "${estimatedMem} GB")
        : null
    } 

    time {
    def fnaFile = ["plant", "invertebrate", "vertebrate_mammalian", "vertebrate_other", "bacteria", "archaea", "bacteria/bact_sink"]
        .collect { category -> file("${params.path_reference_dbs}/${category}/${input_dirs}") }
        .findResult { dir ->
            dir.exists() ? dir.listFiles()?.find { it.name.endsWith('.fna.gz') } : null
        }

    if (!fnaFile) {
        throw new RuntimeException("No .fna.gz file found in any category under ${input_dirs}")
    }

    def sizeGB = fnaFile.size() / 1e9
    def maxTimeH = params.time_max?.replaceAll(/[^\d]/, '')?.toInteger() ?: 24
    def estimatedTime = Math.max(3, Math.min((double)(sizeGB * 0.5), (double)maxTimeH)).round(0)

    return params.executor == 'slurm'
        ? (params.time ?: "${estimatedTime}h")
        : null
    }

    errorStrategy { 'ignore' }

    container "${ workflow.containerEngine == 'singularity'
        'https://depot.galaxyproject.org/singularity/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:2c4c4e771c5f7d6e311c74234a98ccf71669d6fb-0'}"
    
    publishDir "${params.out_dir}/logs/", mode: 'move', overwrite: true, pattern: '*.log'
    publishDir "${params.out_dir}/bams/", mode: 'move', overwrite: true, pattern: '*.bam'

    output:
        path("*.bam")                           , emit: mapped_bam
        path("*.log")                           , emit: mapped_log
        val(assembly_number)                    , emit: assembly_ids

    when:
        def hasBt2Files = ["plant", "invertebrate", "vertebrate_mammalian", "vertebrate_other", "bacteria", "archaea", "bacteria/bact_sink"].any { category ->
            def files = file("${params.path_reference_dbs}/${category}/${input_dirs}").listFiles()
            files && files.any { it.name.endsWith('.bt2') || it.name.endsWith('.bt2l') }
        }

        def bt2File = ["plant", "invertebrate", "vertebrate_mammalian", "vertebrate_other", "bacteria", "archaea", "bacteria/bact_sink"]
            .collect { category -> 
                file("${params.path_reference_dbs}/${category}/${input_dirs}").listFiles()
            }
            .find { it } // skip null if directory doesn't exist
            ?.find { it.name.endsWith('.rev.1.bt2') || it.name.endsWith('.rev.1.bt2l') }
        
        def bamFile = file("${params.out_dir}/bams/${assembly_number}_mapped.bam")
        
        return hasBt2Files && (!bamFile.exists() || bamFile.lastModified() < bt2File.lastModified())

        //return hasBt2Files && !file("${params.out_dir}/bams/${assembly_number}_mapped.bam").exists()

        
    //--np set to 0 because we filer out ambiguous characters anyway
    //-N sets the number of miss-matches in the seed
    //--ignore-quals so the alignment score is not impacted by the qualty if the read.
    // -D 1000 -R 1 -L 22 -i S,0,0.50 
    //(time bowtie2 --rg-id \${RG_TAG} --mm -p "${params.threads}" --no-unal -q \${fastq} -x ${input_dirs}/*fna.gz -N 0 --very-sensitive-local --ignore-quals --np 0 --rdg 999,999 --rfg 999,999 > \${assbly_id}_\${RG_TAG_file}_mapped.sam)2> \${assbly_id}_\${RG_TAG_file}.log

    script:

    """

    ref_seq=\$(echo ${input_dirs}/*fna.gz)
    assembly_id=\$(echo \${ref_seq} | sed 's/\\/.*//g')
    
    (time bowtie2 -p ${task.cpus} --no-unal -q ${input_fastqs} -x ${input_dirs}/*fna.gz -N 0 --sensitive --ignore-quals --np 0 --rdg 999,999 --rfg 999,999 --score-min L,-1.2,-1.2 | samtools view -b - > \${assembly_id}_mapped.bam) 2> \${assembly_id}.log
        
    """

}