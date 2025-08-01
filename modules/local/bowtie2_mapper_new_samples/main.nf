process BOWTIE2_MAPPER {
    executor 'local'
    //errorStrategy 'retry'
    maxForks params.maxForks_SGE
    executor 'sge'
    penv 'smp'
    maxRetries 2
    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    cpus { task.attempt == 1 ? params.threads : params.threads+8 }

    //memory '16 GB'
    //clusterOptions '-l class=*'

    //container "${ workflow.containerEngine == 'singularity'
    //    'https://depot.galaxyproject.org/singularity/bowtie2:2.5.3--py38he00c5e5_0'}"

    container "${ workflow.containerEngine == 'singularity'
        'https://depot.galaxyproject.org/singularity/mulled-v2-c742dccc9d8fabfcff2af0d8d6799dbc711366cf:2c4c4e771c5f7d6e311c74234a98ccf71669d6fb-0'}"

    input:
        tuple val(assembly_number), path(input_dirs), path(mapped_bam), path(input_fastqs)
    
    
    publishDir "${params.out_dir}/logs/", mode: 'move', overwrite: true, pattern: '*.log'
    publishDir "${params.out_dir}/bams/", mode: 'move', overwrite: true, pattern: '*.bam'

    output:
    path "*_mapped.bam", includeInputs: true                        , emit: mapped_bam
    path("*.log")                               , emit: mapped_log
    val(assembly_number)                        , emit: assembly_ids
    

    when:
        def hasBt2Files = ["plant", "invertebrate", "vertebrate_mammalian", "vertebrate_other"].any { category ->
            def files = file("${params.path_reference_dbs}/${category}/${input_dirs}").listFiles()
            files && files.any { it.name.endsWith('.bt2') || it.name.endsWith('.bt2l') }
        }
        return hasBt2Files && file("${params.out_dir}/bams/${assembly_number}_mapped.bam").exists()

        
    //--np set to 0 because we filer out ambiguous characters anyway
    //-N sets the number of miss-matches in the seed
    //--ignore-quals so the alignment score is not impacted by the qualty if the read.
    // -D 1000 -R 1 -L 22 -i S,0,0.50 
    //(time bowtie2 --rg-id \${RG_TAG} --mm -p "${params.threads}" --no-unal -q \${fastq} -x ${input_dirs}/*fna.gz -N 0 --very-sensitive-local --ignore-quals --np 0 --rdg 999,999 --rfg 999,999 > \${assbly_id}_\${RG_TAG_file}_mapped.sam)2> \${assbly_id}_\${RG_TAG_file}.log
    // (time bowtie2 -p ${task.cpus} --no-unal -q ${input_fastqs} -x ${input_dirs}/*fna.gz -N 0 --sensitive --ignore-quals --np 0 --rdg 999,999 --rfg 999,999 | samtools view -b - > \${assembly_id}_mapped.bam) 2> \${assembly_id}.log
    
    script:

    """
    ref_seq=\$(echo ${input_dirs}/*fna.gz)
    assembly_id=\$(echo \${ref_seq} | sed 's/\\/.*//g')
    (time bowtie2 -p ${task.cpus} --no-unal -q ${input_fastqs} -x ${input_dirs}/*fna.gz -N 0 --sensitive --ignore-quals --np 0 --rdg 999,999 --rfg 999,999 | samtools view -b - > \${assembly_id}_new_samples_mapped.bam) 2> \${assembly_id}.log
    samtools merge -@ ${task.cpus} \${assembly_id}_merged.bam ${mapped_bam} \${assembly_id}_new_samples_mapped.bam
    mv \${assembly_id}_merged.bam \${assembly_id}_mapped.bam
    rm \${assembly_id}_new_samples_mapped.bam
    
    """

}