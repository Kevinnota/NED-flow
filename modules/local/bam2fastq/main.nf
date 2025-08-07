process BAM2FASTQ {
    maxForks params.maxForks
    
    container "${ workflow.containerEngine == 'singularity'
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_1' }"

    errorStrategy { 'ignore' }

    input:
        tuple path(bam), val(name)

    output:
        path("*.fastq.gz")      , emit: fastq
        path "versions.yml"     , emit: versions
    
    when:
        !file("${params.out_dir}/filterd_fastq/${name}_clean.dust.fastq.gz").exists()

    script:

    // -f 4 -F 192 keeps only merged reads, and removed all reads that are not merged
    // #samtools view -f 4 -F 192 ${bam} | awk 'length(\$10) >= ${params.min_read_length}' | samtools fastq -t - | gzip - > \${output_name}
    """
    output_name=\$(echo ${bam} | sed 's/.bam/.fastq.gz/')

    
    samtools view -f 4 -F 192 ${bam} | samtools fastq -t - | gzip - > \${output_name}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       samtools: \$(samtools --version)
    END_VERSIONShttps://coredb.eva.mpg.de/sequencing
    """
}