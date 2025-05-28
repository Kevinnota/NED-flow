process SAMTOOLS_SAM2BAM {
container "${ workflow.containerEngine == 'singularity'
        'https://depot.galaxyproject.org/singularity/samtools:1.9--h91753b0_8' }"

executor = 'local'

input:
    path(mapped_sam)

output:
    path("*.bam")           , emit: bam
    path "versions.yml"     , emit: versions


publishDir "${params.out_dir}/bams", mode: 'move', overwrite: true

script:
"""
    echo ${mapped_sam}
    for file in ${mapped_sam}; do
        echo \${file}
        out_file=\$(echo \$file | sed 's/\\.sam/\\.bam/g')
        
        samtools view -bS \$file > \${out_file}  

        rm "\$(readlink \$file)"
        rm \${file}
           
    done

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
        samtools: \$(samtools --version)
        END_VERSIONS  

"""

}