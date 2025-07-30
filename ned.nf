#!/usr/bin/env nextflow

// load all workflows
include { bam2fastq                            } from './workflows/001_bam2fastq.nf'
include { preprocessing                        } from './workflows/002_preprocess.nf'
include { dustmasker                           } from './workflows/003_dustmasker.nf'
include { mergerA                              } from './workflows/004A_fastq_merger.nf'
include { mergerB                              } from './workflows/004B_fastq_merger.nf'
include { bowtie2_mapperA                      } from './workflows/005A_bowtie2_mapping.nf'
include { bowtie2_mapperB                      } from './workflows/005B_bowtie2_mapping.nf'
include { make_log                             } from './workflows/006_make_log.nf'
include { bowtie2_indexer                      } from './workflows/01_bowtie2_index.nf'
include { bowtie2_indexer_sink                 } from './workflows/02_bowtie2_index_sink.nf'
include { add_id2fastq                         } from './workflows/001B_fastq2id.nf'

// BIG THING TO SOLVE - updated genomes are not deleted and remapped.
    // solved, the bam file is replaced if the index data was after the formation of the bamfile in the bams/ directory
//The colors
red = "\033[0;31m"
white = "\033[0m"
yellow = "\033[0;33m"

// Define some functions

def exit_with_error_msg(error, text){
    println "[Shotgun analysis]: ${red}${error}: ${text}${white}"
    exit 0
}
def get_warn_msg(text){
    return "[Shotgun analysis]: ${yellow}(WARN): ${text}${white}"
}
def get_info_msg(text){
    return "[Shotgun analysis]: ${text}"
}
def exit_missing_required(flag){
    exit_with_error_msg("ArgumentError", "missing required argument ${flag}")
}

// --build this makes the bowtie2 index for reference genomes/ 

workflow{

if (params.build){
    input_dirs = Channel.fromPath(params.path_reference_dbs, type: 'dir')    

    // Bowtie2 index the reference genomes that need to be indexed, and than move them from the Work dir to the dir it supposed to be.
    if (params.refgenomes){
        bowtie2_indexer(input_dirs)
    }

    if (params.sink){
        bowtie2_indexer_sink(input_dirs)
    } 
    }


// --preprocessing, does the library triming, identical PCR duplicate removal, and quality filtering

if (params.preprocessing) {
    // reading in bams as a channel
    if (params.bams || params.bams_tsv) {
        
        if (params.bams) {
        bams = Channel.fromPath(params.bams)
        }
        
        if (params.bams_tsv) {
        bams = Channel
            .fromPath(params.bams_tsv)
            .splitCsv(header: true, sep: '\t')
            .map { row -> 
                def lib = row['lib'].replace('.', '_')
                def lane = row['lane']
                def run_id = row['run_id']
                file("/mnt/ngs_data/${run_id}/results/final/s_${lane}_${lib}*.bam")
            }
        }
    
    // from bam to fastq + removes all unmerged reads + minimal read filter

        bam2fastq(bams)
        fastq = bam2fastq.out.fastq
    }
    
    if (params.input_fastq) {
        fastq = Channel.fromPath(params.input_fastq)
       }
    
    if (params.input_fastq_tsv){
        input = Channel
            .fromPath(params.input_fastq_tsv)
            .splitCsv(header: true, sep: '\t')
            .map { row -> 
                def lib = row['lib_name']
                def file_path = file(row['file_path'])
                return [lib, file_path]
        }
        add_id2fastq(input)
        input.view()
        fastq = add_id2fastq.out.id2fastq
        fastq.view()
    }
    
    // preprocessing fastp - remove low complexity (repetitive nucleotides) and low quality reads + dedup
    preprocessing(fastq)
    fastp_fastq = preprocessing.out.fastp_fastq
    
    // dustmasking
    dustmasker(fastp_fastq)
    
}
    
// --mapping, this does the bowtie2 mapping.
if (params.mapping) {
    
    // merging fastq files
    input_fastqs = Channel.fromPath(params.fastq_file, type: 'dir')
    mergerA(input_fastqs)
    merged_fastq = mergerA.out.merged_fastq

    // bowtie2 mapping
        input_dirs = Channel.empty()
        // make index channels

        // this is to set a specific groups
        plant = false
        invertebrate = false
        vertebrate_mammalian = false
        vertebrate_other = false


        if (params.all) {
            //input = Channel.fromPath("${params.path_reference_dbs}/*/GCA_*", type: 'dir')
            //  .filter { !['bacteria', 'archaea', 'fungi'].contains(it.getParent().getName()) }
            //input_dirs = input_dirs.concat(input)
            plant = true
            invertebrate = true
            vertebrate_mammalian = true
            vertebrate_other = true
        }

        if (params.vertebrate) {
            //input = Channel.fromPath("${params.path_reference_dbs}/*/GCA_*", type: 'dir')
            //  .filter { !['bacteria', 'archaea', 'fungi'].contains(it.getParent().getName()) }
            //input_dirs = input_dirs.concat(input)
            vertebrate_mammalian = true
            vertebrate_other = true
        }


        if ( params.plant || plant ){
            input = Channel.fromPath("${params.path_reference_dbs}/plant/GCA_*", type: 'dir')
            input_dirs = input_dirs.concat(input)
        }
        if (params.invertebrate || invertebrate ){
            input = Channel.fromPath("${params.path_reference_dbs}/invertebrate/GCA_*", type: 'dir')
            input_dirs = input_dirs.concat(input)
        }
        if ( params.vertebrate_mammalian || vertebrate_mammalian ){
            input = Channel.fromPath("${params.path_reference_dbs}/vertebrate_mammalian/GCA_*", type: 'dir')
            input_dirs = input_dirs.concat(input)
        }
        if ( params.vertebrate_other || vertebrate_other ){
            input = Channel.fromPath("${params.path_reference_dbs}/vertebrate_other/GCA_*", type: 'dir')
            input_dirs = input_dirs.concat(input)
        }
        if ( params.bacteria ){
            input = Channel.fromPath("${params.path_reference_dbs}/bacteria/GCA_*", type: 'dir')
            input_dirs = input_dirs.concat(input)
        }
	if (params.bact_sink){
            input = Channel.fromPath("${params.path_reference_dbs}/bacteria/bact_sink/*concat", type: 'dir')
            input_dirs = input_dirs.concat(input)
        }

    input_index_assembly = input_dirs.map{
            [it, it.baseName]
        }
    input = merged_fastq.combine(input_index_assembly)
    
    bowtie2_mapperA(input)
    assembly_ids = bowtie2_mapperA.out.assembly_ids.collect()
    
    // write sample that have been maps and the database version
    sample_list = mergerA.out.sample_list
    make_log(sample_list, assembly_ids)

}


}