.. _ned_mapper-page:

NED mapping
==========

The idea of NED-flow is that all reads are mapped indipendent to each reference genomes/assemblies in the database. This will create one bam file per mapping. To reduce the number of files, all reads are concatinated into one fastq file before mapping. This is all done automatically, but NED-flow does require the reads to have the library_id, or some kind of identifier at the end of the read header. This is used by downstream classifier for corretly assigning a read to a sample. This is also automatically done, but requires an input file with fastq_name and library_identifier. 
