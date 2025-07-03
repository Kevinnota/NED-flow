# NED-flow

NED-flow is a pipeline developed for using <ins>**N**</ins>uclear reference genomes to classify reads obtained from shotgun-sequenced ancient samples <ins>**E**</ins>nvironmental <ins>**D**</ins>NA using a Next<ins>**flow**</ins> backbone. The tool contains tools for database management & building and for preprocessing reads, mapping and taxonomic classifying. The pipeline was designed specifically for nuclear reference genomes/assemblies obtained from GenBank.

In short, NED-flow is downloading all reference genomes deposited on NCBI. Ned-ref-manager.py is written in a way that it will download all reference genomes the first time – and only updates the existing reference genomes subsequently. Only new assemblies and updated assemblies are downloaded. The make-ref-sink.py is a tool that will concatenate bacteria/fungi/archaea genomes into sinks that are used for mapping and removing potential microorganism contamination in reference genomes. The ned.nf --build is a Nextflow tool that will index only the genomes that need to be indexed. It can do this both locally and on an SGE/SLURM cluster.

> [!WARNING]
> Downloading all plant, vertebrate, and invertebrate genomes available at present requires <ins>**~23T**</ins> of free disk space. Before starting to download and build the database, make sure there is enough space available. 

## Installation

NED-flow requires Nextflow version >=22.10. Besides that it requires the instalation of a list of python libraries. All required libraries can be found in the ned-py-install.txt file and can be installed with pip.

One by one.
```
pip intsall <package_name>
```
or all at the same time.

```
pip install -r ned-py-install.txt
```

## Database management


### ned-ref-manager.py
```
~/NED-flow/ned-ref-manager.py -db plant
```

### make-ref-sink.py 

### ned.nf --build

## preproccessing and mapping
### ned.nf --preproccessing

### ned.nf --mapping

## Taxonomic classifying
### ned-classifier.py
### Output table info