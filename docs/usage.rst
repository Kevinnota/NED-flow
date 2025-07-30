Usage
=====

This section explains how to install, configure, and run NED-flow.

Installation
------------

To install NED-flow, clone the repository and set up the environment:

.. code-block:: bash

   git clone https://github.com/kevinnota/NED-flow.git
   cd NED-flow
   conda env create -f environment.yml
   conda activate nedflow

If you're using it on a cluster (SGE or SLURM), make sure Nextflow is properly configured.

Running the pipeline
--------------------

To build the reference index (e.g., for the first run):

.. code-block:: bash

   nextflow run ned.nf --build

To run the full pipeline on your data:

.. code-block:: bash

   nextflow run ned.nf --reads data/*.fastq.gz --db refdb/ --out results/

Arguments
---------

``--build``
    Builds or updates the reference database using GenBank.

``--reads``
    Path to input FASTQ files (can be gzipped).

``--db``
    Path to the reference genome directory.

``--out``
    Output directory where results will be stored.

``--profile``
    Specify execution environment (e.g., `local`, `sge`, `slurm`).

Updating the reference database
-------------------------------

NED-flow uses `ned-ref-manager.py` to download and update nuclear reference genomes.

.. code-block:: bash

   python scripts/ned-ref-manager.py --update

This will only download **new** or **updated** assemblies from GenBank.

Removing microbial contamination
--------------------------------

Use `make-ref-sink.py` to concatenate microbial reference genomes into "sink" genomes:

.. code-block:: bash

   python scripts/make-ref-sink.py --input microbe_refs/ --output sinks/

This helps map and filter out possible contamination in your target nuclear reference database.

Tips
----

- Always use the same reference DB version for consistent results.
- Use SLURM or SGE profiles for large datasets or clusters.
- Review logs in `logs/` to debug issues.
