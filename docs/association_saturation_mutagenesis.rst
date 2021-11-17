.. _Association saturation mutagenesis:

==========================================
Association saturation mutagenesis
==========================================

This workflow is about assocation variant calls with barcodes. 
Variants are introduced by an error-prone PCR. 
The workflow takes the sequencing of the region, with barcodes in index read and the reference sequence and maps the reads to the reference, calls variants and associates them with the corresponding barcode.

Input files
===============

Fastq Files
-----------
3 Fastq files from library association sequencing
--Reference sequencing, 1 forward and 1 reverse read
--Barcode sequence, 1 read covering the barcode

Reference file
---------------
Fasta file of  the referencesequence describing the mutated sequence

Example file:

.. code-block:: text

  >TERT
  GATCTGCGATCTAAGTAAGCCCAGGACCGCGCTTCCCACGTGGCGGAGGGACTGGGGACCCGGGCACCCGTCCTGCCCCT
  TCACCTTCCAGCTCCGCCTCCTCCGCGCGGACCCCGCCCCGTCCCGACCCCTCCCGGGTCCCCGGCCCAGCCCCCTCCGG
  GCCCTCCCAGCCCCTCCCCTTCCTTTCCGCGGCCCCGCCCTCTCCTCGCGGCGCGAGTTTCAGGCAGCGCTGCGTCCTGC
  TGCGCACGTGGGAAGCCCTGGCCCCGGCCACCCCCGCGAAAGCTTGCATGCCCTGCAGG

association_saturationMutagenesis.nf
============================

Options
---------------

With :code:`--help` or :code:`--h` you can see the help message.

Usage:

**Mandatory arguments:**
  --fastq-insert                Full path to library association fastq for insert
  --fastq-insertPE              Full path to library association fastq for read2
  --fastq-bc                    Full path to library association fastq for bc
  --design                      Full path to fasta of reference sequence (only one reference sequence)
  --name                        Name of the association. Files will be named after this.

**Optional:**
  --bc-length                   Barcode length (default 15)
  --clipping-penalty            bwa mem clipping penalty (default 80)
  --min-ireads                  minimum number gapped reads for indel candidates (default: 3)
  --split                       Split up the fastq read into chunks with max limit of reads (default: 2000000)
  --outdir                      The output directory where the results will be saved and what will be used as a prefix (default outs)
  --h, --help                   Print the help message

Processes
-------------

Processes run by nextflow in the Association Saturation Mutagenesis Utility.

clean_design
  Removes any illegal characters (defined by Piccard) in the reference file.

create_BWA_ref
  Creates a BWA reference based on the design file

get_name
  Recieves the name written in the header of the reference fasta file

create_BAM
  Merges the forward and reverse reads and stores them in BAM format. This is run for each chunk.

collect_chunks
  combine the chunks (merge bams)

PE_mapping
  Map the reads to the reference.

get_count
  Get the barcodes and count for each barcode. Filter them

extract_reads
  Create a BAM file for each barcode with the corresponding reads in it. This create a large number of files!

call_variants:
  Call variants for each barcode seperately.

combine_variants:
  Combine barcode/variant calls to one final output file.



Output
==========

The output can be found in the folder defined by the option :code:`--outdir` and the subfolder definedby option :code:`--name`. It is structured in folders of the condition as

Files
-------------


design_rmIllegalChars.fa
  Reference file with remove illegal characters
${datasetID}.variants.txt.gz
  Barcode to variant association. Named by the header in the reference.
counts_${datasetID}.filtered.tsv.gz
  Filtered barcode counts. Named by the header in the reference.
counts_${datasetID}.tsv.gz
  All barcode counts. Named by the header in the reference.
