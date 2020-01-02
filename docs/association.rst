=====================
Association
=====================

.. todo::
  Document association workflow
  
  
conda activate MPRAflow
nextflow run association.nf --fastq-insert "${fastq_prefix}_R1_001.fastq.gz" --design "ordered_candidate_sequences.fa" --fastq-bc "${fastq_prefix}_R2_001.fastq.gz"

