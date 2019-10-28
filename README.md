[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.30.1-brightgreen.svg)](https://www.nextflow.io/)

# MPRA_nextflow

#### This pipeline processes sequencing data (FASTA files) from Massively Parallel Reporter Assays (MPRA) to create count tables for candidate sequences tested in the experiment. 

### This package contains two utilites:

### ASSOCIATION:
##### This utility takes in library association sequencing data (FASTQ) and a design file (FASTA) to assign enhancers to barcodes. Functionality includes filtering for quality and coverage of barcodes. This utility must be run before the COUNT utility. 

### COUNT:
##### This utility processes sequence data (FASTQ) of barcodes from the DNA and RNA fractions of the MPRA experiment and outputs count tables labeled with the candidate sequence and a label provided in the design file. This utility can process multiple replicates (batches), and conditions in a parallelized manner combining the results into a single count matrix compatible with MPRAnalyze. 


## Installation

### Required packages

- Samtools
- BWA and/or Bowtie
- Picard
- BCFtools
- conda
- shedskin

### Clone repository 

```bash
git clone https://github.com/shendurelab/MPRAflow.git
```

### Set up conda environment:
This pipeline uses python2.7 and python3.6 and is set up to run on a linux system.   
Two .yml files are provided to create the appropriate environments: mpra.yml py36.yml.
After installing conda on your system create environment by running code below.

```bash
cd ~/MPRAflow/conf
conda env create -f mpraflow_py27.yml
conda env create -f mpraflow_py36.yml
```
##### ****WARNING: currently certain required packages are not available in conda 4.7 as they have discontinued the 'free' channel (https://www.anaconda.com/why-we-removed-the-free-channel-in-conda-4-7/). If you are running 4.7 please follow the instructions below:
```bash
conda config --set restore_free_channel true
conda env create -f mpra.yml
conda config --set restore_free_channel false
conda env create -f py36.yml
```

## Running the pipeline

#### Steps to run the pipeline

This pipeline comes with a `nextflow.config` file to run on SGE systems, allowing each process to be run as a separate 'qsub' command. Please use a submit script for steps 2 and 3. For help messages run:

   ```bash
   conda activate py36
   nextflow run count.nf --help
   nextflow run association.nf --help
   ```

1. Create an 'expermient' csv in the format below:
 
   ```
   condition,batch,dna,rna,name
   cell1,1,DNA FASTQ prefix, RNA FASTQ prefix, desired name
   cell1,2,DNA FASTQ prefix, RNA FASTQ prefix, desired name
   cell2,1,DNA FASTQ prefix, RNA FASTQ prefix, desired name
   cell2,2,DNA FASTQ prefix, RNA FASTQ prefix, desired name
   ```
   an example file is here: `dna_rna_count/experiment.csv`

2. Create a 'label' tsv in the format below:
 
   ```
   insert1_name	insert1_label
   insert2_name insert2_label
   ```
   The insert names must exactly match the names in the design FASTA file
    
3. Run Assocaition

   ```bash 
   cd ~/MPRAflow
   conda activate mpraflow_py36
   nextflow run association.nf --fastq_insert ${fastq}_R1_001.fastq.gz --design pilot_library_noprimer.fa" --fastq_bc ${fastq}_R2_001.fastq.gz" --condaloc '~/miniconda3/bin/activate'
   ```

4. Run Count

   ```bash 
   cd ~/MPRAflow
   conda activate mpraflow_py36
   nextflow run count.nf --dir bulk_FASTQ_directory --e experiment.csv --design pilot_library_noprimer.fa --association output_filtered_coords_to_barcodes.p --condaloc '~/miniconda3/bin/activate'
   ```
   




![MPRA_nextflow](https://github.com/shendurelab/MPRAflow/blob/master/MPRA_nextflow.png)
