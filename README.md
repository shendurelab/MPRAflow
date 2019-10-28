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
git clone https://github.com/JJZhao123/MPRA_nextflow
```

### Set up conda environment:
This pipeline uses python2.7 and python3.6 and is set up to run on a linux system.   
Two .yml files are provided to create the appropriate environments: mpra.yml py36.yml.
After installing conda on your system create environment by running code below.

```bash
conda env create -f mpra.yml
conda env create -f py36.yml
```
##### ****WARNING: currently certain required packages are not available in conda 4.7 as they have discontinued the 'free' channel (https://www.anaconda.com/why-we-removed-the-free-channel-in-conda-4-7/). If you are running 4.7 please follow the instructions below:
```bash
conda config --set restore_free_channel true
conda env create -f mpra.yml
conda config --set restore_free_channel false
conda env create -f py36.yml
```


### Create config files to point system paths to conda c++ libraries and compilers
#### for each of the environments follow the steps below, where `$CONDA_PREFIX` is the location of your conda environment 
(ex `miniconda3/envs/mpra`)


1. Locate the directory for the conda environment in your terminal window by running in the terminal 
    ```bash
    echo $CONDA_PREFIX
    ```

2. Enter that directory and create these subdirectories and files:

    ```bash 
    cd $CONDA_PREFIX
    mkdir -p ./etc/conda/activate.d
    mkdir -p ./etc/conda/deactivate.d
    touch ./etc/conda/activate.d/env_vars.sh
    touch ./etc/conda/deactivate.d/env_vars.sh
    ```


3. Edit `./etc/conda/activate.d/env_vars.sh` as follows:

    ```bash
    #!/bin/sh
    export LD_LIBRARY_PATH=~/miniconda3/lib/ 
    ```

4. Edit `./etc/conda/deactivate.d/env_vars.sh` as follows:

    ```bash
    #!/bin/sh
    unset LD_LIBRARY_PATH
    ```

### Install Shedskin, a python compiler 
#### documentation: https://shedskin.github.io/

1. Download the .tgz file and run the following commands:

    ```bash
    conda activate mpra
    tar zxvf shedskin-0.9.4.tgz
    python setup.py install 
    shedskin test.py
    make
    ./test
    ```

##### **** if this results in printing hello world!, congrats you are done!! If not you likely are lacking shared object files (*.so). See below.
-----------------------------------------
##### Within mpra_py27 env, we hwill need to modify the pipeline2.0/Makefile to point to the GC library.

    shedskin -L ~/miniconda3/envs/mpra_py27/include test
    sed -i '3s|$| -L ~/miniconda3/envs/mpra_py27/lib|' Makefile

------------------------

##### If you are working on UCSF's Wynton cluster, I've had the admins do this already, but if you are not you might need to install the Boehm garbage collector from this website (http://www.hboehm.info/gc/) and build it yourself. If you are working on a shared system you may not have permissions to edit these files- in that case please contact your administrator. 


#### After completing these steps you should have an environment that is ready for running the bulk processing pipeline. 
    
 

## Running the pipeline

#### Steps to run the pipeline

This pipeline comes with a `nextflow.config` file to run on SGE systems, allowing each process to be run as a separate 'qsub' command. Please use a submit script for steps 2 and 3. For help messages run:

   ```bash
   conda activate py36
   nextflow run association/main_gg.nf --help
   nextflow run dna_rna_count/main_gg.nf --help
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
    
3. Run main_gg.txt in the association folder

   ```bash 
   cd association
   conda activate py36
   nextflow run main_gg.nf --fastq_insert ${fastq}_R1_001.fastq.gz --design pilot_library_noprimer.fa" --fastq_bc ${fastq}_R2_001.fastq.gz" --condaloc '~/miniconda3/bin/activate'
   ```

4. Run main_gg.txt in the dna_rna_count folder

   ```bash 
   cd dna_rna_count
   conda activate py36
   nextflow run main_gg.nf --dir bulk_FASTQ_directory --e experiment.csv --design pilot_library_noprimer.fa --association output_filtered_coords_to_barcodes.p --condaloc '~/miniconda3/bin/activate'
   ```
   




![MPRA_nextflow](https://github.com/shendurelab/MPRAflow/blob/master/MPRA_nextflow.png)
