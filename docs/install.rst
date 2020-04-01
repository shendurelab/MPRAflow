.. _Installation:

=====================
Installation
=====================

Installation should take less than 10 minutes

System Requirements
===================

CentOS Linux 7 or above

Required packages
==================

.. code-block:: bash

  	conda 4.6 or above

Download here: https://docs.conda.io/en/latest/miniconda.html

Clone repository
=================

.. code-block:: bash

    git clone https://github.com/shendurelab/MPRAflow.git

Set up conda environment
========================

This pipeline uses python2.7 and python3.6 with additional R scripts. Three `.yml` files are provided to create the appropriate environments and is completely handled by nextflow. The whole pipeline is set up to run on a Linux system. The general environment with nextflow located in the home directory called `environment.yml`.

Install the the conda environment. The general conda environment is called MPRAflow.

.. code-block:: bash

    cd MPRAflow
    conda env create -n MPRAflow -f environment.yml

If you do not have access to the internet, you have to run the previous command on a node with internet. Afterwards you need to start nextflow too (see Steps to run the pipeline). After creation of the second conda environment by nextflow you can cancel it and start it on your internal node. Be aware that folders must have access on all nodes.

Nextflow has problems using conda 4.7 and highet, because the source activate command is replaced by conda activate. If you get error messages after running you can make a symbolik link of the activate command from you bin folder of the conda or miniconda folder to your MPRAflow environment bin folder. E.g. like:

.. code-block:: bash

    ln -s ~/miniconda3/bin/activate ~/miniconda3/envs/MPRAflow/bin/activate

Quick test
============

.. code-block:: bash

    conda activate MPRAflow
    nextflow run count.nf --help
    nextflow run association.nf --help
    nextflow run saturationMutagenesis.nf --help
