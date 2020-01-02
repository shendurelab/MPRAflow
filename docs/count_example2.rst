.. _Count for Saturation Mutagenesis of the TERT promoter:

.. role:: bash(code)
   :language: bash

=====================================================================
Count for Saturation Mutagenesis of the TERT promoter
=====================================================================


This example runs the count workflow on saturation mutagenesis data of the TERT promoter from `Kircher et al. 2019 <https://doi.org/10.1038/s41467-019-11526-w>`_. It will only do the count part. For unraveling single varaint effects please go to the example: :ref:`Saturation mutagenesis of the TERT promoter`. The same saturation mutagenesis library was used in four different experiments. We will focus only on the experiments in HEK293T and in glioblastoma SF7996 (GBM) cells (two different conditions).


Prerequirements
======================

This example depends on the following data and software:


Installation of MPRAflow
----------------------------------------

Please install conda, the MPRAflow environment and clone the actual MPRAflow master branch. You will find more help under :ref:`Installation`.


Reads
--------

We have two conditions (HEK293T and GBM cells). Each of them has three technical replicates and one replicate exists of forward, reverse and index reads for DNA and RNA. These data has to be downloaded. All data is public available on the short read archive (SRA). We will use the SRA-tools to download the reads.

.. note:: You need at least XX GB disk space to download the data!

.. code-block:: bash

    conda install sra-tools
    mkdir -p Count_TERT/data
    cd Count_TERT/data
    fastq-dump --gzip --split-files SRR8647059 SRR8647060 SRR8647061 SRR8647062 SRR8647063 SRR8647064 SRR8647119 SRR8647120 SRR8647121 SRR8647122 SRR8647123 SRR8647124 SRR8647125 SRR8647126 SRR8647127 SRR8647128 SRR8647129 SRR8647130
    cd ..


.. note:: Please be sure that all files are downloaded completely without errors! Depending on your internet connection this can take a while. If you just want some data to run MPRAflow you can just limit yourself to one condition and/or just one replicate.

With

.. code-block:: bash

    tree data


the folder should look like this:

.. code-block:: text

    tree data



Here is an overview of the files:

.. csv-table:: TERT data
   :header: "Condition", "GEO Accession", "SRA Accession", SRA Runs
   :widths: 40, 10, 10, 20

   "TERT-GBM-DNA-1: TERT-GBM transfection DNA replicate 1", GSM3604284, SRX5444854, "SRR8647059"
   "TERT-GBM-DNA-2: TERT-GBM transfection DNA replicate 2", GSM3604285, SRX5444855, "SRR8647060"
   "TERT-GBM-DNA-3: TERT-GBM transfection DNA replicate 3", GSM3604286, SRX5444856, "SRR8647061"
   "TERT-GBM-RNA-1: TERT-GBM transfection RNA replicate 1", GSM3604287, SRX5444857, "SRR8647062"
   "TERT-GBM-RNA-2: TERT-GBM transfection RNA replicate 2", GSM3604288, SRX5444858, "SRR8647063"
   "TERT-GBM-RNA-3: TERT-GBM transfection RNA replicate 3", GSM3604289, SRX5444859, "SRR8647064"
   "TERT-HEK-DNA-1: TERT-HEK transfection DNA replicate 1", GSM3604302, SRX5444888, "SRR8647119, SRR8647120"
   "TERT-HEK-DNA-2: TERT-HEK transfection DNA replicate 2", GSM3604303, SRX5444889, "SRR8647121, SRR8647122"
   "TERT-HEK-DNA-3: TERT-HEK transfection DNA replicate 3", GSM3604304, SRX5444890, "SRR8647123, SRR8647124"
   "TERT-HEK-RNA-1: TERT-HEK transfection RNA replicate 1", GSM3604305, SRX5444891, "SRR8647125, SRR8647126"
   "TERT-HEK-RNA-2: TERT-HEK transfection RNA replicate 2", GSM3604306, SRX5444892, "SRR8647127, SRR8647128"
   "TERT-HEK-RNA-3: TERT-HEK transfection RNA replicate 3", GSM3604307, SRX5444893, "SRR8647129, SRR8647130"


.. note:: The runs SRR8647120, SRR8647122, SRR8647124, SRR8647126, SRR8647128, and SRR8647130 have two index reads (forward and reverse). In the publication by `Kircher et al. 2019 <https://doi.org/10.1038/s41467-019-11526-w>`_ merging and trimming is used to combine them. For simplification we will discard the reverse index reads: :bash:`rm data/*_4.fastq.gz`

.. code-block:: bash



Also two different sequencing runs where made in condition TERT-HEK. Therefore We have to combine the reads:

.. code-block:: bash

    for i in 1 2 3; do
       zcat {SRR8647119,SRR8647120}_$i.fastq.gz | gzip -c > SRR8647119_SRR8647120_$i.fastq.gz;
       zcat {SRR8647121,SRR8647122}_$i.fastq.gz | gzip -c > SRR8647121_SRR8647122_$i.fastq.gz;
       zcat {SRR8647123,SRR8647124}_$i.fastq.gz | gzip -c > SRR8647123_SRR8647124_$i.fastq.gz;
       zcat {SRR8647125,SRR8647126}_$i.fastq.gz | gzip -c > SRR8647125_SRR8647126_$i.fastq.gz;
       zcat {SRR8647127,SRR8647128}_$i.fastq.gz | gzip -c > SRR8647127_SRR8647128_$i.fastq.gz;
       zcat {SRR8647129,SRR8647130}_$i.fastq.gz | gzip -c > SRR8647129_SRR8647130_$i.fastq.gz;
    done

MPRAflow
=================================

Now we are close to start MPRAflow and count the number of barcodes. But before we need to generate an environment csv file to tell nextflow the conditions, replicates and the corresponding reads.

Create environment.csv
---------------------------

Our environment file looks exactly like this:

.. code-block:: text

  Condition,Replicate,DNA_R1,DNA_R2,DNA_R3,RNA_R1,RNA_R2,RNA_R3
  TERT-GBM, 1, SRR8647059_1.fastq.gz,SRR8647059_2.fastq.gz,SRR8647059_3.fastq.gz, SRR8647062_1.fastq.gz,SRR8647062_2.fastq.gz,SRR8647062_3.fastq.gz
  TERT-GBM, 2, SRR8647060_1.fastq.gz,SRR8647060_2.fastq.gz,SRR8647060_3.fastq.gz, SRR8647063_1.fastq.gz,SRR8647063_2.fastq.gz,SRR8647063_3.fastq.gz
  TERT-GBM, 3, SRR8647061_1.fastq.gz,SRR8647061_2.fastq.gz,SRR8647061_3.fastq.gz, SRR8647064_1.fastq.gz,SRR8647064_2.fastq.gz,SRR8647064_3.fastq.gz
  TERT-HEK, 1,
  TERT-HEK, 2,
  TERT-HEK, 3,

Save it into the :code:`Count_TERT/data` folder under :code:`environment.csv`

Run nextflow
------------------------------


Results
-----------------



.. todo::
  add count satMut example
