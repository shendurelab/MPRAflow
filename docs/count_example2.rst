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

.. note:: You need 16 GB disk space to download the data!

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

    data
    ├── SRR8647059_1.fastq.gz
    ├── SRR8647059_2.fastq.gz
    ├── SRR8647059_3.fastq.gz
    ├── SRR8647060_1.fastq.gz
    ├── SRR8647060_2.fastq.gz
    ├── SRR8647060_3.fastq.gz
    ├── SRR8647061_1.fastq.gz
    ├── SRR8647061_2.fastq.gz
    ├── SRR8647061_3.fastq.gz
    ├── SRR8647062_1.fastq.gz
    ├── SRR8647062_2.fastq.gz
    ├── SRR8647062_3.fastq.gz
    ├── SRR8647063_1.fastq.gz
    ├── SRR8647063_2.fastq.gz
    ├── SRR8647063_3.fastq.gz
    ├── SRR8647064_1.fastq.gz
    ├── SRR8647064_2.fastq.gz
    ├── SRR8647064_3.fastq.gz
    ├── SRR8647119_1.fastq.gz
    ├── SRR8647119_2.fastq.gz
    ├── SRR8647119_3.fastq.gz
    ├── SRR8647120_1.fastq.gz
    ├── SRR8647120_2.fastq.gz
    ├── SRR8647120_3.fastq.gz
    ├── SRR8647120_4.fastq.gz
    ├── SRR8647121_1.fastq.gz
    ├── SRR8647121_2.fastq.gz
    ├── SRR8647121_3.fastq.gz
    ├── SRR8647122_1.fastq.gz
    ├── SRR8647122_2.fastq.gz
    ├── SRR8647122_3.fastq.gz
    ├── SRR8647122_4.fastq.gz
    ├── SRR8647123_1.fastq.gz
    ├── SRR8647123_2.fastq.gz
    ├── SRR8647123_3.fastq.gz
    ├── SRR8647124_1.fastq.gz
    ├── SRR8647124_2.fastq.gz
    ├── SRR8647124_3.fastq.gz
    ├── SRR8647124_4.fastq.gz
    ├── SRR8647125_1.fastq.gz
    ├── SRR8647125_2.fastq.gz
    ├── SRR8647125_3.fastq.gz
    ├── SRR8647126_1.fastq.gz
    ├── SRR8647126_2.fastq.gz
    ├── SRR8647126_3.fastq.gz
    ├── SRR8647126_4.fastq.gz
    ├── SRR8647127_1.fastq.gz
    ├── SRR8647127_2.fastq.gz
    ├── SRR8647127_3.fastq.gz
    ├── SRR8647128_1.fastq.gz
    ├── SRR8647128_2.fastq.gz
    ├── SRR8647128_3.fastq.gz
    ├── SRR8647128_4.fastq.gz
    ├── SRR8647129_1.fastq.gz
    ├── SRR8647129_2.fastq.gz
    ├── SRR8647129_3.fastq.gz
    ├── SRR8647130_1.fastq.gz
    ├── SRR8647130_2.fastq.gz
    ├── SRR8647130_3.fastq.gz
    └── SRR8647130_4.fastq.gz

    0 directories, 60 files



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

Also two different sequencing runs where made in condition TERT-HEK using the same library. Therefore, we have to process both runs together. So we will combine the reads. But in this csase bot sequencing runs have different read lengths (20 and 50 bp). Different read length will make false assumption when merging paired-end reads. Therefore we cut them down to 20 bp.

.. code-block:: bash

    for i in 1 2 3; do
       zcat data/{SRR8647119,SRR8647120}_$i.fastq.gz | cut -c 1-20 | gzip -c > data/SRR8647119_SRR8647120_$i.fastq.gz;
       zcat data/{SRR8647121,SRR8647122}_$i.fastq.gz | cut -c 1-20 | gzip -c > data/SRR8647121_SRR8647122_$i.fastq.gz;
       zcat data/{SRR8647123,SRR8647124}_$i.fastq.gz | cut -c 1-20 | gzip -c > data/SRR8647123_SRR8647124_$i.fastq.gz;
       zcat data/{SRR8647125,SRR8647126}_$i.fastq.gz | cut -c 1-20 | gzip -c > data/SRR8647125_SRR8647126_$i.fastq.gz;
       zcat data/{SRR8647127,SRR8647128}_$i.fastq.gz | cut -c 1-20 | gzip -c > data/SRR8647127_SRR8647128_$i.fastq.gz;
       zcat data/{SRR8647129,SRR8647130}_$i.fastq.gz | cut -c 1-20 | gzip -c > data/SRR8647129_SRR8647130_$i.fastq.gz;
    done

.. note:: If you combine multiple sequence runs (e.g. you need more reads) you have to combine the reads before. Otherwise barcodes with the same UMI can be count twice. But it is important that all read lengths are the same. The easiest workaround is to cut them down to the minimum length. If you have a different library (but with the same barcodes) you should run the count utility with both runs separately using different conditions. Later you have to combine the final counts of both conditions.

MPRAflow
=================================

Now we are close to start MPRAflow and count the number of barcodes. But before we need to generate an :code:`environment.csv` file to tell nextflow the conditions, replicates and the corresponding reads.

Create experiment.csv
---------------------------

Our experiment file looks exactly like this:

.. code-block:: text

    Condition,Replicate,DNA_BC_F,DNA_UMI,DNA_BC_R,RNA_BC_F,RNA_UMI,RNA_BC_R
    TERT-GBM,1,SRR8647059_1.fastq.gz,SRR8647059_3.fastq.gz,SRR8647059_2.fastq.gz,SRR8647062_1.fastq.gz,SRR8647062_3.fastq.gz,SRR8647062_2.fastq.gz
    TERT-GBM,2,SRR8647060_1.fastq.gz,SRR8647060_3.fastq.gz,SRR8647060_2.fastq.gz,SRR8647063_1.fastq.gz,SRR8647063_3.fastq.gz,SRR8647063_2.fastq.gz
    TERT-GBM,3,SRR8647061_1.fastq.gz,SRR8647061_3.fastq.gz,SRR8647061_2.fastq.gz,SRR8647064_1.fastq.gz,SRR8647064_3.fastq.gz,SRR8647064_2.fastq.gz
    TERT-HEK,1,SRR8647119_SRR8647120_1.fastq.gz,SRR8647119_SRR8647120_3.fastq.gz,SRR8647119_SRR8647120_2.fastq.gz,SRR8647125_SRR8647126_1.fastq.gz,SRR8647125_SRR8647126_3.fastq.gz,SRR8647125_SRR8647126_2.fastq.gz
    TERT-HEK,2,SRR8647121_SRR8647122_1.fastq.gz,SRR8647121_SRR8647122_3.fastq.gz,SRR8647121_SRR8647122_2.fastq.gz,SRR8647127_SRR8647128_1.fastq.gz,SRR8647127_SRR8647128_3.fastq.gz,SRR8647127_SRR8647128_2.fastq.gz
    TERT-HEK,3,SRR8647123_SRR8647124_1.fastq.gz,SRR8647123_SRR8647124_3.fastq.gz,SRR8647123_SRR8647124_2.fastq.gz,SRR8647129_SRR8647130_1.fastq.gz,SRR8647129_SRR8647130_3.fastq.gz,SRR8647129_SRR8647130_2.fastq.gz


Save it into the :code:`Count_TERT/data` folder under :code:`experiment.csv`.

Run nextflow
------------------------------

Now we have everything at hand to run the count MPRAflow pipeline. Therefore we have to be in the cloned MPRAflow folder. But we will change the working and output directory to the :code:`Count_TERT` folder. For the TERT example the barcode length is 20 bp and the UMI length 10 bp. The MPRAflow count command is:


.. code-block:: bash

    cd <path/to/MPRAflow>/MPRAflow
    conda activate MPRAflow
    nextflow run -resume -w <path/to/TERT>/Count_TERT/work  count.nf --experiment-file "<path/to/TERT>/Count_TERT/data/experiment.csv" --dir "<path/to/TERT>/Count_TERT/data" --outdir "<path/to/TERT>/Count_TERT/output" --bc-length 20 --umi-length 10

.. note:: Please check your :code:`conf/cluster.config` file if it is correctly configured (e.g. with your SGE cluster commands).

If everything works fine the following 5 processes will run: :code:`create_BAM (make idx)` :code:`raw_counts`, :code:`filter_counts`, :code:`final_counts`, :code:`dna_rna_merge_counts`.

.. code-block:: text

    [fe/d8ac14] process > create_BAM (make idx)    [100%] 12 of 12 ✔
    [7d/b56129] process > raw_counts (12)          [100%] 12 of 12 ✔
    [06/2c938d] process > filter_counts (12)       [100%] 12 of 12 ✔
    [2d/ce1afe] process > final_counts (12)        [100%] 12 of 12 ✔
    [68/df8db0] process > dna_rna_merge_counts (6) [100%] 6 of 6 ✔
    Completed at: 09-Jan-2020 15:38:32
    Duration    : 3h 45m 17s
    CPU hours   : 21.8
    Succeeded   : 54


Results
-----------------

All needed output files will be in the :code:`Count_TERT/output` folder. In this tutorial we are only interested in the counts per barcode, because we can use these outputs in the :ref:`Saturation mutagenesis of the TERT promoter` tutorial.

.. code-block:: bash

    cd <path/to/TERT>/Count_TERT
    tree output -P "*[123]_counts.tsv.gz"

.. code-block:: text

    output
    ├── TERT-GBM
    │   ├── 1
    │   │   └── TERT-GBM_1_counts.tsv.gz
    │   ├── 2
    │   │   └── TERT-GBM_2_counts.tsv.gz
    │   └── 3
    │       └── TERT-GBM_3_counts.tsv.gz
    └── TERT-HEK
        ├── 1
        │   └── TERT-HEK_1_counts.tsv.gz
        ├── 2
        │   └── TERT-HEK_2_counts.tsv.gz
        └── 3
            └── TERT-HEK_3_counts.tsv.gz

    8 directories, 6 files

The count files are tab separated and contain the barcode, the number number of unique UMI DNA counts and the umber of unique RNA counts. E.g. this is an example count file:

.. code-block:: bash

     zcat output/TERT-GBM/1/TERT-GBM_1_counts.tsv.gz | head

.. code-block:: text

    AAAAAAAAAAAAAAA 2       76
    AAAAAAAAAAAAAAC 1       2
    AAAAAAAAAAAAAAT 1       7
    AAAAAAAAAAAAAGA 1       5
    AAAAAAAAAAAAATA 2       7
    AAAAAAAAAAAAATG 2       6
    AAAAAAAAAAAAATT 1       2
    AAAAAAAAAAAATGA 3       1
    AAAAAAAAAAAATGG 1       3
    AAAAAAAAAAAATTA 1       2
