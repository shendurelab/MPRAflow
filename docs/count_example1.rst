.. _Basic count workflow:

.. role:: bash(code)
      :language: bash

=====================
Basic Count workflow
=====================

This example runs the count workflow on 5'/5' WT MRPA data in the HEPG2 cell line from `Klein J., Agarwal, V., Keith, A., et al. 2019 <https://www.biorxiv.org/content/10.1101/576405v1.full.pdf>`_.

Prerequirements
======================

This example depends on the following data and software:


Installation of MPRAflow
----------------------------------------

Please install conda, the MPRAflow environment and clone the actual MPRAflow master branch. You will find more help under :ref:`Installation`.

Producing an association pickle
------------------------------------
This workflow requires a python dictionary of candidate regulatory sequence (CRS) mapped to their barcodes in a pickle format. For this example the file can be generated using :ref:`Association example`.


Reads
----------

There is one condition (HEPG2) with three technical replicates. Each replicate contains a forward (barcode-forward), reverse (barcode-reverse), and index (unique molecular identifier) read for DNA and RNA. These data must be downloaded. All data is publically available on the short read archive (SRA). We will use SRA-toolkit to obtain the data.

.. note:: You need 9 GB disk space to download the data!

.. code-block:: bash

    conda install sra-tools
    mkdir -p Count_Basic/data
    cd Count_Basic/data
    fastq-dump --gzip --split-files SRR10800881 SRR10800882 SRR10800883 SRR10800884 SRR10800885 SRR10800886
    cd ..

For large files and unstable internet connection we reccommend the comand `prefetch` from SRA tools before running `fastq-dump`. This command is much smarter in warnings when something went wrong.

conda install sra-tools
cd Count_Basic/data
prefetch SRR10800881 SRR10800882 SRR10800883 SRR10800884 SRR10800885 SRR10800886
fastq-dump --gzip --split-files SRR10800986
cd ..



.. note:: Please be sure that all files are downloaded completely without errors! Depending on your internet connection this can take a while. If you just want some data to run MPRAflow you can just limit yourself to one condition and/or just one replicate.

With

.. code-block:: bash

    tree data


the folder should look like this:

.. code-block:: text

    data

Here is an overview of the files:

.. csv-table:: HEPG2 data
   :header: "Condition", "GEO Accession", "SRA Accession", SRA Runs
   :widths: 40, 10, 10, 20

   "HEPG2-DNA-1: HEPG2 DNA replicate 1", GSM4237863, SRX7474781, "SRR10800881"
   "HEPG2-RNA-1: HEPG2 RNA replicate 1", GSM4237864, SRX7474782, "SRR10800882"
   "HEPG2-DNA-2: HEPG2 DNA replicate 2", GSM4237865, SRX7474783, "SRR10800883"
   "HEPG2-RNA-2: HEPG2 RNA replicate 2", GSM4237866, SRX7474784, "SRR10800884"
   "HEPG2-DNA-3: HEPG2 DNA replicate 3", GSM4237867, SRX7474785, "SRR10800885"
   "HEPG2-RNA-3: HEPG2 RNA replicate 3", GSM4237868, SRX7474786, "SRR10800886"



MPRAflow
=================================

Now we are close to starting MPRAflow and count the number of barcodes. But before we need to generate an environment csv file to tell nextflow the conditions, replicates and the corresponding reads.

Create experiment.csv
---------------------------

Our experiment file looks exactly like this:

.. code-block:: text

    Condition,Replicate,DNA_BC_F,DNA_UMI,DNA_BC_R,RNA_BC_F,RNA_UMI,RNA_BC_R
    HEPG2,1,SRR10800881_1.fastq.gz,SRR10800881_3.fastq.gz,SRR10800881_2.fastq.gz,SRR10800882_1.fastq.gz,SRR10800882_3.fastq.gz,SRR10800882_2.fastq.gz
    HEPG2,2,SRR10800883_1.fastq.gz,SRR10800883_3.fastq.gz,SRR10800883_2.fastq.gz,SRR10800884_1.fastq.gz,SRR10800884_3.fastq.gz,SRR10800884_2.fastq.gz
    HEPG2,3,SRR10800885_1.fastq.gz,SRR10800885_3.fastq.gz,SRR10800885_2.fastq.gz,SRR10800886_1.fastq.gz,SRR10800886_3.fastq.gz,SRR10800886_2.fastq.gz

Save it into the :code:`Count_Basic/data` folder under :code:`experiment.csv`.

Run nextflow
------------------------------

Now we have everything at hand to run the count MPRAflow pripeline. Therefore we have to be in the cloned MPRAflow folder. But we will change the working and output directory to the :code:`Count_Basic` folder. The MPRAflow count command is:


.. code-block:: bash

    cd <path/to/MPRAflow>/MPRAflow
    conda activate MPRAflow
    nextflow run count.nf -w <path/to/Basic>/Count_Basic/work --experiment-file "<path/to/Basic>/Count_Basic/data/experiment.csv" --dir "<path/to/Basic>/Count_Basic/data" --outdir "<path/to/Basic>/Count_Basic/output" --design "<path/to/design/fasta/design.fa" --association "<path/to/association/pickle/SRR10800986_filtered_coords_to_barcodes.pickle"

.. note:: Please check your :code:`conf/cluster.config` file if it is correctly configured (e.g. with your SGE cluster commands).

If everything works fine the following 5 processes will run: :code:`create_BAM (make idx)` :code:`raw_counts`, :code:`filter_counts`, :code:`final_counts`, :code:`dna_rna_merge_counts`, :code:`calc_correlations`, :code:`make_master_tables`.


Results
-----------------

All output files will be in the :code:`Count_Basic/output` folder.

We expect the program to output the following status when complete:

.. code-block:: text

    start analysis
    executor >  sge (32)
    [23/09474b] process > create_BAM (make idx)    [100%] 6 of 6 ✔
    [0f/4ee034] process > raw_counts (6)           [100%] 6 of 6 ✔
    [01/6ac02f] process > filter_counts (6)        [100%] 6 of 6 ✔
    [4f/b23748] process > final_counts (6)         [100%] 6 of 6 ✔
    [86/4ded79] process > dna_rna_merge_counts (3) [100%] 3 of 3 ✔
    [29/0813f8] process > dna_rna_merge (3)        [100%] 3 of 3 ✔
    [1d/4e7d56] process > calc_correlations (1)    [100%] 1 of 1 ✔
    [9c/4714cb] process > make_master_tables (1)   [100%] 1 of 1 ✔
    Completed at: 07-Jan-2020 04:29:07
    Duration    : 11h 28m 5s
    CPU hours   : 41.5
    Succeeded   : 32
