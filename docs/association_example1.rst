.. _Association example:

.. role:: bash(code)
   :language: bash

============================
Basic association workflow
============================

This example runs the association workflow on 5'/5' WT MRPA data in the HEPG2 cell line from `Klein J., Agarwal, V., Keith, A., et al. 2019 <https://www.biorxiv.org/content/10.1101/576405v1.full.pdf>`_.

Prerequirements
======================

This example depends on the following data and software:


Installation of MPRAflow
----------------------------------------

Please install conda, the MPRAflow environment and clone the actual MPRAflow master branch. You will find more help under :ref:`Installation`.

Meta Data
___________

It is necessary to get the ordered oligo array so that each enhancer sequence can be labeled in the analysis and to trim any adaptors still in the sequence, in this case we trim off 15bp from the end of each sequence

.. code-block:: bash

    mkdir -p Assoc_Basic/data
    cd Assoc_Basic/data
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4237nnn/GSM4237954/suppl/GSM4237954_9MPRA_elements.fa.gz

    zcat GSM4237954_9MPRA_elements.fa.gz |awk '{ count+=1; if (count == 1) { print } else { print substr($1,1,171)}; if (count == 2) { count=0 } }' > design.fa

Reads
----------

There is one set of association sequencing for this data, which contains a forward (CRS-forward), reverse (CRS-reverse), and index (barcode) read for DNA and RNA. These data must be downloaded. All data is publically available on the short read archive (SRA). We will use SRA-toolkit to obtain the data.

.. note:: You need 10 GB disk space to download the data!

.. code-block:: bash

    conda install sra-tools
    cd Assoc_Basic/data
    fastq-dump --gzip --split-files SRR10800986
    cd ..

For large files and unstable internet connection we reccommend the comand `prefetch` from SRA tools before running `fastq-dump`. This command is much smarter in warnings when something went wrong.

.. code-block:: bash

    conda install sra-tools
    cd Assoc_Basic/data
    prefetch SRR10800986
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

.. csv-table:: HEPG2 association data
   :header: "Condition", "GEO Accession", "SRA Accession", SRA Runs
   :widths: 40, 10, 10, 20

   "HEPG2-association: HEPG2 library association", GSM4237954, SRX7474872, "SRR10800986"


MPRAflow
=================================

Now we are ready to run MPRAflow and create CRS-barcode mappings.

Run nextflow
------------------------------

Now we have everything at hand to run the count MPRAflow pripeline. Therefore we have to be in the cloned MPRAflow folder. But we will change the working and output directory to the :code:`Assoc_Basic` folder. The MPRAflow count command is:


.. code-block:: bash

    cd <path/to/MPRAflow>/MPRAflow
    conda activate MPRAflow
    nextflow run association.nf -w <path/to/Basic>/Assoc_Basic/work --fastq-insert "<path/to/Basic>/Assoc_Basic/data/SRR10800986_1.fastq.gz" --fastq-insertPE "<path/to/Basic>/Assoc_Basic/data/SRR10800986_3.fastq.gz" --fastq-bc "<path/to/Basic>/Assoc_Basic/data/SRR10800986_2.fastq.gz" --design "<path/to/Basic>/Assoc_Basic/data/design.fa" --name assoc_basic

.. note:: Please check your :code:`conf/cluster.config` file if it is correctly configured (e.g. with your SGE cluster commands).

If everything works fine the following 7 processes will run: :code:`count_bc_nolab` :code:`create_BWA_ref`, :code:`PE_merge`, :code:`align_BWA_PE`, :code:`collect_chunks`, :code:`map_element_barcodes`, :code:`filter_barcodes`.


Results
-----------------

All needed output files will be in the :code:`Assoc_Basic/output` folder.
