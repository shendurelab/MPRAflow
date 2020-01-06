.. _Saturation mutagenesis of the TERT promoter:

============================================
Saturation mutagenesis of the TERT promoter
============================================

This example runs the saturation mutagenesis workflow on saturation mutagenesis data of the TERT promoter from `Kircher et al. 2019 <https://doi.org/10.1038/s41467-019-11526-w>`_.
The same saturation mutagenesis library was used in four different experiments.
We will use the experiments in HEK293T and in glioblastoma SF7996 (GBM) cells in this workflow to see differences between the two cell lines (conditions).

Prerequirements
======================

This example depends on the following data and software:


Installation of MPRAflow
----------------------------------------

Please install conda, the MPRAflow environment and clone the actual MPRAflow master branch. You will find more help under :ref:`Installation`.


Assignment file
----------------------------------------

This file is a tab separated files that assigns variants to barcodes. We will create a new working folder and download the file into it

.. code-block:: bash

    mkdir -p SatMut_TERT/data
    cd SatMut_TERT/data
    wget http
    cd ..


Count tables
----------------

We need the count tables of the count workflow. Please go to the :ref:`Count for Saturation Mutagenesis of the TERT promoter` and run it first. Afterwards copy the count tables into the data folder or use symbolic links:

.. code-block:: bash

    ln -s ../Count_TERT/output/TERT-GBM/*/TERT-GBM_{1,2,3}_counts.tsv.gz data/
    ln -s ../Count_TERT/output/TERT-HEK/*/TERT-HEK_{1,2,3}_counts.tsv.gz data/


Now the data folder should have the following files:

.. code-block:: bash

    tree data

.. code-block:: text

    data
    ├── TERT-GBM_1_counts.tsv.gz
    ├── TERT-GBM_2_counts.tsv.gz
    ├── TERT-GBM_3_counts.tsv.gz
    ├── TERT-HEK_1_counts.tsv.gz
    ├── TERT-HEK_2_counts.tsv.gz
    ├── TERT-HEK_3_counts.tsv.gz
    └── TERT.variants.txt.gz

    0 directories, 7 files


MPRAflow
=================================

Now we are close to start MPRAflow and find out individual variant effects. But before we need to generate an :code:`environment.csv` file to tell nextflow the conditions, replicates and the count files.

Create environment.csv
---------------------------

Our experiment file looks exactly like this:

.. code-block:: text

    Condition,Replicate,COUNTS
    TERT-GBM,1,TERT-GBM_1_counts.tsv.gz
    TERT-GBM,2,TERT-GBM_2_counts.tsv.gz
    TERT-GBM,3,TERT-GBM_3_counts.tsv.gz
    TERT-HEK,1,TERT-HEK_1_counts.tsv.gz
    TERT-HEK,2,TERT-HEK_2_counts.tsv.gz
    TERT-HEK,3,TERT-HEK_3_counts.tsv.gz

Save it into the :code:`SatMut_TERT/data` folder under :code:`experiment.csv`.


Run nextflow
------------------------------
Now we have everything at hand to run the saturation mutagenesis MPRAflow pripeline. Therefore we have to be in the cloned MPRAflow folder. But we will change the working and output directory to the :code:`SatMut_TERT` folder. The MPRAflow saturation mutagenesis command is:


.. code-block:: bash

    cd <path/to/MPRAflow>/MPRAflow
    conda activate MPRAflow
    nextflow run -resume -w <path/to/TERT>/SatMut_TERT/work  saturationMutagenesis.nf --experiment-file "<path/to/TERT>/SatMut_TERT/data/experiment.csv" --assignment "<path/to/TERT>/SatMut_TERT/data/TERT.variants.txt.gz" --dir "<path/to/TERT>/SatMut_TERT/data" --outdir "<path/to/TERT>/SatMut_TERT/output"

.. note:: Please check your :code:`nextflow.config` file if it is correctly configured (e.g. with your SGE cluster commands).

If everything works fine the following 5 processes will run: :code:`create_BAM (make idx)` :code:`raw_counts`, :code:`filter_counts`, :code:`final_counts`, :code:`dna_rna_merge_counts`.

..code-block:: text

    [49/53495c] process > create_BAM (make idx)    [100%] 12 of 12 ✔
    [92/f2a68d] process > raw_counts (12)          [100%] 12 of 12 ✔
    [af/398836] process > filter_counts (12)       [100%] 12 of 12 ✔
    [63/fb29b6] process > final_counts (12)        [100%] 12 of 12 ✔
    [75/f412e8] process > dna_rna_merge_counts (5) [100%] 6 of 6 ✔
    Completed at: 03-Jan-2020 19:55:10
    Duration    : 6h 16m 17s
    CPU hours   : 34.6
    Succeeded   : 54


Results
-----------------




.. todo::
  Add sat mut example
