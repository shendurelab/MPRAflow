.. _Saturation mutagenesis of the TERT promoter:

==================================================
Saturation mutagenesis of the TERT promoter
==================================================

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
    wget http https://github.com/shendurelab/MPRAflow/raw/master/examples/saturationMutagenesis/TERT.variants.txt.gz
    cd ..

    It is also possible to get using the workflow :ref:`Association for Saturation mutagenesis of TERT example`.

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
Now we have everything at hand to run the saturation mutagenesis MPRAflow pipeline. Therefore we have to be in the cloned MPRAflow folder. But we will change the working and output directory to the :code:`SatMut_TERT` folder. The MPRAflow saturation mutagenesis command is:


.. code-block:: bash

    cd <path/to/MPRAflow>/MPRAflow
    conda activate MPRAflow
    nextflow run -resume -w <path/to/TERT>/SatMut_TERT/work  saturationMutagenesis.nf --experiment-file "<path/to/TERT>/SatMut_TERT/data/experiment.csv" --assignment "<path/to/TERT>/SatMut_TERT/data/TERT.variants.txt.gz" --dir "<path/to/TERT>/SatMut_TERT/data" --outdir "<path/to/TERT>/SatMut_TERT/output"

.. note:: Please check your :code:`conf/cluster.config` file if it is correctly configured (e.g. with your SGE cluster commands).

If everything works fine the following 11 processes will run: :code:`calc_assign_variantMatrix` :code:`calc_assign_variantMatrixWith1bpDel`, :code:`fitModel`, :code:`summarizeVariantMatrix`, :code:`statsWithCoefficient`, :code:`plotCorrelation`, :code:`plotStatsWithCoefficient`, :code:`fitModelCombined`, :code:`combinedStats`, :code:`statsWithCoefficientCombined`, and :code:`plotStatsWithCoefficientCombined`.

.. code-block:: text

    [3c/835d00] process > calc_assign_variantMatrix (1)           [100%] 6 of 6 ✔
    [7a/887135] process > calc_assign_variantMatrixWith1bpDel (1) [100%] 6 of 6 ✔
    [ca/a90b00] process > fitModel (8)                            [100%] 12 of 12 ✔
    [67/3a3e8a] process > summarizeVariantMatrix (12)             [100%] 12 of 12 ✔
    [56/846670] process > statsWithCoefficient (12)               [100%] 12 of 12 ✔
    [74/466bfb] process > plotCorrelation (1)                     [100%] 12 of 12 ✔
    [a5/baf1ef] process > plotStatsWithCoefficient (12)           [100%] 12 of 12 ✔
    [ac/d38378] process > fitModelCombined (3)                    [100%] 4 of 4 ✔
    [0b/600d8b] process > combinedStats (2)                       [100%] 4 of 4 ✔
    [32/80f6a6] process > statsWithCoefficientCombined (2)        [100%] 4 of 4 ✔
    [2f/817e76] process > plotStatsWithCoefficientCombined (1)    [100%] 4 of 4 ✔
    Completed at: 07-Jan-2020 11:31:00
    Duration    : 22m 41s
    CPU hours   : 1.0
    Succeeded   : 88


Results
-----------------

All needed output files will be in the :code:`SatMut_TERT/output` folder.