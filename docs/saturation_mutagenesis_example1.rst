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

We need the count tables of the count workflow. Please go to the :ref:`Count for Saturation Mutagenesis of the TERT promoter` and run it first. Afterwards copy the count tables into the data folder:

.. code-block:: bash

    cp <path/to/xyz/> data/


Now the data folder should have the following files:

.. code-block:: bash

    tree data

.. code-block:: text

    _build/
    ├── doctrees
    │   ├── association.doctree
    │   ├── association_example1.doctree
    │   ├── authors.doctree


MPRAflow
=================================

Create environment.csv
---------------------------


Run nextflow
------------------------------


Results
-----------------



.. todo::
  Add sat mut example
