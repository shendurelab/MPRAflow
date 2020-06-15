.. _faq:

==========================
Frequently Asked Questions
==========================

If you have more question please write us a ticket on `github <https://github.com/shendurelab/MPRAflow/issues>`_.

MPRAflow is not able to create a Conda environment
    If you get a message like::

        Caused by: java.lang.IllegalStateException: Failed to create Conda environment
        command: conda env create --prefix /home/user/MPRAflow/work/conda/mpraflow_py27-a6601743cee3b1029d4f3c810b7ebf02 --file /home/user/MPRAflow/conf/mpraflow_py27.yml`

    Try to run conda separately using::

        conda env create --prefix /home/user/MPRAflow/work/conda/mpraflow_py27-a6601743cee3b1029d4f3c810b7ebf02 --file /home/user/MPRAflow/conf/mpraflow_py27.yml

    Afterwards try MPRAflow again. Please be sure that you are connected to the internet!



Can I use STARR-seq with MPRAflow?
    No. For more details have a look at this `comment <https://github.com/shendurelab/MPRAflow/issues/27#issuecomment-636515565>`_.
