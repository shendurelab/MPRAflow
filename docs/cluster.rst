===============================
Running MPRAflow on HPC Cluster
===============================

Nextflow gives us the opportunity to run MPRAflow in a cluster environment. Right now we split up processes into two main groups: `longtime` and `shorttime`. We can define different job setting for both groups. As you can imagine from the names `longtime` defines processes that takes a while when running. Sometimes several days. `shortime` defines processes that are quicker and are usually done in several minutes.

To enable the submission to your cluster you have to edit the `conf/cluster.config`, allowing each process to be run as a separate `qsub`, `sbatch` or similar command. The config contains example code for SGE, LSF, and SLURM architectures. Please remove the `\\` for the architecture you would like to use and place `\\` in front of any architectures not currently in use. A `\\` in front of all of them runs the pipeline on your local machine (default). If you run MPRAflow on a cluster system make sure be that you export all environment variables. E.g. this can be done with the `-V` option by SGE.

  .. note:: Please consult your cluster's wiki page for cluster specific commands and change clusterOptions = to reflect these specifications. Additionally, for large libraries, more memory can be specified in this location.
