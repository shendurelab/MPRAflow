# MPRAflow Changelog

## development

## v2.2

No workflow changes. Only a few fixes and some restructuring of configs. Using nextflow version 20.01 now!

### global changes

* nextflow version 20.01 is needed because of multiMap() function
* introducing new config file `conf/global.config` with global variables like the min. required nextflow version and the actual MPRAflow version.
* moving cluster config to a sepArate file: `conf/cluster.config`. Try to adapt the times to the sort and longtime labels. Modify SLURM queue to SLUM not SGE options.
* improved documentation

### saturationMutagenesis.nf

* Bugfix of default out dir. It was not set to `params.outdir = "outs"` so it tries to create a folder `null`. Now in `params.outdir`
* removing default `params.version` and `params.nf_required_version`. Now in `conf/global.config`
* Catching cases when barcode/p-value filtering produces 0 variants
* Change update depricated fork method. Now works with nextflow 20.01

### count.nf

* removing default `params.version`, `params.nf_required_version` and `params.outdir`. Now in `conf/global.config`.

### association.nf

* removing default `params.version`, `params.nf_required_version` and `params.outdir`. Now in `conf/global.config`.


## v2.1

Initial MPRAflow version for publication.
