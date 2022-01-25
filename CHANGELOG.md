# MPRAflow Changelog

## v2.3.3

### association.nf

* Bugfix in `src/nf_ori_map_barcodes.py` script (see issue #41). Now >= mapq instead of > mapq is used.

## v2.3.2

### count.nf

* Bugfix merge_all.py script (see issue #55)

## v2.3.1

### association.nf

* Bugfix empty design file (see issue #45)

## v2.3

### global changes

* Correcting typos in documentation
* adding new process label `highmem` to `conf/cluster.config`

### association_saturationMutagenesis.nf

New association saturation mutagenesis workflow. This workflow is about assocation variant calls with barcodes. Variants are introduced by an error-prone PCR. The workflow takes the sequencing of the region, with barcodes in index read and the reference sequence and maps the reads to the reference, calls variants and associates them with the corresponding barcode. it is a pre-step of `saturationMutagenesis.nf`.

### count.nf

* using BC threshold input for `plot_perInsertCounts_correlation.R` instead of hard-coded. Modify process `calc_correlations` to use the new input theshold.

### association.nf

* Remove "windows" characters in design fasta


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
