#!/usr/bin/env nextflow

params.version=2.1
/*
========================================================================================
                         MPRAflow
========================================================================================
MPRA Analysis Pipeline. Started 2019-07-29.
Saturation Mutagenesis Utility
#### Homepage / Documentation
https://github.com/shendurelab/MPRAflow
#### Authors
Max Schubach <max.schubach@bihealth.de>
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
    shendurelab/MPRAflow v${params.version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run saturationMutagenesis.nf
    Mandatory arguments:
      --dir                         Directory of count files (must be surrounded with quotes)
      --assignment                  Variant assignment file
      --e, --experiment-file        Experiment csv file

    Options:
      --outdir                      The output directory where the results will be saved (default outs)
      --thresh                      Minimum number of observed barcodes to retain variant (default 10)
      --pvalue                      pValue cutoff for significant different cvariant effects. For variant effect plots only. (default 1e-5)
    """.stripIndent()
}

/*
* SET UP CONFIGURATION VARIABLES
*/

// Show help message
if (params.containsKey('h') || params.containsKey('help')){
    helpMessage()
    exit 0
}


// Configurable variables
//defaults
results_path = params.outdir
params.nf_required_version="19.10"
params.thresh = 10
params.pvalue = 1e-5

// Validate Inputs

// experiment file saved in params.experiment_file
if (params.containsKey('e')){
    params.experiment_file=file(params.e)
} else if (params.containsKey("experiment-file")) {
    params.experiment_file=file(params["experiment-file"])
} else {
    exit 1, "Experiment file not specified with --e or --experiment-file"
}
if( !params.experiment_file.exists()) exit 1, "Experiment file ${params.experiment_file} does not exist"


// Association file in params.association_file
if ( params.containsKey("assignment")){
    params.assignment_file=file(params.assignment)
    if( !params.assignment_file.exists() ) exit 1, "Assignment file ${params.assignment_file} does not exist"
} else {
    exit 1, "Assignment file not specified with --assignment"
}

// Create FASTQ channels
Channel.fromPath(params.experiment_file).splitCsv(header: true).flatMap{
  row -> [[row.Condition, row.Replicate,file([params.dir,"/",row.COUNTS].join())]]
}.into{count_channel; count_1bpDel_channel}

// Header log info
log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'
MPRAflow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'shendurelab/MPRAflow'
summary['Pipeline Version'] = params.version

summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Working dir']      = workflow.workDir
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
summary['Experiment File']  = params.experiment_file
summary['BC threshold']     = params.thresh

log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

println 'start analysis'

/*
* STEP 1: Calculate assignment variant Matrix with and witout indels
*/

process 'calc_assign_variantMatrix' {
  publishDir "$params.outdir/$cond/$rep", mode:'copy'
  label 'shorttime'

  conda 'conf/mpraflow_py27.yml'

  input:
    file(assignments) from params.assignment_file
    tuple val(cond), val(rep), file(counts) from count_channel
  output:
    tuple val(cond), val(rep), val("VarMatrix"), file("${cond}_${rep}.VarMatrix.tsv.gz") into variantMatrix, variantMatrixCombined
  shell:
    """
    python ${"$baseDir"}/src/satMut/extractVariantMatrix.py -a ${assignments} ${counts} | \
    gzip -c > ${cond}_${rep}.VarMatrix.tsv.gz
    """
}

process 'calc_assign_variantMatrixWith1bpDel' {
  publishDir "$params.outdir/$cond/$rep", mode:'copy'
  label 'shorttime'

  conda 'conf/mpraflow_py27.yml'

  input:
    file(assignments) from params.assignment_file
    tuple val(cond), val(rep), file(counts) from count_1bpDel_channel
  output:
    tuple val(cond), val(rep), val("VarMatrix1bpDel"), file("${cond}_${rep}.VarMatrix1bpDel.tsv.gz") into variantMatrix1bpDel, variantMatrix1bpDelCombined
  shell:
    """
    python ${"$baseDir"}/src/satMut/extractVariantMatrixHandleIndels.py -a ${assignments} ${counts} | \
    gzip -c > ${cond}_${rep}.VarMatrix1bpDel.tsv.gz
    """
}

/*
* STEP 2: FIT model
*/

process 'fitModel' {
  publishDir "$params.outdir/$cond/$rep", mode:'copy'
  label 'longtime'

  conda 'conf/mpraflow_r.yml'

    variantMatrixBoth = variantMatrix.concat(variantMatrix1bpDel)

  input:
    tuple val(cond), val(rep), val(type), file(variantMatrix) from variantMatrixBoth
  output:
    tuple val(cond), val(rep), val(type), file(variantMatrix), file("${cond}_${rep}.${type}.ModelCoefficients.txt") into variantMatrixModelCoefficients
  shell:
    """
    Rscript ${"$baseDir"}/src/satMut/fitModel.R $variantMatrix ${cond}_${rep}.${type}.ModelCoefficients.txt
    """
}


process 'summarizeVariantMatrix' {
  publishDir "$params.outdir/$cond/$rep", mode:'copy'
  label 'shorttime'

  conda 'conf/mpraflow_py27.yml'

  input:
    tuple val(cond), val(rep), val(type), file(variantMatrix), file(modelCoefFile) from variantMatrixModelCoefficients
  output:
    tuple val(cond), val(rep), val(type), file(modelCoefFile), file("${variantMatrix}.stats") into variantMatrixModelCoefficientsStats, variantMatrixModelCoefficientsStatsForCombined
  shell:
    """
    python ${"$baseDir"}/src/satMut/summarizeVariantMatrix.py $variantMatrix
    """
}


process 'statsWithCoefficient' {
  publishDir "$params.outdir/$cond/$rep", mode:'copy'
  label 'shorttime'

  input:
    tuple val(cond), val(rep), val(type), file(modelCoefFile), file(statsFile) from variantMatrixModelCoefficientsStats
  output:
    tuple val(cond), val(rep), val(type), file("${cond}_${rep}.${type}.ModelCoefficientsVsStats.txt") into variantMatrixModelCoefficientsVsStats, variantMatrixModelCoefficientsVsStats2
  shell:
    """
    (
      echo -e "Position\\tBarcodes\\tDNA\\tRNA\\tCoefficient\\tStdError\\ttValue\\tpValue";
      join -t"\$(echo -e '\\t')" <(sort $statsFile) <( sed s/^X//g $modelCoefFile | sort );
    ) > ${cond}_${rep}.${type}.ModelCoefficientsVsStats.txt;
    """
}

List combinationsOf(List list, int r) {
    assert (0..<list.size()).contains(r) // validate input
    def combs = [] as Set
    list.eachPermutation {
        combs << it.subList(0, r).sort { a, b -> a <=> b }
    }
    combs as List
}

process 'plotCorrelation' {
  publishDir "$params.outdir/$cond", mode:'copy'
  label 'shorttime'

  conda 'conf/mpraflow_r.yml'

  result = variantMatrixModelCoefficientsVsStats2.groupTuple(by: [0,2]).map{n -> [n[0],n[2],combinationsOf(n[1],2),combinationsOf(n[3],2)]}.transpose(by:[2,3]).map{n -> [n[0],n[1],n[2][0],n[2][1],n[3][0],n[3][1]]}
  input:
    tuple val(cond), val(type), val(rep1), val(rep2), file(results1), file(results2) from result
    val(thresh) from params.thresh
  output:
    tuple val(cond), val(type), file("${cond}.${type}.${rep1}_${rep2}.correlation.png") into correlations
  shell:
    """
    Rscript ${"$baseDir"}/src/satMut/plotCorrelation.R $results1 $results2 ${cond}.${type}.Replicate_${rep1} ${cond}.${type}.Replicate_${rep2} $thresh "${cond}.${type}.${rep1}_${rep2}.correlation.png"
    """
}

process 'plotStatsWithCoefficient' {
  publishDir "$params.outdir/$cond/$rep", mode:'copy'
  label 'shorttime'

  conda 'conf/mpraflow_r.yml'

  input:
    tuple val(cond), val(rep), val(type), file(results) from variantMatrixModelCoefficientsVsStats
    val(threshhold) from params.thresh
    val(pvalue) from params.pvalue
  output:
    tuple val(cond), val(rep), val(type), file("${cond}.${type}.saturationMutagenesis.png") into satMutPlot
  shell:
    """
    Rscript ${"$baseDir"}/src/satMut/plotElements.R $results ${cond}_${rep}.${type} $threshhold $pvalue ${cond}.${type}.saturationMutagenesis.png

    """
}


process 'fitModelCombined' {
  publishDir "$params.outdir/$cond", mode:'copy'
  label 'longtime'

  conda 'conf/mpraflow_r.yml'

  result = variantMatrixCombined.concat(variantMatrix1bpDelCombined).groupTuple(by: [0,2]).fork{i ->
                            cond: i[0]
                            type: i[2]
                            replicate: i[1].join(" ")
                            files: i[3]
                          }

  input:
    val(cond) from result.cond
    val(type) from result.type
    val(replicates) from result.replicate
    file(variantMatrix) from result.files
  output:
    tuple val(cond), val(type), file("${cond}.Combined.${type}.ModelCoefficients.txt") into variantMatrixModelCoefficientsCombined
  shell:
    """
    Rscript ${"$baseDir"}/src/satMut/fitModelCombined.R ${cond}.Combined.${type}.ModelCoefficients.txt $variantMatrix $replicates
    """
}


process 'combinedStats' {
  publishDir "$params.outdir/$cond", mode:'copy'
  label 'shorttime'

  conda 'conf/mpraflow_r.yml'

  result = variantMatrixModelCoefficientsStatsForCombined.groupTuple(by: [0,2]).fork{i ->
                            cond: i[0]
                            type: i[2]
                            replicate: i[1].join(" ")
                            files: i[4]
                          }

  input:
    val(cond) from result.cond
    val(type) from result.type
    val(replicates) from result.replicate
    file(stats) from result.files
  output:
    tuple val(cond), val(type), file("${cond}.Combined.${type}.stats") into variantMatrixModelCoefficientsStatsCombined
  shell:
    """
    Rscript ${"$baseDir"}/src/satMut/combineStats.R ${cond}.Combined.${type}.stats $stats

    """
}

process 'statsWithCoefficientCombined' {
  publishDir "$params.outdir/$cond", mode:'copy'
  label 'shorttime'

  result = variantMatrixModelCoefficientsStatsCombined.concat(variantMatrixModelCoefficientsCombined).groupTuple(by: [0,1]).fork{i ->
                            cond: i[0]
                            type: i[1]
                            files: i[2]
                          }

  input:
    val(cond) from result.cond
    val(type) from result.type
    tuple file(stat), file(model) from result.files
  output:
    tuple val(cond), val(type), file("${cond}.Combined.${type}.ModelCoefficientsVsStats.txt") into variantMatrixModelCoefficientsVsStatsCombined
  shell:
    """
    (
      echo -e "Position\\tBarcodes\\tDNA\\tRNA\\tCoefficient\\tStdError\\ttValue\\tpValue";
      join -t"\$(echo -e '\\t')" <(sort $stat | sort) <( sed s/^X//g $model | sort );
    ) > ${cond}.Combined.${type}.ModelCoefficientsVsStats.txt;
    """
}


process 'plotStatsWithCoefficientCombined' {
  publishDir "$params.outdir/$cond", mode:'copy'
  label 'shorttime'

  conda 'conf/mpraflow_r.yml'

  input:
    tuple val(cond), val(type), file(results) from variantMatrixModelCoefficientsVsStatsCombined
    val(thresh) from params.thresh
    val(pvalue) from params.pvalue
  output:
    tuple val(cond), val(type), file("${cond}.Combined.${type}.saturationMutagenesis.png") into combinedSatMutPlot
  shell:
    """
    Rscript ${"$baseDir"}/src/satMut/plotElements.R $results ${cond}.Combined.${type} $thresh $pvalue ${cond}.Combined.${type}.saturationMutagenesis.png
    """
}
