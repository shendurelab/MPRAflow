#!/usr/bin/env nextflow

/*
========================================================================================
                         MPRAflow
========================================================================================
MPRA Analysis Pipeline. Started 2019-07-29.
Library Association package 

#### Homepage / Documentation
https://github.com/shendurelab/MPRAflow
#### Authors
Gracie Gordon <gracie.gordon@ucsf.edu>
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     shendurelab/MPRAflow v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run MPRA-nextflow -profile singularity,test

    Mandatory arguments:
      --fastq_insert                Path to library association fastqs for insert (must be surrounded with quotes)
      --fastq_bc                    Path to library association fastq for bc (must be surrounded with quotes)
      --design                      fasta of ordered oligo sequences
      --out                         prefix for outputs (sample name)
      --condaloc                    location of conda instilation activate file ex: ~/miniconda3/bin/activate (must be surrounded with quotes)
      --labels                      tsv with the oligo pool fasta and a group label (ex: positive_control) if no labels desired simply add NA instead of a label in this file

    Options:
      --aligner                     alignment bwa-mem (1) bowtie2 (2) (currently only supports bwa-mem)
      --min_cov                     minimum coverage of bc to count it (default 2)
      --min_frac                    minimum fraction of bc map to single insert (default 0.5)
      --mapq                        map quality (default 30)
      --baseq                       base quality (default 30)
      --cigar                       require exact match ex: 200M (default none) 
      --outdir                      The output directory where the results will be saved (default outs)
      
    Extras:
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/py36.yaml"
params.email = false
params.plaintext_email = false
multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

//defaults
params.aligner="1"
params.min_cov="2"
params.min_frac="0.5"
params.baseq="30"
params.mapq="30"
params.cigar="n"
params.outdir="outs"
params.nf_required_version="19.07.0"
params.out="output"
//params.condaloc='/netapp/home/ggordon/tools/miniconda3/bin/activate'

// Validate inputs
if ( params.fastq_insert ){
    fastq_insert = file(params.fastq_insert)
    if( !fastq_insert.exists() ) exit 1, "Fastq insert file not found: ${params.fastq_insert}"
}

if ( params.fastq_bc ){
    fastq_bc = file(params.fastq_bc)
    if( !fastq_bc.exists() ) exit 1, "Fastq barcode file not found: ${params.fastq_bc}"
}

if ( params.design ){
    design = file(params.design)
    if( !design.exists() ) exit 1, "Fasta oligo design file not found: ${params.design}"
}

/*
if ( params.out ){
    out= params.out
    if( !out.exists() ) exit 1, "prefix not specified: ${params.out}"
}
*/

if (params.condaloc){
    condaloc=file(params.condaloc)
    if (!condaloc.exists()) exit 1, "conda location not specified ${params.condaloc}"
}

if (params.labels){
    labels=file(params.labels)
    if (!labels.exists()) exit 1, "label file not specified ${params.labels}"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}



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
summary['Pipeline Name']    = 'AhituvLab/MPRA_nextflow'
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Fastq insert']     = params.fastq_insert
summary['Fastq barcode']    = params.fastq_bc
summary['design fasta']     = params.design
summary['alignment']        = params.aligner
summary['minimum BC cov']   = params.min_cov
summary['map quality']      = params.mapq
summary['base quality']     = params.baseq
summary['cigar string']     = params.cigar
summary['output id']        = params.out

//summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Working dir']      = workflow.workDir
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
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


/*
 * Parse software version numbers
 */
//conda??
/*
process get_software_versions {

    input:
    file(trimmomatic_jar_path)

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    module load sra/2.8.0
    module load igvtools/2.3.75
    module load fastqc/0.11.5
    module load bedtools/2.25.0
    module load bowtie/2.2.9
    module load samtools/1.8

    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    java -XX:ParallelGCThreads=1 -jar ${trimmomatic_jar_path} -version &> v_trimmomatic.txt
    multiqc --version > v_multiqc.txt
    samtools --version &> v_samtools.txt
    bowtie2 --version &> v_bowtie2.txt
    fastq-dump --version &> v_fastq-dump.txt
    bedtools --version &> v_bedtools.txt
    igvtools &> v_igv-tools.txt
    macs2 --version &>  v_macs2.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}
*/

/*
* STEP 1: Align
* Process 1A: create BWA reference
*/

process 'create_BWA_ref' {
    tag "make ref"
    label 'shorttime'
    input:
    file(design) from design 
    file(params.condaloc) 

    output:
    file "${design}.fai" into reference_fai
    file "${design}.bwt" into reference_bwt
    file "${design}.sa" into reference_sa
    file "${design}.pac" into reference_pac
    file "${design}.ann" into reference_ann
    file "${design}.amb" into reference_amb
    file "${design}.dict" into reference_dict
    script:
    """
    #!/bin/bash
    source ${params.condaloc} mpraflow_py36

    bwa index -a bwtsw ${design}
    samtools faidx ${design}
    picard CreateSequenceDictionary REFERENCE=${design} OUTPUT=${design}".dict"

    """
}


/*
* 
* Process 1B: align with BWA
*/


process 'align_BWA' {
    tag "align"
    label 'longtime'
    publishDir params.outdir, mode:'copy'    

    input:
    file(design) from design
    file(fastq_insert) from fastq_insert
    file(params.out)
    file(reference_fai) from reference_fai
    file reference_bwt from reference_bwt
    file reference_sa from reference_sa
    file reference_pac from reference_pac
    file reference_ann from reference_ann
    file reference_amb from reference_amb
    file reference_dict from reference_dict

    output:
    file "${params.out}.bam" into bam
    file "${params.out}.sorted.bam" into s_bam
    file 'count_bam.txt' into bam_ch

    script:
    """
    #!/bin/bash
    source ${params.condaloc} mpraflow_py36
    
    bwa index $fastq_insert
    bwa mem $design $fastq_insert > ${params.out}.sam
    samtools view -S -b ${params.out}.sam > ${params.out}.bam
    echo 'bam made'
    samtools view ${params.out}.bam | head

    #sort bam
    samtools sort ${params.out}.bam -o ${params.out}.sorted.bam
    samtools view ${params.out}.sorted.bam | head

    samtools view ${params.out}".bam" | head
    samtools view ${params.out}".bam" | wc -l > count_bam.txt
    """
}



/*
* STEP 2 pre: count fastq and bam length
*/

process 'count_bc' {
    tag 'count'
    label 'shorttime'
    publishDir params.outdir, mode:'copy'

    input:
    file(fastq_bc) from fastq_bc

    output:
    file 'count_fastq.txt' into bc_ch

    """
    #!/bin/bash
    source ${params.condaloc} mpraflow_py36
    zcat $fastq_bc | wc -l 
    zcat $fastq_bc | wc -l  > count_fastq.txt

    """

}

/*
* STEP 2: Assign barcodes to enhancer sequences
*/

process 'map_enhancer_barcodes' {
    tag "assign"
    label "shorttime"
    publishDir params.outdir, mode:'copy'

    input:
    params.mapq
    params.baseq
    params.cigar
    file(fastq_bc) from fastq_bc 
    file count_fastq from bc_ch
    file count_bam from bam_ch     
    file bam from bam

    output:
    file "${params.out}_coords_to_barcodes.pickle" into map_ch
    file "${params.out}_barcodes_per_candidate.feather"
    file "${params.out}_barcodes_per_candidate-no_repeats-no_jackpots.feather" into count_table_ch
    file "${params.out}_barcodes_per_candidate-no_repeats.feather"
    file "${params.out}_barcodes_per_candidate-no_jackpots.feather"
    file "${params.out}_barcode_counts.pickle"

    script:

    file1=file(count_bam)
    println file1

    """
    #!/bin/bash 
    source ${params.condaloc} mpraflow_py36
    echo "test assign inputs"
    echo ${params.mapq}
    echo ${params.baseq}
    echo ${params.baseq}
    echo $fastq_bc
    zcat $fastq_bc | head
    
    echo ${count_fastq} 
    echo ${count_bam}    
    cat ${count_fastq}     
    cat ${count_bam} 
    
    python ${"$PWD"}/src/nf_ori_map_barcodes.py ${"$PWD"} ${fastq_bc} ${count_fastq} $bam ${count_bam} ${params.out} ${params.mapq} ${params.baseq} ${params.cigar}
    """
    

}

/*
* STEP 3: Filter barcodes for minimum coverage and unique mapping
*/

process 'filter_barcodes' {
    tag "$filter"
    label "shorttime"
    publishDir params.outdir, mode:'copy'
    input:
        params.min_cov
        params.out
        file(map) from map_ch
        file(table) from count_table_ch
    
    output:
         file "${params.out}_filtered_coords_to_barcodes.p"
         file "${params.out}_original_counts.png"
         file "original_count_summary.txt"
         file "${params.out}_filtered_counts.png"
         file "filtered_count_summary.txt" 
 
    script:
    """
    #!/bin/bash 
    source ${params.condaloc} mpraflow_py36
    python ${"$PWD"}/src/nf_filter_barcodes.py ${params.out} ${map} ${table} ${params.min_cov} ${params.min_frac} ${params.labels}
    """


}





/*
 * Completion e-mail notification
 */
/*
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[NCBI-Hackathons/ATACFlow] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[NCBI-Hackathons/ATACFlow] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", attach1: "$baseDir/results/Documentation/pipeline_report.html", attach2: "$baseDir/results/pipeline_info/NCBI-Hackathons/ATACFlow_report.html", attach3: "$baseDir/results/pipeline_info/NCBI-Hackathons/ATACFlow_timeline.html" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[NCBI-Hackathons/ATACFlow] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[NCBI-Hackathons/ATACFlow] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[NCBI-Hackathons/ATACFlow] Pipeline Complete"

}

*/
