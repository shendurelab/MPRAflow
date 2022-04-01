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
Max Schubach <max.schubach@bihealth.de>
Sean Whalen <sean.whalen@gladstone.ucsf.edu>
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
      --fastq-insert                Full path to library association fastq for insert (must be surrounded with quotes)
      --fastq-bc                    Full path to library association fastq for bc (must be surrounded with quotes)
      --design                      Full path to fasta of ordered oligo sequences (must be surrounded with quotes)
      --name                        Name of the association. Files will be named after this.

    Options:
      --fastq-insertPE              Full path to library association fastq for read2 if the library is paired end (must be surrounded with quotes)
      --min-cov                     minimum coverage of bc to count it (default 3)
      --min-frac                    minimum fraction of bc map to single insert (default 0.5)
      --mapq                        map quality (default 30)
      --baseq                       base quality (default 30)
      --cigar                       require exact match ex: 200M (default none)
      --outdir                      The output directory where the results will be saved and what will be used as a prefix (default outs)
      --split                       Number read entries per fastq chunk for faster processing (default: 2000000)
      --labels                      tsv with the oligo pool fasta and a group label (ex: positive_control) if no labels desired a file will be automatically generated
      --variants                    tsv with reference_name variant_positions ref_bases alt_bases, only input for variant analyses workflow

    Extras:
      --h, --help                   Print this help message
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
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
params.email = false
params.plaintext_email = false
output_docs = file("$baseDir/docs/output.md")

//defaults
params.min_cov="3"
params.min_frac="0.5"
params.baseq="30"
params.mapq="30"
params.cigar="n"
params.split=2000000


// Validate inputs
if ( !params.containsKey("name") ){
    exit 1, "Pleas especify a name of this workflow using --name"
}

if ( params.containsKey("fastq-insert") ){
    params.fastq_insert_file = file(params['fastq-insert'])
    if( !params.fastq_insert_file.exists() ) exit 1, "Fastq insert file not found: ${params.fastq_insert_file}"
} else {
    exit 1, "Fastq insert file not specified with --fastq-insert"
}

if(params.containsKey("fastq-insertPE")){
    params.fastq_insertPE_file = file(params['fastq-insertPE'])
    if( !params.fastq_insertPE_file.exists() ) exit 1, "Fastq paired-end insert file not found: ${params.fastq_insertPE_file}"
} else {
  params.fastq_insertPE_file = null
}

// Fastq barcode file in params.fastq_bc_file
if ( params.containsKey("fastq-bc")){
    params.fastq_bc_file = file(params['fastq-bc'])
    if( !params.fastq_bc_file.exists() ) exit 1, "Fastq barcode file not found: ${params.fastq_bc_file}"
} else {
    exit 1, "Fastq barcode file not specified with --fastq-bc"
}

// design file saved in params.design_file
if ( params.containsKey("design")){
    params.design_file=file(params.design)
    if( !params.design_file.exists() ) exit 1, "Design file ${params.design} does not exist"
} else {
    exit 1, "Design file not specified with --design"
}

// label file saved in params.label_file
if (params.containsKey("labels")){
    params.label_file=file(params.labels)
    if (!params.label_file.exists()) exit 1, "Label file ${params.label_file} does not exist"
} else {
    params.label_file=null
}

// variants file saved in params.variants_file
if (params.containsKey("variants")){
    params.variants_file=file(params.variants)
    if (!params.variants_file.exists()) exit 1, "variants file ${params.variants_file} does not exist"
} else {
    params.variants_file=null
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
summary['Pipeline Name']    = 'MPRAflow'
summary['Pipeline Version'] = params.version
summary['Fastq insert']     = params.fastq_insert_file
summary['fastq paired']     = params.fastq_insertPE_file
summary['Fastq barcode']    = params.fastq_bc_file
summary['design fasta']     = params.design_file
summary['minimum BC cov']   = params.min_cov
summary['map quality']      = params.mapq
summary['base quality']     = params.baseq
summary['cigar string']     = params.cigar
summary['min % mapped']     = params.min_frac
summary['Output dir']       = params.outdir
summary['Run name'] = params.name
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
if(params.variants_file!=null) summary['Variants file'] = params.variants_file
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['base directory']   = "$baseDir"
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile

//summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
//summary['Max CPUs']         = params.max_cpus
//summary['Max Time']         = params.max_time


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
* count fastq and bam length remove the illegal regex characters
* and make label file if missing
* contributions: Gracie Gordon & Max Schubach
*/

if (params.label_file != null) {
    process 'count_bc' {
        tag 'count'
        label 'shorttime'
        publishDir "${params.outdir}/${params.name}", mode:'copy'

        input:
            file(fastq_bc) from params.fastq_bc_file
            file(design) from params.design_file
            file(labels) from params.label_file
        output:
            file 'count_fastq.txt' into bc_ch
            file "label_rmIllegalChars.txt" into fixed_label
            file "design_rmIllegalChars.fa" into fixed_design
        shell:
            """
            #!/bin/bash
            ## Get rid of illegal regex characters
            awk '{gsub(/\\[/,"_")}1' $labels > t_new_label.txt
            awk '{gsub(/\\]/,"_")}1' t_new_label.txt > label_rmIllegalChars.txt


            awk '{gsub(/\\[/,"_")}1' $design > t_new_design.txt
            awk '{gsub(/\\]/,"_")}1' t_new_design.txt > design_rmIllegalChars.fa

            zcat $fastq_bc | wc -l  > count_fastq.txt
            """
    }
}

/*
* count fastq and bam length remove the illegal regex characters
* and make design file
* contributions: Gracie Gordon & Max Schubach
*/

if (params.label_file == null) {
    process 'count_bc_nolab' {
        tag 'count'
        label 'shorttime'
        publishDir "${params.outdir}/${params.name}", mode:'copy'

        input:
            file(fastq_bc) from params.fastq_bc_file
            file(design) from params.design_file
        output:
            file 'count_fastq.txt' into bc_ch
            file "label_rmIllegalChars.txt" into fixed_label
            file "design_rmIllegalChars.fa" into fixed_design
        shell:
            """
            #!/bin/bash
            #CREATE LABEL FILE and remove illegal regex characters
            awk -F'\t' 'BEGIN {OFS = FS} NR%2==1 {print substr(\$1,2,length(\$1)),"na"}' $design > labels.txt
            awk '{gsub(/\\[/,"_")}1' labels.txt > t_new_label.txt
            awk '{gsub(/\\]/,"_")}1' t_new_label.txt > label_rmIllegalChars.txt


            awk '{gsub(/\\[/,"_")}1' $design | \\
            awk '{gsub(/\\]/,"_")}1' | \\
            sed 's/\\r//g' > design_rmIllegalChars.fa

            zcat $fastq_bc | wc -l  > count_fastq.txt
            """
    }
}


/*
* STEP 1: Align
* Process 1A: create BWA reference
* contributions: Gracie Gordon
*/

process 'create_BWA_ref' {
    tag "make ref"
    label 'shorttime'

    conda 'conf/mpraflow_py36.yml'

    input:
        file(design) from fixed_design
        file(label) from fixed_label
    output:
        file "${design}.fai" into reference_fai
        file "${design}.bwt" into reference_bwt
        file "${design}.sa" into reference_sa
        file "${design}.pac" into reference_pac
        file "${design}.ann" into reference_ann
        file "${design}.amb" into reference_amb
        file "${design}.dict" into reference_dict
    shell:
        """
        #!/bin/bash
        bwa index -a bwtsw $design
        samtools faidx $design
        picard CreateSequenceDictionary REFERENCE=$design OUTPUT=$design".dict"
        """
}

/*
*CHUNKING FASTQ
*/

Channel
    .fromPath(params.fastq_insert_file)
    .splitFastq( by: params.split, file: true )
    .set{ R1_ch }

if (params.fastq_insertPE_file != null) {
    Channel
        .fromPath(params.fastq_insertPE_file)
        .splitFastq( by: params.split, file: true )
        .set{ R3_ch }
}



/*
*Process 1B: merge Paired end reads
* contributions: Gracie Gordon
*/
if (params.fastq_insertPE_file != null) {
    process 'PE_merge' {
        tag 'merge'
        label 'shorttime'

        conda 'conf/mpraflow_py36.yml'

        input:
            file(fastq_insert) from R1_ch
            file(fastq_insertPE) from R3_ch
        output:
            file "*merged.fastqjoin" into mergedPE
        shell:
            """
            fastq-join $fastq_insert $fastq_insertPE -o ${fastq_insert}_merged.fastq
            """
    }
}


/*
* Process 1C: align with BWA
* contributions: Gracie Gordon
*/

//paired ends
if (params.fastq_insertPE_file != null) {
    process 'align_BWA_PE' {
        tag "align"
        label 'longtime'

        conda 'conf/mpraflow_py36.yml'

        input:
            file(design) from fixed_design
            file(chunk) from mergedPE
            val(name) from params.name
            file(reference_fai) from reference_fai
            file reference_bwt from reference_bwt
            file reference_sa from reference_sa
            file reference_pac from reference_pac
            file reference_ann from reference_ann
            file reference_amb from reference_amb
            file reference_dict from reference_dict
        output:
            file "${name}.${chunk}.sorted.bam" into s_bam
            file '*count_bam.txt' into bam_ch
        shell:
            """
            bwa mem $design $chunk | \
            samtools sort - -o ${name}.${chunk}.sorted.bam

            echo 'bam made'

            samtools view ${name}.${chunk}.sorted.bam | head

            samtools view ${name}.${chunk}.sorted.bam | \
            wc -l > ${chunk}.count_bam.txt
            """
    }
} else {
    //single end
    process 'align_BWA_S' {
        tag "align"
        label 'longtime'

        conda 'conf/mpraflow_py36.yml'

        input:
            file(design) from fixed_design
            file(chunk) from R1_ch
            val(name) from params.name
            file(reference_fai) from reference_fai
            file reference_bwt from reference_bwt
            file reference_sa from reference_sa
            file reference_pac from reference_pac
            file reference_ann from reference_ann
            file reference_amb from reference_amb
            file reference_dict from reference_dict
        output:
            file "${name}.${chunk}.sorted.bam" into s_bam
            file '*count_bam.txt' into bam_ch
        shell:
            """
            bwa mem $design $chunk | samtools sort - -o ${name}.${chunk}.sorted.bam
            echo 'bam made'
            samtools view ${name}.${chunk}.sorted.bam | head
            samtools view ${name}.${chunk}.sorted.bam | wc -l > ${chunk}_count_bam.txt
            """
    }
}

/*
*COLLCT FASTQ CHUNCKS
*/

process 'collect_chunks'{
    label 'shorttime'

    conda 'conf/mpraflow_py36.yml'

    input:
        file sbam_listFiles from s_bam.collect()
        file count_bamFiles from bam_ch.collect()
    output:
        file 's_merged.bam' into s_merge
        file 'count_merged.txt' into ch_merge
    script:
        count_bam = count_bamFiles.join(' ')
        sbam_list = sbam_listFiles.join(' ')
    shell:
        """
        #collect sorted bams into one file
        samtools merge all.bam $sbam_list
        samtools sort all.bam -o s_merged.bam

        #collect bam counts into one file

        samtools view s_merged.bam | wc -l > count_merged.txt
        """
}




/*
* Assign barcodes to element sequences
* contributions: Sean Whalen
*/

if (params.variants_file == null){
    process 'map_element_barcodes' {
        tag "assign"
        label "shorttime"
        publishDir "${params.outdir}/${params.name}", mode:'copy'
    
        conda 'conf/mpraflow_py36.yml'
    
        input:
            val(name) from params.name
            val(mapq) from params.mapq
            val(baseq) from params.baseq
            val(cigar) from params.cigar
            file(fastq_bc) from params.fastq_bc_file
            file count_fastq from bc_ch
            file count_bam from ch_merge
            file bam from s_merge
        output:
            file "${name}_coords_to_barcodes.pickle" into map_ch
            file "${name}_barcodes_per_candidate-no_repeats-no_jackpots.feather" into count_table_ch
            file "${name}_barcode_counts.pickle"
        shell:
            """
            echo "test assign inputs"
            echo ${mapq}
            echo ${baseq}
            echo $fastq_bc
            zcat $fastq_bc | head
    
            echo ${count_fastq}
            echo ${count_bam}
            cat ${count_fastq}
            cat ${count_bam}
    
            python ${"$baseDir"}/src/nf_ori_map_barcodes.py ${"$baseDir"} ${fastq_bc} ${count_fastq} \
            $bam ${count_bam} ${name} ${mapq} ${baseq} ${cigar}
            """
    }
} else {
    process 'map_element_barcodes_mut' {
        tag "assign"
        label "shorttime"
        publishDir "${params.outdir}/${params.name}", mode:'copy'
    
        conda 'conf/mpraflow_py36.yml'
    
        input:
            val(name) from params.name
            val(mapq) from params.mapq
            val(baseq) from params.baseq
            val(cigar) from params.cigar
            file(fastq_bc) from params.fastq_bc_file
            file count_fastq from bc_ch
            file count_bam from ch_merge
            file bam from s_merge
            file(variants) from params.variants_file
        output:
            file "${name}_coords_to_barcodes.pickle" into map_ch
            file "${name}_barcodes_per_candidate-no_repeats-no_jackpots.feather" into count_table_ch
            file "${name}_barcode_counts.pickle"
            
        shell:
            """
            echo "test assign inputs"
            echo ${mapq}
            echo ${baseq}
            echo $fastq_bc
            zcat $fastq_bc | head
    
            echo ${count_fastq}
            echo ${count_bam}
            cat ${count_fastq}
            cat ${count_bam}
    
            python ${"$baseDir"}/src/nf_ori_map_barcodes.py ${"$baseDir"} ${fastq_bc} ${count_fastq} \
            $bam ${count_bam} ${name} ${mapq} ${baseq} ${cigar} --mutations ${variants}
            """
    }
}


/*
* Filter barcodes for minimum coverage and unique mapping
* contributions: Gracie Gordon
*/
if (params.variants_file == null) {
    process 'filter_barcodes' {
        tag "$filter"
        label "shorttime"
        publishDir "${params.outdir}/${params.name}", mode:'copy'
    
        conda 'conf/mpraflow_py36.yml'
    
        input:
            val(min_cov) from params.min_cov
            val(min_frac) from params.min_frac
            val(out) from params.name
            file(map) from map_ch
            file(table) from count_table_ch
            file(label) from fixed_label
        output:
            file "${out}_filtered_coords_to_barcodes.pickle"
            file "${out}_original_counts.png"
            file "original_count_summary.txt"
            file "${out}_filtered_counts.png"
            file "filtered_count_summary.txt"
    
        shell:
            """
            python ${"$baseDir"}/src/nf_filter_barcodes.py ${out} ${map} ${table} \
            ${min_cov} ${min_frac} $label
            """
    }
} else {
    process 'filter_barcodes_mut' {
        tag "$filter"
        label "shorttime"
        publishDir "${params.outdir}/${params.name}", mode:'copy'
    
        conda 'conf/mpraflow_py36.yml'
    
        input:
            val(min_cov) from params.min_cov
            val(min_frac) from params.min_frac
            val(out) from params.name
            file(map) from map_ch
            file(table) from count_table_ch
            file(label) from fixed_label
            file(design) from fixed_design
            file(variants) from params.variants_file
        output:
            file "${out}_filtered_coords_to_barcodes.pickle"
            file "${out}_original_counts.png"
            file "original_count_summary.txt"
            file "${out}_filtered_counts.png"
            file "filtered_count_summary.txt"
            file "label_mutExpand.txt"
            file "design_mutExpand.fa"
    
        shell:
            """
            python ${"$baseDir"}/src/expand_variants.py ${design} ${label} ${variants} design_mutExpand.fa label_mutExpand.txt
            python ${"$baseDir"}/src/nf_filter_barcodes.py ${out} ${map} ${table} \
            ${min_cov} ${min_frac} label_mutExpand.txt
            """
    }
}



/*
 * Completion e-mail notification
 */
/*
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[MPRAflow] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[MPRAflow] FAILED: $workflow.runName"
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
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", attach1: "$baseDir/results/Documentation/pipeline_report.html", attach2: "$baseDir/results/pipeline_info//MPRAflow_report.html", attach3: "$baseDir/results/pipeline_info//MPRAflow_timeline.html" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[/MPRAflow] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[/MPRAflow] Sent summary e-mail to $params.email (mail)"
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

    log.info "[/MPRAflow] Pipeline Complete"

}

*/
