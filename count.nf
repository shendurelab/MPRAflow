#!/usr/bin/env nextflow

/*
========================================================================================
                         MPRAflow
========================================================================================
MPRA Analysis Pipeline. Started 2019-07-29.
Count Utility

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

    nextflow run main.nf

    Mandatory arguments:
      --dir                         fasta directory (must be surrounded with quotes)
      --association                 pickle dictionary from library association process
      --design                      fasta of ordered insert sequences
      --s                           barcode length (default 26)          
      --l                           sample index length (default 10)
      --m                           UMI length (default 16)
      --e                           experiment csv file
      --out                         sample id to use as prefix (default output)
      --condaloc                    location of conda instilation activate file ex: ~/miniconda3/bin/activate (must be surrounded with quotes)
      --labels                      tsv with the oligo pool fasta and a group label (ex: positive_control) if no labels desired simply add NA instead of a label in this file

    Options:
      --outdir                      The output directory where the results will be saved (default outs)
      --merge_intersect             Only retain barcodes in RNA and DNA fraction (TRUE/FALSE, default: FALSE)
      --mpranalyze                  Only generate MPRAnalyze outputs (True:1, False:0 default 0)
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
//params.multiqc_config = "$baseDir/conf/py36.yaml"
params.email = false
params.plaintext_email = false
//multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

//defaults
params.outdir="outs"
results_path = params.outdir
params.nf_required_version="19.07.0"
params.out="output"
//params.condaloc='/netapp/home/ggordon/tools/miniconda3/bin/activate'
params.sample_idx="GATCCGGTTG"
params.s='26'
params.l='10'
params.m='16'
params.merge_intersect="FALSE"
params.mpranalyze=0

//params.dir="bulk_dna"
//params.e='/wynton/group/ye/ggordon/MPRA_nextflow/barcode_matching_scripts/toy_experiment.csv'
//params.design='/wynton/group/ye/ggordon/MPRA_nextflow/lib_assoc_scripts/pilot_library_noprimer.fa'
//params.association='/wynton/group/ye/ggordon/MPRA_nextflow/lib_assoc_scripts/outs/output_filtered_coords_to_barcodes.p'

// Validate Inputs

/*
if ( params.out ){
    out= params.out
    if( !out.exists() ) exit 1, "prefix not specified: ${params.out}"
}
*/

if (params.e){
    env=file(params.e)
    if( !env.exists() ) exit 1, "environment file not specified ${params.e}"
}

if ( params.design ){
    design=file(params.design)
    if( !design.exists() ) exit 1, "design file not specified ${params.design}"
}

if ( params.association ){
    assoc=file(params.association)
    if( !assoc.exists() ) exit 1, "association pickle not specified ${params.association}"
}

if ( params.condaloc ){
    cloc=file(params.condaloc)
    if( !cloc.exists() ) exit 1, "conda location not provided ${params.condaloc}"
}

if (params.labels){
    labels=file(params.labels)
    if (!labels.exists()) exit 1, "label file not specified ${labels}"
}

// Create FASTQ channels
r1=params.dir+'/*R1*.fastq.gz'
//r2=params.dir+'/*R2*.fastq.gz'
//r3=params.dir+'/*R3*.fastq.gz'
// define channelsi
R1_fastq = Channel
    .fromPath( r1 )
    .map { file -> tuple(file.simpleName, file) }

//if(params.m !=0){
//    println 'test'
//    R2_fastq = Channel
//        .fromPath( r2 )
//        .map { file -> tuple(file.simpleName, file) }
//}

//R3_fastq = Channel
//    .fromPath( r3 )
//    .map { file -> tuple(file.simpleName, file) }


//R1_fastq.subscribe { println "value: $it" }
//R3_fastq.subscribe { println "value: $it" }

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
summary['Pipeline Name']    = 'shendurelab/MPRAflow'
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
summary['output id']        = params.out

//summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
//summary['Container Engine'] = workflow.containerEngine
//if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Working dir']      = workflow.workDir
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
summary['Experiment File']  = params.e
summary['design file']      = params.design
summary['r1']               = r1
//summary['r2']               = r2
//summary['r3']               = r3
summary['m']                = params.m
summary['s']                = params.s

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


println 'start analysis'

/*
* STEP 1: Create BAM files 
*/

//if UMI
if(params.m !=0){
    process 'create_BAM' {
        tag "make idx"
 
        input:
        file(params.condaloc)
        val(params.sample_idx)
        set datasetID, file(r1_fastq) from R1_fastq
    
    
        output:
        set datasetID, file("${datasetID}_index.lst") into idx_list
        set datasetID, file("${datasetID}.bam") into clean_bam
        """
        #!/bin/bash
        source $params.condaloc mpraflow_py27
       
        echo "sample idx"
        echo $params.sample_idx
      
        echo $r1_fastq

        new_var2=\$(echo $r1_fastq | awk -F"_R1_" '{print \$1"_R2_"\$2}')
        echo \$new_var2
        new_2=$params.dir"/"\$new_var2
        echo \$new_2

        new_var3=\$(echo $r1_fastq | awk -F"_R1_" '{print \$1"_R3_"\$2}')
        echo \$new_var3
        new_3=$params.dir"/"\$new_var3
        echo \$new_3
     
        echo $params.sample_idx'        '${datasetID} >> ${datasetID}_index.lst 
       
        #paste <( zcat $r1_fastq ) <(zcat \$new_3 ) <(zcat \$new_2) | awk 'BEGIN{ counter=0 }{ counter+=1; if (counter == 2) { print \$1\$2 } else { if (counter==4) { print \$1\$2; counter=0 } else { print \$1 }}}' | python ${"$PWD"}/src/FastQ2BAM.py -p -s $params.s -m $params.m | python ${"$PWD"}/src/MergeTrimReadsBAM.py -p --mergeoverlap > ${datasetID}".bam" 
      
        #paste <( zcat $r1_fastq ) <(zcat \$new_3 ) <(zcat \$new_2) | awk 'BEGIN{ counter=0 }{ counter+=1; if (counter == 2) { print \$1"GATCCGGTTG"\$2\$3 } else { if (counter==4) { print \$1"IIIIIIIIII"\$2\$3; counter=0 } else { print \$1 }}' | python ${"$PWD"}/src/SplitFastQdoubleIndexBAM.py -p -s $params.s -l $params.l -m $params.m -i ${datasetID}_index.lst | python ${"$PWD"}/src/MergeTrimReadsBAM.py -p --mergeoverlap > ${datasetID}".bam" 
        paste <( zcat $r1_fastq ) <(zcat \$new_3 ) <(zcat \$new_2) | awk 'BEGIN{ counter=0 }{ counter+=1; if (counter == 2) { print \$1"GATCCGGTTG"\$2\$3 } else { if (counter==4) { print \$1"IIIIIIIIII"\$2\$3; counter=0 } else { print \$1 }}}' > ${datasetID}_combined.fastq 
        python ${"$PWD"}/src/SplitFastQdoubleIndexBAM.py -s $params.s -l $params.l -m $params.m -i ${datasetID}_index.lst --outprefix ${datasetID}_demultiplex --remove --summary ${datasetID}_combined.fastq 
        python ${"$PWD"}/src/MergeTrimReadsBAM.py --mergeoverlap --outprefix ${datasetID} ${datasetID}_demultiplex.bam


        
        #paste <( zcat $r1_fastq ) <(zcat \$new_3 ) <(zcat \$new_2) | awk 'BEGIN{ counter=0 }{ counter+=1; if (counter == 2) { print \$1\$2 } else { if (counter==4) { print \$1\$2; counter=0 } else { print \$1 }}}' | python ${"$PWD"}/src/SplitFastQdoubleIndexBAM.py -p -s $params.s -l $params.l -m $params.m -i ${datasetID}_index.lst | python ${"$PWD"}/src/MergeTrimReadsBAM.py -p --mergeoverlap > ${datasetID}".bam"
 
        """   
    }   
}   

//if no UMI
if(params.m==0){
    process 'create_BAM_noUMI' {
        tag "make idx"
    
        input:
        file(params.condaloc)
        val(params.sample_idx)
        set datasetID, file(r1_fastq) from R1_fastq 
    
        output:
        set datasetID, file("${datasetID}_index.lst") into idx_list
        set datasetID, file("${datasetID}.bam") into clean_bam   
 
        """
        #!/bin/bash
        source $params.condaloc mpraflow_py27
        echo "sample idx"
        echo $params.sample_idx
        
        echo $r1_fastq
        
        new_var=\$(echo $r1_fastq | awk -F"_R1_" '{print \$1"_R3_"\$2}')
        echo \$new_var
        new_3=$params.dir"/"\$new_var
        echo \$new_3

        echo $params.sample_idx'        '${datasetID} >> ${datasetID}_index.lst

        #paste <( zcat $r1_fastq ) <(zcat \$new_3 ) | awk 'BEGIN{ counter=0 }{ counter+=1; if (counter == 2) { print \$1\$2 } else { if (counter==4) { print \$1\$2; counter=0 } else { print \$1 }}}' | python ${"$PWD"}/src/FastQ2BAM.py -p -s $params.s | python ${"$PWD"}/src/MergeTrimReadsBAM.py -p --mergeoverlap > ${datasetID}".bam" 
       
        paste <( zcat $r1_fastq ) <(zcat \$new_3 ) | awk 'BEGIN{ counter=0 }{ counter+=1; if (counter == 2) { print \$1"GATCCGGTTG"\$2 } else { if (counter==4) { print \$1"IIIIIIIIII"\$2; counter=0 } else { print \$1 }}}' > ${datasetID}_combined.fastq
        python ${"$PWD"}/src/SplitFastQdoubleIndexBAM.py -p -s $params.s -l $params.l -m $params.m -i ${datasetID}_index.lst --outprefix ${datasetID}_demultiplex --remove --summary ${datasetID}_combined.fastq 
        python ${"$PWD"}/src/MergeTrimReadsBAM.py --mergeoverlap --outprefix ${datasetID} ${datasetID}_demultiplex.bam
        
        #paste <( zcat $r1_fastq ) <(zcat \$new_3 ) | awk 'BEGIN{ counter=0 }{ counter+=1; if (counter == 2) { print \$1\$2 } else { if (counter==4) { print \$1\$2; counter=0 } else { print \$1 }}}' | python ${"$PWD"}/src/SplitFastQdoubleIndexBAM.py -p -s $params.s -l $params.l -m $params.m -i ${datasetID}_index.lst | python ${"$PWD"}/src/MergeTrimReadsBAM.py -p --mergeoverlap > ${datasetID}".bam" 


        """   
    }
}


/*
* STEP 2: create raw counts
*/

process 'raw_counts'{

    publishDir "$params.outdir/$datasetID"

    input:
    set datasetID, file(bam) from clean_bam

    output:
    set datasetID, file("${datasetID}_raw_counts.tsv.gz") into raw_ct

    script:
    if(params.m==0)
        """
        #!/bin/bash
        source $params.condaloc mpraflow_py36
       
        samtools view -F -r $bam | awk '{print \$10}' | sort | gzip -c > ${datasetID}_raw_counts.tsv.gz 

        #samtools view -F -r $bam | awk '{print \$10}' | sort | uniq -c |  awk 'BEGIN{ OFS="\t" }{ print \$2,\$1 }' | gzip -c > ${datasetID}_raw_counts.tsv.gz  

        """

    else if(params.m!=0)
        """
        #!/bin/bash
        source $params.condaloc mpraflow_py36

        samtools view -F -r $bam | awk 'BEGIN{ OFS= "\t" }{ for (i=12; i<=NF; i++) { if (\$i ~ /^XJ:Z:/) print \$10,substr(\$i,6,16) }}' | sort | uniq -c | awk 'BEGIN{ OFS="\t" }{ print \$2,\$3,\$1 }' | awk '{if(\$2~"GGGGGGGGGGGGGGG" || \$2~"NNNNNN"); else{print}}' | gzip -c > ${datasetID}_raw_counts.tsv.gz
        """ 

}

/*
* STEP 3: Filter counts for correct barcode length
*/

process 'filter_counts'{
    publishDir "$params.outdir/$datasetID"

    input:
    set datasetID, file(rc) from raw_ct 

    output:
    set datasetID, file("${datasetID}_filtered_counts.tsv.gz") into filter_ct

    """
    #!/bin/bash
    source $params.condaloc mpraflow_py36
    
    zcat $rc | grep -v "N" | awk 'BEGIN{ OFS="\t" }{ if (length(\$1) == 15) { print } }' | gzip -c > ${datasetID}_filtered_counts.tsv.gz

    """

}

/*
* STEP 4: Record overrepresended UMIs and final count table
*/

process 'final_counts'{

    publishDir "$params.outdir/$datasetID"

    input:
    set datasetID, file(fc) from filter_ct

    output:
    set datasetID, file("${datasetID}_counts.tsv") into final_count

    script:
    if(params.m==0)
        """
        #!/bin/bash
        source $params.condaloc mpraflow_py36
        
        zcat $fc | awk '{print \$1}' | uniq -c > ${datasetID}_counts.tsv
        
        """
    else if(params.m!=0)
        """
        
        #!/bin/bash
        source $params.condaloc mpraflow_py36
       
        for i in $fc; do echo \$(basename \$i); zcat \$i | cut -f 2 | sort | uniq -c | sort -nr | head; echo; done > ${params.outdir}/${datasetID}/${datasetID}_freqUMIs.txt
        zcat $fc | awk '{print \$1}' | uniq -c > ${datasetID}_counts.tsv

        """

}

/*
* STEP 5: MPRAnalyze input generation (if option selected)
*/

//MPRAnalyze option
if(params.mpranalyze != 0){
    /*
    * STEP 5: Merge each DNA and RNA file
    */
    
    process 'dna_rna_mpranalyze_merge'{
        publishDir "$params.outdir/", mode:'copy'
    
        input:
        file(clist) from final_count.collect()
    
        output:
        file "*tmp.csv" into merged_ch
    
        """
        #!/bin/bash
        source $params.condaloc mpraflow_py36
        #
        python ${"$PWD"}/src/merge_counts.py ${params.e} ${"$PWD"}/${params.outdir}/
        
        """
    }
    
    
    /*
    * STEP 6: Merge all DNA/RNA counts into one big file
    */
    
    process 'final_merge'{
        publishDir "$params.outdir/"
    
        input:
        file(pairlist) from merged_ch.collect()
    
        output:
        file "*.csv" into merged_out
    
        """
        #!/bin/bash
        source $params.condaloc mpraflow_py36
        python ${"$PWD"}/src/merge_all.py ${params.e} ${"$PWD"}/${params.outdir}/ ${params.out} ${params.design}
    
        """
    }
    
    
    /*
    * STEP 7: Add label to outfile
    */
    
    process 'final_label'{
    
        publishDir "$params.outdir/", mode:'copy'
    
        input:
        file(table) from merged_out
    
        output:
        file "*_final_labeled_counts.txt" into labeled_out
    
        """
        #!/bin/bash
        source $params.condaloc mpraflow_py36
        python ${"$PWD"}/src/label_final_count_mat.py $table ${params.association} ${params.out}"_final_labeled_counts.txt"  ${params.design}
        """
    }
    
    /*
    * STEP 8: Generate inputs
    */
    
    process 'gen_mpranalyze'{
        publishDir "$params.outdir/", mode:'copy'
        
        input:
        file(t) from labeled_out
    
        """
        #!/bin/bash
        source $params.condaloc mpraflow_py36
        python ${"$PWD"}/src/mpranalyze_compiler.py $t ${"$PWD"}/$params.outdir/
       
        """
    }

}


/*
* STEP 5: Merge each DNA and RNA file label with sequence and enhancer and normalize
*/
//merge and normalize
if(params.mpranalyze == 0){

    process 'dna_rna_merge'{
        publishDir "$params.outdir/", mode:'copy'
       
     
        input:
        file(clist) from final_count.collect()
    
        output:
        file "*.tsv" into merged_ch
    
        """
        #!/bin/bash
        source $params.condaloc mpraflow_py36
        
        #run this in parallel not the most elegant solution, but seems to work
        iter='a'
        itera='b'
        head -1 ${params.e} > tmp.header.txt
        sed 1d ${params.e} | while read d; do itera=\$itera\$iter; echo \${itera}; cat tmp.header.txt > tmp.file_\${itera}.txt; echo \$d >> tmp.file_\${itera}.txt; python ${"$PWD"}/src/merge_label.py tmp.file_\${itera}.txt ${"$PWD"}/${params.outdir}/ ${params.association} ${params.design} ${params.merge_intersect} & done 
        #python ${"$PWD"}/bin/merge_label.py ${params.e} ${"$PWD"}/${params.outdir}/	${params.association} ${params.design} ${params.merge_intersect}
        sleep 2m
        wait
        echo 'all jobs are done!'
     
        """
    
    }
    
    
    
    
    
    /*
    * STEP 6: Calculate correlations between Replicates
    */
    process 'calc_correlations'{
        publishDir "$params.outdir/", mode:'copy'
    
        input:
        file(pairlist) from merged_ch.collect()
    
        output:
        file "*.png"
    
        """
        #!/bin/bash
        source $params.condaloc mpraflow_py36
        Rscript ${"$PWD"}/src/plot_perInsertCounts_correlation.R ${params.e} ${"$PWD"}/${params.outdir}/ ${"$PWD"}/${params.outdir}/${params.out} ${params.labels}
    
        """
    }

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
