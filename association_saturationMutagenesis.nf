#!/usr/bin/env nextflow

/*
========================================================================================
                         MPRAflow
========================================================================================
MPRA Analysis Pipeline. Started 2019-07-29.
Saturation Mutagenesis Association Utility

#### Homepage / Documentation
https://github.com/shendurelab/MPRAflow
#### Authors
Max Schubach <max.schubach@bihealth.de>
Martin Kircher <martin.kircher@bihealth.de>
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
      --fastq-insert                Full path to library association fastq for insert
      --fastq-insertPE              Full path to library association fastq for read2
      --fastq-bc                    Full path to library association fastq for bc
      --design                      Full path to fasta of reference sequence (only one reference sequence)
      --name                        Name of the association. Files will be named after this.

    Options:
      
      --bc-length                   Barcode length (default 15)
      --clipping-penalty            bwa mem clipping penalty (default 80)
      --min-ireads                  minimum number gapped reads for indel candidates (default: 3)
      --outdir                      The output directory where the results will be saved and what will be used as a prefix (default outs)

    Extras:
      --h, --help                   Print this help message
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

//defaults
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
  exit 1, "Fastq insert file not specified with --fastq-insertPE"
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
// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// BC length
if (params.containsKey("bc-length")){
    params.bc_length = params["bc-length"]
} else {
    params.bc_length = 15
}

// BC length
if (params.containsKey("clipping-penalty")){
    params.bc_length = params["clipping-penalty"]
} else {
    params.clipping_penalty = 80
}

// min-ireads
if (params.containsKey("min-ireads")){
    params.min_ireads = params["min-ireads"]
} else {
    params.min_ireads = 3
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
summary['BC length']        = params.bc_length
summary['clipping_penalty'] = params.clipping_penalty
summary['min ireads']       = params.min_ireads
summary['Output dir']       = params.outdir
summary['Run name']         = params.name
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['base directory']   = "$baseDir"
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile

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
* remove the illegal regex characters from reference
* contributions: GMax Schubach
*/

process 'clean_design' {
    tag 'count'
    label 'shorttime'
    publishDir "${params.outdir}/${params.name}", mode:'copy'

    input:
        file(design) from params.design_file
    output:
        file "design_rmIllegalChars.fa" into fixed_design,fixed_design2
    shell:
        """
        cat $design | awk '{gsub(/\\[/,"_")}1' | awk '{gsub(/\\]/,"_")}1' > design_rmIllegalChars.fa
        """
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
    output:
        file "${design}.fai" into reference_fai,reference_fai2
        file "${design}.bwt" into reference_bwt
        file "${design}.sa" into reference_sa
        file "${design}.pac" into reference_pac
        file "${design}.ann" into reference_ann
        file "${design}.amb" into reference_amb
        file "${design}.dict" into reference_dict
    shell:
        """
        bwa index -a bwtsw $design
        samtools faidx $design
        picard CreateSequenceDictionary REFERENCE=$design OUTPUT=$design".dict"
        """
}

process 'get_name' {
    label 'shorttime'

    conda 'conf/mpraflow_py36.yml'

    input:
        file(design_fai) from reference_fai
    output:
         env ELEMENT into element,element2,element3,element4
    shell:
        """
        ELEMENT=\$(cat $design_fai | head -n 1 | cut -f 1)
        """
}

Channel
    .fromPath(params.fastq_insert_file)
    .splitFastq( by: params.split, file: true, compress: true)
    .set{ FWD_ch }
Channel
    .fromPath(params.fastq_insertPE_file)
    .splitFastq( by: params.split, file: true, compress: true)
    .set{ REV_ch }
Channel
    .fromPath(params.fastq_bc_file)
    .splitFastq( by: params.split, file: true, compress: true)
    .set{ BC_ch }

/*
* Merge reads
* contributions: Max Schubach
*/
process 'create_BAM' {
    label 'longtime'

    conda 'conf/mpraflow_py27.yml'

    input:
        file(fw_fastq) from FWD_ch
        file(rev_fastq) from REV_ch
        file(fastq_bc_file) from BC_ch
        val datasetID from element
        val bc_length from params.bc_length
    output:
        file "${datasetID}.*.bam" into clean_bam
    shell:
        """
        fwd_length=`zcat $fw_fastq | head -2 | tail -1 | wc -c`;
        fwd_length=\$(expr \$((\$fwd_length-1)));
        rev_start=\$(expr \$((\$fwd_length+$bc_length+1)));
        echo \$rev_start

        paste <( zcat $fw_fastq ) <( zcat $fastq_bc_file | cut -c 1-$bc_length ) <( zcat $rev_fastq ) | \
        awk '{
            if (NR % 4 == 2 || NR % 4 == 0) {
                print \$1\$2\$3
            } else {
                print \$1
            }
        }' | \
        python ${"$baseDir"}/src/FastQ2doubleIndexBAM.py -p -s \$rev_start -l $bc_length -m 0 --RG ${datasetID} |
        python ${"$baseDir"}/src/MergeTrimReadsBAM.py --FirstReadChimeraFilter '' --adapterFirstRead '' --adapterSecondRead '' -p --mergeoverlap \
        > ${datasetID}.${fw_fastq}.bam
    
        """
}

/*
*COLLCT FASTQ CHUNCKS
*/

process 'collect_chunks'{
    label 'shorttime'

    conda 'conf/mpraflow_py36.yml'

    input:
        file bams from clean_bam.collect()
        val datasetID from element2
    output:
        tuple val(datasetID),file("${datasetID}.bam") into bams_merged
    script:
        bam_list = bams.join(' ')
    shell:
        """
        #collect sorted bams into one file
        samtools merge ${datasetID}.bam $bam_list
        """
}


/*
* Mapp reads
* contributions: Max Schubach
*/
process 'PE_mapping' {
    tag 'align'
    label 'longtime'

    conda 'conf/mpraflow_py36.yml'

    input:
        tuple val(datasetID),file(bam) from bams_merged
        file(design) from fixed_design
        file(reference_fai) from reference_fai
        file reference_bwt from reference_bwt
        file reference_sa from reference_sa
        file reference_pac from reference_pac
        file reference_ann from reference_ann
        file reference_amb from reference_amb
        file reference_dict from reference_dict
        val clipping from params.clipping_penalty
    output:
        tuple val(datasetID),file("aligned_${datasetID}.bam"),file("aligned_${datasetID}.bam.bai") into aligned_bam,aligned_bam2
        file "aligned_${datasetID}.bam_stats" into aligned_bam_stats
    shell:
        """
        (
            samtools view -H $bam | grep '^@RG';
            bwa mem -L $clipping -M -C $design <(
                samtools view -F 513 $bam | \
                awk 'BEGIN{ FS="\\t"; OFS="\\n" }{ 
                    split(\$0,a,"\\t");
                    helper = ""; 
                    for (i=12; i <= length(a); i++) { helper = helper""a[i]"\\t"}; 
                    sub("\\t\$","",helper); 
                    print "@"\$1" "helper,\$10,"+",\$11;
                }'
            ); 
            bwa mem -L $clipping -M -C $design <(
                samtools view -f 64 $bam | \
                awk 'BEGIN{ FS="\\t"; OFS="\\n" }{
                    split(\$0,a,"\\t"); 
                    helper = ""; 
                    for (i=12; i <= length(a); i++) { helper = helper""a[i]"\\t"};
                    sub("\\t\$","",helper); print "@"\$1" "helper,\$10,"+",\$11 
                }'
            ) <(
                samtools view -f 128 $bam | \
                awk 'BEGIN{ FS="\\t"; OFS="\\n" }{
                    split(\$0,a,"\\t");
                    helper = ""; 
                    for (i=12; i <= length(a); i++) { helper = helper""a[i]"\\t"}; 
                    sub("\\t\$","",helper); print "@"\$1" "helper,\$10,"+",\$11 
                }'
            ) | grep -v "^@" \
        ) | samtools view -Su - | \
        samtools sort - > aligned_${datasetID}.bam;

        samtools index aligned_${datasetID}.bam;

        samtools flagstat aligned_${datasetID}.bam > aligned_${datasetID}.bam_stats
        """
}



/*
* Process 1C: align with BWA
* contributions: Max Schubach
*/
process 'get_count' {
    label 'shorttime'

    conda 'conf/mpraflow_py27.yml'

    input:
        tuple val(datasetID),file(bam),file(bai) from aligned_bam
        file fixed_design from fixed_design
        val bc_length from params.bc_length
    output:
        tuple val(datasetID),file("counts_${datasetID}.filtered.tsv.gz") into filtered_counts
        file "counts_${datasetID}.tsv.gz" into counts
    shell:
        """
        python ${"$baseDir"}/src/satMut/extractBarcodesBAM.py --BAMfield -f $bc_length $bam | \
        sort -k1,1 -k2,2 | uniq -c | \
        awk 'BEGIN{OFS="\\t"}{ print \$2,\$3,\$1}' | gzip -c > counts_${datasetID}.tsv.gz;

        zcat counts_${datasetID}.tsv.gz | \
        python ${"$baseDir"}/src/satMut/filterAssignmentBarcodes.py -f $fixed_design \
        -r ${datasetID}:20:-20 | \
        gzip -c > counts_${datasetID}.filtered.tsv.gz
        """
}

process 'extract_reads' {
    label 'highmem'

    conda 'conf/mpraflow_py27.yml'

    input:
        tuple val(datasetID),file(counts) from filtered_counts
        tuple val(datasetID_bam),file(bam),file(bai) from aligned_bam2
        file fixed_design from fixed_design
        val bc_length from params.bc_length
    output:
        file("reads/*/*.bam") into reads
    shell:
        """
        python ${"$baseDir"}/src/satMut/extractReadsAssignmentSimple.py --BAMfield -f $bc_length -a $counts -m 10000 -o reads $bam;
        """
}

reads
    .flatten()
    .map { file -> tuple(file.name[0..4], file) }
    .groupTuple()
    .set { grouped_reads }



process 'call_variants' {
    label 'shorttime'

    conda 'conf/mpraflow_py27.yml'

    input:
        tuple val(prefix),file(read_bam) from grouped_reads
        file(design_fai) from reference_fai2
        file(design) from fixed_design2
        val datasetID from element3
        val(m) from params.min_ireads
    output:
        file("variants_${prefix}.txt.gz") into prefix_variants
    script:
        read_bam_list = read_bam.collect{"$it"}.join(' ')
    shell:
        """
        region=\$(grep -i $datasetID $design_fai | awk '{ print \$1":20-"\$2-20 }');
        for i in $read_bam_list; do
            bcftools mpileup -A -m $m --fasta-ref $design \$i | \
            bcftools call -c -f GQ | \
            python ${"$baseDir"}/src/satMut/extractVariants.py -r \$region;
        done | gzip -c > variants_${prefix}.txt.gz
        """
}


process 'combine_variants' {
    publishDir "${params.outdir}/${params.name}", mode:'copy'
    label 'shorttime'

    conda 'conf/mpraflow_py36.yml'

    input:
        file(variants) from prefix_variants.collect()
        val element from element4
    output:
        file("${element}_variants.txt") into final_variants
    script:
        variant_list = variants.collect{"$it"}.join(' ')
    shell:
        """
        zcat $variant_list > ${element}_variants.txt
        """
}
