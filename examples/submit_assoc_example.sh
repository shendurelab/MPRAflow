#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o .                        #-- output directory (fill in)
#$ -e .                        #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=1G                  #-- submits on nodes with enough free memory (required)
##$ -l arch=linux-x64               #-- SGE resources (CPU type)
##$ -l netapp=1G,scratch=1G         #-- SGE resources (home and scratch disks)
#$ -l h_rt=24:00:00                #-- runtime limit (see above; this requests 24 hours)
##$ -t 1-10                        #-- remove first '#' to specify the number of
                                   #-- tasks if desired (see Tips section)

source ~/miniconda3/bin/activate mpraflow_py36


dir='$FASTQ_PATH'
pref=$x

echo $pref

nextflow run association.nf --fastq_insert ${dir}${pref}"_R1_001.fastq.gz" --fastq_insertPE ${dir}${pref}"_R3_001.fastq.gz" --design "example_design.fa" --fastq_bc ${dir}${pref}"_R2_001.fastq.gz"

