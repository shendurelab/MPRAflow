#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Merge/Adapter trim reads stored in BAM

:Author: Martin Kircher
:Contact: martin.kircher@bihealth.de
:Date: *19.04.2018
:Type: tool
:Input: BAM
:Output: BAM

"""

import sys, os
import math
import random

import pysam
from optparse import OptionParser,OptionGroup
import string

from MergeTrimReads import set_adapter_sequences, set_options, set_keys, process_SR, process_PE
table = string.maketrans('TGCA','ACGT') # COMPLEMENT DNA

maxadapter_comp = 30
min_length = 5
notified = set()

def revcompl(seq):
  global table
  seq=seq.translate(table)
  return seq[::-1]

# ADD FLAG TO QUALITY FILTER REAG TAG
def add_quality_flag(tags,qtype):
  qc_tag_helper = []
  foundzq = False
  if tags != None:
    for tag,value in tags:
      if tag == "ZQ": 
        foundzq = True
        qc_tag = list(set(list(value+qtype)))
        qc_tag.sort()
        qc_tag_helper.append((tag,"".join(qc_tag)))
      else: qc_tag_helper.append((tag,value))
  if not foundzq: qc_tag_helper.append(("ZQ",qtype))
  return qc_tag_helper

# CREATE NEW BAM READ FROM ORIGINAL READS AND MERGED SEQUENCE/QUALITY
def create_new_read(new_seq,new_qual,read1,read2):
  global notified
  new = pysam.AlignedRead()
  if not read1.qname.startswith("M_"): new.qname = "M_"+read1.qname
  else: new.qname = read1.qname
  new.seq = new_seq
  new.qual = new_qual
  new.is_unmapped = True
  new.pos = -1
  new.mpos = -1

  new.is_qcfail = read1.is_qcfail and read2.is_qcfail
  if read1.tags != None: htags = dict(read1.tags)
  else: htags = {}
  if (len(new_seq) < min_length):
    new.is_qcfail = True
    if "ZQ" in htags: htags["ZQ"]+="L"
    else: htags["ZQ"]="L"
  stags = set()
  new_tags = []
  if read2.tags != None:
    for tag,value in read2.tags:
      stags.add(tag)
      if tag == "NM" or tag == "MD": continue
      elif tag in htags and value != htags[tag]: # NEW TAG DIFF VALUE
        if tag == "ZQ":
          qc_tag = list(set(list(value+htags[tag])))
          qc_tag.sort()
          new_tags.append((tag,"".join(qc_tag)))
        else:
          if tag not in notified:
            sys.stderr.write("Do not know how to combine %s BAM tags. Information of one of the reads will get lost during merging.\n"%tag)
            notified.add(tag)
      elif tag in htags and value == htags[tag]: # SAME TAG AND VALUE
        new_tags.append((tag,value))
      else: # NEW TAG
        new_tags.append((tag,value))
  for tag,value in htags.iteritems():
    if tag not in stags: new_tags.append((tag,value))
  new.tags = new_tags
  return new


options_adapter_F="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG"
options_adapter_S="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT"
options_adapter_chimera="ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA"

parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from and write it to PIPE",default=False,action="store_true")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'MergeTrim').",default="MergeTrim")
parser.add_option("", "--SAM", dest="SAM", help="Output SAM not BAM.",default=False,action="store_true")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")

group = OptionGroup(parser, "Paired End merging/Single Read trimming  options")
group.add_option("--mergeoverlap",dest="mergeoverlap",help="Merge PE reads of molecules longer than read length that show a minimum overlap (default False)",default=False,action="store_true")
group.add_option("--onlyoverlap",dest="onlyoverlap",help="Only sequence overlapping between two reads (default False)",default=False,action="store_true")
group.add_option("-f", "--adapterFirstRead", dest="adapter_F", help="Adapter that is observed after the forward read (def. Multiplex: %s)"%options_adapter_F[:maxadapter_comp],default=options_adapter_F)
group.add_option("-s", "--adapterSecondRead", dest="adapter_S", help="Adapter that is observed after the reverse read (def. Multiplex: %s)"%options_adapter_S[:maxadapter_comp],default=options_adapter_S)
group.add_option("-c", "--FirstReadChimeraFilter", dest="adapter_chimera", help="If the forward read looks like this sequence, the cluster is filtered out. Provide several sequences separated by comma.(def. Multiplex: %s)"%(",".join(map(lambda x:x[:maxadapter_comp],options_adapter_chimera.split(',')))),default=options_adapter_chimera)
group.add_option("-k","--key",dest="key",help="Key sequence with which each sequence starts. Comma separate for forward and reverse reads. (def. '')",default="")
group.add_option("-i","--allowMissing",dest="allowMissing",help="Allow one base in one key to be missing or wrong. (def. False)",default=False,action="store_true")
group.add_option("-t","--trimCutoff",dest="trimCutoff",help="Lowest number of adapter bases to be observed for Single Read trimming (default 1)",default=1,type="int")
group.add_option("--insertOne",dest="insertOne",help="For reads failing regular merging (not overlapping read ends!), check whether they can be merged by inserting 1 base (default off)",default=False,action="store_true")
group.add_option("--noSecondTrim",dest="noSecondTrim",help="Do not trim single end reads starting in T_ (def. False)",default=False,action="store_true")
group.add_option("--minoverlap",dest="minoverlap",help="Minimum overlap of merged PE reads (default 11)",default=11,type="int")
parser.add_option_group(group)
(options, args) = parser.parse_args()

if options.onlyoverlap: options.mergeoverlap = True

set_adapter_sequences(options.adapter_F, options.adapter_S, options.adapter_chimera, maxadapter_comp)
set_options(options.trimCutoff, options.allowMissing, options.mergeoverlap, options.onlyoverlap, options.minoverlap-1)
if not set_keys(options.key):
  sys.exit()

if options.outprefix == "":
  sys.stderr.write("Outprefix can not be empty!\n")
  sys.exit()

if options.outdir != None and not os.path.isdir(options.outdir):
  sys.stderr.write("Output folder does not exist!\n")
  sys.exit()
elif options.outdir == None:
  options.outdir = ""
else:
  options.outdir = options.outdir.rstrip('/')+'/'

################################
## PROCESS FILES
################################

## CREATE OUTPUT FILE(S)/STREAM
fileflags = 'wb'
if options.SAM: fileflags = 'w'

files = args
if options.pipe: files = [None]

outfile = None
count_all = 0
count_fkey = 0
count_merged = 0
count_merged_overlap = 0
count_trimmed = 0
count_nothing = 0
count_chimera = 0

for filename in files:
  if filename == None:
    infile = pysam.Samfile( "-", 'rb' )
  else:
    infile = pysam.Samfile( filename, 'rb' )

  cheader = infile.header
  if ('HD' in cheader)  and ('SO' in cheader['HD']) and (cheader['HD']['SO'] != 'queryname'):
    sys.stderr.write("Input BAM has wrong sorting order. Skipping input...\n")
    continue

  if outfile == None:
    if options.verbose: sys.stderr.write("Creating output files/streams...\n")

    if options.pipe:
      if len(cheader) == 0:
        cheader['HD'] = {'VN': '1.0'}
      outfile = pysam.Samfile( "-", fileflags, header = cheader)
      if options.verbose: sys.stderr.write("BAM/SAM output on stdout...\n")
    else:
      outfilename = options.outdir+options.outprefix+".bam"
      if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
      outfile = pysam.Samfile( outfilename, fileflags, header = cheader)

  bread = pysam.AlignedRead()
  for read in infile:
    if read.is_paired and bread.qname == "":
      bread = read
      continue
    elif read.is_paired and bread.qname != "":
      if read.qname == bread.qname:
        count_all += 1
        # SAVE SEQUENCES IN MERGING VARIABLES
        if bread.is_read1:
          #sys.stderr.write("%s %s\n"%(bread.is_reverse,read.is_reverse))
          #sys.stderr.write("+ %s\n  %s\n"%(bread.seq,read.seq))
          if bread.is_reverse:
            read1 = revcompl(bread.seq)
            qual1 = bread.qual[::-1]
          else:
            read1 = bread.seq
            qual1 = bread.qual
          if read.is_reverse:
            read2 = revcompl(read.seq)
            qual2 = read.qual[::-1]
          else:
            read2 = read.seq
            qual2 = read.qual
          #sys.stderr.write("- %s\n  %s\n"%(read1,read2))
            
        else:
          #sys.stderr.write("%s %s\n"%(read.is_reverse,bread.is_reverse))
          #sys.stderr.write("+ %s\n  %s\n"%(read.seq,bread.seq))
          if read.is_reverse:
            read1 = revcompl(read.seq)
            qual1 = read.qual[::-1]
          else:
            read1 = read.seq
            qual1 = read.qual
          if bread.is_reverse:
            read2 = revcompl(bread.seq)
            qual2 = bread.qual[::-1]
          else:
            read2 = bread.seq
            qual2 = bread.qual
          #sys.stderr.write("- %s\n  %s\n"%(read1,read2))
            
        # IF WE HAVE NO QUALITY SCORES, ASSUME EQUAL QUALITY SCORES OF 15
        if qual1 == "*": qual1 = len(read1)*"0"
        if qual2 == "*": qual2 = len(read2)*"0"

        flag,newseq,newqual = process_PE(read1,qual1,read2,qual2)

        if flag != '':
          bread.is_qcfail = True
          bread.tags = add_quality_flag(bread.tags,flag)
          read.is_qcfail = True
          read.tags = add_quality_flag(read.tags,flag)
          if flag == 'K': count_fkey += 1
          elif flag == 'D': count_chimera += 1

        if newseq != '':
          if flag == '': 
            if len(newseq) > max(len(read1),len(read2)): count_merged_overlap += 1
            else: count_merged += 1
          new = create_new_read(newseq,newqual,bread,read)
          outfile.write(new)
        else:
          if options.insertOne and flag == '':
            success = False

            # SET OPTIONS TO DEACTIVATE OVERLAP MERGING
            set_options(options.trimCutoff, options.allowMissing, False, False,options.minoverlap-1)

            # CHECK WHETHER IT ACTUALLY COULD BE AN INDEL PROBLEM
            segment = False
            if len(read1) <= len(read2):
              cutpoint = len(read1)/2
              nread1 = read1[cutpoint:]
              nqual1 = qual1[cutpoint:]
              flag_,newseq_,newqual_ = process_PE(nread1,nqual1,read2,qual2)
              #sys.stderr.write("1a: %s\n    %s\n"%(nread1,revcompl(read2)))
              if newseq_ != "": segment = True
              else:
                nread1 = read1[:cutpoint]
                nqual1 = qual1[:cutpoint]
                flag_,newseq_,newqual_ = process_PE(nread1,nqual1,read2,qual2)
                #sys.stderr.write("4a: %s\n    %s\n"%(nread1,revcompl(read2)))
                if newseq_ != "": segment = True
            else:
              cutpoint = len(read2)/2
              nread2 = read2[cutpoint:]
              nqual2 = qual2[cutpoint:]
              flag_,newseq_,newqual_ = process_PE(read1,qual1,nread2,nqual2)
              #sys.stderr.write("1b: %s\n    %s\n"%(read1,revcompl(nread2)))
              if newseq_ != "": segment = True
              else:
                nread2 = read2[:cutpoint]
                nqual2 = qual2[:cutpoint]
                #sys.stderr.write("4b: %s\n    %s\n"%(read1,revcompl(nread2)))
                flag_,newseq_,newqual_ = process_PE(read1,qual1,nread2,nqual2)
                if newseq_ != "": segment = True
              
            if segment:
              for i in range(len(read1)): # CHECK INSERTION IN FORWARD READ
                nread1 = read1[:i]+"I"+read1[i:]
                nqual1 = qual1[:i]+"!"+qual1[i:]
                flag_,newseq_,newqual_ = process_PE(nread1,nqual1,read2,qual2)
                
                if flag_ != '':
                  bread.is_qcfail = True
                  bread.tags = add_quality_flag(bread.tags,flag_)
                  read.is_qcfail = True
                  read.tags = add_quality_flag(read.tags,flag_)
                  if flag_ == 'K': count_fkey += 1
                  elif flag_ == 'D': count_chimera += 1
                  outfile.write(bread)
                  outfile.write(read)
                  success = True
                  break
                elif newseq_ != '':
                  #sys.stderr.write("Index of insert in forward read: %d of %d\n"%(i,len(read1)))
                  if flag_ == '': 
                    if len(newseq_) > max(len(read1)+1,len(read2)): count_merged_overlap += 1
                    else: count_merged += 1
                  new = create_new_read(newseq_,newqual_,bread,read)
                  outfile.write(new)
                  success = True
                  break

              if not success:
                for i in range(len(read2)): # CHECK INSERTION IN REVERSE READ
                  nread2 = read2[:i]+"I"+read2[i:]
                  nqual2 = qual2[:i]+"!"+qual2[i:]
                  flag_,newseq_,newqual_ = process_PE(read1,qual1,nread2,nqual2)
                  
                  if flag_ != '':
                    bread.is_qcfail = True
                    bread.tags = add_quality_flag(bread.tags,flag_)
                    read.is_qcfail = True
                    read.tags = add_quality_flag(read.tags,flag_)
                    if flag_ == 'K': count_fkey += 1
                    elif flag_ == 'D': count_chimera += 1
                    outfile.write(bread)
                    outfile.write(read)
                    success = True
                    break
                  elif newseq_ != '':
                    #sys.stderr.write("Index of insert in reverse read: %d of %d\n"%(i,len(read2)))
                    if flag_ == '': 
                      if len(newseq_) > max(len(read1),len(read2)+1): count_merged_overlap += 1
                      else: count_merged += 1
                    new = create_new_read(newseq_,newqual_,bread,read)
                    outfile.write(new)
                    success = True
                    break
              
            if not success:
              if flag == '': count_nothing += 1
              outfile.write(bread)
              outfile.write(read)
            # RESET OPTIONS
            set_options(options.trimCutoff, options.allowMissing, options.mergeoverlap, options.onlyoverlap, options.minoverlap-1)
          else:
            if flag == '': count_nothing += 1
            outfile.write(bread)
            outfile.write(read)
        bread = pysam.AlignedRead()
      else:
        if options.verbose: 
          sys.stderr.write('Warning: Skipping read from incomplete pair\n')
        bread = read
    else: # Single End read
      count_all += 1
      if options.noSecondTrim and read.qname.startswith("T_"):
        count_nothing += 1
        outfile.write(read)
        continue
      
      read1 = read.seq
      qual1 = read.qual
      if qual1 == "*": qual1 = len(read1)*"0"
      flag,newseq,newqual = process_SR(read1,qual1)

      if flag != '':
        read.is_qcfail = True
        read.tags = add_quality_flag(read.tags,flag)
        if flag == 'K': count_fkey += 1
        elif flag == 'D': count_chimera += 1

      if newseq != '':
        if not read.qname.startswith("T_"): read.qname = "T_"+read.qname
        read.seq = newseq
        read.qual = newqual
        read.is_unmapped = True
        read.pos = -1
        read.mpos = -1
        read.cigar = None
        read.mapq = 0
        read.tid = 0
        read.mrnm = 0
        new_tags = []
        for tag,value in read.tags:
          if tag != "NM" and tag != "MD": new_tags.append((tag,value))
        read.tags = new_tags
        if flag == '': count_trimmed += 1
        outfile.write(read)
      else:
        if flag == '': count_nothing += 1
        outfile.write(read)

    if options.verbose and (count_all % 10000 == 0) and (count_all > 0):
        sys.stderr.write("<==> Total %d; Merged (trimming) %d; Merged (overlap) %d; Kept PE/SR %d; Trimmed SR %d; Adapter dimers/chimeras %d; Failed Key %d\n"%(count_all,count_merged,count_merged_overlap,count_nothing,count_trimmed,count_chimera,count_fkey))

  if not options.pipe:
    infile.close()

sys.stderr.write("Total %d; Merged (trimming) %d; Merged (overlap) %d; Kept PE/SR %d; Trimmed SR %d; Adapter dimers/chimeras %d; Failed Key %d\n"%(count_all,count_merged,count_merged_overlap,count_nothing,count_trimmed,count_chimera,count_fkey))
if outfile != None:
  outfile.close()
