#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
Extract the index sequence from the middle and end of an Illumina run. Separates reads for Paired End runs.

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *06.12.2019
:Type: tool
:Input: fasta, fastq
:Output: BAM

"""

import sys
import os
import string
table = string.maketrans('.','N')

import pysam
from collections import defaultdict
from library import read_fastq
from optparse import OptionParser,OptionGroup

ireadlength1 = None
ireadlength2 = None

INF_QUALITY = 200

def read_sequence_file(infile,sec_read_start=None):
  global options
  global ireadlength1,ireadlength2
  if ireadlength2 != 0:
    ireadl2 = -ireadlength2
  else:
    ireadl2 = None
  # seqid,oindseq1,oindqual1,oindseq2,oindqual2,r1,q1,r2,q2
  if options.nextera and ireadl2 != None:
    for seqid, seq, qual in read_fastq(infile):
      seq = seq.translate(table)
      if qual != None and options.qualityoffset != 33: qual = "".join(map(lambda x:chr(ord(x)-options.qualityoffset+33),qual))
      if sec_read_start == None: 
        if qual != None:
          yield seqid,seq[-(ireadlength1+ireadlength2):-ireadlength2],qual[-(ireadlength1+ireadlength2):-ireadlength2],seq[-ireadlength2:],qual[-ireadlength2:],seq[:-(ireadlength1+ireadlength2)],qual[:-(ireadlength1+ireadlength2)],None,None
        else:
          yield seqid,seq[-(ireadlength1+ireadlength2):-ireadlength2],None,seq[-ireadlength2:],None,seq[:-(ireadlength1+ireadlength2)],None,None,None
      else:
        if qual != None:
          yield seqid,seq[(sec_read_start-ireadlength1-ireadlength2):(sec_read_start-ireadlength2)],qual[(sec_read_start-ireadlength1-ireadlength2):(sec_read_start-ireadlength2)],seq[(sec_read_start-ireadlength2):(sec_read_start)],qual[(sec_read_start-ireadlength2):(sec_read_start)],seq[:(sec_read_start-ireadlength1-ireadlength2)],qual[:(sec_read_start-ireadlength1-ireadlength2)],seq[sec_read_start:],qual[sec_read_start:]
        else:
          yield seqid,seq[(sec_read_start-ireadlength1-ireadlength2):(sec_read_start-ireadlength2)],None,seq[(sec_read_start-ireadlength2):(sec_read_start)],None,seq[:(sec_read_start-ireadlength1-ireadlength2)],None,seq[sec_read_start:],None

  else:
    for seqid, seq, qual in read_fastq(infile):
      seq = seq.translate(table)
      if qual != None and options.qualityoffset != 33: qual = "".join(map(lambda x:chr(ord(x)-options.qualityoffset+33),qual))
      if sec_read_start == None:
        if qual != None:
          yield seqid,seq[-ireadlength1:],qual[-ireadlength1:],"",None,seq[:-ireadlength1],qual[:-ireadlength1],None,None
        else:
          yield seqid,seq[-ireadlength1:],None,"",None,seq[:-ireadlength1],None,None,None
      else:
        if qual != None:
          if ireadl2 == None:
            yield seqid,seq[(sec_read_start-ireadlength1):sec_read_start],qual[(sec_read_start-ireadlength1):sec_read_start],"",None,seq[:(sec_read_start-ireadlength1)],qual[:(sec_read_start-ireadlength1)],seq[sec_read_start:],qual[sec_read_start:]
          else:
            yield seqid,seq[(sec_read_start-ireadlength1):sec_read_start],qual[(sec_read_start-ireadlength1):sec_read_start],seq[ireadl2:],qual[ireadl2:],seq[:(sec_read_start-ireadlength1)],qual[:(sec_read_start-ireadlength1)],seq[sec_read_start:ireadl2],qual[sec_read_start:ireadl2]
        else:
          if ireadl2 == None:
            yield seqid,seq[(sec_read_start-ireadlength1):sec_read_start],None,"",None,seq[:(sec_read_start-ireadlength1)],None,seq[sec_read_start:],None
          else:
            yield seqid,seq[(sec_read_start-ireadlength1):sec_read_start],None,seq[ireadl2:],None,seq[:(sec_read_start-ireadlength1)],None,seq[sec_read_start:ireadl2],None

  raise StopIteration


parser = OptionParser("%prog [options] seq_files")
parser.add_option("-p","--PIPE",dest="pipe",help="Read from and write to PIPE",default=False,action="store_true")
parser.add_option("-s","--start",dest="start",help="First base of the second read (default None)",type='int')
parser.add_option("-l","--ireadlength",dest="ireadlength1",help="Length of index read (in case it is longer than the indexes in the index file provided)",type="int",default=0)
parser.add_option("-m","--2nd_ireadlength",dest="ireadlength2",help="Length of a second index read (in case it is longer than the indexes in the index file provided)",type="int",default=0)
parser.add_option("--qualityoffset",dest="qualityoffset",help="Offset of quality scores in FastQ input file (default 33)",type="int",default=33)
parser.add_option("--nextera",dest="nextera",help="Order of reads like Nextera",default=False, action="store_true")
parser.add_option("-r","--RG",dest="readgroup",help="Set read group for all reads",default="")

group = OptionGroup(parser, "Quality filtering of index read")
group.add_option("-q","--quality",dest="quality",help="Minimum quality score of bases in both indexes (default 0:off)",type='int',default=0)
group.add_option("--qualityN",dest="qualityN",help="Consider subset of bases for quality score filter (e.g. NyyyyyN removes first and last base)",default="")
group.add_option("--2nd_qualityN",dest="qualityN_2nd",help="Consider subset of bases for quality score filter (e.g. NyyyyyN removes first and last base)",default="")
parser.add_option_group(group)

group = OptionGroup(parser, "Output options")
group.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
group.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'demultiplex').",default="demultiplex")
group.add_option("", "--SAM", dest="SAM", help="Output SAM not BAM.",default=False,action="store_true")
group.add_option("--summary",dest="summary",help="Show summary for reads writen to separate files",default=False,action="store_true")
group.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")
parser.add_option_group(group)

(options, args) = parser.parse_args()

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

if ((options.start != None) and (options.start > 0)):
  options.start-=1
else:
  options.start=None

ireadlength1 = options.ireadlength1
ireadlength2 = options.ireadlength2

options.qualityN = options.qualityN.upper()
options.qualityN_2nd = options.qualityN_2nd.upper()

#sys.stderr.write("%d %d %d %s %s\n"%(ireadlength1,ireadlength2,options.start,isdoubleIndex,options.index))

## CREATE OUTPUT FILE(S)/STREAM
fileflags = 'wb'
if options.SAM: fileflags = 'w'
SAMheader = { 'HD': {'VN': '1.4','SO':'queryname'}, 'SQ': [{'LN': 0, 'SN': '*'}] }
if options.readgroup != "":
  SAMheader['RG'] = [{'ID': options.readgroup, "PL":"Illumina", "LB": options.readgroup, "SM": options.readgroup }]

if options.verbose: sys.stderr.write("Creating output file/stream...\n")
outfiles = {}

if options.pipe:
  if options.verbose: sys.stderr.write("BAM/SAM output on stdout...\n")
  outfiles[None] = [pysam.Samfile( "-", fileflags, header = SAMheader),0,None]
else:
    outfilename = options.outdir+options.outprefix+".bam"
    if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
    outfiles[None] = [pysam.Samfile( outfilename , fileflags, header = SAMheader),0,None]

files = args
if options.pipe: files = [None]

for filename in files:
  if filename == None:
    infile = sys.stdin
  else:
    infile = open(filename)

  for seqid,oindseq1,oindqual1,oindseq2,oindqual2,r1,q1,r2,q2 in read_sequence_file(infile,sec_read_start=options.start):
    #print seqid,oindseq1,oindqual1,oindseq2,oindqual2,r1,q1,r2,q2
    is_qcfail = False
    cind = None
    tags = []

    minqual1 = INF_QUALITY
    if options.quality > 0 and ireadlength1 > 0: minqual1 = get_min_qual(oindqual1,options.qualityN)
    minqual2 = INF_QUALITY
    if options.quality > 0 and ireadlength2 > 0: minqual2 = get_min_qual(oindqual2,options.qualityN_2nd)
    if ireadlength2 == 0: indseq2 = None

    if (ireadlength1 > 0) and (ireadlength2 > 0):
      if (minqual1 < options.quality) or (minqual2 < options.quality): is_qcfail = True
    elif (ireadlength1 > 0):
      if (minqual1 < options.quality): is_qcfail = True
    elif (ireadlength2 > 0):
      if (minqual2 < options.quality): is_qcfail = True

    if r2 != None and len(r2) == 0:
      r2 = None
      q2 = None

    if options.readgroup != "": tags.append(("RG",options.readgroup))

    if (len(oindseq1) != 0):
      tags.append(("XI",oindseq1))
      tags.append(("YI",oindqual1))
    if (len(oindseq2) != 0):
      tags.append(("XJ",oindseq2))
      tags.append(("YJ",oindqual2))

    if is_qcfail: tags.append(("ZQ","I"))
    forward = pysam.AlignedRead()
    forward.qname = seqid
    forward.seq = r1
    if q1 != None: forward.qual = q1
    else: forward.qual = "*"
    forward.is_unmapped = True
    forward.is_qcfail = is_qcfail
    forward.tags = tags
    forward.pos = -1
    forward.mpos = -1
    if r2 != None:
      forward.is_read1 = True
      forward.is_paired = True
      reverse = pysam.AlignedRead()
      reverse.qname = seqid
      reverse.is_read2 = True
      reverse.is_paired = True
      reverse.seq = r2
      if q2 != None: reverse.qual = q2
      else: reverse.qual = "*"
      reverse.pos = -1
      reverse.mpos = -1
      reverse.is_qcfail = is_qcfail
      reverse.is_unmapped = True
      forward.mate_is_unmapped = True
      reverse.mate_is_unmapped = True
      reverse.tags = tags
      outfiles[None][1]+=1
      outfiles[None][0].write(forward)
      outfiles[None][0].write(reverse)
    else:
      outfiles[None][1]+=1
      outfiles[None][0].write(forward)

  # REMOVE EMPTY FILES AND PRINT SUMMARY
  summary = []
  closed = False
  for tag,value in outfiles.iteritems():
    if value[2] != None and not closed: 
      value[0].close()
      close = True
    if value[1] == 0 and value[2] != None:
      try:
        if options.verbose: sys.stderr.write("Removing: %s\n"%(value[2]))
        os.remove(value[2])
      except: sys.stderr.write("Error: Removing output files %s %s\n"(value[2]))
    elif options.summary:
      sys.stderr.write("Summary: Wrote %10d sequences\n"%value[1])

  if filename != None:
    infile.close()
