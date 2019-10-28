#!/usr/bin/env python
# -*- coding: ASCII -*-

"""
Convert fast(a|q) reads into BAM

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *19.09.2011
:Type: tool
:Input: fasta, fastq
:Output: BAM

"""

import sys
import os
sys.path.append("/home/sbsuser/lib/python")
sys.path.append("/mnt/solexa/bin/")

import pysam
from library import read_fastq
from optparse import OptionParser,OptionGroup

parser = OptionParser("%prog [options] seq_files")
parser.add_option("-p","--PIPE",dest="pipe",help="Read from and write to PIPE",default=False,action="store_true")
parser.add_option("-s","--start",dest="start",help="First base of the second read (default None)",type='int')
parser.add_option("--qualityoffset",dest="qualityoffset",help="Offset of quality scores in FastQ input file (default 33)",type="int",default=33)
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'outbam').",default="outbam")
parser.add_option("", "--SAM", dest="SAM", help="Output SAM not BAM.",default=False,action="store_true")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")
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
  if options.verbose: sys.stderr.write("Second read starts at cycle: %d\n"%(options.start+1))
else:
  options.start=None

## CREATE OUTPUT FILE(S)/STREAM
fileflags = 'wb'
if options.SAM: fileflags = 'w'
SAMheader = { 'HD': {'VN': '1.4','SO':'queryname'}, 'SQ': [{'LN': 0, 'SN': '*'}] }

if options.verbose: sys.stderr.write("Creating output file/stream...\n")
outfile = None
if options.pipe:
  outfile = [pysam.Samfile( "-", fileflags, header = SAMheader),0]
  if options.verbose: sys.stderr.write("BAM/SAM output on stdout...\n")
else:
  outfilename = options.outdir+options.outprefix+".bam"
  if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
  outfile = [pysam.Samfile( outfilename , fileflags, header = SAMheader),0]

files = args
if options.pipe: files = [None]

for filename in files:
  if filename == None:
    infile = sys.stdin
  else:
    infile = open(filename)
    # ACTUAL INDEX IDENTIFICATION AND READ SORTING
  for seqid, seq, qual in read_fastq(infile):
    seqid = seqid.split()[0]
    seq2,qual2 = None,None
    if qual != None and options.qualityoffset != 33: qual = "".join(map(lambda x:chr(ord(x)-options.qualityoffset+33),qual))
    if options.start != None:
      seq2 = seq[options.start:]
      seq = seq[:options.start]
      if qual != None:
        qual2 = qual[options.start:]
        qual = qual[:options.start]

    forward = pysam.AlignedRead()
    forward.qname = seqid
    forward.seq = seq
    if qual != None: forward.qual = qual
    else: forward.qual = "*"
    forward.is_unmapped = True
    forward.pos = -1
    forward.mpos = -1
    if seq2 != None:
      forward.is_read1 = True
      forward.is_paired = True
      reverse = pysam.AlignedRead()
      reverse.qname = seqid
      reverse.is_read2 = True
      reverse.is_paired = True
      reverse.seq = seq2
      if qual2 != None: reverse.qual = qual2
      else: reverse.qual = "*"
      reverse.pos = -1
      reverse.mpos = -1
      reverse.is_unmapped = True
      forward.mate_is_unmapped = True
      reverse.mate_is_unmapped = True
      outfile[1]+=1
      outfile[0].write(forward)
      outfile[0].write(reverse)
    else:
      outfile[1]+=1
      outfile[0].write(forward)

  # CREATE SUMMARY AND CLEAN UP NOT USED FILES
  if options.verbose: sys.stderr.write("Converted %d sequences.\n"%outfile[1])

  if filename != None:
    infile.close()
