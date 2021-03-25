#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *10.07.2016
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict
import pysam

def count_softclip_cigar(cigarlist):
  #M     BAM_CMATCH      0
  #I     BAM_CINS        1
  #D     BAM_CDEL        2
  #N     BAM_CREF_SKIP   3
  #S     BAM_CSOFT_CLIP  4
  #H     BAM_CHARD_CLIP  5
  #P     BAM_CPAD        6
  #=     BAM_CEQUAL      7
  #X     BAM_CDIFF       8
  alnLength = 0
  clipLength = 0
  for (operation,length) in cigarlist:
    if operation in [0,2,8,7]: alnLength+=length
    elif operation == 4 or operation == 5: clipLength+=length
  return clipLength,alnLength

tagName = "XI"
def getTagSequence(tags):
  global tagName
  global options
  tag = None
  if options.fixTags:
    ntags = []
    for tag in tags:
      if type(tag[1]) != type(""):
        ntags.append(tag)
      elif " " not in tag[1]:
        ntags.append(tag)
      else:
        fields = tag[1].split()
        ntags.append((tag[0],fields[0]))
        for field in fields[1:]:
          name = field[:2]
          ctype = field[3]
          value = field[5:]
          tags.append((name,value))
          
  for key,value in tags:
    if key == tagName: 
      tag = value[:options.first]
      if options.check != "" and value[options.first:(options.first+len(options.check))] != options.check: 
        return None
      #sys.stderr.write("%s %s\n"%(tag,value[options.first:(options.first+len(options.check))]))
      break
  return tag

parser = OptionParser("%prog [options]")
parser.add_option("--BAMfield", dest="BAMfield", help="Expect barcode in BAM XI field rather than sequence name",default=False,action="store_true")
parser.add_option("-f","--first", dest="first", help="Use first N bases of barcode (default 20)",default=20,type="int")
parser.add_option("-c","--check", dest="check", help="Check that bases following the barcode match this sequence (default '')",default="")
parser.add_option("-s","--switch", dest="switch", help="Tag is stored in XJ rather than XI",default=False,action="store_true")
parser.add_option("--fixTags", dest="fixTags", help="Fix SAM tag field on-the-fly (default off)",default=False,action="store_true")
(options, args) = parser.parse_args()

if options.switch: tagName = "XJ"

for filename in args:
  if os.path.exists(filename):
    infile = pysam.Samfile(filename, "rb")
    for read in infile:
      if read.is_secondary: continue
      
      if read.is_paired:
        if read.is_read1 and (not read.is_unmapped) and (read.mapq > 0) and (not read.mate_is_unmapped) and (read.tid == read.rnext) and (abs(read.tlen) < 10000): 
          # Extract barcode sequence
          tag = None
          if options.BAMfield:
            tag = getTagSequence(read.tags)
            if tag == None: continue
          else:
            tag = read.qname.split("#")[-1][:options.first]  
          # Extract alignment information
          readID,chrom,pos,cigar,isize = read.qname,infile.getrname(read.tid),min(read.pnext,read.pos),read.cigar,abs(read.tlen)
          clipLength,alnLength = count_softclip_cigar(cigar)
          if (clipLength > 5) or (isize == 0): continue
          if read.pnext > read.pos+alnLength:
            sys.stdout.write("%s\t%s:%d-%d,%s:%d-%d\n"%(tag,chrom,read.pos,read.pos+alnLength,chrom,read.pnext,pos+isize))
          else:
            sys.stdout.write("%s\t%s:%d-%d\n"%(tag,chrom,pos,pos+isize))
      else:
        if (not read.is_unmapped) and (read.mapq > 0):
          tag = None
          if options.BAMfield:
            tag = getTagSequence(read.tags)
            if tag == None: continue
          else:
            tag = read.qname.split("#")[-1][:options.first]
          readID,chrom,pos,cigar,isize = read.qname,infile.getrname(read.tid),read.pos,read.cigar,len(read.seq)
          clipLength,alnLength = count_softclip_cigar(cigar)
          if (clipLength > 5) or (isize == 0): continue
          sys.stdout.write("%s\t%s:%d-%d\n"%(tag,chrom,pos,pos+alnLength))
