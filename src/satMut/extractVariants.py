#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *11.09.2015
"""

import sys, os
from optparse import OptionParser
import gzip

parser = OptionParser("%prog [options]")
parser.add_option("-r","--region", dest="region", help="Genomic region to keep (def '')",default="")
(options, args) = parser.parse_args()

chrom,start,end = None,None,None
if len(options.region) != 0:
  try:
    chrom,region = options.region.split(":")
    start,end = map(int,region.split("-"))
  except:
    sys.stderr.write("Parsing region (%s) failed\n"%(options.region))
    chrom,start,end = None,None,None

# Usage:
# module load bcftools/latest
# /net/shendure/vol1/home/mkircher/bin/samtools-1.2/samtools mpileup -q 10 -A -m 3 -R -f reference/reference.fa -u test.bam | bcftools call -c -f GQ | /net/shendure/vol1/home/mkircher/bin/regulatory_assays/extractVariants.py

for line in sys.stdin:
  if line.startswith("#"): continue
  else:
    fields = line.rstrip().split("\t")
    if chrom != None and (fields[0] != chrom or (int(fields[1]) < start) or (int(fields[1]) > end)): continue
    call = dict(zip(fields[8].split(":"),fields[-1].split(":")))
    info = dict(map(lambda x: (x,True) if len(x.split("=")) != 2 else tuple(x.split("=")),fields[7].split(";")))
    if ('GT' in call) and call['GT'] != "0/0" and call['GT'] != "0|0":
      if "DP4" in info:
        counts = map(lambda x:int(x),info["DP4"].split(","))
        if (counts[0]+counts[1]) < (counts[2]+counts[3]):
          print "%s:%s:%s>%s"%(fields[0],fields[1],fields[3],fields[4].split(",")[0])
      else:
        if ('GQ' in call):
          if int(call['GQ']) >= 30:
            allele = fields[4].split(",")[int(filter(lambda x: x!="0",call['GT'].replace("|","/").split("/"))[0])-1]
            print "%s:%s:%s>%s"%(fields[0],fields[1],fields[3],allele)
        else:
          print "%s:%s:%s>%s"%(fields[0],fields[1],fields[3],fields[4].split(",")[0])
