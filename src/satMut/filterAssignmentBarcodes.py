#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *11.07.2016
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict
import pysam

fRegion = None
rsize = 0
fStart = 1
fEnd = None

def checkMinCoverage(barcode,regions,counts):
  global rsize,fStart,fEnd,fRegion
  global options
  region = [0]*rsize
  conflict = 0
  for (elem,obs) in regions:
    rname,rrange = elem.split(":")
    if rname == fRegion:
      rstart,rend = map(int,rrange.split("-"))
      rstart = rstart-1-fStart
      rend = min(rend-1-fStart,rsize)
      for i in range(rstart,rend):
        region[i]+=obs
    else:
      conflict+=1
  minCoverage = min(region[options.indel:-options.indel])
  if (minCoverage >= options.coverage) and (conflict/(1 if minCoverage == 0 else minCoverage) <= options.contamination):
    sys.stdout.write("%s\t%s\t%d\t%d\n"%(barcode,",".join(map(lambda (x,y):"%dx%s"%(y,x),regions)),counts,minCoverage))
  

parser = OptionParser("%prog [options]")
parser.add_option("-c","--coverage", dest="coverage", help="Minimum coverage (default 3)",default=3,type="int")
parser.add_option("-m","--contamination", dest="contamination", help="Maximum contamination from other contigs (default 0.2)",default=0.2,type="float")
parser.add_option("-i","--indel", dest="indel", help="Maximum deleted bases at the edges (default 3)",default=3,type="int")
parser.add_option("-f","--fasta", dest="reference", help="Fasta index reference genome (default /net/shendure/vol1/home/mkircher/regulatory_tests/sat_mutagenesis/new_assignments/reference/reference.fa)",default="/net/shendure/vol1/home/mkircher/regulatory_tests/sat_mutagenesis/new_assignments/reference/reference.fa")
parser.add_option("-r","--region", dest="region", help="Region to limit (default CHR:START:END, '')",default='')
(options, args) = parser.parse_args()

if options.region != "":
  try:
    fRegion,fStart,fEnd = options.region.split(":")[0],int(options.region.split(":")[1]),int(options.region.split(":")[2])
  except:
    fRegion = None
    fStart = 1
    fEnd = None

genome = pysam.Fastafile(options.reference)
try:
  seq = genome.fetch(fRegion)
  if fEnd == None: fEnd = len(seq)
  elif fEnd < 0: fEnd = len(seq)+fEnd
  rsize = fEnd-fStart
  #print fRegion,fStart,fEnd,fEnd-fStart
except:
  sys.stderr.write("Invalid region defined!\n")
  sys.exit()

fStart-=1
fEnd-=1

lBarcode = None
lregions = []
lcounts = 0
for line in sys.stdin:
  fields = line.rstrip().split("\t")
  if len(fields) < 3: continue
  if lBarcode != fields[0]:
    if (lcounts >= options.coverage) and (lBarcode != None):
      checkMinCoverage(lBarcode,lregions,lcounts)
    lBarcode = fields[0]
    lregions = []
    lcounts = 0
  val = int(fields[2])
  for elem in fields[1].split(","):
    lregions.append((elem,val))
  lcounts+=val

if (lcounts >= options.coverage) and (lBarcode != None):
  checkMinCoverage(lBarcode,lregions,lcounts)