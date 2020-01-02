#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: martin.kircher@bihealth.de
:Date: *26.02.2018
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict
import gzip

parser = OptionParser("%prog [options]")
#parser.add_option("-m","--minimum", dest="minimum", help="Minimum barcode count considered (def 1.0)",default=1.0,type='float')
#parser.add_option("-f","--fraction", dest="fraction", help="Minimum fraction explained by insert (def 0.5)",default=0.5,type='float')
(options, args) = parser.parse_args()

for filename in args:
  if os.path.exists(filename):
    infile = gzip.open(filename)
    header = []
    totals = defaultdict(int)
    DNA = defaultdict(int)
    RNA = defaultdict(int)
    for line in infile:
      if line.startswith("#"):
        header = line[1:].split()[1:]
      else:
        fields = map(int,line.split()[1:])
        for enum,(cid,count) in enumerate(zip(header,fields)):
          totals[cid]+=count
          if enum > 1 and count > 0:
            DNA[cid]+=fields[0]
            RNA[cid]+=fields[1]
    infile.close()
    
    outfile = open(filename+".stats",'w')
    for key,value in sorted(totals.iteritems()):
      outfile.write("%s\t%d\t%d\t%d\n"%(key,value,DNA[key],RNA[key]))
    outfile.close()
