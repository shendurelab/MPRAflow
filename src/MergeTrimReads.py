#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Merge/Adapter trim reads stored in BAM

:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *02.01.2012
:Type: module for C cross-compilation

"""

import sys, os
import math
import random
import string

table = string.maketrans('TGCA','ACGT') # COMPLEMENT DNA


#######################################
# ARBITRARY PARAMETERS AND DEFAULTS
#######################################

offset = 33
cutoff_merge_trim = 0.80
cutoff_merge_seqs_early = 0.95
cutoff_merge_seqs = 0.90
options_min_overlap_seqs = 10 # corresponds n+1 bases overlap
max_prob_same_channel = 0.5
max_prob_N = 0.25
maxadapter_comp = 30
min_length = 5

options_adapter_F=["AGATCGGAAGAGCACACGTCTGAACTCCAGTCACIIIIIIIATCTCGTATGCCGTCTTCTGCTTG"]
options_adapter_S=["AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT"]
options_adapter_chimera="ACACTCTTTCCCTACACGTCTGAACTCCAG,ACACTCTTTCCCACACGTCTGAACTCCAGT,ACACTCTTTCCCTACACACGTCTGAACTCC,CTCTTTCCCTACACGTCTGAACTCCAGTCA,GAAGAGCACACGTCTGAACTCCAGTCACII,GAGCACACGTCTGAACTCCAGTCACIIIII,GATCGGAAGAGCACACGTCTGAACTCCAGT,AGATCGGAAGAGCACACGTCTGAACTCCAG,AGAGCACACGTCTGAACTCCAGTCACIIII,ACACGTCTGAACTCCAGTCACIIIIIIIAT,GTGCACACGTCTGAACTCCAGTCACIIIII,AGCACACGTCTGAACTCCAGTCACIIIIII,CGTATGCCGTCTTCTGCTTGAAAAAAAAAA"
adapter_chimeras = (options_adapter_chimera).split(",")
options_trimCutoff = 1
keys = ('','')
len_key1 = len(keys[0])
len_key2 = len(keys[1])
handle_key = (len_key1 > 0) or (len_key2 > 0)
options_allowMissing = False
options_mergeoverlap = False
options_onlyoverlap = False

#######################################
# HELPER FUNCTIONS
#######################################

# CALCULATE THE EDIT DISTANCE BETWEEN TWO READS
def edits(seq1,seq2):
  lmin = min(len(seq1),len(seq2))
  lmax = max(len(seq1),len(seq2))
  dist = lmax-lmin
  for pos in range(lmin):
    if (seq1[pos] != seq2[pos]): dist+=1
  return dist


##CALCULATE QUALITY SCORE CORRECTED IDENTITY BETWEEN TWO READS 
##FIRST READ DEFINES MAX LENGTH OF COMPARISION
def quality_ident(seq1,qual1,seq2,qual2=[],maxcomp=-1):
  global max_prob_same_channel
  global max_prob_N
  ident = 0.0
  res = 0.0

  if len(qual2) != 0 and len(seq2) < len(seq1):
    hseq1 = seq1
    hqual1 = qual1
    qual1 = qual2
    seq1 = seq2
    qual2 = hqual1
    seq2 = hseq1

  if (maxcomp < 0) or (maxcomp > len(seq1)):
    maxcomp=len(seq1)

  if (len(seq1) >= maxcomp) and (len(seq2) >= maxcomp) and (len(qual1) >= maxcomp):
    if len(qual2) == 0: # IN CASE WE DO NOT HAVE A SECOND QUALITY SCORE
      for pos in range(maxcomp):
        if (seq1[pos] != seq2[pos]) and (seq1[pos] != "I") and (seq2[pos] != "I"):
          ident += 1.0-qual1[pos]
        elif (seq1[pos] == "I") or (seq2[pos] == "I"):
          maxcomp-=1
        else:
          ident += 1.0
      #for pos in range(maxcomp):
        #if (seq1[pos] != seq2[pos]):
          #if (seq1[pos] == "A" and seq2[pos] == "C") or  (seq1[pos] == "C" and seq2[pos] == "A") or (seq1[pos] == "G" and seq2[pos] == "T") or (seq1[pos] == "T" and seq2[pos] == "G"):
            #ident += 1.0-min(max_prob_same_channel,qual1[pos])
          #elif (seq1[pos] == "I") or (seq2[pos] == "I"):
            #ident += 1.0
          #elif (seq1[pos] == "N") or (seq2[pos] == "N"):
            #ident += 1.0-min(max_prob_N,qual1[pos])
          #else:
            #ident += 1.0-qual1[pos]
        #else:
          #ident += 1.0
    else:
      for pos in range(maxcomp):
        if (seq1[pos] != seq2[pos]) and (seq1[pos] != "I") and (seq2[pos] != "I"):
          val = qual1[pos] if qual1[pos] < qual2[pos] else qual2[pos]
          ident += 1.0-val
        elif (seq1[pos] == "I") or (seq2[pos] == "I"):
          maxcomp-=1
        else:
          ident += 1.0
  else:
    sys.stderr.write("Quality_ident: Call with invalid arguments!\n")
    print (len(seq1) >= maxcomp),(len(seq2) >= maxcomp),(len(qual1) >= maxcomp),len(seq2) , maxcomp
    sys.exit()
    return 0.0

  if maxcomp > 0: 
    res = ident/float(maxcomp)
  return res

# CONVERTE QUALITY STRING TO LOG SCORES (INTS)
def convert_quality_logprob(qualstring):
  global offset
  return map(lambda x:ord(x)-offset,qualstring)

# CONVERTE LOG SCORES TO QUALITY STRING
def convert_logprob_quality(probs):
  global offset
  return "".join(map(lambda x:chr(max(33,min(126,x+offset))),probs))

# CONVERTE LOG SCORES TO PROBABLILITIES FOR IDENTITY CALCULATION
def convert_logprob_prob(lprobs):
  res = []
  for prob in lprobs:
    val = 1.0-math.pow(10.0,prob/-10.0)
    if val > max_prob_N: res.append(val)
    else: res.append(max_prob_N)
  return res

# REVERSE COMPLEMENT A DNA SEQUENCE
def revcompl(seq):
  global table
  return "".join(seq.translate(table)[::-1])

# WARNING: NUMERIC INSTABLE CALCULATION WITH REAL NUMBERS
def cons_base_prob(base1,base2,prob1,prob2):
  aprob1 = math.log10(1.0-prob1)-math.log10(3.0)
  lprob1 = math.log10(prob1)
  lprob2 = math.log10(prob2)
  aprob2 = math.log10(1.0-prob2)-math.log10(3.0)
  bases = [base1,base2]
  if base1 == base2: bases = [base1]

  total_prob = 0.0
  for char_elem in 'ACGT':
    help = 0.0
    if base1 == char_elem: help+=lprob1
    else: help+=aprob1
    if base2 == char_elem: help+=lprob2
    else: help+=aprob2
    total_prob+=(10**help)
  total_prob = math.log10(total_prob)
  call_base = 'N'
  call_quality = -1
  for char_call in bases:
    thelp = 0.0
    for char_elem in 'ACGT':
      if char_elem != char_call:
        help = 0.0
        if base1 == char_elem: help+=lprob1
        else: help+=aprob1
        if base2 == char_elem: help+=lprob2
        else: help+=aprob2
        thelp+=(10**help)
    thelp = math.log10(thelp)
    val = int(round(-10.0*(thelp-total_prob)))
    hqual = 60 if val > 60 else val
    #sys.stderr.write("\t".join([base1,base2,prob1,prob2,lprob1,lprob2,aprob1,aprob2,help,total_prob,hqual])+"\n")
    if (hqual > call_quality) or ((hqual == call_quality) and (random.random() >=0.5)):
      call_base = char_call
      call_quality = hqual
  #sys.stderr.write("\t".join([call_base,call_quality])+"\n")
  return call_base,call_quality


# OUT SOURCED TEST FOR MERGING TRIMMED READS, ROBUST TO DIFFERENT LENGTH
def check_merge(read1,qual1,pqual1,read2,qual2,pqual2):
  global cutoff_merge_seqs
  global options_onlyoverlap
  new_seq = ""
  new_qual = []
  lread1 = len(read1)
  lread2 = len(read2)
  if (lread1 > 0) and (lread2 > 0):
    oident = quality_ident(read2,qual2,read1,qual1)
    if oident > cutoff_merge_seqs:
      spos = -1
      new_lseq = list(read1)
      new_qual = pqual1
      for pos in range(lread1,lread2):
        new_lseq.append(read2[pos])
        new_qual.append(pqual2[pos])
      for pos in range(lread2):
        if spos == -1 and read2[pos] != "I" and read1[pos] != "I": spos = pos
        if (new_lseq[pos] == "I") or ((new_lseq[pos] == "N") and (read2[pos] != "N") and (read2[pos] != "I")):
          new_lseq[pos] = read2[pos]
          new_qual[pos] = pqual2[pos]
        elif (pos < lread1) and (read1[pos] != "N") and (read2[pos] != "N") and (read1[pos] != "I") and (read2[pos] != "I"):
          nbase,nqual = cons_base_prob(read1[pos],read2[pos],qual1[pos],qual2[pos])
          new_qual[pos] = nqual
          new_lseq[pos] = nbase
      if spos == -1 : spos = 0
      if options_onlyoverlap:
        new_seq = "".join(new_lseq[spos:lread1+1])
        new_qual = new_qual[spos:lread1+1]
        #if read2.startswith("I") or read1.startswith("I"): 
          #sys.stderr.write("OnlyOverlap mode %d %d!\n"%(spos,lread1+1))
      else:
        #if read2.startswith("I") or read1.startswith("I"): 
          #sys.stderr.write("Not in OnlyOverlap mode %d %d!\n"%(spos,lread1+1))
        new_seq = "".join(new_lseq)
  return new_seq,convert_logprob_quality(new_qual)


#######################################
# 'INTERFACE' FUNCTIONS
#######################################

def set_options(trimcutoff=1,allowMissing=False,mergeoverlap=False,onlyoverlap=False, min_overlap_seqs=10):
  global options_trimCutoff
  global options_allowMissing
  global options_mergeoverlap
  global options_onlyoverlap
  global options_min_overlap_seqs
  options_trimCutoff = trimcutoff
  options_allowMissing = allowMissing
  options_mergeoverlap = mergeoverlap
  options_onlyoverlap = onlyoverlap
  options_min_overlap_seqs = min_overlap_seqs

def set_adapter_sequences(forward='',reverse='',chimera='',max_comp=30):
  global options_adapter_F, options_adapter_S
  global options_adapter_chimera
  global adapter_chimeras
  global maxadapter_comp
  options_adapter_F = (forward.upper()).split(',')
  options_adapter_S = reverse.upper().split(',')
  options_adapter_chimera = chimera.upper()
  adapter_chimeras = (options_adapter_chimera).split(",")
  for adapter_F in options_adapter_F:
    if adapter_F not in adapter_chimeras: adapter_chimeras.append(adapter_F)
  helpremove = []
  for ind in range(len(adapter_chimeras)):
    adapter_chimeras[ind]=adapter_chimeras[ind].strip()
    if (len(adapter_chimeras[ind]) > 0) and (len(adapter_chimeras[ind]) < maxadapter_comp+1):
      #print "Extending Chimera sequence by Is (sequence shorter than maxadapter_comp)"
      while (len(adapter_chimeras[ind]) < (maxadapter_comp+1)):
        adapter_chimeras[ind] += "I"
    elif (len(adapter_chimeras[ind]) == 0):
      helpremove.append(ind)
  for ind in helpremove[::-1]:
    adapter_chimeras.pop(ind)

  for i in range(len(options_adapter_F)):
    if (len(options_adapter_F[i]) < maxadapter_comp):
      #print "Extending first adapter by Is (Adapter shorter than read)"
      while (len(options_adapter_F[i]) < maxadapter_comp):
        options_adapter_F[i] += "I"

  for i in range(len(options_adapter_S)):
    if (len(options_adapter_S[i]) < maxadapter_comp):
      #print "Extending second adapter by Is (Adapter shorter than read)"
      while (len(options_adapter_S[i]) < maxadapter_comp):
        options_adapter_S[i] += "I"

def set_keys(key_text):
  global keys
  global len_key1
  global len_key2
  global handle_key
  key_text = key_text.upper()
  if key_text.count(",") == 1:
    keys = (key_text.split(',')[0],key_text.split(',')[1])
  elif key_text.count(",") == 0:
    keys = (key_text,key_text)
  else:
    sys.stderr.write("Unexpected number of keys specified.\n")
    return False
  len_key1 = len(keys[0])
  len_key2 = len(keys[1])
  handle_key = (len_key1 > 0) or (len_key2 > 0)
  return True


def process_PE(read1,qual1,read2,qual2):
  if len(read1) != len(qual1) or len(read2) != len(qual2): 
    sys.stderr.write("Error: Length of quality score and read strings do not match!\n")
    return '','',''
  
  # CHECK KEY, IF KEY SEQUENCE SPECIFIED AND MATCHES REMOVE KEY FROM BOTH ENDS
  if handle_key and len(read1) > 0:
    if ((read1[:len_key1] == keys[0]) and (read2[:len_key2] == keys[1])) or (options_allowMissing and (edits(read1[:len_key1],keys[0]) == 1) and (read2[:len_key2] == keys[1])) or (options_allowMissing and (read1[:len_key1] == keys[0]) and (edits(read2[:len_key2],keys[1]) == 1)):
      read1 = read1[len_key1:]
      qual1 = qual1[len_key1:]
      read2 = read2[len_key2:]
      qual2 = qual2[len_key2:]
    elif options_allowMissing and (read1[:len_key1-1] == keys[0][1:]) and (read2[:len_key2] == keys[1]):
      read1 = read1[len_key1-1:]
      qual1 = qual1[len_key1-1:]
      read2 = read2[len_key2:]
      qual2 = qual2[len_key2:]
    elif options_allowMissing and (read1[:len_key1] == keys[0]) and (read2[:len_key2-1] == keys[1][1:]):
      read1 = read1[len_key1:]
      qual1 = qual1[len_key1:]
      read2 = read2[len_key2-1:]
      qual2 = qual2[len_key2-1:]
    else:
      return 'K','',''

  # CHECK ADAPTER CHIMERAS
  if len(adapter_chimeras) > 0:
    lqual1 = convert_quality_logprob(qual1)
    pqual1 = convert_logprob_prob(lqual1)
    lread1 = len(read1)
    for chimera in adapter_chimeras:
      if quality_ident(read1,pqual1,chimera,maxcomp=maxadapter_comp) > cutoff_merge_trim:
        read1 = ""
        read2 = ""
        break
      elif (quality_ident(read1,pqual1,chimera[1:],maxcomp=maxadapter_comp) > cutoff_merge_trim): # CONSIDER LOSING FIRST BASE...
        read1 = ""
        read2 = ""
        break
    if read1 == "" or read2 == "":
      return 'D','',''

  lread1 = len(read1)
  lread2 = len(read2)
  ## CONVERTE QUALITY STRINGS TO PROBABILIITIES...
  lqual1 = convert_quality_logprob(qual1)
  pqual1 = convert_logprob_prob(lqual1)
  lqual2 = convert_quality_logprob(qual2)
  pqual2 = convert_logprob_prob(lqual2)

  rread2 = revcompl(read2)
  rqual2 = list(pqual2)
  rqual2.reverse()
  rlqual2 = list(lqual2)
  rlqual2.reverse()
  mlength = min(lread1,lread2)
  clength = max(lread1,lread2)

  ## CHECK FOR MERGING WITH ADAPTERS ...
  have_merged = False
  for start in range(mlength+1)[::-1]: ## COMMON OVERLAP
    cval1 = 0.0
    for i in range(len(options_adapter_F)):
      hcval = quality_ident(read1[start:],pqual1[start:],options_adapter_F[i],maxcomp=maxadapter_comp)
      if hcval > cval1: cval1 = hcval
    cval2 = 0.0
    for i in range(len(options_adapter_S)):
      hcval = quality_ident(read2[start:],pqual2[start:],options_adapter_S[i],maxcomp=maxadapter_comp)
      if hcval > cval2: cval2 = hcval

    new_seq,new_qual = check_merge(read1[:start],pqual1[:start],lqual1[:start],rread2[max(0,lread2-start):],rqual2[max(0,lread2-start):],rlqual2[max(0,lread2-start):])
    if (new_seq != "") and (cval1 > cutoff_merge_trim or cval2 > cutoff_merge_trim):
      have_merged = True
      return '',new_seq,new_qual

  if not have_merged:
    for start in range(mlength+1,clength+1): ## SEQUENCE LEFT FOR ASYMMETRIC READ LENGTH
      cval1 = 0.0
      for i in range(len(options_adapter_F)):
        hcval = quality_ident(read1[start:],pqual1[start:],options_adapter_F[i],maxcomp=maxadapter_comp)
        if hcval > cval1: cval1 = hcval
      cval2 = 0.0
      for i in range(len(options_adapter_S)):
        hcval = quality_ident(read2[start:],pqual2[start:],options_adapter_S[i],maxcomp=maxadapter_comp)
        if hcval > cval2: cval2 = hcval

      help_rread2 = ""
      help_lqual2 = []
      help_qual2 = []
      for ipos in range(max(0,start-lread2)):
        help_rread2 += "I"
        help_qual2.append(0.0)
        help_lqual2.append(0)
      new_seq,new_qual = check_merge(read1[:start],pqual1[:start],lqual1[:start],help_rread2+rread2[max(0,lread2-start):],help_qual2+rqual2[max(0,lread2-start):],help_lqual2+rlqual2[max(0,lread2-start):])
      if (new_seq != "") and ((cval1 > cutoff_merge_trim or cval2 > cutoff_merge_trim) or ((read1[start:]==read2[start:]) and (len(read1[start:]) == 0))):
        have_merged = True
        if len(new_seq) < min_length:
          return 'D','',''
        else:
          return '',new_seq,new_qual
        break

  if not have_merged and options_mergeoverlap:
    # SO FAR NOTHING. INSERT MIGHT BE LONGER THAN READ LENGTH. TRY TO MERGE SEQUENCES WITH OVERLAP
    max_value1,max_pos1 = -1,-1
    for start in range(lread1-options_min_overlap_seqs):
      if lread1-start <= lread2:
        cval = quality_ident(read1[start:],pqual1[start:],rread2,rqual2)
        if cval > cutoff_merge_seqs_early:
          max_value1,max_pos1=cval,start
          break
        elif (max_value1 == -1) or (cval > max_value1):
          max_value1,max_pos1=cval,start

    # ABOVE CUTOFF?
    if (max_value1 > cutoff_merge_seqs):
      if options_onlyoverlap:
        new_lseq = list(read1[max_pos1:mlength])
        new_lqual = lqual1[max_pos1:mlength]
        for pos in range(mlength-max_pos1):
          new_lseq[pos],new_lqual[pos] = cons_base_prob(new_lseq[pos],read1[max_pos1+pos],rqual2[pos],pqual1[max_pos1+pos])
      else:
        new_lseq = list(read1[:max_pos1]+rread2)
        new_lqual = lqual1[:max_pos1]+rlqual2
        for pos in range(max_pos1,mlength):
          #new_lseq[pos],new_lqual[pos] = cons_base_prob(new_lseq[pos],read1[pos],rqual2[pos],pqual1[pos])
          new_lseq[pos],new_lqual[pos] = cons_base_prob(new_lseq[pos],read1[pos],rqual2[pos-max_pos1],pqual1[pos])
      return '',"".join(new_lseq),convert_logprob_quality(new_lqual)
  return '','',''


def overlap_reads(read1,qual1,read2,qual2):
  if len(read1) != len(qual1) or len(read2) != len(qual2): 
    sys.stderr.write("Error: Length of quality score and read strings do not match!\n")
    return '','',''
  
  lread1 = len(read1)
  lread2 = len(read2)
  ## CONVERTE QUALITY STRINGS TO PROBABILIITIES...
  lqual1 = convert_quality_logprob(qual1)
  pqual1 = convert_logprob_prob(lqual1)
  lqual2 = convert_quality_logprob(qual2)
  pqual2 = convert_logprob_prob(lqual2)

  rread2 = revcompl(read2)
  rqual2 = list(pqual2)
  rqual2.reverse()
  rlqual2 = list(lqual2)
  rlqual2.reverse()
  mlength = min(lread1,lread2)

  max_value1,max_pos1 = -1,-1
  for start in range(mlength-options_min_overlap_seqs):
    if lread2-start <= lread1:
      cval = quality_ident(read1,pqual1,rread2[start:],rqual2[start:])
      if cval > cutoff_merge_seqs_early:
        max_value1,max_pos1=cval,start
        break
      elif (max_value1 == -1) or (cval > max_value1):
        max_value1,max_pos1=cval,start

  # ABOVE CUTOFF?
  if (max_value1 > cutoff_merge_seqs):
    new_lseq = list(rread2[max_pos1:mlength])
    new_lqual = rlqual2[max_pos1:mlength]
    for pos in range(mlength-max_pos1):
      new_lseq[pos],new_lqual[pos] = cons_base_prob(new_lseq[pos],read1[pos],rqual2[max_pos1+pos],pqual1[pos])
    return '',"".join(new_lseq),convert_logprob_quality(new_lqual)
  else:
    return '','',''


def consensus_reads(read1,qual1,read2,qual2):
  if len(read1) != len(qual1) or len(read2) != len(qual2) or len(read1) != len(read2): 
    sys.stderr.write("Error: Length of quality score and read strings do not match!\n")
    return '','',''
  
  ## CONVERTE QUALITY STRINGS TO PROBABILIITIES...
  lqual1 = convert_quality_logprob(qual1)
  pqual1 = convert_logprob_prob(lqual1)
  lqual2 = convert_quality_logprob(qual2)
  pqual2 = convert_logprob_prob(lqual2)

  cval = quality_ident(read1,pqual1,read2,pqual2)
  if cval > 0.75:
    new_lseq = list(read2)
    new_lqual = lqual2
    for pos in range(len(read1)):
      new_lseq[pos],new_lqual[pos] = cons_base_prob(new_lseq[pos],read1[pos],pqual2[pos],pqual1[pos])
    return '',"".join(new_lseq),convert_logprob_quality(new_lqual)
  else:
    return '','',''


def process_SR(read1,qual1):
  if len(read1) != len(qual1): 
    sys.stderr.write("Error: Length of quality score and read strings do not match!\n")
    return '','',''

  if handle_key:
    key_handled = False
    if (read1[:len_key1] == keys[0]) or (options_allowMissing and (edits(read1[:len_key1],keys[0]) == 1)):
      read1 = read1[len_key1:]
      qual1 = qual1[len_key1:]
      key_handled = True
    elif options_allowMissing and (read1[:len_key1-1] == keys[0][1:]):
      read1 = read1[len_key1-1:]
      qual1 = qual1[len_key1-1:]
      key_handled = True
    if not key_handled:
      return 'K','',''

  # CONVERTE QUALITY STRING TO PROBABILIITIES...
  lqual1 = convert_quality_logprob(qual1)
  pqual1 = convert_logprob_prob(lqual1)
  lread1 = len(read1)

  # CHECK WHETHER WE (STILL) HAVE A SEQUENCE
  if len(adapter_chimeras) > 0:
    for chimera in adapter_chimeras:
      if (quality_ident(read1,pqual1,chimera,maxcomp=maxadapter_comp) > cutoff_merge_trim):
        read1 = ""
        break
      elif (quality_ident(read1,pqual1,chimera[1:],maxcomp=maxadapter_comp-1) > cutoff_merge_trim): # CONSIDER LOSING FIRST BASE...
        read1 = ""
        break
    if read1 == "": 
      return 'D','',''

  adapter_pos = -2
  # CHECK FOR ADAPTER OF READ11...
  max_value,max_pos = -1,-1
  for start in range(lread1):
    cval = 0.0
    for i in range(len(options_adapter_F)):
      hcval = quality_ident(read1[start:],pqual1[start:],options_adapter_F[i],maxcomp=maxadapter_comp)
      if hcval > cval: cval = hcval

    if (cval > cutoff_merge_seqs_early):
      max_value,max_pos=cval,start
      break
    elif ((cval > cutoff_merge_trim) and (cval > max_value) and (lread1-start > min_length)) or ((lread1-start <= min_length) and (max_value == -1)):
      max_value,max_pos=cval,start

  if (max_value > cutoff_merge_trim) and ((lread1-max_pos) >= options_trimCutoff): # ADAPTER FOUND
    adapter_pos = max_pos
    read1 = read1[:adapter_pos]
    if (len(read1) < min_length):
      return 'D','',''
    else:
      return '',read1,qual1[:adapter_pos]
  return '','',''

if __name__ == '__main__':
  set_adapter_sequences(",".join(options_adapter_F),",".join(options_adapter_S),options_adapter_chimera,maxadapter_comp)
  set_options(options_trimCutoff,options_allowMissing,True,False)
  set_keys(',')
  print "Read cosensus",consensus_reads('ACGTACGTACGT','000000000000','ACGTACGTCCGT','000000000000')
  print "SR",process_SR('ACGTACGTACGT','000000000000')
  print "PE fail",process_PE('ACGTACGTACGT','000000000000','TGCATGCATGCA','000000000000')
  print "PE no adapter",process_PE('ACGTACGTACGT','000000000000','ACGTACGTACGT','000000000000')
  print "PE adapter",process_PE('ATAAACATATGGCAAACATGGTTCTAGATCGGAAGAGCACACGT','00000000000000000000000000000000000000000000','AGAACCATGTTTGCCATATGTTTATAGATCGGAAGAGCGTCGTG','00000000000000000000000000000000000000000000')
  print "PE overlap consensus", process_PE('ATAAACATATGGCAAACATGGTTCTAAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAATCACTAACACATTCATTTTCCAGACACCCTTCACACTACCGTCGGATCGTGCGTGTACCTCTGAATCTCGTATGCCGTCTTCT',  '55???BBBDDDDDDDDFFFFF>EEHBF?EFFFFG@FHHHHHHHHHHHH@GHHHHHHHHHHHHHHF?FFHEGDG=EHG?FFGCFGBGHHHFGCFCFHHHHHBD.?C--@=<+<CECE=;46=@@+=,=,5=5,AB,,55A;*5,C=:;18?##',  'TCTGGAAAATGAATGTGTTAGTGATTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTAGAACCATGTTTGCCATATGTTTATCTTCAGCTTCCCGATTACGGATCTCGTATGTGTAGATCTCGGTGGTCGCCG',  '+?BDDDDDFFFCCBFEFHFFF>EF@GGHHHHFHFHHHFFHHHHHHFFHHHHHGHHHHHHHHHHHHHHH-A-CFF,C//AECFD?EEFDFGH/CAFHF/ACFHDFFEEEHDH<DDD)7?;?+@7@-7,66B?D?,+@@@######')
  print "Overlap read consensus", overlap_reads('ATAAACATATGGCAAACATGGTTCTAAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAATCACTAACACATTCATTTTCCAGACACCCTTCACACTACCGTCGGATCGTGCGTGTACCTCTGAATCTCGTATGCCGTCTTCT',  '55???BBBDDDDDDDDFFFFF>EEHBF?EFFFFG@FHHHHHHHHHHHH@GHHHHHHHHHHHHHHF?FFHEGDG=EHG?FFGCFGBGHHHFGCFCFHHHHHBD.?C--@=<+<CECE=;46=@@+=,=,5=5,AB,,55A;*5,C=:;18?##',  'TCTGGAAAATGAATGTGTTAGTGATTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTAGAACCATGTTTGCCATATGTTTATCTTCAGCTTCCCGATTACGGATCTCGTATGTGTAGATCTCGGTGGTCGCCG',  '+?BDDDDDFFFCCBFEFHFFF>EF@GGHHHHFHFHHHFFHHHHHHFFHHHHHGHHHHHHHHHHHHHHH-A-CFF,C//AECFD?EEFDFGH/CAFHF/ACFHDFFEEEHDH<DDD)7?;?+@7@-7,66B?D?,+@@@######')
  print
  set_adapter_sequences("ATTGCGTGAACCGAAGCTCATCAAGATCTG","GTTGATCCGGTCCTAGGCAGTGAAGATCTC","",maxadapter_comp)
  print 
  print process_SR("TTTACGGCTCATTGCGTGAACCGAAGCTCATCAAGATCTGGCCTCGGCGGCCAAGCTTAGTCGCCTATACGGTGATGGGTGCATNNTTCATNNNNATNNNCTCCCCGTTTAATCCCATATCTCGTATGCCGTCTTCTGCTTGAAAAAAAA","=======++5<5@9@@>CCEEE>CCCCCEEFFGFFFFFFDFFEFEEEECCDDCE+ACEDD5CEDEEECDE@D=9E+CDD@;@DE##113D9####00###21*28@@2<E?E=E=EEE</;(6?EEE<E<<E;6<EE#############")
  set_adapter_sequences("GTTGATCCGGTCCTAGGCAGTGAAGATCTC","ATTGCGTGAACCGAAGCTCATCAAGATCTG","",maxadapter_comp)
  print process_SR("AGCCGTAAAAGTTGATCCGGTCCTAGGCAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAACAAGGAAAATAAAAAAAAAAACCCAGCCAACTAACAAAAAACAAACAAAACAAAAACGAAACCAAAACCAAAAAATAC","????????BDB<9<BBFFBCCC7ECCACEECFGFCDBCBFFFEEDC?>7C7>EEFGH=CDEGGHHHDE##################################################################################")
  print 
  print process_PE("TTTACGGCTCATTGCGTGAACCGAAGCTCATCAAGATCTGGCCTCGGCGGCCAAGCTTAGTCGCCTATACGGTGATGGGTGCATNNTTCATNNNNATNNNCTCCCCGTTTAATCCCATATCTCGTATGCCGTCTTCTGCTTGAAAAAAAA","=======++5<5@9@@>CCEEE>CCCCCEEFFGFFFFFFDFFEFEEEECCDDCE+ACEDD5CEDEEECDE@D=9E+CDD@;@DE##113D9####00###21*28@@2<E?E=E=EEE</;(6?EEE<E<<E;6<EE#############","AGCCGTAAAAGTTGATCCGGTCCTAGGCAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAACAAGGAAAATAAAAAAAAAAACCCAGCCAACTAACAAAAAACAAACAAAACAAAAACGAAACCAAAACCAAAAAATAC","????????BDB<9<BBFFBCCC7ECCACEECFGFCDBCBFFFEEDC?>7C7>EEFGH=CDEGGHHHDE##################################################################################")