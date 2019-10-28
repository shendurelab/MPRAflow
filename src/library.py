import sys
import os
import math
import subprocess

"""
Just some frequently used functions stripped from the scripts

:Author: Michael Siebauer
:Contact: michael_siebauer@eva.mpg.de
:Date: 02.07.2010
:Type: library

"""

def is_fastq(filehandle):
    line = filehandle.readline()
    filehandle.seek(0)
    return line.startswith("@")


def get_filter_IDS(options):
  if (options.filter != None) and os.path.exists(options.filter):
      filterIDs = set()
      infile = open(options.filter)
      for line in infile:
        filterIDs.add(line.strip())
      infile.close()
      if len(filterIDs) > 0:
            sys.stderr.write("Read %i IDs for filter\n"%(len(filterIDs)))
            return filterIDs
  return None

def move_partially(sourcef,targetf, is_mock=False):
  print "Starting step-by-step moving because "+targetf+" already exists."
  folders_to_iterate = [(sourcef,targetf)]
  while len(folders_to_iterate) > 0:
    csource,ctarget = folders_to_iterate.pop()
    files_to_move = []
    folders_to_move = []
    for elem in os.listdir(csource):
      if os.path.exists(csource+elem):
        files_to_move.append(csource+elem)
      elif os.path.isdir(csource+elem) and not(os.path.isdir(ctarget+elem)):
        folders_to_move.append(csource+elem)
      elif os.path.isdir(csource+elem) and os.path.isdir(ctarget+elem):
        folders_to_iterate.append((csource+elem+"/",ctarget+elem+"/"))
    if len(files_to_move)>1:
      print "Moving",len(files_to_move),"files to",ctarget
      if is_mock:
        print "mv "+" ".join(files_to_move)+" "+ctarget
      else:
        move = subprocess.Popen("mv "+" ".join(files_to_move)+" "+ctarget,shell=True)
        move.wait()
    if len(folders_to_move)>1:
      print "Moving",len(folders_to_move),"folders to",ctarget
      if is_mock:
        print "mv "+" ".join(folders_to_move)+" "+"/".join(ctarget.split("/")[:-1])+"/"
      else:
        move = subprocess.Popen("mv "+" ".join(folders_to_move)+" "+"/".join(ctarget.split("/")[:-1])+"/",shell=True)
        move.wait()
    if len(os.listdir(csource)) == 0 and not(is_mock):
      move = subprocess.Popen("rmdir "+csource,shell=True)
      move.wait()



def make_output_filename(inputname, options, suffix):
    if suffix:
        if "." in inputname:
            inputname = ".".join(inputname.split(".")[:-1])+suffix
        else:
            inputname = inputname + suffix

    if ("outdir" in vars(options) and  options.outdir and os.path.isdir(options.outdir)):
          inputname = inputname.split("/")[-1]
          inputname = options.outdir.rstrip("/")+"/"+ inputname
    return inputname
        

def read_fastq(filehandle):
    """ Reads fastq and fasta entries from filehandle
        supports multiline fasta and fastq
    """
    count = 0
    seq = ""            
    id = ""             
    qual = None
    is_fasta = False
    for line in filehandle:
        line = line.rstrip("\n\r")
        if ( (count == 0) and (line.startswith("@") or line.startswith(">"))): # Read identifier
            id = line[1:]
            count+=1
            is_fasta = line.startswith(">")

        elif count == 1:        # read sequence
            seq = line
            count+=1
        
        elif count == 2:  # multiple case: a) quality identifier (fastq), b) more sequence (fastq,fasta), c) next sequence identifier (fasta)
            if line.startswith("+"):  # case a)
                id2 = line[1:]
                if ( (len(id2) > 0) and (id2 != id)): # optional sanity check
                    sys.stderr.write("[NOTE] sequences identifier does not match quality identifier: " + id + " " + id2 + "\n")
                count+=1

            elif is_fasta and ( line.startswith(">") or line.startswith("@")) : # case c)
                yield id, seq, None
                id, seq, qual = None, None, None
                count = 1
                id = line[1:]
                is_fasta = line.startswith(">")

            else: # case b)
                seq = seq + line
        
        elif count == 3:
            if qual == None:
                qual = line
            else:
                qual = qual + line
                
            if (len(qual) > len(seq)): # Another sanity check
                sys.stderr.write("[NOTE] sequence and quality line differ in length\n")
            if (len(qual) >= len(seq)):
                count = 0
                yield id, seq, qual
                id, seq, qual = None, None, None
        else:
          sys.stderr.write("Unexpected line:" + str(line.strip()) + "\n")
          count = 0

    if id and seq:
        yield id, seq, qual

    raise StopIteration

def write_fastq(filehandle, id, seq, qual, compress=False):
    if (not id.startswith("@")):
        filehandle.write("@")
    filehandle.write(id.rstrip("\n")+"\n")
    filehandle.write(seq.rstrip("\n") + "\n")
    filehandle.write("+")
    if (not compress):
        if id.startswith("@"):
            filehandle.write(id[1:].rstrip("\n"))
        else:
            filehandle.write(id.rstrip("\n"))
    filehandle.write("\n")
    filehandle.write(qual.rstrip("\n")+"\n")

def write_fasta(filehandle, id, seq):
    if (not id.startswith(">")):
        filehandle.write(">")
    filehandle.write(id.rstrip("\n")+"\n")
    max_size = 80
    seq = "".join(seq.split('\n'))
    pos = 0
    while pos < len(seq):
      filehandle.write(seq[pos:pos+max_size]+"\n")
      pos += max_size

def write_out(filehandle, id, seq, qual):
    if qual:
        write_fastq(filehandle, id, seq, qual)
    else:
        write_fasta(filehandle, id, seq)



def parse_range_string(options, silent=False, write_to=sys.stdout):
    start = 0
    end = None
    try:
      fields = options.range.strip().split("-")
      start = max(int(fields[0])-1,0)
      
      if fields[1].upper() != "MAX":
        end = int(fields[1])
        if end < start: end = None
      if (not silent):
          write_to.write("Set range parameter to:\n")
          if end <> None:
            write_to.write("Start: " + str(start+1) + " End: " + str(end) + "\n")
          else:
            write_to.write("Start: " + str(start+1) + " End: MAX\n")
      return (start, end)
    except:
      if (not silent):
        write_to.write("Couldn't evaluate range parameter, falling back to default: 1-MAX\n")
      start = 0
      end = None
      return (start, end)

def parse_rangestr(rangestr):
  res = []
  fields = rangestr.split(',')
  for elem in fields:
    if "-" in elem:
      se = elem.split('-')
      if len(se) == 2:
        try:
          start = int(se[0])
          end = int(se[1])
          res.extend(range(start,end+1))
        except: return None
      else: return None
    else:
      try:
        pos = int(elem)
        res.append(pos)
      except: return None
  if len(res) == 0: return None
  else:
    res = list(set(res))
    res.sort()
    return res


def is_complex_comp(seq,cutoff):
  counts = []
  for x in 'ACGT':
    counts.append(seq.count(x))
  total = float(sum(counts))
  return (max(counts) <= cutoff*total)

def is_complex_entropy(seq,cutoff):
  return (entropy(seq) >= cutoff)

def entropy(seq):
  counts = []
  for x in 'ACGT':
    counts.append(seq.count(x))
  total = float(sum(counts))
  entropy = 0
  for elem in counts:
    if (total > 0) and (elem/total <> 0):
      entropy -= elem/total*math.log(elem/total,2)
  return entropy
