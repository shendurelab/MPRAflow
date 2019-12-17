#Author Gracie Gordon 2019

import pickle as p
import pandas as pd
import numpy as np
import sys
from Bio import SeqIO

#starting position of barcode in sequence 
bc_pos=int(sys.argv[1])
#barcode length
bc_len=int(sys.argv[2])
#record suffix length for barcode (e.g. CRS_name:001 would be 4)
bc_name_len=int(sys.argv[3])
#file names
input_fa=sys.argv[4]
#'GSM2221119_liverEnhancer_design.fa'
out_p=sys.argv[5]
#'inoue2017_library_new.p'

bc_list=[]
record_list=[]
out_dict={}

for record in SeqIO.parse(input_fa, "fasta"):
    bc_list.append(str(record.seq[bc_pos:bc_pos+bc_len]))
    record_list.append(record.id[:-bc_name_len])

data=pd.DataFrame({'id': record_list, 'bc':bc_list})

for insert in set(record_list):
	tmp=data.loc[data['id'] == insert]
	b=list(tmp['bc'])
	out_dict[insert]=b

p.dump( out_dict, open( out_p, "wb" ) )
