#merge tables

import sys
import pandas as pd
import numpy as np
import dask.dataframe as dd

cond=sys.argv[1]
rep=sys.argv[2]
#BCs=pd.DataFrame(pd.read_csv(mutual_BCS.txt))


#HepG2      1  RZ_H1
o=sys.argv[2]
d='/'
c='_counts.tsv'

name=[cond,rep].join("_")
name_dna=[name,"DNA"].join("_")
name_rna=[name,"RNA"].join("_")
dnaloc=o+name_dna+d+name_dna+c
rnaloc=o+name_rna+d+name_rna+c
print(dnaloc)
print(rnaloc)

dna=pd.DataFrame(pd.read_csv(dnaloc,delim_whitespace=True,names=['dna_count','BC']))
rna=pd.DataFrame(pd.read_csv(rnaloc,delim_whitespace=True,names=['rna_count','BC']))

print(dna.head())
print(rna.head())

dk_dna=dd.from_pandas(dna,npartitions=1)
print(dk_dna.head())
dk_rna=dd.from_pandas(rna,npartitions=1)

out=dd.merge(dk_dna,dk_rna, on=['BC'],how='outer')

out=out[sorted(out.columns)]
print(out.head())

def name(i):
    return (label_name+'_tmp')
out.to_csv('*.csv', index=False,name_function=name)
