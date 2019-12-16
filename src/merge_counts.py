#Author: Gracie Gordon 2019
#merge tables

import sys
import pandas as pd
import numpy as np
import dask.dataframe as dd

exp=pd.DataFrame(pd.read_csv(sys.argv[1]))
#BCs=pd.DataFrame(pd.read_csv(mutual_BCS.txt))

print(exp)


#HepG2      1  RZ_H1
o=sys.argv[2]
d='/'
c='_counts.tsv'
for index, row in exp.iterrows():
    print(row['condition'], row['replicate'],row['dna'],row['rna'],row['name'])
    dnaloc=o+row['dna']+d+row['dna']+c
    rnaloc=o+row['rna']+d+row['rna']+c
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
        return (str(row['name'])+'_tmp')
        #return (str(row['name'])
    out.to_csv('*.csv', index=False,name_function=name)



