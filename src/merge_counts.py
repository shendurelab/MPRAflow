#Author: Gracie Gordon 2019
#merge tables

import sys
import pandas as pd
import numpy as np
import dask.dataframe as dd

if (sys.argv[1]=="DNA"):
    dnaloc=sys.argv[2]
    rnaloc=sys.argv[3]
else:
    dnaloc=sys.argv[3]
    rnaloc=sys.argv[2]

outfile=sys.argv[4]
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

out.to_csv([outfile], index=False)
