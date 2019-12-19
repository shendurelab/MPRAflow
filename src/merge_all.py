#Author Gracie Gordon 2019
#merge tables

import sys
import pandas as pd
import numpy as np
import dask.dataframe as dd


cond=sys.argv[1]
outfile=sys.argv[2]

dk_full_df= None

replicates=int((len(sys.argv)-3)/2)
for i in range(3,(len(sys.argv)-replicates-1)):
    file=sys.argv[i]
    rep=sys.argv[i+replicates]

    #DNA 1 (condition A, replicate 1)
    colnames=["Barcode", "DNA %s (condition %s, replicate %s)" % (rep,cond,rep),
                         "RNA %s (condition %s, replicate %s)" % (rep,cond,rep)]
    cur=pd.DataFrame(pd.read_csv(file,header='infer'))
    print(cur.head())
    cur.columns=colnames
    cur_dk=dd.from_pandas(cur,npartitions=1)
    print(cur.head())

    if (dk_full_df is not None):

        tmp=dd.merge(dk_full_df,cur_dk, on=["Barcode"],how='outer')
        dk_full_df=tmp
    else:
        dk_full_df=cur_dk


    print(dk_full_df.head())


dk_full_df=dk_full_df[sorted(dk_full_df.columns)]
print(dk_full_df.head())


dk_full_df.to_csv([outfile], index=False)
