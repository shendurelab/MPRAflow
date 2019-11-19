#merge tables

import sys
import pandas as pd
import numpy as np
import dask.dataframe as dd

exp=pd.DataFrame(pd.read_csv(sys.argv[1]))
#BCs=pd.DataFrame(pd.read_csv(mutual_BCS.txt))

print(exp)


#HepG2      1  RZ_H1
#out directory
o=sys.argv[2]
outfile=sys.argv[3]
d='/'
c='_counts.tsv'
end='_tmp.csv'
dk_full_df= None 
iter=1
for index, row in exp.iterrows():
    print(row['condition'], row['replicate'],row['dna'],row['rna'],row['name'])
    #dnaloc=o+row['dna']+d+row['dna']+c
    #rnaloc=o+row['rna']+d+row['rna']+c
    df_loc=o+row['name']+end
    print(df_loc)
    

    #DNA 1 (condition A, replicate 1)
    colnames=["Barcode", "DNA "+str(iter)+" (condition "+str(row['condition'])+', replicate '+str(row['replicate'])+")","RNA "+str(iter)+" (condition "+str(row['condition'])+', replicate '+str(row['replicate'])+")"]
    cur=pd.DataFrame(pd.read_csv(df_loc,header='infer'))
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

    iter+=1
   

dk_full_df=dk_full_df[sorted(dk_full_df.columns)]
print(dk_full_df.head())
    
def name(i):
    return str(outfile)

dk_full_df.to_csv('*.csv', index=False,name_function=name)



