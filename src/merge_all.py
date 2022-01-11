#Author Gracie Gordon 2019
#merge tables

import sys
import pandas as pd
import numpy as np
import dask.dataframe as dd

import click

# options
@click.command()
@click.option('--condition',
              required=True, 
              type=str,
              help='Name of the condition.')
@click.option('--counts',
              'counts',
              required=True,
              nargs=2,
              multiple=True,
              type=(str, click.Path(exists=True, readable=True)),
              help='Replicate name and Count file. Can be used multiple times')
@click.option('--output',
              'output_file',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')
def cli(condition, counts, output_file):

    dk_full_df= None
    
    for replicate_count in counts:
        
        rep=replicate_count[0]
        file=replicate_count[1]
        
        #DNA 1 (condition A, replicate 1)
        colnames=["Barcode", "DNA %s (condition %s, replicate %s)" % (rep,condition,rep),
                             "RNA %s (condition %s, replicate %s)" % (rep,condition,rep)]
        cur=pd.DataFrame(pd.read_csv(file, sep='\t', header=None))
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


    dk_full_df.compute().to_csv(output_file, index=False)
    
    
if __name__ == '__main__':
    cli()