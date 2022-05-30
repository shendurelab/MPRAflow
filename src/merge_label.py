#Author: Gracie Gordon 2019
#merge tables

import sys
import pandas as pd
import numpy as np
import dask.dataframe as dd

import math
import pickle

from Bio import SeqIO

import click

# options
@click.command()
@click.option('--control-type',
              'control_type',
              required=True, 
              type=click.Choice(['DNA', 'RNA'], case_sensitive=False),
              help='control is DNA or RNA.')
@click.option('--control',
              'control_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='Control file counts')
@click.option('--experiment',
              'experiment_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='experiment file counts')
@click.option('--coord',
              'coord_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='experiment file counts')
@click.option('--design',
              'design_file',
              required=True,
              type=click.Path(exists=True, readable=True),
              help='experiment file counts')
@click.option('--merge-intersect',
              'merge_inter',
              required=True, 
              type=click.Choice(['FALSE', 'TRUE'], case_sensitive=False),
              help='Merge intersections')
@click.option('--bc-length',
              'bc_length',
              required=True, 
              type=int,
              help='Length of BC')
@click.option('--output',
              'outfile',
              required=True,
              type=click.Path(writable=True),
              help='Output file.')
def cli(control_type, control_file, experiment_file, coord_file, design_file, merge_inter, bc_length, outfile):

    #read in files
    if (control_type=="DNA"):
        dnaloc=control_file
        rnaloc=experiment_file
    else:
        dnaloc=experiment_file
        rnaloc=control_file


    #process fastq
    design=open(design_file)
    fasta_dict = {rec.id : rec.seq for rec in SeqIO.parse(design, "fasta")}

    #counts=pd.read_csv(data,header='infer',sep=',')

    assoc=pickle.load( open(coord_file,'rb'))

    BC_key = {}
    for k,v in assoc.items():
        for x in v:
            BC_key.setdefault(x,k)


    #dna rna pair, merge, and label
    #get files

    print(dnaloc)
    print(rnaloc)

    #get count dfs
    dna=pd.DataFrame(pd.read_csv(dnaloc,delim_whitespace=True,names=['dna_count','Barcode']))
    rna=pd.DataFrame(pd.read_csv(rnaloc,delim_whitespace=True,names=['rna_count','Barcode']))
    print('original')
    print(dna.head())
    print(rna.head())

    #convert to dask df for merging
    dk_dna=dd.from_pandas(dna,npartitions=1)
    print(dk_dna.head())
    dk_rna=dd.from_pandas(rna,npartitions=1)

    if(merge_inter.upper()=="TRUE"):
        print('merge')
        out=dd.merge(dk_dna,dk_rna, on=['Barcode'])
        print(out.head())
    else:
        out=dd.merge(dk_dna,dk_rna, on=['Barcode'],how='outer')
        print('merge')
        print(out.head())
        #get rid of NAs in dna and set rna NAs to zero
        #out.fillna(0)
        #out.isna()


    out=out[sorted(out.columns)]
    out=out.fillna(0)
    print('sorted')
    print(out.head())


    #back to pandas for labeling
    counts=out.compute()
    print('back2 pandas')
    print(counts.head())
    #fill in labels from dictionary
    label=[]
    for i in counts.Barcode:
        try:
                label.append(BC_key[i])
                #print(BC_key[i])
        except:
                label.append('no_BC')

    #counts['label']=label
    seqs=[]
    for l in label:
        #print(l)
        #print(seqs)
        try:
            #print('sequence')
            #print(fasta_dict[l])
            seqs.append(str(fasta_dict[l]).upper())
        except:
            seqs.append('NA')
    counts.insert(0,'Sequence',seqs)
    counts.insert(0, 'Label', label)
    #print(counts)

    mask=(counts['Barcode'].str.len() == bc_length)
    #print(mask)
    counts[mask]
    counts_filtered_t = counts[mask]
    #counts_filtered=counts_filtered_t.sort_values(by=['log2'])


    #res <- as.data.frame(t(sapply(unique(data$name),FUN=function(x) { sel <- which(data$name == x); c(((sum(data$X[sel])+1)/(length(sel)+1))/sum(data$X)*10^6,((sum(data$Y[sel])+1)/(length(sel)+1))/sum(data$Y)*10^6,length(sel)) } )))
    #res=''
    #res=pd.DataFrame()
    #normalize inserts
    #for i in set(counts_filtered_t.Label):
    #    sel=counts_filtered_t.loc[counts_filtered_t['Label']==i]
    #    #print(sel)
    #
    #    #new formula
    #    #dna=(sum(sel.dna_count)+1)/((len(sel.dna_count)+1))/(sum(counts_filtered_t.dna_count)/(10**6))
    #    #rna=(sum(sel.rna_count)+1)/((len(sel.rna_count)+1))/(sum(counts_filtered_t.rna_count)/(10**6))
    #
    #    #copied formula
    #    dna=(sum(sel.dna_count)+1)/((len(sel.dna_count)+1))/sum(counts_filtered_t.dna_count)*10**6
    #    rna=(sum(sel.rna_count)+1)/((len(sel.rna_count)+1))/sum(counts_filtered_t.rna_count)*10**6
    #    #rna=((sum(counts_filtered_t.rna_count[sel])+1)/(length(sel)+1))/sum(counts_filtered_t.rna)*10^6
    #    res_temp=(pd.DataFrame([dna,rna,len(sel.dna_count)]))
    #    res_t=res_temp.transpose()
    #    res_t.rename(index={0:str(i)},inplace=True)
    #
    #    if 'res' in locals():
    #        res=pd.concat([res,res_t])
    #    else:
    #        res=res_t
    #print(i)

    #normalize inserts
    # way more efficient, same result
    res = counts_filtered_t.groupby("Label", sort=False).agg({'rna_count':[('rna_sum','sum')],
                                                            'dna_count':[("dna_sum",sum)],
                                                            'Label':[("n_obs_bc",np.count_nonzero)]})
    res.columns = ['rna_sum','dna_sum','n_obs_bc']

    dna_total = sum(counts_filtered_t.dna_count)
    rna_total = sum(counts_filtered_t.rna_count)
    res.insert(0, 'dna_count',(res.dna_sum+1) / (res.n_obs_bc+1) / dna_total * 10**6)
    res.insert(1, 'rna_count',(res.rna_sum+1) / (res.n_obs_bc+1) / rna_total * 10**6)

    print('test')
    print(res.head())
    res = res[['dna_count','rna_count','n_obs_bc']]
    res.index.name = 'name'

    res.reset_index(inplace=True)
    res.insert(3, 'ratio',res.rna_count/res.dna_count)
    res.insert(4, 'log2',np.log2(res.ratio))
    print('merged')
    print(res.head())

    counts_filtered=dd.from_pandas(res,npartitions=1)
    #counts_filtered.fillna(0)
    print(counts_filtered.head())

    print(outfile)

    counts_filtered.to_csv([outfile], index=False,sep='\t')

    del res

    # this script processes the RNA and DNA counts and assigns the enhancer tag
    #outputs dataframe

    #CMD: python label_count_mat.py test.merged.H2.tsv ../lib_assoc_scripts/mp_assoc_original/bc_info_mp/Gracie_mp_filtered_coords_to_barcodes.pickle test.log2.fold.txt

if __name__ == '__main__':
    cli()