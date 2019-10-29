# this script processes the RNA and DNA counts and assigns the enhancer tag
#outputs dataframe

#CMD: python label_count_mat.py test.merged.H2.tsv ../lib_assoc_scripts/mp_assoc_original/bc_info_mp/Gracie_mp_filtered_coords_to_barcodes.pickle test.log2.fold.txt 

import pickle
import numpy
import pandas
import sys

from Bio import SeqIO


#outfile name
name="final_labeled_counts.txt"


#read in files
data=sys.argv[1]
coord_file=sys.argv[2]
#rna_f=sys.argv[2]
#assoc_f=sys.argv[2]
out_f=sys.argv[3]
design_file=sys.argv[4]

design=open(design_file)

fasta_dict = {rec.id : rec.seq for rec in SeqIO.parse(design, "fasta")}

#fasta_dict=SeqIO.to_dict(SeqIO.parse(design, "fasta"))


#print(fasta_dict['HepG2_negatives_TFs_C:SLEA_hg18:chr2:210861483-210861650|54:V_AHRARNT_02:GGGGATCGCGTGCCAGCCC;76:V_COUPTF_Q6:CCCCCTGACCTTTGCCCCCTGCC;102:V_HNF3ALPHA_Q6:TGTTTGCTTTG:001'])

#DNA=pandas.read_csv(dna_f,header=None,sep='\s+')
#RNA=pandas.read_csv(rna_f,header=None,sep='\s+')
counts=pandas.read_csv(data,header='infer',sep=',')

assoc=pickle.load( open(coord_file,'rb'))
#assoc=pickle.load( open('coords_to_barcodes.pickle','rb'))

#print(DNA.head())
#print(RNA.head())


#calc log2(RNA/DNA)

#counts.columns=['BC','countDNA','countRNA']
#DNA.columns=['countDNA','BC']
#RNA.columns=['countRNA','BC']
#result=DNA.merge(RNA, left_on='BC', right_on='BC')

#counts["log2"]=numpy.log2(counts['countRNA']/counts['countDNA'])

#print(result)
#flip dict so that BCs are the Keys
BC_key = {}
for k,v in assoc.items():
    for x in v:
         BC_key.setdefault(x,k)

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
        print(fasta_dict[l])
        seqs.append(str(fasta_dict[l]).upper())
    except:
        seqs.append('NA')
counts.insert(0,'Sequence',seqs)
counts.insert(0, 'label', label)
print(counts)

#filter results
#mask=pandas.Series(result['BC'].str.len() == 10,name='bools')

mask=(counts['Barcode'].str.len() == 15)
print(mask)
counts[mask]
counts_filtered_t = counts[mask]
#counts_filtered=counts_filtered_t.sort_values(by=['log2'])

counts_filtered=counts_filtered_t

print(counts_filtered)
#result.filtered = result.loc[mask]
counts_filtered.to_csv(out_f,header=True,index=False,sep="\t") 



