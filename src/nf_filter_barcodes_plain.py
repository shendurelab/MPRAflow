import sys
import pickle
import pandas as pd
from collections import Counter
from collections import defaultdict
import pyarrow
from matplotlib import pyplot as plt
import seaborn
import re

pref=sys.argv[1]
enh_bc=sys.argv[2]
count_table=sys.argv[3]


#enh_bc=pref+'_coords_to_barcodes.pickle'
#count_table=pref+'_barcodes_per_candidate-no_repeats-no_jackpots.feather'

cov_thresh=sys.argv[4]
label_file=sys.argv[5]
print(enh_bc)
print(count_table)

#function to get rid of lowly covered BCs and get rid of BCs that map to multiple other enhancers
def covered_no_perm(min_threshold, mydict):
    out_dict={}
    ##get rid of BCs that show up less than 3x
    for k,v in mydict.items():
        bc_counts=(Counter(v))
        filtered=({x : bc_counts[x] for x in bc_counts if bc_counts[x] >= int(min_threshold) and bc_counts[x]<= int(50)})
        new_BCs=(list(filtered.keys()))
        mydict[k]=new_BCs

    #get rid of BCs that are premiscuous 
    all_bcs=(mydict.values())
    flat_BC = [item for sublist in all_bcs for item in sublist]
    counts=Counter(flat_BC)
    permis=({x : counts[x] for x in counts if counts[x] > 1})
    dups=(list(permis.keys()))
    print(dups)

    ##remove duplicates and get rid of poly A/G/T/C (7 in a row maybe change to more)
    for k,v in mydict.items():
        new_v=[]
        for i in v:
            if i not in dups:
                m=0
                for match in re.finditer(r"(\w)\1\1\1\1\1\1", i):   
                    m=i[match.start():match.end()]
                if(m==0):
                    new_v.append(i)
        out_dict[k]=new_v
    return(out_dict)



#takes in a pandas series count dataframe and plots the coverage coloring by type of control sequence
def create_lib_plots(c_obj, plot_name, sum_name,lab_opt="FALSE"):
    #change counter to dataframe
    df=c_obj
    #df = pd.DataFrame.from_dict(c_obj,orient='index').reset_index()
    print(df)
    #add label column for 'class'
    if lab_opt =="TRUE":
        label=[]
        for i in df.coord:
            tmp=(str.split(str(i),'_C:')[0])
            label.append(str(str.split(str(tmp),'_r')[0]))
        df['label']=label
        df=df.sort_values('label')
        #plot
        #a4_dims = (11.7, 8.27)
        plt.figure(figsize=(16, 6))
        p=seaborn.barplot(data=df,x="coord", y="n_barcodes",hue="label",dodge=False)
        print(plot_name)
        p.set(xticklabels=[])
        fig=p.get_figure()
        fig.savefig(plot_name)
    else:
        plt.figure(figsize=(16, 6))
        p=seaborn.barplot(data=df,x="coord", y="n_barcodes",dodge=False)
        p.set(xticklabels=[])
        fig=p.get_figure()
        fig.savefig(plot_name)

    print(df.describe)
    #summary=pd.DataFrame(df.describe)
    #summary.to_csv(sum_name, index=False,)
    df.describe().to_csv(sum_name)

#plot original data
create_lib_plots(pd.read_feather(count_table),pref+"_original_counts.png","original_count_summary.txt","TRUE")
#filter and plot new data
dic=pickle.load(open(enh_bc,'rb'))
filtered_dict=covered_no_perm(cov_thresh, dic)

#new counts
filtered_counts=(pd.Series(filtered_dict, name = 'n_barcodes').str.len().rename_axis('coord').reset_index())
create_lib_plots(filtered_counts, pref+"_filtered_counts.png","filtered_count_summary.txt","FALSE")
pickle.dump(filtered_dict,open(pref+"_filtered_coords_to_barcodes.p","wb"))











