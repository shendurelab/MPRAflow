#Adapted from Vikram Agarwal by Gracie Gordon

args <- commandArgs(trailingOnly = TRUE)


#exp=args[1]
exp=read.table(args[1], sep=',', header=T)
folder=args[2]
thresh=args[3]
#out=args[3]

print(folder)

nrm_reps=c()
all_reps=c()
##MAKE MASTER TABLE
for (i in 1:nrow(exp)){
   rep=exp[i,2]
   cond=exp[i,1]
   name=exp[i,5]
   file=(paste0(folder,name,"-byInsert.tsv"))
   
   tab=as.data.frame(read.table(file,header=TRUE))
   
   filter_tab=tab[tab$n_obs_bc >= thresh,]

   n_inserts=(dim(filter_tab)[1])
   
   cond_col=as.data.frame(rep(cond,n_inserts))
   rep_col=as.data.frame(rep(rep,n_inserts))
   colnames(cond_col)='condition'
   colnames(rep_col)='replicate'

   pref=as.data.frame(cbind(cond_col,rep_col))
   lab_tab=as.data.frame(cbind(pref,filter_tab)) 
  
   ##NORMALISZE PER REPLCATE TO COMBINE IN NEXT FUNCTION 
   nrm_reps1=lab_tab
   nrm_reps1$ratio=(nrm_reps1$ratio)/(median(nrm_reps1$ratio))
   nrm_reps1$log2=round(log2(nrm_reps1$ratio),8)
   nrm_reps=rbind(nrm_reps,nrm_reps1)
   all_reps=rbind(all_reps,lab_tab)
}

write.table(all_reps,file=paste0("allreps.tsv"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE )

uniq_names=unique(all_reps$name)
##MAKE AVERAGED ACROSS REPLICATES
all_avg=c()
for (cond in unique(exp$cond)){
   for (insert in uniq_names){
       filt_insert=all_reps[(nrm_reps$name == insert & nrm_reps$cond == cond),]
       m_count=mean(filt_insert$n_obs_bc)
       m_ratio=(mean(filt_insert$ratio))
       m_log2=(log2(mean(filt_insert$ratio)))
       line=c(cond,insert,m_ratio,m_log2,m_count)
       all_avg=rbind(all_avg,line)      

   } 

}
colnames(all_avg)=c('condition','name','mean_ratio','mean_log2','mean_n_obs_bc')


write.table(all_avg,file=paste0("average_allreps.tsv"),quote=FALSE,sep='\t',row.names=FALSE,col.names=TRUE )

print('done')
