
#adapted from Vikram Agarwal by Gracie Gordon

cpath <- grep('conda', .libPaths(), value=TRUE, ignore.case=TRUE)
.libPaths(cpath)


library(ggplot2)
library(dplyr)

stringsAsFactors = FALSE
args <- commandArgs(trailingOnly = TRUE)


# condition
cond=args[1]
# labels
if (file.exists(args[2])) {
  label_f=as.data.frame(read.table(args[2],sep='\t',header=T,stringsAsFactors = F))
  colnames(label_f)=c('name','label')
  useLabels=TRUE
} else {
  useLabels=FALSE
}

# threshold
thresh=strtoi(args[3])

# replicates and count files
num_replicates=((length(args)-3)/2)
files=args[4:(4+num_replicates-1)]
replicates=args[(4+num_replicates):length(args)]

data <- data.frame(File=files,Replicate=replicates)
data$Condition <- cond

print(data)

# pairwise comparison only if more than one replicate
if(data %>% nrow >1){

  # make pairwise combinations
  selected <- combn(data$Replicate,2)
  print('sel')
  print(selected)


  print('reps')
  for(i in seq(1,dim(selected)[2])){
    print(selected[,i])
    r1=selected[1,i]
    r2=selected[2,i]
    data1 <- read.table(as.character((data %>% filter(Replicate == r1))$File),as.is=T,sep="\t",header=T,stringsAsFactors = F)
    data2 <- read.table(as.character((data %>% filter(Replicate == r2))$File),as.is=T,sep="\t",header=T,stringsAsFactors = F)

    #FIXME what is this no_BC name?? we should docuent it
    data1<-data1 %>% filter(n_obs_bc > thresh, name != 'no_BC')
    data2<-data2 %>% filter(n_obs_bc > thresh, name != 'no_BC')

    res <- data1 %>% inner_join(data2,by=c('name'))
    if (useLabels){
      res <- res %>% inner_join(label_f, by=c('name'))
    } else {
      res$label = args[2]
    }

    dna_p <- ggplot(res, aes(log2(dna_count.x), log2(dna_count.y))) +
                geom_point(aes(colour = label), show.legend = TRUE) +
                xlim(-5,5) + ylim(-5,5) +
                xlab(sprintf(paste("log2 Normalized DNA count per insert,\n replicate", r1))) +
                ylab(sprintf(paste("log2 Normalized DNA count per insert,\n replicate", r2))) +
                geom_text(x=0, y=4.5,label=sprintf("   r = %.2f", cor(log2(res$dna_count.x),log2(res$dna_count.y),method="pearson")),size=10) +
                geom_text(x=0, y=4, label=sprintf("rho = %.2f", cor(res$dna_count.x,res$dna_count.y,method="spearman")),size=10) +
                geom_abline(intercept = 0, slope = 1) +
                theme_classic(base_size = 30)
    ggsave(sprintf("%s_%s_%s_DNA_pairwise.png",cond,as.character(r1),as.character(r2)),width=15,height=15)

    rna_p <- ggplot(res, aes(log2(rna_count.x), log2(rna_count.y))) +
                geom_point(aes(colour = label), show.legend = TRUE) +
                xlim(-5,5) + ylim(-5,5) +
                xlab(sprintf(paste("log2 Normalized RNA count per insert,\n replicate", r1))) +
                ylab(sprintf(paste("log2 Normalized RNA count per insert,\n replicate", r2))) +
                geom_text(x=0, y=4.5,label=sprintf("   r = %.2f", cor(log2(res$rna_count.x),log2(res$rna_count.y),method="pearson")),size=10) +
                geom_text(x=0, y=4, label=sprintf("rho = %.2f", cor(res$rna_count.x,res$rna_count.y,method="spearman")),size=10) +
                geom_abline(intercept = 0, slope = 1) +
                theme_classic(base_size = 30)
    ggsave(sprintf("%s_%s_%s_RNA_pairwise.png",cond,as.character(r1),as.character(r2)),plot=rna_p,width=15,height=15)

    ratio_p <- ggplot(res, aes(log2(ratio.x), log2(ratio.y))) +
                  geom_point(aes(colour = label), show.legend = TRUE) +
                  xlim(-5,5) + ylim(-5,5) +
                  xlab(sprintf(paste("log2 RNA/DNA per insert,\n replicate", r1))) +
                  ylab(sprintf(paste("log2 RNA/DNA per insert,\n replicate", r2))) +
                  geom_text(x=0, y=4.5,label=sprintf("   r = %.2f", cor(log2(res$ratio.x),log2(res$ratio.y),method="pearson")),size=10) +
                  geom_text(x=0, y=4, label=sprintf("rho = %.2f", cor(res$ratio.x,res$ratio.y,method="spearman")),size=10) +
                  geom_abline(intercept = 0, slope = 1) +
                  theme_classic(base_size = 30)
    ggsave(sprintf("%s_%s_%s_Ratio_pairwise.png",cond,as.character(r1),as.character(r2)),plot=ratio_p,width=15,height=15)



    outs = as.character(sprintf("%s vs %s:  RNA: %.5f DNA: %.5f ratio: %.5f NormSymmetry: %d",r1,r2,
        cor(res$rna_count.x,res$rna_count.y,method="spearman"),
        cor(res$dna_count.x,res$dna_count.y,method="spearman"),
        cor(res$ratio.x,res$ratio.y,method="spearman"),
        abs(length(which((res$ratio.x-res$ratio.y)>0))-length(which((res$ratio.x-res$ratio.y)<0)))
        + abs(length(which((res$ratio.x-res$ratio.y)>0))-length(which((res$ratio.x-res$ratio.y)<0)))
        + abs(length(which((res$ratio.x-res$ratio.y)>0))-length(which((res$ratio.x-res$ratio.y)<0)))
    ))

    write(outs,file=sprintf("%s_%s_%s_correlation.txt",cond,as.character(r1),as.character(r2)),append=TRUE)

  }
}


print('hist')
nrows=ceiling(num_replicates/2)
png(sprintf("%s_barcodesPerInsert.png",cond),width=300,height=244,type="cairo")
par(mfrow=c(nrows,2), mar=c(2.5,3,2.5,3))

for(n in 1:(data%>%nrow)){
    data1 <- read.table(as.character(data[n,]$File),as.is=T,sep="\t",header=T,stringsAsFactors = F) %>% filter(name != 'no_BC')

    print(head(data1))
    print(head(data1$n_obs_b))
    #hist(as.numeric(data1$n_obs_bc), main=paste("replicate", n, sep=' '), xlab=NULL, las=1)
    hist(as.numeric(data1$n_obs_bc),breaks=100,xlim=c(0,300), main=paste("replicate", data[n,]$Replicate, sep=' '), xlab=NULL, las=1)
    abline(v=median(as.numeric(data1$n_obs_bc)), lwd=6, col='red')
}

#
print('boxplot')
iter=0
all=''
dlog=''
dname=''

png(sprintf("%s_barcodesPerInsert_box.png",cond),width=800,height=800,type="cairo")

for(n in 1:(data%>%nrow)){
    data1=read.table(as.character(data[n,]$File),as.is=T,sep="\t",header=T,stringsAsFactors = F)  %>% filter(name != 'no_BC')

    dlog=c(dlog,data1$log2)
    dname=c(dname,data1$name)

 }
 #all=as.data.frame(cbind(((dname)),as.numeric(dlog)))
 all=as.data.frame(cbind((as.character(dname)),as.character(dlog)))
 colnames(all)=c('name','log2')
 print(head(all))

 if (useLabels){
   all <- all %>% inner_join(label_f, by=c('name'))
 } else {
   all$label = args[2]
 }

 print('merged')
 print(head(all))
 all$name <- factor(all$name)
 all$log2 <- as.numeric(as.character(all$log2))
 all$label <- as.factor(all$label)
 #all=as.data.frame((all))
 all=all[-1,]
 all=all[order(all$log2),]
 #all$name=as.factor(all$name)
 #all$log2=as.numeric(all$log2)

 print(head(all))
 print(str(all))
 #colnames(dlog)=c('name','log2')
 #colnames(all) <- as.character(unlist(all[1,]))
 #all=all[-1,]
 bymedian=with(all,reorder(name,-log2,median,order=TRUE))
 #fac <- with(boxData, reorder(temp, value, median, order = TRUE))
 all$name=factor(all$name,levels=levels(bymedian))
 #boxData$temp <- factor(boxData$temp, levels = levels(fac))
 #all$name=with(all,reorder(name,-log2,median))
 #print(head(all$bymedian))
 png(sprintf("%s_all_barcodesPerInsert_box.png",cond),width=800,height=800,type="cairo")
 bp<-ggplot(all, aes(x=name, y=log2, color=label)) +geom_boxplot()+xlab('insert')+ylab('log2 fold change')+
 theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_text(size=15),axis.title.y= element_text(size=15),axis.text.y=element_text(size=15),legend.text = element_text(size=15))
 print(bp)
 #boxplot(all$log2~bymedian, main='activity',xlab='enhancer', ylab='log2 fold change',las=2)
 #boxplot(all$log2~all$name,main='activity',xlab='enhancer', ylab='log2 fold change',las=2)
 dev.off()
 png(sprintf("%s_group_barcodesPerInsert_box.png",cond),width=800,height=800,type="cairo")
 bp<-ggplot(all, aes(x=label, y=log2, fill=label)) +geom_violin()+geom_boxplot(width=0.1,fill='white')+xlab('insert')+ylab('log2 fold change')+
 theme(axis.text.x = element_text(angle = 90, hjust=1,size=15),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_text(size=15),axis.title.y= element_text(size=15),axis.text.y=element_text(size=15),legend.text = element_text(size=15))
 print(bp)
 dev.off()
