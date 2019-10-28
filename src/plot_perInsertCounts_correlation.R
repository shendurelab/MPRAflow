#adapted from Vikram Agarwal
list.of.packages <- c("ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)

stringsAsFactors = FALSE
args <- commandArgs(trailingOnly = TRUE)


#exp=args[1]
exp=read.table(args[1], sep=',', header=T)
folder=args[2]
out=args[3]
label_f=as.data.frame(read.table(args[4],sep='\t',header=T))
colnames(label_f)=c('name','label')
head(label_f)
maxRep <- 1000
print(exp)
o_list=c()
#for (cond in unique(exp$condition)){
for (cond in unique(exp$condition)){
    print(cond)
    c=exp[which(exp$condition==cond),]
    #c=exp
    nReps=(dim(c)[1])
    colnames(c)=colnames(exp)
    #selected <- combn(1:nReps,2)
    print(c)
    if(dim(c)[1]>1){
        selected <- combn(c$name,2)
        print('sel')
        print(selected)
        #png(sprintf("%s_all_pairwise.png",cond),width=800,height=800,type="cairo")
        print('tests') 
        print(dim(selected))
        print(dim(selected)[1])
        print(dim(selected)[2])

        #apply(selected,2,function(x) {
        print('reps')
        print(dim(selected))
        for(i in seq(1,dim(selected)[2])){
            print(selected[,i])
        
            #print(x)
            r1=as.character(selected[1,i])
            r2=as.character(selected[2,i])
            print(r1)
            print(r2)
            data1 <- read.table(sprintf("%s/%s-byInsert.tsv",folder,r1),as.is=T,sep="\t",header=T)
            data2 <- read.table(sprintf("%s/%s-byInsert.tsv",folder,r2),as.is=T,sep="\t",header=T)
             
            thresh=10
            data1=data1[data1$n_obs_bc > thresh,]
            data2=data2[data2$n_obs_bc > thresh,]
            

            data1_f=data1[data1$name != 'no_BC',]
            data2_f=data2[data2$name != 'no_BC',]
            #print(head(data1_f))
            #print(head(data2_f))


            label_t=label_f
            print('merge')
            res <- merge(data1_f,data2_f,by="name")
            res1 <- merge(res,label_t, by="name")
            print(head(res1))
        
            selSamples <- c(sprintf('%s',r1),sprintf('%s',r2))
            
            png(sprintf("%s_%s_%s_DNA_pairwise.png",cond,r1,r2),width=800,height=800,type="cairo")
            dna_p<-ggplot(res1, aes(log2(dna_count.x), log2(dna_count.y)))+
            geom_point(aes(colour = label), show.legend = TRUE)+xlim(-5,5)+ylim(-5,5)+xlab(sprintf(paste("log2 Normalized DNA count per insert,\n replicate", r1)))+
            ylab(sprintf(paste("log2 Normalized DNA count per insert,\n replicate", r2)))+
            geom_text(x=0, y=4.5,label=sprintf("   r = %.2f", cor(log2(res1$dna_count.x),log2(res1$dna_count.y),method="pearson")),size=10)+
            geom_text(x=0, y=4, label=sprintf("rho = %.2f", cor(res1$dna_count.x,res1$dna_count.y,method="spearman")),size=10)+geom_abline(intercept = 0, slope = 1)+ theme_classic(base_size = 30)
            print(dna_p)
            dev.off()

            png(sprintf("%s_%s_%s_RNA_pairwise.png",cond,r1,r2),width=800,height=800,type="cairo")
            rna_p<-ggplot(res1, aes(log2(rna_count.x), log2(rna_count.y)))+
            geom_point(aes(colour = label), show.legend = TRUE)+xlim(-5,5)+ylim(-5,5)+xlab(sprintf(paste("log2 Normalized RNA count per insert,\n replicate", r1)))+
            ylab(sprintf(paste("log2 Normalized RNA count per insert,\n replicate", r2)))+
            geom_text(x=0, y=4.5,label=sprintf("   r = %.2f", cor(log2(res1$rna_count.x),log2(res1$rna_count.y),method="pearson")),size=10)+
            geom_text(x=0, y=4, label=sprintf("rho = %.2f", cor(res1$rna_count.x,res1$rna_count.y,method="spearman")),size=10)+geom_abline(intercept = 0, slope = 1)+ theme_classic(base_size = 30)
            print(rna_p)
            dev.off()

            png(sprintf("%s_%s_%s_Ratio_pairwise.png",cond,r1,r2),width=800,height=800,type="cairo")
            ratio_p<-ggplot(res1, aes(log2(ratio.x), log2(ratio.y)))+
            geom_point(aes(colour = label), show.legend = TRUE)+xlim(-5,5)+ylim(-5,5)+xlab(sprintf(paste("log2 RNA/DNA per insert,\n replicate", r1)))+
            ylab(sprintf(paste("log2 RNA/DNA per insert,\n replicate", r2)))+
            geom_text(x=0, y=4.5,label=sprintf("   r = %.2f", cor(log2(res1$ratio.x),log2(res1$ratio.y),method="pearson")),size=10)+
            geom_text(x=0, y=4, label=sprintf("rho = %.2f", cor(res1$ratio.x,res1$ratio.y,method="spearman")),size=10)+geom_abline(intercept = 0, slope = 1)+ theme_classic(base_size = 30)
            print(ratio_p)
            dev.off()

            #plot(log2(res1$dna_count.x),log2(res1$dna_count.y), pch=19, cex=0.3, col = "darkblue", las=1, cex.axis=1.5, las = 2, bty="n", xlab = "", ylab = "", xlim=c(-2,3), ylim=c(-2,3))
            #mtext(text = sprintf(paste("log2 Normalized DNA count per insert,\n replicate", r1)), side = 1, line = 6, cex=1.1)
            #mtext(text = sprintf(paste("log2 Normalized DNA count per insert,\n replicate", r2)), side = 2, line = 3.5, cex=1.1)
            #abline(0,1,lty=2)
            #legend(-2.3, 3, cex=1.7, bg="white", bty="n", legend = c(sprintf("   r = %.2f", cor(log2(res1$dna_count.x),log2(res1$dna_count.y),method="pearson")),
            #               sprintf("rho = %.2f", cor(res1$dna_count.x,res1$dna_count.y,method="spearman"))))
            
            #plot(log2(res1$rna_count.x),log2(res1$rna_count.y), pch=19, cex=0.3, col = "darkgreen", las=1, cex.axis=1.5, las = 2, bty="n", xlab = "", ylab = "", xlim=c(-2,3), ylim=c(-2,3)) 
            #mtext(text = sprintf(paste("log2 Normalized RNA count per insert,\n replicate", r1)), side = 1, line = 6, cex=1.1)
            #mtext(text = sprintf(paste("log2 Normalized RNA count per insert,\n replicate", r2)), side = 2, line = 3.5, cex=1.1)
            #abline(0,1,lty=2)
            #legend(-2.3, 3, cex=1.7, bg="white", bty="n", legend = c(sprintf("   r = %.2f", cor(log2(res1$rna_count.x),log2(res1$rna_count.y),method="pearson")),
            #               sprintf("rho = %.2f", cor(res1$rna_count.x,res1$rna_count.y,method="spearman"))))
        
            #plot(log2(res1$ratio.x),log2(res1$ratio.y), pch=19, cex=0.3, col = "darkred", las=1, cex.axis=1.5, xlim=c(-2,4), ylim=c(-2,4), las = 2, bty="n", xlab = "", ylab = "")
            #mtext(text = sprintf(paste("log2 RNA/DNA per insert,\n replicate",r1)), side = 1, line = 6, cex=1.1)
            #mtext(text = sprintf(paste("log2 RNA/DNA per insert,\n replicate", r2)), side = 2, line = 3, cex=1.1)
            #abline(0,1,lty=2)
            #legend(-2.3, 4, cex=1.7, bg="white", bty="n", legend = c(sprintf("   r = %.2f", cor(log2(res1$ratio.x),log2(res1$ratio.y),method="pearson")),
            #               sprintf("rho = %.2f", cor(res1$ratio.x,res1$ratio.y,method="spearman"))))

       
            #sink(paste0(out,'.txt'))
            outs=as.character(sprintf("%s vs %s:  RNA: %.5f DNA: %.5f ratio: %.5f NormSymmetry: %d",r1,r2,
                  cor(res1$rna_count.x,res1$rna_count.y,method="spearman"),
                  cor(res1$dna_count.x,res1$dna_count.y,method="spearman"),
                  cor(res1$ratio.x,res1$ratio.y,method="spearman"),
                  abs(length(which((res1$ratio.x-res1$ratio.y)>0))-length(which((res1$ratio.x-res1$ratio.y)<0)))
                  +abs(length(which((res1$ratio.x-res1$ratio.y)>0))-length(which((res1$ratio.x-res1$ratio.y)<0)))
                  +abs(length(which((res1$ratio.x-res1$ratio.y)>0))-length(which((res1$ratio.x-res1$ratio.y)<0)))
            ))
            #print(outs)
            #o_list=c(o_list,outs)

            write(outs,file=paste0(out,'.txt'),append=TRUE)
                  
        }
    }

}

print('hist')
for (cond in unique(exp$condition)){
    c=exp[which(exp$condition==cond),]
    nReps=(dim(c)[1])
    nrows=ceiling(nReps/2)
    png(sprintf("%s_barcodesPerInsert.png",cond),width=300,height=244,type="cairo")
    par(mfrow=c(nrows,2), mar=c(2.5,3,2.5,3))
    
    for(n in c$name){
        data1_f <- read.table(sprintf("%s/%s-byInsert.tsv",folder,n),as.is=T,sep="\t",header=T)
        data1=data1_f[data1_f$name != 'no_BC',]

        print(head(data1))
        print(head(data1$n_obs_b))    
        #hist(as.numeric(data1$n_obs_bc), main=paste("replicate", n, sep=' '), xlab=NULL, las=1)
        hist(as.numeric(data1$n_obs_bc),breaks=100,xlim=c(0,300), main=paste("replicate", n, sep=' '), xlab=NULL, las=1)
        abline(v=median(as.numeric(data1$n_obs_bc)), lwd=6, col='red')
    }
    dev.off()  

} 

print('boxplot')
for (cond in unique(exp$condition)){
    print(cond)
    iter=0
    all=''
    dlog=''
    dname=''
    c=exp[which(exp$condition==cond),]
    png(sprintf("%s_barcodesPerInsert_box.png",cond),width=800,height=800,type="cairo")
    
    for(n in c$name){
        data=read.table(sprintf("%s/%s-byInsert.tsv",folder,n),as.is=T,sep="\t",header=T)
        data1=data[data$name != 'no_BC',]
        
        dlog=c(dlog,data1$log2)
        dname=c(dname,data1$name)
        #d_log=as.data.frame(cbind(data1$name,as.double(data1$log2)))
        #colnames(d_log)=c('name','log2')
        #colnames(d_log)=n
        
        #if(iter !=0){
        #    all=rbind(all,d_log)
            #all=cbind(all, d_log)        
        #    iter=iter+1
        #}
        #else{
            #temp=as.data.frame(data1$name)
            #colnames(temp)='name'
            #all=cbind(temp,d_log)
        #    all=d_log
        #    iter=iter+1
        #}
     }
     #all=as.data.frame(cbind(((dname)),as.numeric(dlog)))
     all1=as.data.frame(cbind((as.character(dname)),as.character(dlog)))
     print(head(label_f))
     label_t=label_f
     colnames(all1)=c('name','log2')
     print(head(all1))
     all=merge(all1,label_t, by="name")
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

}






