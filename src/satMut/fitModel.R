args <- commandArgs(trailingOnly = TRUE)

data <- read.table(args[1],as.is=T,header=T,sep="\t",comment.char="~")

model <- glm(as.formula(sprintf("log2(RNA) ~ log2(DNA) + %s",paste(colnames(data)[4:dim(data)[2]],collapse=" + "))),data=data)
model$data <- NULL

# save(model,file=sprintf("nobackup/%s-VarMatrix.Model.R",cname))
write.table(summary(model)$coefficients,file=args[2],row.names=T,col.names=F,quote=F,sep="\t")
