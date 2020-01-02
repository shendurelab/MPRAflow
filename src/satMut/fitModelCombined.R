library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

replicates=((length(args)-1)/2)
l <- c(0)
reps <- c()
for (i in 2:(length(args)-replicates)){

  rep=args[i+replicates]
  reps <- c(reps,rep)

  file=args[i]
  print(file)
  print(rep)
  input <- read.table(file,as.is=T,header=T,sep="\t",comment.char="~")
  if (i == 2){
    data <- input
  } else {
    input <- input[,colnames(input) %in% colnames(data)]
    data <- data[,colnames(data) %in% colnames(input)]
    data <- rbind(data,input)
  }
  l <- c(l,dim(data)[1])
}

# Create indicator variables
for (i in 1:length(reps)) {
  name <- paste0("isRep",reps[i])
  data[[name]] <- 0
  start <- l[i]+1
  end <- l[i+1]

  data$isRep1[start:end] <- 1
}

linComb = ""
for (i in 1:length(reps)) {
  name <- paste0("isRep",reps[i])
  data[[name]] <- data[[name]]*log2(data$DNA)
  if (i == 1) {
    linComb <- name
  } else {
    linComb = paste(linComb,name, sep=" + ")
  }
}

print("Start glm")
model <- glm(as.formula(sprintf("log2(RNA) ~ %s + %s",linComb,paste(colnames(data)[4:(dim(data)[2]-3)],collapse=" + "))),data=data)
model$data <- NULL
write.table(summary(model)$coefficients,file=args[1],row.names=T,col.names=F,quote=F,sep="\t")
