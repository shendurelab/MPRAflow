library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

for (i in 2:length(args)){

  file=args[i]

  input <- read.table(file,as.is=T,header=F,sep="\t",comment.char="~")
  colnames(input) <- c("Variant","BC","DNA","RNA")
  if (i == 2){
    data <- input
  } else {
    data <- data %>% inner_join(input, by=c("Variant")) %>%
          mutate(BC = BC.x+BC.y,DNA =DNA.x+DNA.y,RNA=RNA.x+RNA.y) %>%
          select(Variant,BC,DNA,RNA)
  }
}

write.table(data,file=args[1],row.names=F,col.names=F,quote=F,sep="\t")
