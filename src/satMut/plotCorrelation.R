library(dplyr)
library(tidyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

## styles
size.text = 20;
size.title = 25;
size.line = 1;
size.geom_line = 2;
size.geom_point = 4
standard_style <- theme_bw() + theme( plot.title = element_text(size = size.title, face="bold",hjust = 0.5),
                  panel.grid.major = element_blank() , panel.grid.minor = element_blank(), panel.border = element_blank(),
                  axis.text = element_text(colour = "black",size=size.text), axis.title = element_text(colour = "black",size=size.title), axis.ticks = element_line(colour = "black", size=1), axis.line.y = element_line(color="black", size = size.line), axis.line = element_line(colour = "black", size=size.line),
                  legend.key =  element_blank(), legend.text = element_text(size=size.title), legend.position=c(0.8, 0.1), legend.key.size = unit(2, 'lines'), legend.title=element_blank(),
                  legend.background = element_rect(size=1, linetype="solid",
                                                    colour ="black"))+ theme(axis.line.x = element_line(color="black", size = size.line))

colours=c("#333399","#339966")

readData <- function(file, name, threshold) {
  data <- read.table(file,as.is=T,header=T,sep="\t") %>%
            separate('Position', c("Position","Ref","Alt"), sep = "([\\_\\.])") %>%
            mutate(Alt = if_else(Alt == 'd', '-', Alt)) %>%
            mutate(Type=if_else(Alt=="-","1 bp deletions", "SNVs")) %>%
            filter(Barcodes >= threshold) %>%
            select(Position, Ref, Alt , Coefficient, Type)
  data$name <- name
  return(data)
}

getPlot <- function(file1,file2, name1,name2, threshold) {
  data <- readData(file1,name1,threshold) %>% inner_join(readData(file2,name2,threshold),by=c("Position", "Ref", "Alt", "Type"))

  correlation <- round(cor.test(data$Coefficient.x,data$Coefficient.y)$estimate,3)


  p <- ggplot(data, aes(x=data$Coefficient.x,y=data$Coefficient.y,colour=Type)) + geom_abline(intercept=0, colour="gray") +
      geom_point(size=3)  +
      ggtitle(paste(name1,"vs.",name2,"(rho_p",paste0(correlation, ")" ))) + labs(x = paste0("Log2 variant effect\n",name1), y= paste0(name2,"\nLog2 variant effect")) +
      standard_style + scale_colour_manual(values = colours) + guides(colour = guide_legend(override.aes = list(size=5)))
  return(p)
}


file1 <- args[1]
file2 <- args[2]
name1 <- args[3]
name2 <- args[4]
threshold <- as.double(args[5])

output <- args[6]

p <- getPlot(file1,file2,name1,name2,threshold)


ggsave(output,plot=p,height=20,width=20)
