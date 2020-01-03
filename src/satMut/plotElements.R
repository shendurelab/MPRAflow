library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
args <- commandArgs(trailingOnly = TRUE)


defaultColours=c("A"="#0f9447","C"="#235c99","T"="#d42638","G"="#f5b328","-"="#cccccc",
                         "Significant"="#005500","Not significant"="red")


standard_SatMut_region_style <- function() {
  ## styles
  size.text = 20;
  size.title = 25;
  size.line = 1;
  size.geom_line = 2;
  size.geom_point = 4
  standard_style <- theme_bw() +
    theme(plot.title = element_text(size = size.title, face="bold",hjust = 0.5),
          panel.grid.major = element_blank() , panel.grid.minor = element_blank(), panel.border = element_blank(),
          axis.text = element_text(colour = "black",size=size.text), axis.title = element_text(colour = "black",size=size.title), axis.ticks = element_line(colour = "black", size=1), axis.line.y = element_line(color="black", size = size.line), axis.line = element_line(colour = "black", size=size.line),
          legend.key =  element_blank(), legend.text = element_text(size=size.text),
          legend.position="top", legend.box.just = "left",  legend.background = element_rect(fill = "transparent", colour = "transparent"), legend.margin = margin(0, 0, 0, 0),
          legend.key.size = unit(2, 'lines'), legend.title=element_text(size=size.text))  +
    theme(axis.line.x = element_line(color="black", size = size.line))
}


modify.filterdata <- function(data, barcodes=10, threshold=1e-5, deletions=TRUE, range=NULL) {
  data <- data %>% select(Position,Ref,Alt,Barcodes, Coefficient, pValue) %>%
    filter(Barcodes >= barcodes) %>%
    mutate(significance=if_else(pValue<threshold,"Significant", "Not significant")) %>%
    mutate(printpos=if_else(Alt=="A",as.double(Position)-0.4,if_else(Alt=="T",as.double(Position)-0.2, if_else(Alt=="G",as.double(Position)+0.0,if_else(Alt=="C",as.double(Position)+0.2,as.double(Position)+0.4)))))
  if (!deletions) {
    data <- data %>% filter(Alt != "-")
  }
  if (!is.null(range)) {
    data <- data %>% filter(Position >= range[1] & Position <= range[2])
  }
  return(data)
}

getPlot <- function(data,name) {
  colours <- defaultColours

  refs <- data$Ref %>% unique()
  aesRefsValues <- c(if_else("A" %in% refs,15,c()),if_else("C" %in% refs,16,c()),if_else("G" %in% refs,17,c()),if_else("T" %in% refs,18,c()))
  aesRefsShape <- c(if_else("A" %in% refs,0,c()),if_else("C" %in% refs,1,c()),if_else("G" %in% refs,2,c()),if_else("T" %in% refs,5,c()))
  aesRefsValues <- aesRefsValues[!is.na(aesRefsValues)]
  aesRefsShape <- aesRefsShape[!is.na(aesRefsShape)]

  alts <- data$Alt %>% unique()
  sigs <- data$significance %>% unique()
  aesSize<-c(rep(8,length(alts)), rep(1,length(sigs)))
  altBreaks<-c(as.character(alts), sigs)

  data <- data %>% select(printpos,Coefficient,significance,Alt,Ref) %>% dplyr::rename(Position=printpos)
  p <- ggplot() +
    geom_linerange(data = data, aes(x=Position,ymin=0,ymax=Coefficient, colour=significance), size=1, show.legend = TRUE) +
    geom_point(data= data, aes(x=Position,y=Coefficient,colour=Alt, shape=Ref), size=3, show.legend = TRUE) +
    scale_shape_manual("Reference:",values=aesRefsValues, guide=guide_legend(override.aes = list(size=7, linetype=0, shape=aesRefsShape,nrow=2))) + #
    scale_colour_manual("Alternative:", values = colours, breaks=altBreaks, labels=altBreaks,
                        guide= guide_legend(override.aes = list(size=aesSize, linetype=1, shape=32))) +
    ggtitle(name) + labs(x = paste0("Position"), y= "Log2 variant effect") + standard_SatMut_region_style()
  return(p)
}

file <- args[1]
name <- args[2]
barcodes <- args[3]
significance <- args[4]
output <- args[5]

input <- read.table(args[1],as.is=T,header=T,sep="\t") %>%
          separate('Position', c("Position","Ref","Alt"), sep = "([\\_\\.])") %>%
          mutate(Alt = if_else(Alt == 'd', '-', Alt))
input <- modify.filterdata(input, barcodes=barcodes, threshold=significance, deletions=TRUE, range=NULL)

p <- getPlot(input, name)


ggsave(filename=output,plot=p, width=20, height=10)
