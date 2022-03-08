library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
require(cowplot)

input_data <- read.table(file = 'NTD_DMS_scores.tsv', sep = '\t', header = TRUE)
mutant <- c('WT','N17Q','N61Q','N74Q','N122Q','N149Q','N165Q','N234Q','N282Q')
scitrans_data <- c(100,53.9147,61.9767,85.9303,66.9767,63.3333,49.2248,82.7261,54.3928)
input_data[input_data$mut=='WT',]$Exp_score
dms_data <- c(0.7430105,input_data[input_data$mut=='N17Q',]$Exp_score,input_data[input_data$mut=='N61Q',]$Exp_score,input_data[input_data$mut=='N74Q',]$Exp_score,input_data[input_data$mut=='N122Q',]$Exp_score,input_data[input_data$mut=='N149Q',]$Exp_score,input_data[input_data$mut=='N165Q',]$Exp_score,input_data[input_data$mut=='N234Q',]$Exp_score,input_data[input_data$mut=='N282Q',]$Exp_score)
results <-data.frame(mutant, scitrans_data, dms_data)

plot_dms_vs_pseudovirus <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=dms_data, y=scitrans_data)) +
    geom_point(size=0.5,pch=16, alpha=0.5) +
    #geom_violin(size=0.5) +
    #geom_boxplot(width=0.3, size = 0.3, color="black", outlier.shape=NA, alpha=0.5) +
    #geom_beeswarm(size=0.3, pch=16, alpha=0.5) +
    theme_cowplot(12) +
    theme(plot.title=element_blank(),
          plot.background = element_rect(fill = "white"),
          axis.title=element_text(size=textsize,face="bold"),
          axis.text=element_text(size=textsize,face="bold"),
          legend.key.size=unit(0.1,'in'),
          legend.spacing.x=unit(0.03, 'in'),
          legend.title=element_blank(),
          legend.text=element_text(size=textsize-1,face="bold"),
          legend.position='right') +
    labs(x=bquote(bold(paste('NTD DMS Expression Score'))),y=bquote(bold(paste('Pseudovirus Production'))))
  ggsave(graphname, p, height=2, width=2)
}