#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
library(qualpalr)
library(ggbeeswarm)
require(cowplot)

plot_RSA_vs_param <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=RSA, y=expr_avg)) +
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
         labs(x=bquote(bold(paste('Relative solvent accessibility'))),y=bquote(bold(paste('Mean expression score'))))
  ggsave(graphname, p, height=2, width=2)
  }

df <- read_tsv('result/RBD_exp_RSA.tsv')
plot_RSA_vs_param(df, 'graph/RBD_exp_vs_RSA.png')
print ((paste("Cor between expression score and RSA:", cor(df$expr_avg, df$RSA, method='spearman'))))
