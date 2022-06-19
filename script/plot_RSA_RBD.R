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
  p <- ggplot(df,aes(x=RSA, y=mean_exp_score)) +
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
    scale_y_continuous(breaks=c(0.25, 0.5, 0.75, 1, 1.25))+
    labs(x=bquote(bold(paste('Relative solvent accessibility'))),y=bquote(bold(paste('Mutational tolerability'))))
  ggsave(graphname, p, height=2, width=2)
}

resi_classification <- function(RSA_monomer, delta_RSA){
  if (RSA_monomer < 0.05){return ('Buried')}
  else if (delta_RSA > 0.5*RSA_monomer){return ('Interface')}
  else (return ('Exposed'))
}

df  <- read_tsv("new_RBD_RSA.tsv")
plot_RSA_vs_param(df, 'Exp_vs_RSA_RBD.png')
print ((paste("Cor between expression score and RSA:", cor(df$mean_exp_score, df$RSA, method='spearman'))))