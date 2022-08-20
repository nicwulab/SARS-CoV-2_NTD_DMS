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
library(tidyquant)
require(cowplot)

plot_mean_exp_score_by_rep <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=mean_exp_score_rep1,y=mean_exp_score_rep2)) +
    geom_point(size=0.5,pch=16, alpha=0.3) +
    scale_x_continuous(breaks = seq(0.5, 1.5, by = 0.25),limits = c(0.48,1.6))+
    scale_y_continuous(breaks = seq(0.5, 1.5, by = 0.25))+
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
    labs(x=bquote(bold(paste('Mutation tolerability (replicate 1)'))),y=bquote(bold(paste('Mutation tolerability (replicate 2)'))))
  ggsave(graphname, p, height=2.2, width=2.2)
}

df_mean<- read_tsv('result/NTD_DMS_scores_by_resi_and_replicates.tsv')%>%
  filter(count >=6)
print(nrow(df))
plot_mean_exp_score_by_rep(df_mean, 'graph/QC_replicate_mean_exp.png')
print (paste("Correlation between replicates:", cor(df_mean$mean_exp_score_rep1, df_mean$mean_exp_score_rep2, method='spearman')))