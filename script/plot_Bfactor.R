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

plot_dist_vs_param <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=Bfactor,y=mean_exp_score)) +
         geom_point(size=0.5,pch=16, alpha=0.5) +
         #geom_smooth(method = "loess", span=1.5) +
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
         labs(x=bquote(bold(paste('B factor'))),y=bquote(bold(paste('Mean expression score'))))
  ggsave(graphname, p, height=2, width=2)
  }

exp  <- read_tsv('result/NTD_DMS_scores_by_resi.tsv')
Bfac <- read_tsv('result/NTD_Bfactor.tsv')
df <- merge(x=exp, y=Bfac, by='pos') %>%
        filter(count >= 6)
print (paste("Cor between expression score and B factor:", cor(df$mean_exp_score, df$Bfactor, method='spearman')))
plot_dist_vs_param(df, 'graph/Exp_vs_Bfactor.png')
