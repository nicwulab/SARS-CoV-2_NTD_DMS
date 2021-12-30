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
  p <- ggplot(df,aes(x=RSA_trimer, y=mean_exp_score)) +
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

resi_classification <- function(RSA_monomer, delta_RSA){
  if (RSA_monomer < 0.05){return ('Buried')}
  else if (delta_RSA > 0.5*RSA_monomer){return ('Interface')}
  else (return ('Exposed'))
  }

exp <- read_tsv('result/NTD_DMS_scores_by_resi.tsv')
rsa <- read_tsv('result/NTD_RSA.tsv')
df  <- merge(x=exp, y=rsa, by = "pos", all=TRUE) %>%
         filter(count >= 6) %>%
         mutate(resi_type=mapply(resi_classification, RSA_monomer, delta_RSA)) %>%
         mutate(resi_type=factor(resi_type, levels=c('Buried','Interface','Exposed')))
plot_RSA_vs_param(df, 'graph/Exp_vs_RSA.png')
print ((paste("Cor between expression score and RSA:", cor(df$mean_exp_score, df$RSA_trimer, method='spearman'))))
print (filter(df, resi_type == 'Buried' & mean_exp_score > 1))
