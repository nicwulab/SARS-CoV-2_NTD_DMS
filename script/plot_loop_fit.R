#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(plyr)
library(dplyr)
library(gridExtra)
library(qualpalr)
library(sinaplot)
library(ggforce)
require(cowplot)

plot_by_class <- function(df, graphname, ylab){
  textsize <- 7
  p <- ggplot(df,aes(x=mut_class, y=Exp_score, group=mut_class)) +
         geom_violin(width=1, color="black") +
         geom_sina(pch=16, size=0.1,method="counts", bin_limit=0.4, scale="width", maxwidth=0.5, color='black', alpha=0.2) +
         geom_boxplot(width=0.3, color="black", outlier.shape=NA, alpha=0) + 
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title.x=element_blank(),
               axis.title.y=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
         ylab(ylab)
  ggsave(graphname, p, height=2, width=2)
  }

extract_sites <- function(site_info){
  residues <- c()
  for (resi_range in site_info$residues){
    resi_range <- as.numeric(str_split(resi_range, '-')[[1]])
    residues <- c(residues, seq(resi_range[1], resi_range[2], 1))
    }
  return (residues)
  }

site_info <- read_tsv('data/site_info.tsv')
loop_resi <- extract_sites(site_info)

df <- read_tsv('result/NTD_DMS_scores.tsv') %>%
        filter(Input_freq >= 0.000075) %>%
        filter(mut_class!='WT') %>%
        mutate(pos=as.numeric(str_sub(mut,2,-2))) %>%
        mutate(mut_class=ifelse(mut_class=='missense' & pos %in% loop_resi, 'loop', mut_class))
plot_by_class(df, 'graph/Exp_by_site.png', 'Expression score')
df <- select(df, mut, mut_class, Exp_score)
print (arrange(filter(df, mut_class=='loop'), Exp_score))
