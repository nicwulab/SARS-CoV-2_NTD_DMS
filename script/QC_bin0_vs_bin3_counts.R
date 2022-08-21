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
library(stringr)
require(cowplot)

plot_bin_0_bin_3_ratio_rep_1 <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=bin0_ratio_rep1,y=bin3_ratio_rep1)) +
    geom_point(size=0.5,pch=16, alpha=0.3) +
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
    labs(x=bquote(bold(paste('Bin 0 ratio (replicate 1)'))),y=bquote(bold(paste('Bin 3 ratio (replicate 1)'))))
  ggsave(graphname, p, height=2.2, width=2.2)
}

plot_bin_0_bin_3_ratio_rep_2 <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=bin0_ratio_rep2,y=bin3_ratio_rep2)) +
    geom_point(size=0.5,pch=16, alpha=0.3) +
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
    labs(x=bquote(bold(paste('Bin 0 ratio (replicate 2)'))),y=bquote(bold(paste('Bin 3 ratio (replicate 2)'))))
  ggsave(graphname, p, height=2.2, width=2.2)
}

df <- read_tsv('result/NTD_DMS_scores.tsv') %>%
  filter(!grepl('silent',mut_class)) %>%
  filter(!grepl('WT',mut_class)) %>%
  filter(avg_total_freq > 0.000075)

plot_bin_0_bin_3_ratio_rep_1(df, "graph/bin0_bin_3_ratio_rep_1_all.png")
plot_bin_0_bin_3_ratio_rep_2(df, "graph/bin0_bin_3_ratio_rep_2_all.png")

df_high_in_rep_1 <- df %>% filter(bin0_bin3_ratio_rep1>0.5)
df_high_in_rep_2 <- df %>% filter(bin0_bin3_ratio_rep2>0.5)

plot_bin_0_bin_3_ratio_rep_1(df_high_in_rep_1, "graph/bin0_bin_3_ratio_rep_1.png")
plot_bin_0_bin_3_ratio_rep_2(df_high_in_rep_2, "graph/bin0_bin_3_ratio_rep_2.png")