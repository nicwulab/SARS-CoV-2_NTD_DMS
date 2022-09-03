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
    scale_x_continuous(limit=c(0,0.6),breaks=c(0, 0.2, 0.4, 0.6))+
    scale_y_continuous(limit=c(0,0.6),breaks=c(0, 0.2, 0.4, 0.6))+
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
    scale_x_continuous(limit=c(0,0.6),breaks=c(0, 0.2, 0.4, 0.6))+
    labs(x=bquote(bold(paste('Bin 0 ratio (replicate 2)'))),y=bquote(bold(paste('Bin 3 ratio (replicate 2)'))))
  ggsave(graphname, p, height=2.2, width=2.2)
}

plot_bin_0_by_replicates <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=bin0_ratio_rep1,y=bin0_ratio_rep2)) +
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
    scale_x_continuous(limit=c(0,0.6),breaks=c(0, 0.2, 0.4, 0.6))+
    labs(x=bquote(bold(paste('Bin 0 ratio (replicate 1)'))),y=bquote(bold(paste('Bin 0 ratio (replicate 2)'))))
  ggsave(graphname, p, height=2.2, width=2.2)
}

plot_bin_3_by_replicates <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=bin3_ratio_rep1,y=bin3_ratio_rep2)) +
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
         scale_x_continuous(limit=c(0,0.6),breaks=c(0, 0.2, 0.4, 0.6))+
         scale_y_continuous(limit=c(0,0.8),breaks=c(0, 0.2, 0.4, 0.6, 0.8))+
    labs(x=bquote(bold(paste('Bin 3 ratio (replicate 1)'))),y=bquote(bold(paste('Bin 3 ratio (replicate 2)'))))
  ggsave(graphname, p, height=2.2, width=2.2)
}

df <- read_tsv('result/NTD_DMS_scores.tsv') %>%
  filter(!grepl('silent',mut_class)) %>%
  filter(!grepl('WT',mut_class)) %>%
  filter(avg_total_freq > 0.000075) %>%
  mutate(bimodality = NA)

for (i in 1 : nrow(df)){
  bin03_count <- c(df[i,]$Expr_bin0_rep1, df[i,]$Expr_bin0_rep2)
  min_bin03_count <- min(bin03_count)
  bin12_count <- c(df[i,]$Expr_bin1_rep1, df[i,]$Expr_bin2_rep2)
  max_bin12_count <- max(bin12_count)
  if (min_bin03_count>max_bin12_count){
    df[i,]$bimodality <- TRUE
  } else {
    df[i,]$bimodality <- FALSE
  }
}
print(table(df$bimodality))

plot_bin_0_bin_3_ratio_rep_1(df, "graph/bin0_bin_3_ratio_rep_1_all.png")
plot_bin_0_bin_3_ratio_rep_2(df, "graph/bin0_bin_3_ratio_rep_2_all.png")

plot_bin_0_by_replicates(df, "graph/bin0_by_replicates.png")
print(cor(df$bin0_ratio_rep1, df$bin0_ratio_rep2, method= c("spearman")))
plot_bin_3_by_replicates(df, "graph/bin3_by_replicates.png")
print(cor(df$bin3_ratio_rep1, df$bin3_ratio_rep2, method= c("spearman")))