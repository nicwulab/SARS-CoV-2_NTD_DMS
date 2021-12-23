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
library(viridis)
library(qualpalr)
require(cowplot)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
  }

plot_hist <- function(t, main_title){
  textsize <- 7
  p <- ggplot(t,aes(x=fit)) +
         geom_histogram(binwidth=0.1) +
         theme_cowplot(12) +
         theme(plot.title=element_text(size=textsize,face="bold", hjust=0.5),
               axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_blank(),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
     ggtitle(main_title) +
     scale_fill_manual(values=c('black'),drop=FALSE) +
     labs(y=expression(bold('count')),x=expression(bold('fitness'))) +
     coord_cartesian(xlim=c(-2.5,1))
  return (p)
  }

plot_replicate_cor <- function(df, graphname, param){
  print (paste('correlation for:', graphname, cor(df$rep1, df$rep2)))
  textsize <- 7
  df$density <- get_density(df$rep1, df$rep2, n = 100)
  p <- ggplot(df,aes(x=rep1, y=rep2, color=density)) +
         geom_hex(bins = 70) +
         scale_fill_continuous(type = "viridis") +
         theme_cowplot(12) +
         theme(plot.title=element_blank(),
               plot.background = element_rect(fill = "white"),
               axis.title=element_text(size=textsize,face="bold"),
               axis.text=element_text(size=textsize,face="bold"),
               legend.key.size=unit(0.1,'in'),
               legend.spacing.x=unit(0.03, 'in'),
               legend.title=element_text(size=textsize,face="bold"),
               legend.text=element_text(size=textsize,face="bold"),
               legend.position='right') +
         labs(x=bquote(bold(paste(.(param),' (replicate 1)'))),y=bquote(bold(paste(.(param),' (replicate 2)'))))
  ggsave(graphname, p, height=2, width=2.5)
  }

mut_classification <- function(mut_class, resi){
  if (mut_class=='missense' & resi=='W64'){return ('W64X')}
  else(return (mut_class))
  }

plot_exp_vs_fus <- function(df, graphname){
  print (paste('correlation for:', graphname, cor(df$Exp_score, df$Fus_score)))
  textsize <- 7
  palette <- qualpal(n = 3, list(h = c(0, 360), s = c(0.4, 0.6), l = c(0.5, 0.85)))$hex
  df <- df %>%
          filter(mut_class != 'WT') %>%
          mutate(resi=str_sub(mut,1,-2)) #%>%
          #mutate(mut_class=mapply(mut_classification, mut_class, resi)) %>%
          #mutate(mut_class=factor(mut_class, levels=c('missense','nonsense','silent','W64X')))
  p <- ggplot(df,aes(x=Exp_score, y=Fus_score, color=mut_class)) +
         geom_point(size=0.1,pch=16, alpha=0.5) +
         scale_color_manual(values=palette) +
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
         guides(color = guide_legend(override.aes = list(size = 1))) +
         labs(x=bquote(bold(paste('Expression score'))),y=bquote(bold(paste('Fusion score'))))
  ggsave(graphname, p, height=2, width=2.5)
  }

plot_by_class <- function(df, graphname, ylab){
  df <- df %>%
          filter(mut_class != 'WT')
  textsize <- 7
  p <- ggplot(df,aes(x=mut_class, y=score, group=mut_class)) +
         geom_violin(width=1, color="black") +
         geom_boxplot(width=0.3, color="black", outlier.shape=NA) + 
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

df <- read_tsv('result/NTD_DMS_scores.tsv') %>%
        filter(Input_freq >= 0.000075)
print (nrow(df))
df_exp <- df %>%
            rename(rep1=Exp_score_rep1) %>%
            rename(rep2=Exp_score_rep2) %>%
            rename(score=Exp_score)
df_fus <- df %>%
            rename(rep1=Fus_score_rep1) %>%
            rename(rep2=Fus_score_rep2) %>%
            rename(score=Fus_score)
plot_replicate_cor(df_exp, 'graph/QC_replicate_exp.png', "Expression score")
plot_replicate_cor(df_fus, 'graph/QC_replicate_fus.png', "Fusion score")
plot_exp_vs_fus(df, 'graph/Exp_vs_fus.png')
plot_by_class(df_exp, 'graph/Exp_by_class.png', 'Expression score')
plot_by_class(df_fus, 'graph/Fus_by_class.png', 'Fusion score')
