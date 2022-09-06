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
library(plyr)
library(ggforce)
require(cowplot)

plot_position_type <- function(data_table, graphname, ylab,violin){
  summary <- ddply(data_table,c("type"), summarise, mean = mean(Exp_score), sd = sd(Exp_score))
  textsize <- 7
  if (violin=='yes'){
    p <-  data_table %>%
      arrange(Exp_score) %>%   
      mutate(name = factor(type, levels=c("Circulating mutations","Other missense mutations"))) %>% 
      ggplot() +
      geom_violin(aes(x=name,y=Exp_score),size=0.3, scale='area', width=0.8)
    ggsave(graphname, p, width=2, height=2.5, dpi=600)
  }
  else{
    p <-  ggplot()
  }
  p <- p +
    geom_sina(data=data_table,aes(x=type,y=Exp_score),
              pch=16, size=0.1,method="counts", bin_limit=0.4, scale="width", maxwidth=0.5, color='gray40', alpha=0.3) +
    geom_boxplot(data=data_table,aes(x=type,y=Exp_score),width=0.2, outlier.shape=NA, size=0.3, color='black',alpha=0.2) +
    theme_cowplot(12) +
    theme(plot.background = element_rect(fill = "white"),
          axis.text=element_text(size=textsize,face="bold",colour = 'black'),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
          axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
          axis.title=element_text(size=textsize,face="bold"),
          legend.title=element_blank(),
          legend.key.height = unit(0.15, 'in'),
          legend.key.width = unit(0.05, 'in'),
          legend.text=element_text(size=textsize,face="bold"),
          legend.position='right') +
    guides(colour = guide_legend(override.aes = list(size=0.5))) +
    xlab("") +
    ylab(ylab)+ scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    theme(axis.text.x = element_text(angle = 90))
  ggsave(graphname, p, width=2, height=2.2, dpi=600)
}

plot_dist_vs_param <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=Exp_score_rep1,y=Exp_score_rep2)) +
    geom_point(size=0.5,pch=16, alpha=0.5) +
    #geom_point(data=not_b, aes(x=`dist_to_RBD/S2`,y=mean_exp_score), color="red", size=0.5,pch=16, alpha=0.7)+
    #geom_point(data=kinda_b, aes(x=`dist_to_RBD/S2`,y=mean_exp_score), color="blue", size=0.5,pch=16, alpha=0.7)+
    #geom_point(data=b1, aes(x=`dist_to_RBD/S2`,y=mean_exp_score), color="green", size=0.5,pch=16, alpha=0.7)+
    #geom_point(data=b2, aes(x=`dist_to_RBD/S2`,y=mean_exp_score), color="orange", size=0.5,pch=16, alpha=0.7)+
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
          scale_y_continuous(limit=c(0,2),breaks=c(0, 0.5, 1, 1.5, 2))+
          scale_x_continuous(limit=c(0,2),breaks=c(0, 0.5, 1, 1.5, 2))+
    labs(x=bquote(bold(paste('Expression score (replicate 1)'))),y=bquote(bold(paste('Expression score (replicate 2)'))))
  ggsave(graphname, p, height=2, width=2)
}

plot_indel_vs_other <- function(data_table, graphname, ylab,violin){
  summary <- ddply(data_table,c("type"), summarise, mean = mean(mean_exp_score), sd = sd(mean_exp_score))
  textsize <- 7
  if (violin=='yes'){
    p <-  data_table %>%
      arrange(mean_exp_score) %>%   
      mutate(name = factor(type, levels=c("In/del sites","Other sites"))) %>% 
      ggplot() +
      geom_violin(aes(x=name,y=mean_exp_score),size=0.3, scale='area', width=0.8)
    ggsave(graphname, p, width=2, height=2.2, dpi=600)
  }
  else{
    p <-  ggplot()
  }
  p <- p +
    geom_sina(data=data_table,aes(x=type,y=mean_exp_score),
              pch=16, size=0.1,method="counts", bin_limit=0.4, scale="width", maxwidth=0.5, color='gray40', alpha=0.5) +
    geom_boxplot(data=data_table,aes(x=type,y=mean_exp_score),width=0.2, outlier.shape=NA, size=0.3,color='black', alpha=0.2) +
    theme_cowplot(12) +
    theme(plot.background = element_rect(fill = "white"),
          axis.text=element_text(size=textsize,face="bold",colour = 'black'),
          axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
          axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
          axis.title=element_text(size=textsize,face="bold"),
          legend.title=element_blank(),
          legend.key.height = unit(0.15, 'in'),
          legend.key.width = unit(0.05, 'in'),
          legend.text=element_text(size=textsize,face="bold"),
          legend.position='right') +
    guides(colour = guide_legend(override.aes = list(size=0.5))) +
    scale_y_continuous(limit=c(0.5,1.5),breaks=c(0.5, 0.75, 1, 1.25, 1.5))+
    xlab("") +
    ylab(ylab)+ scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    theme(axis.text.x = element_text(angle = 90))
  ggsave(graphname, p, width=2, height=2.2, dpi=600)
}

df  <- read_tsv('result/NTD_DMS_scores.tsv')%>%
  filter(avg_total_freq > 0.000075) %>%
  filter(mut_class == "missense") %>%
  mutate(pos=str_sub(mut,2,-2))

print(nrow(df))
NTD_circulating_mutation<- c('L18F', 'T19R','T19I', 'T20N', 'P26S', 'Q52R','A67V', 'G75V', 'T76I', 'D80A', 'T95I', 'D138Y', 'G142D', 'K147E', 'W152R', 'W152C', 'E154K', 'F157L', 'R158G', 'R190S', 'I210V', 'V213G', 'D215G', 'D253G', 'G257S')

NTD_mutation_sites <- c(18,19,20,26,52,67,75,76,80,95,138,142,147,152,154,157,158,190,210,213,215,253,257) #52,80,95,138,157,158,190,210,213,215

print(length(NTD_circulating_mutation))



df <- df %>% mutate (type = "Other missense mutations")
df["type"][df$mut %in% NTD_circulating_mutation,] <- "Circulating mutations"

print (table(df$type))

plot_position_type(df, "graph/ntd_circulating_mutations_vs_other_mutations.png", "Expression scores","yes")
types <- c("Ciculating mutations","Sites of mutation")


df_mut <- df %>% filter(type == "Circulating mutations")
plot_dist_vs_param(df_mut, "graph/ntd_circultating_mutation_correlation.png")

print(cor(df_mut$Exp_score_rep1, df_mut$Exp_score_rep2))

print(t.test(filter(df,type == "Other missense mutations")$Exp_score,filter(df,type == "Circulating mutations")$Exp_score)) 

df_mean <- read_tsv('result/NTD_DMS_scores_by_resi.tsv')%>% filter(count>=6)

NTD_deletion_sites <- c(24,25,26,27,69,70,142,143,144,145,156,157,211,212,214,241,242,243,246,247,248,249,250,251,252)#27,211,212,214,241,242,243

print(length(NTD_deletion_sites))

df_mean <- df_mean %>% mutate(type = "Other sites")
df_mean["type"][df_mean$pos %in% NTD_deletion_sites,] <- "In/del sites"

print(table(df_mean$type))

plot_indel_vs_other(df_mean, "graph/ntd_indel_sites_vs_other_sites.png", "Mutational tolerability","yes")

print(t.test(filter(df_mean,type=="In/del sites")$mean_exp_score,filter(df_mean,type=="Other sites")$mean_exp_score)) 

