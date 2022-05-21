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
  summary <- ddply(data_table,c("type"), summarise, mean = mean(mean_exp_score), sd = sd(mean_exp_score))
  textsize <- 7
  if (violin=='yes'){
    p <-  data_table %>%
      arrange(mean_exp_score) %>%   
      mutate(name = factor(type, levels=c("N1","N2","N3","N4","N5","Other residues"))) %>% 
      ggplot() +
      geom_violin(aes(x=name,y=mean_exp_score),size=0.3, scale='area', width=1)
    ggsave(graphname, p, width=2, height=2.2, dpi=600)
  }
  else{
    p <-  ggplot()
  }
  p <- p +
    geom_boxplot(data=data_table,aes(x=type,y=mean_exp_score),width=0.2, outlier.shape=NA, size=0.3) +
    geom_sina(data=data_table,aes(x=type,y=mean_exp_score),
              pch=16, size=0.1,method="counts", bin_limit=0.4, scale="width", maxwidth=0.5, color='gray40', alpha=0.5) +
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
    ylab(ylab)+ scale_x_discrete(labels = function(x) str_wrap(x, width = 10))
  ggsave(graphname, p, width=2, height=2.2, dpi=600)
}

df  <- read_tsv('result/NTD_DMS_scores_by_resi.tsv')%>%
  filter(count >=6)

N1 <- seq(14,26)
N2 <- seq(67,79)
N3 <- seq(141,156)
N4 <- seq(177,186)
N5 <- seq(246,260)



df <- df %>% mutate (type = "Other residues")
df["type"][df$pos %in% N1,] <- 'N1'
df["type"][df$pos %in% N2,] <- 'N2'
df["type"][df$pos %in% N3,] <- 'N3'
df["type"][df$pos %in% N4,] <- 'N4'
df["type"][df$pos %in% N5,] <- 'N5'
print (head(df))

plot_position_type(df, "graph/ntd_loop_vs_other_residues.png", "Mutational tolerability","yes")
types <- c('N1','N2','N3','N4','N5','Other residues')

for (pos_type1 in types){
  p_value <- t.test(filter(df,type=="Other residues")$mean_exp_score, filter(df,type==pos_type1)$mean_exp_score)$p.value
  #p_value <- t.test(filter(fit_table,type==pos_type1)$RSA_tetramer, filter(fit_table,type==pos_type2)$RSA_tetramer)$p.value
  print (paste('Other residues vs', pos_type1, ':', p_value))
}
