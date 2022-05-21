#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(plyr)
library(gridExtra)
library(ggforce)
library(qualpalr)
library(tidyquant)
require(cowplot)

plot_dist_vs_param <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=`dist_to_RBD/S2`,y=mean_exp_score)) +
    geom_point(size=0.5,pch=16, alpha=0.3) +
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
    labs(x=bquote(bold(paste('Distance to RBD/S2 (',ring(A),')',sep=''))),y=bquote(bold(paste('Mutational tolerability'))))
  ggsave(graphname, p, height=2, width=2)
}

plot_position_type <- function(data_table, graphname, ylab,violin){
  summary <- ddply(data_table,c("type"), summarise, mean = mean(mean_exp_score), sd = sd(mean_exp_score))
  textsize <- 7
  if (violin=='yes'){
    p <-  data_table %>%
      arrange(mean_exp_score) %>%   
      mutate(name = factor(type, levels=c("16-26","64-74","142-152","176-187","244-261","Other residues"))) %>% 
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

hot_spot_1 <- seq(16,26)
hot_spot_2 <- seq(64,74)
hot_spot_3 <- seq(142,152)
hot_spot_4 <- seq(176,187)
hot_spot_5 <- seq(244,261)



df <- df %>% mutate (type = "Other residues")
df["type"][df$pos %in% hot_spot_1,] <- '16-26'
df["type"][df$pos %in% hot_spot_2,] <- '64-74'
df["type"][df$pos %in% hot_spot_3,] <- '142-152'
df["type"][df$pos %in% hot_spot_4,] <- '176-187'
df["type"][df$pos %in% hot_spot_5,] <- '244-261'

print (head(df))

plot_position_type(df, "graph/hot_spots_vs_other_residues.png", "Mutational tolerability","yes")
types <- c('16-26','64-74','142-152','176-187','244-261','Other residues')

for (pos_type1 in types){
  p_value <- t.test(filter(df,type=="Other residues")$mean_exp_score, filter(df,type==pos_type1)$mean_exp_score)$p.value
  #p_value <- t.test(filter(fit_table,type==pos_type1)$RSA_tetramer, filter(fit_table,type==pos_type2)$RSA_tetramer)$p.value
  print (paste('Other residues vs', pos_type1, ':', p_value))
}
