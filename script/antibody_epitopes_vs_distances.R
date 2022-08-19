#R code
library(ggplot2)
library(ggforce)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(plyr)
library(gridExtra)
library(qualpalr)
library(tidyquant)
require(cowplot)

plot_dist_vs_param <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=`dist_to_RBD/S2`,y=mean_exp_score)) +
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
    scale_x_continuous(limit=c(0,53),breaks=c(0, 10, 20, 30, 40, 50))+
    labs(x=bquote(bold(paste('Distance to RBD/S2 (Å)'))),y=bquote(bold(paste('Mutational tolerability'))))
  ggsave(graphname, p, height=2, width=2)
}

plot_position_type <- function(data_table, graphname, ylab,violin){
  summary <- ddply(data_table,c("type"), summarise, mean = mean(`dist_to_RBD/S2`), sd = sd(`dist_to_RBD/S2`))
  textsize <- 7
  if (violin=='yes'){
    p <-  data_table %>%
      arrange(`dist_to_RBD/S2`) %>%   
      mutate(name = factor(type, levels=c("C1520","C1717","C1791","Antigenic supersite"))) %>% 
      ggplot() +
      geom_violin(aes(x=name,y=`dist_to_RBD/S2`),size=0.3, scale='area', width=1)
    ggsave(graphname, p, width=2, height=2.2, dpi=600)
  }
  else{
    p <-  ggplot()
  }
  p <- p +
    geom_boxplot(data=data_table,aes(x=type,y=`dist_to_RBD/S2`),width=0.2, outlier.shape=NA, size=0.3) +
    geom_sina(data=data_table,aes(x=type,y=`dist_to_RBD/S2`),
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

df <- read_tsv('result/Dist_NTD_to_RBD_S2.tsv') 
group_1_2_aka_NTD_supersite_resi <- c(seq(14,20),seq(140,158),seq(245,264))
group_3_resi <- c(seq(27,32),seq(57,60),seq(210,218),seq(286,303))
group_4_resi <- c(97,98,99,100,101,102,122,149,152,153,154,155,156,157,158,seq(178,188))
group_5_resi <- c(14,15,16,17,18,19,20,21,22,23,24,25,26,61,122,123,124,125,167,168,169,234)

Antigenic_supersite <- df[df$pos %in% group_1_2_aka_NTD_supersite_resi,] %>% mutate(type = 'Antigenic supersite')
c1717 <- df[df$pos %in% group_3_resi,] %>% mutate(type = 'C1717')
c1520 <- df[df$pos %in% group_4_resi,] %>% mutate(type = 'C1520')
c1791 <- df[df$pos %in% group_5_resi,] %>% mutate(type = 'C1791')
combined_data_set <- full_join(full_join(Antigenic_supersite,c1717), full_join(c1520,c1791))

#plot_dist_vs_param(df, 'Exp_vs_dist_by_anitbody_broadness.png')
#plot_dist_vs_param(df, 'Exp_vs_dist_NTD_supersite_anitbody.png')
#plot_dist_vs_param(df, 'Exp_vs_dist_C1717.png')
#plot_dist_vs_param(df, 'Exp_vs_dist_C1520.png')
#plot_dist_vs_param(df, 'Exp_vs_dist_C1791.png')

plot_position_type(combined_data_set, "graph/antibody_epitopes_vs_distance.png", "Distance to RBD/S2 (Å)","yes")
antibody_types <- c('Antigenic supersite','C1717','C1520','C1791')

  for (pos_type1 in antibody_types){
    p_value <- t.test(filter(combined_data_set,type=="Antigenic supersite")$`dist_to_RBD/S2`, filter(combined_data_set,type==pos_type1)$`dist_to_RBD/S2`)$p.value
    #p_value <- t.test(filter(fit_table,type==pos_type1)$RSA_tetramer, filter(fit_table,type==pos_type2)$RSA_tetramer)$p.value
    print (paste('NTD supersite vs', pos_type1, ':', p_value))
  }