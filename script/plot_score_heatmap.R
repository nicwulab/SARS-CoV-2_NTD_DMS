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

plot_score_heatmap <- function(fitness_table, WTresibox, start_resi, end_resi, legend_title){
  textsize <- 5
  fitness_table <- fitness_table %>%
                     filter(Pos >= start_resi & Pos <= end_resi)
  WTresibox     <- WTresibox %>%
                     filter(Pos >= start_resi & Pos <= end_resi) %>%
                     mutate(x=x-min(x)+1)
  p <-  ggplot() +
          geom_tile(data=fitness_table,aes(x=resi,y=aa,fill=parameter)) +
          scale_fill_gradientn(colours=c("blue","blue","white","white","white","red","red"),
                limits=c(-0.6,2.1),
                values=rescale(c(-0.6,0,0.8,1,1.2,2,2.1)),
                breaks=c(0,1,2),
                labels=c('0','1','2'),
                guide="colorbar",
                na.value="grey") +
          theme_cowplot(12) +
          theme(plot.background = element_rect(fill = "white"),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=7,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1)) +
          guides(fill = guide_colorbar(title.theme=element_text(size=7,face="bold",colour='black',hjust=0.5),
                                       label.theme=element_text(size=7,face="bold",colour='black'),
                                       frame.colour="black",
                                       frame.linewidth = 1,
                                       ticks = TRUE,
                                       ticks.colour = "black",
                                       barwidth = 0.5, barheight = 6, title=legend_title)) +
          geom_point(data=WTresibox, aes(x=x, y=y), color='black', size=0.2) +
          xlab("") +
          ylab("amino acid")
  }

aa_level <- rev(c('E','D','R','K','H','Q','N','S','T','P','G','C','A','V','I','L','M','F','Y','W','_'))

start_pos <- 14
end_pos   <- 301

df <- read_tsv('result/NTD_DMS_scores.tsv') %>%
  filter(!grepl('silent',mut_class)) %>%
  filter(!grepl('WT',mut_class)) %>%
  mutate(resi=str_sub(mut,1,-2))
residues <- unique(df$resi)

df <- df %>%
  filter(!grepl('silent',mut_class)) %>%
  filter(!grepl('WT',mut_class)) %>%
  filter(Input_freq >= 0.000075) %>%
  mutate(resi=str_sub(mut,1,-2)) %>%
  mutate(aa=str_sub(mut,-1,-1)) %>%
  filter(aa %in% aa_level) %>%
  mutate(aa=factor(aa,levels=aa_level)) %>%
  mutate(resi=factor(resi,levels=unique(residues))) %>%
  complete(resi, aa) %>%
  mutate(Pos=str_sub(resi,2,-1)) %>%
  mutate(Pos=as.numeric(as.character(Pos))) %>%
  arrange(Pos) %>%
  mutate(parameter=Exp_score) %>% #adjust to select the parameter of interest
  mutate(parameter=case_when(str_sub(resi,1,1)==aa ~ 1, TRUE ~ parameter)) %>% #Set WT parameter (usually 0 or 1)
  mutate(Mutation=paste(resi,aa,sep='')) %>%
  select(Mutation, Input_freq, resi, Pos, aa, parameter) #Variable name for input freq may need to chance

WTresibox  <- df %>%
  select(resi,Pos) %>%
  unique() %>%
  mutate(WT_resi=str_sub(resi,1,1)) %>%
  mutate(x=seq(1,end_pos-start_pos+1)) %>%
  mutate(y=match(WT_resi,aa_level)) %>%
  select(resi,WT_resi,Pos,x, y)

print (range(df$parameter,na.rm=T))

legend_title <- "Exp Score"
p1 <- plot_score_heatmap(df, WTresibox, 14, 85, legend_title)
p2 <- plot_score_heatmap(df, WTresibox, 86, 157, legend_title)
p3 <- plot_score_heatmap(df, WTresibox, 158, 229, legend_title)
p4 <- plot_score_heatmap(df, WTresibox, 230, 301, legend_title)
p <- grid.arrange(p1, p2, p3, p4, nrow=4)
ggsave('graph/NTD_exp_heatmap.png',p,width=5.5, height=7, dpi=300)
