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
require(cowplot)

frequency_data <- read.csv(file = 'result/residue_freq.csv', header = TRUE)
scores_by_resi <- read_tsv('result/NTD_DMS_scores_by_resi.tsv')
frequency_data$mean_exp_score <- scores_by_resi$mean_exp_score
filtered_frequency_data <- frequency_data %>% filter(count >= 6)
plot_freq_vs_score <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=alignment_frequency, y=mean_exp_score)) +
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
    labs(x=bquote(bold(paste('Sequence conservation'))),y=bquote(bold(paste('Mutational tolerability'))))
  ggsave(graphname, p, height=2, width=2)
}
plot_freq_vs_score(filtered_frequency_data,"graph/dms_expression_vs_seq_conservation.png")
print(cor(filtered_frequency_data$alignment_frequency, filtered_frequency_data$mean_exp_score,method="spearman", use="complete.obs"))
