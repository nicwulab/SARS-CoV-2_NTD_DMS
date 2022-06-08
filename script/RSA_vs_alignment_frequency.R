library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
require(cowplot)

frequency_data <- read.csv(file = 'result/residue_freq.csv', header = TRUE)
RSA <- read_tsv("result/NTD_RSA.tsv")
RSA_trimer <-RSA$RSA_trimer
freq <- frequency_data$alignment_frequency
count <- frequency_data$count
data_set <- data.frame(RSA_trimer, freq, count)
filtered_data_set <- data_set %>% filter(count >= 6)
plot_RSA_vs_freq <- function(df, graphname){
  textsize <- 7
  p <- ggplot(df,aes(x=freq, y=RSA_trimer)) +
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
    labs(x=bquote(bold(paste('Alignment frequency'))),y=bquote(bold(paste('Relative solvent accesibility'))))
 # sp<-p+scale_y_continuous(limits=c(0, 1))
  ggsave(graphname, height=2, width=2)
}
plot_RSA_vs_freq(filtered_data_set, "graph/RSA_vs_seq_con.png")
print(cor(filtered_data_set$freq, filtered_data_set$RSA_trimer,method= "spearman",use="complete.obs"))
