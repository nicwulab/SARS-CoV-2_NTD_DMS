library(tidyr)
library(dplyr)
library(readr)

rbd_data <- read.csv("data/RBD_DMS.csv")
rbd_data_autofill <- rbd_data %>% fill(Position)
rbd_data_autofill <-rename(rbd_data_autofill, Pos = Position)
rbd_data_autofill <-rename(rbd_data_autofill, WT_Resi = WTaa)

mean_exp_score<- c()
for (i in seq(336, 517)){
  each_resi<-filter(rbd_data_autofill, Pos == i)
  avg <- (sum(each_resi$ACE2.High)+sum(each_resi$ACE2.Low))/20
  mean_exp_score<- c(mean_exp_score, avg)}

pos<-seq(336,517)
exp_scores<-data.frame(pos, mean_exp_score)
write_tsv(exp_scores,"data/RBD_DMS_exp.tsv")
