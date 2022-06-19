library(tidyr)
library(dplyr)
library(readr)

rbd_data<- read.csv("RBD_DMS_data.csv")
rbd_data_autofill<-rbd_data %>% fill(Position..)
rbd_data_autofill<-dplyr::rename(rbd_data_autofill, Pos = Position..)
rbd_data_autofill<-dplyr::rename(rbd_data_autofill, WT_Resi = ï..WT.a.a.)

silent_mutations<-rbd_auto_fill %>% filter (Naive.Freq == "WT")
nonsense_mutations<-rbd_auto_fill %>% filter (Mutation == "*")
silent_average<-(sum(silent_mutations$ACE2.High)+sum(silent_mutations$ACE2.Low))/182
nonsense_average<-(sum(nonsense_mutations$ACE2.High)+sum(nonsense_mutations$ACE2.Low))/182
print(silent_average)
print(nonsense_average)

rbd_data_autofill<-rbd_data_autofill %>% mutate(normalized_exp_score = NA)

for (i in 1:nrow(rbd_data_autofill)){
  sum_of_score<-rbd_data_autofill[i,]$ACE2.High+rbd_data_autofill[i,]$ACE2.Low
  rbd_data_autofill[i,]$normalized_exp_score<-(sum_of_score-nonsense_average)/(silent_average-nonsense_average)
}

mean_exp_score<- c()

for (i in seq(336, 517)){
  each_resi<-filter(rbd_data_autofill, Pos == i)%>% filter(Mutation != "*") %>% filter(Naive.Freq != "WT")
  avg<-sum(each_resi$normalized_exp_score)/19
  mean_exp_score<- c(mean_exp_score, avg)}

pos<-seq(336,517)
exp_scores<-data.frame(pos, mean_exp_score)
combined_data<-inner_join(exp_scores,rsa_data,by="pos")
write_tsv(combined_data,"new_RBD_RSA.tsv")