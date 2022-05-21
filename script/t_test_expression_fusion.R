wildtype_expression<-c(1,1,1,1)
S50Q_expression<-c(1.351052966, 2.431662411, 1.471383975, 1.535767511)
G232E_expression<-c(1.110252096, 0.978228733, 1.300901838, 1.332339791)
DM_expression<-c(2.412623806, 1.839806478, 2.05341658, 2.094076006)
list_expression<- c("S50Q_expression","G232E_expression","DM_expression")
for (i in list_expression){
  p_value<-t.test(wildtype_expression,i)$p.value
  print(paste('wildtype expression vs', i, ':', p_value))
}

wildtype_fusion_3hpm<-c(1,1,1,1)
S50Q_fusion_3hpm<-c(1.03, 1.16, 1.20, 0.89)
G232E_fusion_3hpm<-c(0.860294118,0.692682927,0.775471698,0.74291498)
DM_fusion_3hpm <-c(0.926470588,1.151219512,1.50754717,1.305668016)