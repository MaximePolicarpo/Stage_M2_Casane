#!/usr/bin/env Rscript

library(ggplot2)
setwd("~/Bureau/Scripts_Maxime/resultats/Run_16_02_2017")


Migration1=read.table(file = "out_simu_stop_1000000_generations_100_repeats_1e-08_mu_500.0_ind_0_h_0_s_1e-06,5e-05,1e-05,0.0005,0.0001_mig_0.0001_matriceIndice.csv", sep=",")
Migration2=read.table(file = "out_simu_stop_1000000_generations_100_repeats_1e-08_mu_500.0_ind_0_h_0_s_1e-06,5e-05,1e-05,0.0005,0.0001_mig_0.0005_matriceIndice.csv", sep=",")
Migration3=read.table(file = "out_simu_stop_1000000_generations_100_repeats_1e-08_mu_500.0_ind_0_h_0_s_1e-06,5e-05,1e-05,0.0005,0.0001_mig_1e-05_matriceIndice.csv", sep=",")
Migration4=read.table(file = "out_simu_stop_1000000_generations_100_repeats_1e-08_mu_500.0_ind_0_h_0_s_1e-06,5e-05,1e-05,0.0005,0.0001_mig_5e-05_matriceIndice.csv", sep=",")
Migration5=read.table(file = "out_simu_stop_1000000_generations_100_repeats_1e-08_mu_500.0_ind_0_h_0_s_1e-06,5e-05,1e-05,0.0005,0.0001_mig_1e-06_matriceIndice.csv", sep=",")


V1=c(1e-04, as.numeric(Migration1[50000,])/100)
V2=c(5e-04, as.numeric(Migration2[50000,])/100)
V3=c(1e-05, as.numeric(Migration3[50000,])/100)
V4=c(5e-05, as.numeric(Migration5[50000,])/100)
V5=c(1e-06, as.numeric(Migration5[50000,])/100)


DF=matrix(NA, ncol=31, nrow=5)
DF=as.data.frame(DF)
DF[1,]=V1
DF[2,]=V2
DF[3,]=V3
DF[4,]=V4
DF[5,]=V5

DF=DF[order(DF[,1]),]

par(xpd=TRUE)
plot(DF[,1], DF[,2], type="p", pch=16, ylim=c(0, 1), ylab=("Probability"), xlab=("Migration rate"), main=("Probabilite de pseudogenisation apres 500 000 generations en fonction du taux de migration"), xlim=c(9e-07, 5.1e-04), cex=1.25, col="red")
points(DF[,1], DF[,3], type="p", pch=17, col="blue")
points(DF[,1], DF[,4], type="p", pch=3, col="green")
points(DF[,1], DF[,5], type="p", pch=4, col="yellow")
points(DF[,1], DF[,6], type="p", pch=5)
points(DF[,1], DF[,7], type="p", pch=8)
points(DF[,1], DF[,8], type="p", pch=19)
points(DF[,1], DF[,9], type="p", pch=13)
points(DF[,1], DF[,10], type="p", pch=23)

lines(DF[,1], DF[,2], type="l", col="red")
lines(DF[,1], DF[,3], type="l", col="blue")
lines(DF[,1], DF[,4], type="l", col="green")
lines(DF[,1], DF[,5], type="l", col="yellow")

legend("right", c("0 pseudogene", "1 pseudogene", "2 pseudogenes", "3 pseudogenes"), pch=c(16, 17, 3, 4), col=c("red", "blue", "green", "yellow"))

#legend("topright", c("0 pseudogene", "1 pseudogene", "2 pseudogenes", "3 pseudogenes", "4 pseudogenes", "5 pseudogenes", "6 pseudogenes", "7 pseudogenes", "8 pseudogenes"), pch=c(16, 17, 3, 4, 5, 8, 19, 13, 23))



#### Meme script mais en mieux #####

#!/usr/bin/env Rscript

#library(ggplot2)
#setwd("~/Bureau/Scripts_Maxime/resultats/Run_16_02_2017")


#matrices=list()
#question=as.integer(readline(prompt="Combien de fichiers matrice allez vous utiliser ?: "))

#Mise de toutes les matrices dans une meme liste

#for (i in seq(1, question, 1)){
#  cat("Veuillez mettre le fichier", i)
#  matrices[[i]]=read.table(file = readline(prompt="Nom du fichier.format(i): "), sep=",")
#}


#On prend la generation souhaitee pour chacune de ces matrices
#generation=as.integer(readline(prompt="Quelle generation voulez vous voir: "))

#i=1
#Vmatrix=matrix(data=NA, nrow=question, ncol=30)
#for (v in seq(1, question, 1)){
#  Vmatrix[v,]=(matrices[[v]][generation,])/100
#}





#DF=matrix(NA, ncol=31, nrow=2)   #ncol -> Nombre de genes +1 et nrow -> Nombre de matrices utilisee 

#i=2
#for(j in matrices) {
#  DF[i,]=j
#  i=i+1
#}

