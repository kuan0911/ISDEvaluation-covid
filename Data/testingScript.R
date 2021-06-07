#### File Information #####################################################################################################################
#File Name: testingScript.R
#Date Created: June 7, 2021
#Author: Li-Hao Kuan
#Email: lihao@ualberta.ca

### General Comments ######################################################################################################################
#These are some example to use the ISDEvaluation project. The code is not organize. Please refer to the documentation in Humza Haider's ISDEvaluation repo.
#This file keeps the code that I use to preprocessing the data. The code is not automatic and need to use manually.
#############################################################################################################################################

### preprocessing - write your own preprocessing code according to the dataset. Do not use this code because this code is not readable for others.###

survivalDataset <- read.csv(file = "Data/covid/discharge_wo_argentina.csv")

survivalDataset$chronic_disease=NULL
survivalDataset$travel_hist_date=NULL
survivalDataset$travel_hist_location=NULL

survivalDataset$chronic_disease_binary=NULL
survivalDataset$latitude=NULL
survivalDataset$longitude=NULL

survivalDataset$GDP_per_capita_country=NULL
survivalDataset$GDP_total_country=NULL

survivalDataset$population_density_city=NULL
survivalDataset$population_density_country=NULL

survivalDataset$population_density_city=as.numeric(gsub(",", "", survivalDataset$population_density_city, fixed = TRUE))
survivalDataset$population_density_country=as.numeric(gsub(",", "", survivalDataset$population_density_country, fixed = TRUE))

survivalDataset$delta = survivalDataset$event
survivalDataset$event = NULL
survivalDataset[survivalDataset=="False"]=0
survivalDataset[survivalDataset=="True"]=1
survivalDataset$delta = as.integer(survivalDataset$delta)

survivalDataset$translon = cos(survivalDataset$latitude)*cos(survivalDataset$longitude)
survivalDataset$translat = cos(survivalDataset$latitude)*sin(survivalDataset$longitude)

new_df <- survivalDataset
new_df$city <- factor(new_df$city, exclude = NULL)
new_df <- model.matrix(~.-1, data = new_df[c("city")],contrasts.arg = list(city = contrasts(new_df$city, contrasts = FALSE)))
survivalDataset = cbind(survivalDataset,new_df)
survivalDataset$city = NULL
survivalDataset$countries = NULL
### End of preprocessing ###################################################

source('analysisMaster.R')

### Running ISD Evaluation ###
ISD = analysisMaster(survivalDataset, CoxKP = F, MTLRModel=F,KaplanMeier = F,RSFModel=F,GBMModel=T,numberOfFolds = 5)

### Plot survival curves ###
plotSurvivalCurves(ISD$survivalCurves$Cox, 1:10)

### Dcalibration Histogram ###
binnames = c("[0.9,1]","[0.8,0.9)","[0.7,0.8)","[0.6,0.7)","[0.5,0.6)","[0.4,0.5)","[0.3,0.4)","[0.2,0.3)","[0.1,0.2)","[0,0.1)")
barplot(100*ISD$DcalHistogram$Cox/sum(ISD$DcalHistogram$Cox),horiz=TRUE,xlab='Percentage in bin',names.arg=binnames,las=1,main='CoxKP',cex.axis=1.6,cex.names=1.6,cex.lab=1.6,cex.main=2)


