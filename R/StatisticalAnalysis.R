
#Load packages and source functions 
library(xlsx)
source('./R/BayesFactor_functions.R')
library(tidyverse)

# Read in dataframe with SSP and BPND sheets
SSPBPNDtidy <- read.csv(file = './DerivedData/AllData.csv')

############################# 
# Frequentist correlations
#############################

### whole striatum
#SSP-SD v.s. STR
cor.SD_STR <- cor.test(SSPBPNDtidy$fslSTR,SSPBPNDtidy$SSP_SD,alternative = "greater")
#SSP-PhTA v.s. STR
cor.PhTA_STR <- cor.test(SSPBPNDtidy$fslSTR,SSPBPNDtidy$SSP_PhTA,alternative = "less")

### LST
#SSP-SD v.s. VST
cor.SD_LST <- cor.test(SSPBPNDtidy$fslLST,SSPBPNDtidy$SSP_SD,alternative = "greater")
#SSP-PhTA v.s. VST
cor.PhTA_LST <- cor.test(SSPBPNDtidy$fslLST,SSPBPNDtidy$SSP_PhTA,alternative = "less")

#Create a table of the correlation tests
corTable <- rbind(broom::tidy(cor.SD_STR), 
                  broom::tidy(cor.SD_LST), 
                  broom::tidy(cor.PhTA_STR),
                  broom::tidy(cor.PhTA_LST))
corTable <- data.frame(Scale = c('SD_STR','SD_LST','PhTA_STR','PhTA_LST'),corTable)

################
# Replication BF
################

library(hypergeo)

#SSP-SD v.s. STR
RepBF.SD.STR <- CorrelationReplicationBF(r.orig = 0.54, n.orig = 23-2, r.rep = cor.SD_STR$estimate, n.rep = (cor.SD_STR$parameter+2) ) #smaller N (23-2) to account for loss of 2 df due to covariates
RepBF.SD.LST <- CorrelationReplicationBF(r.orig = 0.52, n.orig = 23-2, r.rep = cor.SD_LST$estimate, n.rep = (cor.SD_LST$parameter+2) )  #smaller N (23-2) to account for loss of 2 df due to covariates

#SSP-PhTA v.s. VST
RepBF.PhTA.STR <- CorrelationReplicationBF(r.orig = -0.36, n.orig = 23-2, r.rep = cor.PhTA_STR$estimate, n.rep = (cor.PhTA_STR$parameter+2) )  #smaller N (23-2) to account for loss of 2 df due to covariates
RepBF.PhTA.LST <- CorrelationReplicationBF(r.orig = -0.51, n.orig = 23-2, r.rep = cor.PhTA_LST$estimate, n.rep = (cor.PhTA_LST$parameter+2) )  #smaller N (23-2) to account for loss of 2 df due to covariates

#Create rep BF table
repBFTable <- rbind(RepBF.SD.STR,
                    RepBF.SD.LST,
                    RepBF.PhTA.STR,
                    RepBF.PhTA.LST)
repBFTable <- data.frame(Scale = c('SD_STR','STR_LST','PhTA_STR','PhTA_LST'),repBFTable,row.names = NULL)

#Write tables
write.csv(corTable, file = './Results/StatisticalAnalysis/FreqCorrelations.csv',row.names = F)
write.csv(repBFTable, file = './Results/StatisticalAnalysis/repBF.csv',row.names = F)

###########
# Table
##########

#Make MD table for publication
Table1<-data.frame( Scale = c('SocDes','SocDes','PhTA','PhTA'),
            Region = c('STR','LST','STR','LST'),
            Pearson.r.orig = c(0.54,0.52,-0.36,-0.51),
            df.orig = c(19,19,19,19),
            p.value.orig = c(0.012,0.015,0.106,0.019),
            Pearson.r.present = c(round(corTable$estimate,2)),
            df.present = c(round(corTable$parameter,2)),
            p.value.present = c(round(corTable$p.value,2)),
            BF01 = round(repBFTable$BF01,1),
            BF10 = round(repBFTable$BF10,2))

library(knitr)
library(kableExtra)
library(xlsx)
library(dplyr)

Table1 <- Table1 %>%
  select(-Scale)

Table1<-kable(Table1, 
      col.names = c('',
                    "r",
                    'df',
                    "p-value",
                    "r",
                    'df',
                    "p-value",
                    'BF01','BF10'),
      align = 'c',
      format = "latex", 
      booktabs = T,  
      caption = "Correlations between SocDes, PhTA and D1-R BPND scores from the original (7) and current study. Replication BFs denotes how much evidence there is in favor of the original correlation compared to no correlation. Note that the correlation between PhTA and STR was not significant in the original study but have still been included here for completeness. ") %>%
  group_rows("SocDes", 1, 2) %>%
  group_rows("PhTA", 3, 4) %>%
  kable_styling(latex_options = c("hold_position"),
                full_width = F,
                position = "center") %>%
  add_header_above(c(" ","Original study[note]" = 3, "Present study [note]" = 3, "Replication BF" = 2)) %>%
  add_footnote(c("two-sided test","one-sided test in direction of original study"))

write.table(as.character(Table1),file = './Results/StatisticalAnalysis/Table1_latex_notation.txt',row.names = F)

############
# Graphics
############

### Make scatter plots
png(filename = './Results/Graphics/ScatterPlots.png',width = 17.5,height = 17.5*0.75, units ='cm',res = 350 )

  par(xaxs="i",yaxs="i") #Make axis connect
  par(mfrow=c(2,2))
  par(mar=c(4,4,1,1) ) #margins
  
  plot(SSPBPNDtidy$fslSTR,SSPBPNDtidy$SSP_SD_Tscore_male,
       xlab = "Striatum", 
       ylab = "Social Desirability",
       cex = 1.1, pch=21,  bg="lightblue",
       xlim = c(1.2,2.2),
       ylim = c(25,75),
       las=1)
  
  #Fit model and create predictions 
  mod<-lm(SSP_SD_Tscore_male ~ fslSTR,data = SSPBPNDtidy)
  newx<-seq(from = min(SSPBPNDtidy$fslSTR) ,to = max(SSPBPNDtidy$fslSTR),length.out = 100)
  prd<-predict(mod,newdata=data.frame(fslSTR=newx),interval = c("confidence"))
  
  #Add reg line and 95% CI
  lines(newx,prd[,1],col="black",lty=1)
  lines(newx,prd[,2],col="black",lty=2)
  lines(newx,prd[,3],col="black",lty=2)
  
  
  plot(SSPBPNDtidy$fslLST,SSPBPNDtidy$SSP_SD_Tscore_male,     
       xlab = "Limbic Striatum", 
       ylab = "Social Desirability",
       cex = 1.1, pch=21,  bg="lightblue",
       xlim = c(0.8,2.1),
       ylim = c(25,75),
       las=1)
  
  #Fit model and create predictions 
  mod<-lm(SSP_SD_Tscore_male ~ fslLST,data = SSPBPNDtidy)
  newx<-seq(from = min(SSPBPNDtidy$fslLST) ,to = max(SSPBPNDtidy$fslLST),length.out = 100)
  prd<-predict(mod,newdata=data.frame(fslLST=newx),interval = c("confidence"))
  
  #Add reg line and 95% CI
  lines(newx,prd[,1],col="black",lty=1)
  lines(newx,prd[,2],col="black",lty=2)
  lines(newx,prd[,3],col="black",lty=2)
  
  
  plot(SSPBPNDtidy$fslSTR,SSPBPNDtidy$SSP_PhTA_Tscore_male,
       xlab = "Striatum", 
       ylab = "Physical Trait Aggression",
       cex = 1.1, pch=21,  bg="lightblue",
       xlim = c(1.2,2.2),
       ylim = c(25,75),
       las=1)
  
  #Fit model and create predictions 
  mod<-lm(SSP_PhTA_Tscore_male ~ fslSTR,data = SSPBPNDtidy)
  newx<-seq(from = min(SSPBPNDtidy$fslSTR) ,to = max(SSPBPNDtidy$fslSTR),length.out = 100)
  prd<-predict(mod,newdata=data.frame(fslSTR=newx),interval = c("confidence"))
  
  #Add reg line and 95% CI
  lines(newx,prd[,1],col="black",lty=1)
  lines(newx,prd[,2],col="black",lty=2)
  lines(newx,prd[,3],col="black",lty=2)
  
  
  plot(SSPBPNDtidy$fslLST,SSPBPNDtidy$SSP_PhTA_Tscore_male,
       xlab = "Limbic Striatum", 
       ylab = "Physical Trait Aggression",
       cex = 1.1, pch=21,  bg="lightblue",
       xlim = c(0.8,2.1),
       ylim = c(25,75),
       las=1)
  
  #Fit model and create predictions 
  mod<-lm(SSP_PhTA_Tscore_male ~ fslLST,data = SSPBPNDtidy)
  newx<-seq(from = min(SSPBPNDtidy$fslLST) ,to = max(SSPBPNDtidy$fslLST),length.out = 100)
  prd<-predict(mod,newdata=data.frame(fslLST=newx),interval = c("confidence"))
  
  #Add reg line and 95% CI
  lines(newx,prd[,1],col="black",lty=1)
  lines(newx,prd[,2],col="black",lty=2)
  lines(newx,prd[,3],col="black",lty=2)
dev.off()

### Make prior-posterior SDR plots
png(filename = './Results/Graphics/repBFplots.png',width = 25,height = 20, units ='cm',res = 350 )
  par(xaxs="i",yaxs="i") #Make axis connect
  par(mfrow=c(2,2))
  par(mar=c(3.2,4,1,1) ) #marginsrepposteriorplot(r.orig = 0.54, n.orig = 23-2, r.rep = cor.SD_STR$estimate, n.rep = (cor.SD_STR$parameter+2),xlab = '')
  
  #SD
  repposteriorplot(r.orig = 0.54, n.orig = 23-2, r.rep = cor.SD_STR$estimate, n.rep = (cor.SD_STR$parameter+2), xlab = 'r: STR vs. SocDes')
  repposteriorplot(r.orig = 0.52, n.orig = 23-2, r.rep = cor.SD_LST$estimate, n.rep = (cor.SD_LST$parameter+2), xlab = 'r: LST vs. SocDes ')
  
  #PhTA
  repposteriorplot(r.orig = -0.36, n.orig = 23-2, r.rep = cor.PhTA_STR$estimate, n.rep = (cor.PhTA_STR$parameter+2) ,xlab = 'r: STR vs. PhTA')
  repposteriorplot(r.orig = -0.51, n.orig = 23-2, r.rep = cor.PhTA_LST$estimate, n.rep = (cor.PhTA_LST$parameter+2), xlab = 'r: LST vs. PhTA')
dev.off()



