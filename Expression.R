---
title: "Expression"
output: html_notebook
---

###**Libraries used**
```{r Libraries used}
library(ggplot2)
library(ggpubr)
library(edgeR)
library(reshape2)
library(tidyverse)
library(dplyr)
```

###**Prepare clinical data table**
```{r Clinical data table}
setwd("~/Desktop/Data4Analysis-selected/")

Insulin_Phase2<-read.table('Insulin_Phase2.txt',header=T,sep='\t')

load('clinical_data_2018_04_06.RData')
cohort<-read.table('Phase2_BMI.txt',header=T)
cohort$colnames<-paste(cohort$Patnr,cohort$Time,sep='') 
```


###**Add new variables and columns**
```{r}
cohort$Type<-ifelse(cohort$substudie=='KONTROLL','NO','OB')
cohort$newcond<-paste(cohort$Type,cohort$Time,sep='.')

#Replacing f2, h2 w POB subjects
cohort<-cohort %>%
  mutate(newcond = str_replace(newcond, 'OB.f2', 'POB.f0')) %>%
  mutate(newcond = str_replace(newcond, 'OB.h2', 'POB.h0'))

cohort <- cohort %>%
  mutate(Type = case_when(
    newcond == "POB.f0" ~ "POB",
    newcond == "POB.h0" ~ "POB",
    newcond == "NO.f0" ~ "NO",
    newcond == "NO.h0" ~ "NO",
    newcond == "OB.f0" ~ "OB",
    newcond == "OB.h0" ~ "OB"
  ))
```


###**Extract only female subjects**
```{r}
female<- clinical[clinical$Gender=='K',]$Patnr
female.df<- clinical[clinical$Gender=='K',]

cohort.female<-cohort[cohort$Patnr %in% female,] 
Insulin_Phase2.female<-Insulin_Phase2[,colnames(Insulin_Phase2) %in% cohort.female$colnames]

rownames(Insulin_Phase2.female) <- Insulin_Phase2[,'TC']
```


###**Remove subjects with a BMI > 30 after weight-loss**
```{r}
femaleBMI <- female.df[female.df$BMI > 30,]
femaleBMI.2high <- femaleBMI[femaleBMI$Time == 'P2',]$Patnr 
femaleBMI.2high  #Patnr used to construct vector of column indices below
BMI.2high <- as.character(femaleBMI.2high)
cohort.female <- cohort.female[!cohort.female$Patnr %in% BMI.2high,]

ix1 <- which(colnames(Insulin_Phase2.female) %in% c("NO10f0","NO10h0", "NO10f2","NO10h2", "NO16f0","NO16h0","NO16f2","NO16h2","NO17f0","NO17h0","NO17f2","NO17h2", "NO5f0", "NO5h0","NO5f2","NO5h2", "NO53f0","NO53h0","NO53f2","NO53h2", "NO54f0", "NO54h0", "NO54f2", "NO54h2"))
clean1 <- Insulin_Phase2.female[, -ix1] 
```


###**Remove outliers**
```{r}
outliers <- c("NG8","NG10","NK5","NO119","NO144","NO138","NO41","NO66")
cohort.female<-cohort.female[!cohort.female$Patnr %in% outliers, ]

ix2 <- which(colnames(clean1) %in% c("NG8f0", "NG8h0", "NG10f0","NG10h0", "NK5f0","NK5h0", "NO119f0","NO119h0","NO119f2","NO119h2", "NO144f0","NO144h0","NO144f2","NO144h2", "NO138f0","NO138h0","NO138f2","NO138h2", "NO41f0","NO41h0","NO41f2","NO41h2", "NO66f0", "NO66h0","NO66f2","NO66h2"))
clean2 <- clean1[, -ix2]
```


###**Normalize table w cpm, remove top 10 highly expressed TC**
```{r}
CpmTbl <- cpm(clean2)    
RS_CpmTbl <- rowSums(CpmTbl)
  
#Make a dataframe with rowsum cpm expressions and TCs
dfA <- data.frame(TC = Insulin_Phase2[,1], rowSums = c(RS_CpmTbl))

positions <- order(dfA$rowSums, decreasing = T)
dfA.sorted <- dfA[positions,]
top10<-head(dfA.sorted,10)
he_TCnames <- top10$TC
  
#Make a new table w/o corresponding rows from CpmTbl
ix3 <- which(Insulin_Phase2$TC %in% he_TCnames)
CpmTblA <- CpmTbl[-ix3, ]
Counts <- clean2[-ix3, ] 
```


###**Remove rows where <20 % of samples have ≥1 cpm**
```{r}
#Make a 0-matrix with dim of CpmTblA
mat <- matrix(data = 0, nrow = nrow(CpmTblA), ncol = ncol(CpmTblA))

#Make for-loop to assign 1 to the corresponding element in mat if cpm > 1
for (i in 1:nrow(CpmTblA)) {
  for (j in 1:ncol(CpmTblA)) {
    if (CpmTblA[i,j] > 1) {
      mat[i,j] <- 1
    }
  }
}
  
#Make for loop and vector of row-indices where < 20% of samples have ≥1 cpm 
q <- c()
for (i in 1:nrow(CpmTblA)) {
  if (sum(mat[i,])/ncol(mat) < 0.2) {
    q <- c(q, i)
  }
}
CpmTblB <- CpmTblA[-q,]
Counts2 <- Counts[-q,]
```

###**Make boxplot w subgroups**
```{r}
NO.f0<-cohort.female[cohort.female$newcond=='NO.f0',]$colnames
NO.h0<-cohort.female[cohort.female$newcond=='NO.h0',]$colnames
OB.f0<-cohort.female[cohort.female$newcond=='OB.f0',]$colnames
OB.h0<-cohort.female[cohort.female$newcond=='OB.h0',]$colnames
POB.f0<-cohort.female[cohort.female$newcond=='POB.f0',]$colnames
POB.h0<-cohort.female[cohort.female$newcond=='POB.h0',]$colnames
  
Groups <- cohort.female$newcond
log2CPM <- c(log2(CpmTblB[1,NO.f0]), log2(CpmTblB[1, NO.h0]), log2(CpmTblB[1, OB.f0]), log2(CpmTblB[1, OB.h0]), log2(CpmTblB[1, POB.f0]), log2(CpmTblB[1, POB.h0]))
Types <- cohort.female$Type
dfggp <- data.frame(Groups, log2CPM, Types)

ggplot(dfggp, aes(x = Groups, y = log2CPM, color = Types)) + geom_boxplot()
```


