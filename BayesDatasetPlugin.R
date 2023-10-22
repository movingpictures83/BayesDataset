


library(vegan)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(scales)
library(stringi)
library(MASS)
#library(stargazer)
library(reporttools)
library(epitools)
library(gdata)
library(car)
library(plyr)
library(dplyr)
library(data.table)
library(tibble)
library(psych)
library(tidyr)
library(janitor)
library(psych)
library(plotrix)
library(slopegraph)
library(Lock5Data)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(treemap)
library (treemapify)
library(ggraph)
library(igraph)



###The code, uses the following "INPUT" files:

dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

input <- function(inputfile) {
  parameters <<- read.table(inputfile, as.is=T);
  rownames(parameters) <<- parameters[,1];
    pfix <<- prefix()
  if (length(pfix) != 0) {
     pfix <<- paste(pfix, "/", sep="")
  }

Res_Score_Two <<- readRDS(paste(pfix, parameters["score2", 2], sep="/"))
Res_Score_One_Two <<- readRDS(paste(pfix, parameters["score12", 2], sep="/"))
ARO_Counts <<- readRDS(paste(pfix, parameters["counts", 2], sep="/"))
Tx1 <<- read.csv(paste(pfix, parameters["reads", 2], sep="/"), sep="\t", header=T, check.names=F, row.names=1, quote="#")
Inf1 <<- read.csv(paste(pfix, parameters["infection", 2], sep="/"), header=T)

#Res_Score_Two <- readRDS("RDS/Res_Score_Two.rds")
#Res_Score_One_Two <- readRDS("RDS/Res_Score_One_Two.rds")
#ARO_Counts <- readRDS("RDS/ARO_Counts.rds")
#Tx1<-read.csv("inputs/bracken_combined_reads.tsv", sep = "\t",header=T,check.names = F, row.names = 1, quote = "#")
#Inf1<-read.csv("inputs/Infection_Data_For_Bayesian_Model.csv", header=T)
}

run <- function() {}

output <- function(outputfile) {
##################CREATE THE DATASET FOR BAYESIAN MODELLING################################################################################################################################

Rcgc12<-Res_Score_One_Two

Rcgc12<-Rcgc12[,c("Antibiotic","CAMBODIA_POP_POOL", "KENYA_POP_POOL", "UK_POP_POOL")]

Rcgc12$Antibiotic <- gsub('AB_', '', Rcgc12$Antibiotic)

Tot_AMR12<-Rcgc12 %>%
  adorn_totals("row") %>%
  mutate_at(vars(-Antibiotic), list(~replace(., row_number() < n(),
                                       (.[-n()]/.[n()])))) 

Tot_AMR12<- Tot_AMR12 %>%
  filter(Antibiotic=="amikacin" | Antibiotic=="ampicillin" | Antibiotic=="cefotaxime" | Antibiotic=="cefoxitin" | Antibiotic=="cefpodoxime" | 
           Antibiotic=="ceftazidime" | Antibiotic=="ceftriaxone" | Antibiotic=="cefuroxime" | Antibiotic=="chloramphenicol" | Antibiotic=="ciprofloxacin" | 
           Antibiotic=="Gentamicin_All" | Antibiotic=="imipenem" | Antibiotic=="meropenem" | Antibiotic=="nalidixicacid" | Antibiotic=="nitrofurantoin" | Antibiotic=="trimethoprimsulfamethoxazole") %>% 
  droplevels()

long_Tot_AMR12 <- Tot_AMR12 %>% gather(Setting, Rcgc_1_2,2:4)

Rcgc2<-Res_Score_Two

Rcgc2<-Rcgc2[,c("Antibiotic","CAMBODIA_POP_POOL", "KENYA_POP_POOL", "UK_POP_POOL")]

Rcgc2$Antibiotic <- gsub('AB_', '', Rcgc2$Antibiotic)

Tot_AMR2<-Rcgc2 %>%
  adorn_totals("row") %>%
  mutate_at(vars(-Antibiotic), list(~replace(., row_number() < n(),
                                             (.[-n()]/.[n()])))) 

Tot_AMR2<- Tot_AMR2 %>%
  filter(Antibiotic=="amikacin" | Antibiotic=="ampicillin" | Antibiotic=="cefotaxime" | Antibiotic=="cefoxitin" | Antibiotic=="cefpodoxime" | 
           Antibiotic=="ceftazidime" | Antibiotic=="ceftriaxone" | Antibiotic=="cefuroxime" | Antibiotic=="chloramphenicol" | Antibiotic=="ciprofloxacin" | 
           Antibiotic=="Gentamicin_All" | Antibiotic=="imipenem" | Antibiotic=="meropenem" | Antibiotic=="nalidixicacid" | Antibiotic=="nitrofurantoin" | Antibiotic=="trimethoprimsulfamethoxazole") %>% 
  droplevels()

long_Tot_AMR2 <- Tot_AMR2 %>% gather(Setting, Rcgc_2,2:4)

ByMo1<-merge(long_Tot_AMR12,long_Tot_AMR2,by=c("Antibiotic","Setting"),all=TRUE)


ByMo2<-ByMo1 %>% 
  mutate(Rcgc2_zero_replaced = case_when(Rcgc_2==0 ~ Rcgc_1_2, is.na(Rcgc_2)~Rcgc_1_2, TRUE ~ Rcgc_2))

ByMo2$Rcgc_2[is.na(ByMo2$Rcgc_2)] <- 0


##attach(Inf1)
##fix(Inf1)

ByMo3<-merge(ByMo2,Inf1,by=c("Antibiotic","Setting"),all=TRUE)



##attach(Tx1)
##fix(Tx1)

Tx2<-Tx1[,grepl("POOL",colnames(Tx1))]
Tx2 <- add_rownames(Tx2, "Taxonomy")
Tx3<- Tx2 %>%
  filter(grepl("k__Bacteria",Taxonomy)) %>% 
  droplevels()

d3<-filter(Tx3 [,2:7])
d3<-d3 %>% replace(is.na(.), 0)
AllBacteria<-as.data.frame(t(colSums(d3 [,])))


Tx3<-Tx3 %>% replace(is.na(.), 0)
ecol<-Tx3[grep("g__Escherichia;s__Escherichia coli", Tx3$Taxonomy), ]
ecol<-filter(ecol [,2:7])
kle<-Tx3[grep("g__Klebsiella;s__Klebsiella pneumoniae", Tx3$Taxonomy), ]
kle<-filter(kle [,2:7])

eb<-Tx3[grep("g__Enterobacter;", Tx3$Taxonomy), ]
eb<-filter(eb [,2:7])
eb<-as.data.frame(t(colSums(eb [,])))
sal<-Tx3[grep("g__Salmonella;", Tx3$Taxonomy), ]
sal<-filter(sal [,2:7])
sal<-as.data.frame(t(colSums(sal [,])))
f_E<-Tx3[grep("f__Enterobacteriaceae;", Tx3$Taxonomy), ]
f_E<-filter(f_E [,2:7])
f_E<-as.data.frame(t(colSums(f_E [,])))
o_E<-Tx3[grep("o__Enterobacterales;", Tx3$Taxonomy), ]
o_E<-filter(o_E [,2:7])
o_E<-as.data.frame(t(colSums(o_E [,])))

ClinicalGroups<-rbind(ecol,kle,eb,sal)
SumClinicalGroups<-as.data.frame(t(colSums(ClinicalGroups)))
Rtaxcolkleebsal<-SumClinicalGroups/AllBacteria
Rtaxf_E<-f_E/AllBacteria
Rtaxo_E<-o_E/AllBacteria

Rtaxcolkleebsal<-t(Rtaxcolkleebsal)
Rtaxcolkleebsal<-as.data.frame(Rtaxcolkleebsal)
Rtaxcolkleebsal<- Rtaxcolkleebsal %>% rownames_to_column("SampleID1")
colnames(Rtaxcolkleebsal)[2]<-"Rtax_colkleebsal_Bracken"

Rtaxf_E<-t(Rtaxf_E)
Rtaxf_E<-as.data.frame(Rtaxf_E)
Rtaxf_E<- Rtaxf_E %>% rownames_to_column("SampleID2")
colnames(Rtaxf_E)[2]<-"Rtax_f_E_Bracken"

Rtaxo_E<-t(Rtaxo_E)
Rtaxo_E<-as.data.frame(Rtaxo_E)
Rtaxo_E<- Rtaxo_E %>% rownames_to_column("SampleID3")
colnames(Rtaxo_E)[2]<-"Rtax_o_E_Bracken"


dataBracken<-cbind(Rtaxcolkleebsal,Rtaxf_E,Rtaxo_E)

names(dataBracken)[names(dataBracken) == "SampleID1"] <- "Setting"
dataBracken$SampleID2<-NULL
dataBracken$SampleID3<-NULL

dataBracken<- dataBracken %>%
  filter(grepl("POP_POOL",Setting)) %>% 
  droplevels()


ByMo4<-merge(ByMo3,dataBracken,by="Setting", all=TRUE)

ByMo4<-ByMo4[,c(1,2,3,4,5,8,9,10,6,7)]

#ARCHIVE DATA FOR BAYESIAN MODEL

#write.csv(ByMo4, "outputs/Dataset_For_Bayesian_Model.csv",row.names=F)
write.csv(ByMo4, outputfile,row.names=F)
}


#########################################################non-metric multidimensional scaling (NMDS) ############################

