library(steponeR)
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(lubridate)
#importing qPCR data with CT values
p2 <- list.files(path = "qPCR", pattern = ".txt", full.names = TRUE)
p2

#Importing melt temperature files, combining them into one dataframe and turning it into a tibble
tst<-list.files(path = "melt", recursive = TRUE,
                 pattern= ".txt", full.names = TRUE)
mlt<-rbindlist(sapply(tst, fread, simplify = FALSE),
               use.names = TRUE, idcol = "File.Name")
mlt<-as_tibble(mlt)

#renaming Sample to Sample.Name so it matches with other dataframes 
mlt<-mlt%>%
  rename(
    Sample.Name = Sample
  )

#if I wantto look at individiual samples with each Ct, merge melting temps here? 

master<-read.csv(header = TRUE, file ="acer_master.csv")

df <- steponeR(files = p2, 
              delim = "\t", 
              target.ratios = c("SYM.CAM"),
              fluor.norm = list(CAM = 0, SYM = 0), 
              ploidy = list(CAM = 2, SYM = 1),
              extract = list(CAM = .982, SYM = .813),
              copy.number = list(CAM = 1, SYM = 1))

acer <- df$result
#Below code removes melt after and before the melt text files in order to combine dframes
mlt$File.Name<-gsub("_melt","",as.character(mlt$File.Name))#removes "melt" from end of files
mlt$File.Name<-gsub("melt/","", as.character(mlt$File.Name))
meltTemp<-select(mlt, File.Name, Well, Sample.Name,Tm )
#combines the SYM and CAM averages
meltTemp<-meltTemp%>%full_join(acer)

#acer<-acer %>%
#  filter(!Sample.Name== "positive",!Sample.Name == "negative") %>%
  #filter(!Sample.Name == "neg", !Sample.Name == "pos", !Sample.Name == "NEC6")%>%
 # filter(!Sample.Name == "277", !Sample.Name == "279")
#nrow(acer)

#did negative template control amplify(should go back and make all negatives the same name)
neg<- meltTemp %>%
  filter(Sample.Name == "neg")
negative<- meltTemp %>%
  filter(Sample.Name == "negative")%>%
  full_join(neg)
rm(neg)

#Are there two values, if not, filter out
repT<-filter(negative, SYM.reps == 2 |CAM.reps == 2)# only two samples 


#test<-filter(negative,is.na(CAM.CT.mean))



#convets numeric timepoint in master sheet to character in order to combine them across the qpcr data
master<-master%>%mutate_if(is.numeric,as.character)%>%
    full_join(acer)

#select(master,Sample.Name,TimePoint,Cam.CT.mean,Sym.CT.mean,Sym.Cam,Sym.CT.sd,Cam.CT.sd)
final<-select(master,Sample.Name,TimePoint,CAM.CT.mean,SYM.CT.mean,SYM.CAM,SYM.CT.sd,CAM.CT.sd)
final<-transform(final, TimePoint = as.numeric(TimePoint)) #converts Timpoint which was originally a "character" type to numeric type

master<-master %>%
  mutate(Date = as_date(Date,format = "%m/%e/%y"))





#NEC 9 
batch<-master%>%
  filter(File.Name == "Acer_1.12-2.txt")

batch9<- master%>%
  filter(File.Name == "Acer_01.18.txt")%>%
  full_join(batch)

#histogram for Cam NEC 9 
ggplot(batch9, aes(x = Cam.CT.mean, fill = str_detect(Sample.Name, pattern = "NEC")))+
  geom_histogram()+
  scale_x_continuous(breaks = 17:40)
  

#histogram for Sym NEC9   
ggplot(batch9, aes(x = Sym.CT.mean, fill = str_detect(Sample.Name, pattern = "NEC")))+
  geom_histogram()

#histogram for NEC 8
batchq<-master%>%
  filter(File.Name == "Acer_01.12_edit.txt")
batch8<-master%>%
  filter(File.Name == "Acer_1.11_edit.txt")%>%
  full_join(batchq)

ggplot(batch8, aes(x = Cam.CT.mean,fill = str_detect(Sample.Name, pattern = "NEC")))+
  geom_histogram()
ggplot(batch8, aes(x = Sym.CT.mean, fill = str_detect(Sample.Name, pattern = "NEC")))+
  geom_histogram()

ggplot(master, aes(x = Date, y = log10(Sym.Cam), color = Treatment, group = interaction(Date,Treatment)))+
  geom_point()+
  geom_violin()

str(master)

