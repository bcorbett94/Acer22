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
mlt$File.Name<-gsub("_melt","",as.character(mlt$File.Name))
mlt$File.Name<-gsub("melt/","", as.character(mlt$File.Name))

meltTemp<-select(mlt, Filename = File.Name , Well, Sample.Name,Tm)

df$result<-rename(df$result,Filename = File.Name)

#combines the SYM and CAM averages
mTemp<-meltTemp%>%
  full_join(df$unknowns, by = c("Well","Filename","Sample.Name"))%>%
  full_join(df$result, by = c("Filename","Sample.Name"))

#Are there two values, if not, filter out
#did negative template control amplify

negN<-filter(mTemp,str_detect(Sample.Name, 'neg'))%>% 
  filter(SYM.reps == 2 & CAM.reps == 2)#If negatives comeback with both samples amplifying with target temp,remove

mTemp<-filter(mTemp, !str_detect(Filename, '1214'))%>%#removes 1214 because negatives had amplification
  filter(!str_detect(Sample.Name, 'pos'))

#ugly way to remove the lower case versions of these columns, was throwing off the code
mTemp$Cam.CT.mean<-NULL
mTemp$Sym.CT.mean<-NULL
mTemp$Sym.reps<-NULL
mTemp$Cam.reps<-NULL
mTemp$Sym.CT.sd<-NULL

rerunSamples<-filter(mTemp, SYM.reps != 2 | CAM.reps != 2)%>% #samples that did not work 
  filter(!str_detect(Sample.Name, 'neg|NEC|pos'))

mTemp<-filter(mTemp,SYM.reps ==2 & CAM.reps ==2)#Getting rid of Samples that did not have both wells work

remove(meltTemp)#Removes meltTemp after its not needed downstream

#separates CAM and SYM data based on target to identify thresholds
cBatch<-mTemp%>%
  filter(Target.Name == "CAM")

sBatch<-mTemp%>%
  filter(Target.Name == "SYM")
 
NECs<-filter(sBatch, str_detect(Sample.Name, 'NEC'))%>%
  filter(SYM.reps == 2)%>%
  filter(87.2>Tm | Tm> 88.293)

NECc<-filter(cBatch,str_detect(Sample.Name, 'NEC'))%>%
    filter(CAM.reps == 2)%>%
    filter(74.344>Tm | Tm>75.384)#Range for CAM target

# Checking how far the NEC with CT's are from the other samples based on target type
cmean<-mean(cBatch$CAM.CT.mean)

smean<-mean(sBatch$SYM.CT.mean)

NECc<-NECc%>%
  mutate(checkCAM = CT - cmean)
NECF<-NECs%>%
  mutate(checkSYM = CT-smean)%>%
  full_join(NECc)
checkk<-NECF%>%
  filter(checkSYM<8 | checkCAM<8)

mTemp<-mTemp%>%
  full_join(checkk)#Joins the NECs that had lower than 8 cycle difference 


#converts numeric timepoint in master sheet to character in order to combine them across the qpcr data
master<-master%>%mutate_if(is.numeric,as.character)%>%
    full_join(acer)

mTemp<-mTemp%>%
  filter(!Tm<74 & !Tm>89)# Removes non target temperatures

#select(master,Sample.Name,TimePoint,Cam.CT.mean,Sym.CT.mean,Sym.Cam,Sym.CT.sd,Cam.CT.sd)
final<-select(master,Sample.Name,TimePoint,CAM.CT.mean,SYM.CT.mean,SYM.CAM,SYM.CT.sd,CAM.CT.sd)
final<-transform(final, TimePoint = as.numeric(TimePoint)) #converts Timpoint which was originally a "character" type to numeric type

master<-master %>%
  mutate(Date = as_date(Date,format = "%m/%e/%y"))

#Histograms separating CAM and SYM targets based on thresholds checked against QuantStudios
ggplot(cBatch, aes(Tm, fill = str_detect(Sample.Name, pattern ="NEC")))+
  geom_histogram(binwidth = .1)+
  coord_cartesian(xlim = c(72,77))+
  geom_vline(xintercept= c(74.34,75.45))#threshold for CAM target
  
ggplot(sBatch, aes(Tm, fill = str_detect(Sample.Name, pattern ="NEC")))+
  geom_histogram (binwidth = .1)+
  coord_cartesian(xlim = c(84.5,90))+
  geom_vline(xintercept = c(87.2,88.35))#threshold for SYM target


#Histograms for Extraction control batches to determine if rerunning is neccesary
#NEC 9 
batch<-master%>%
  filter(File.Name == "Acer_1.12-2.txt")

batch9<- master%>%
  filter(File.Name == "Acer_01.18.txt")%>%
  full_join(batch)

#histogram for Cam NEC 9 
ggplot(batch9, aes(x = CAM<.CT.mean, fill = str_detect(Sample.Name, pattern = "NEC")))+
  geom_histogram()+
  scale_x_continuous(breaks = 17:40)
  

#histogram for Sym NEC9   
ggplot(batch9, aes(x = SYM.CT.mean, fill = str_detect(Sample.Name, pattern = "NEC")))+
  geom_histogram()

#histogram for NEC 8
batchq<-master%>%
  filter(File.Name == "Acer_01.12_edit.txt")
batch8<-master%>%
  filter(File.Name == "Acer_1.11_edit.txt")%>%
  full_join(batchq)

ggplot(batch8, aes(x = CAM.CT.mean,fill = str_detect(Sample.Name, pattern = "NEC")))+
  geom_histogram()
ggplot(batch8, aes(x = SYM.CT.mean, fill = str_detect(Sample.Name, pattern = "NEC")))+
  geom_histogram()

ggplot(master, aes(x = Date, y = log10(SYM.CAM), color = Treatment, group = interaction(Date,Treatment)))+
  geom_violin()

str(master)

