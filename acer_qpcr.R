library(steponeR)
install.packages("plotly")
library(tidyverse)
library(dplyr)
library(plotly)
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
  filter(SYM.reps == 2 | CAM.reps == 2)#If negatives comeback with both samples amplifying with target temp,remove

mTemp<-filter(mTemp, !str_detect(Filename, '1214'))%>%#removes 1214 because negatives had amplification
  filter(!str_detect(Sample.Name, 'pos'))

#ugly way to remove the lower case versions of these columns, was throwing off the code
mTemp$Cam.CT.mean<-NULL
mTemp$Sym.CT.mean<-NULL
mTemp$Sym.reps<-NULL
mTemp$Cam.reps<-NULL
mTemp$Sym.CT.sd<-NULL
mTemp$Cam.CT.sd<-NULL

rerunSamples<-filter(mTemp, SYM.reps != 2 | CAM.reps != 2)%>% #samples that did not work 
  filter(!str_detect(Sample.Name, 'neg|NEC|pos'))

rerunDupes<-rerunSamples%>%
  count(Sample.Name) %>%
  filter(n>=2)
rerunDupes



batch11<-mTemp%>%
  filter(Filename == "acer_20230215_1.txt")

mTemp<-filter(mTemp,SYM.reps ==2 & CAM.reps ==2)#Dataframe with samples that amplified twice 

dupesamp <- mTemp %>%
  count(Sample.Name) %>%
  filter(n>=2)
dupesamp

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

mTemp<-mTemp%>%
  filter(!Tm>75.384 | !88.3<Tm)%>%
  filter(!Tm<74.344 | !87.2>Tm)# Removes non target temperatures

#finalMtemp<-select(mTemp,Sample.Name,CAM.CT.mean,SYM.CT.mean,SYM.CAM,CAM.CT.mean,SYM.CT.mean,SYM.CAM )

#converts numeric timepoint in master sheet to character in order to combine them across the qpcr data
master2<-master%>%mutate_if(is.numeric,as.character)%>%
    full_join(mTemp)%>%
    filter(!str_detect(Sample.Name, 'NEC|pos|neg'))%>%
    filter(!str_detect(Filename, '1221'))#NEC 7 amplified w/ melt temps, removes this batch

#needing to know number of samples are good by removing duplicates temporarily
dupesamp <- mTemp %>%
  count(Sample.Name) %>%
  filter(n>=2)
dupesamp

#select(master,Sample.Name,TimePoint,Cam.CT.mean,Sym.CT.mean,Sym.Cam,Sym.CT.sd,Cam.CT.sd)
final<-select(master2,Sample.Name,TimePoint,CAM.CT.mean,SYM.CT.mean,SYM.CAM,SYM.CT.sd,CAM.CT.sd)
final<-transform(final, TimePoint = as.numeric(TimePoint)) #converts Timpoint which was originally a "character" type to numeric type

master2<-master2 %>%
  mutate(Date = as_date(Date,format = "%m/%e/%y"))

#Histograms separating CAM and SYM targets based on thresholds checked against QuantStudios
ggplot(cBatch, aes(Tm, fill = str_detect(Sample.Name, pattern ="NEC")))+
  geom_histogram(binwidth = .1)+
  coord_cartesian(xlim = c(72,77))+
  geom_vline(xintercept= c(74.34,75.409))#threshold for CAM target
  
ggplot(sBatch, aes(Tm, fill = str_detect(Sample.Name, pattern ="NEC")))+
  geom_histogram (binwidth = .1)+
  coord_cartesian(xlim = c(84.5,90))+
  geom_vline(xintercept = c(87.2,88.3))#threshold for SYM target


#Data separated in batches to determine how many cycles away samples are from NECs that amplified twice
#columns symchk and camchk reflect differences based on NEC values

batch9<-mTemp%>%
  filter(Filename == "acer_20230118.txt")%>%
  mutate(cam_chk = 32.430-CAM.CT.mean,sym_chk = 28.70-SYM.CT.mean)%>%
  filter(!cam_chk<8 & !sym_chk<8)

bb9<-mTemp%>%
  filter(Filename == "acer_20230118.txt")%>%
  filter(!str_detect(Sample.Name, 'NEC'))%>%
  mutate(cam_chk = 32.43 -CAM.CT.mean, sym_chk = 28.7-SYM.CT.mean)%>%
  filter(cam_chk<8 | sym_chk<8)

#Batch8 filtering the list of redos 
batch8<-mTemp%>%
  filter(Filename == "acer_20230111.txt")%>%
  mutate(cam_chk = 31.79 -CAM.CT.mean, sym_chk = 35.65-SYM.CT.mean)%>%
  filter(!cam_chk<8 & !sym_chk<8)

bb8<-mTemp%>%
  filter(Filename == "acer_20230111.txt")%>%
  filter(!str_detect(Sample.Name, 'NEC'))%>%
  mutate(cam_chk = 31.79 -CAM.CT.mean, sym_chk = 35.65-SYM.CT.mean)%>%
  filter(cam_chk<8 | sym_chk<8)

#Batch10 filtering the list of redos 
#Finds the highest 
batch10<-mTemp%>%
  filter(Filename == "acer_20230119.txt")%>%
  mutate(cam_chk = 31.827 -CAM.CT.mean, sym_chk = 26.762-SYM.CT.mean)%>%
  filter(!cam_chk<8 & !sym_chk<8)


bb10<-mTemp%>%
  filter(Filename == "acer_20230119.txt")%>%
  filter(!str_detect(Sample.Name, 'NEC'))%>%
  mutate(cam_chk = 31.79 -CAM.CT.mean, sym_chk = 35.65-SYM.CT.mean)%>%
  filter(cam_chk<8 | sym_chk<8)

batch11<-mTemp%>%
  filter(Filename == "acer_20230215_1.txt")%>%
  mutate(cam_chk = 31.42 -CAM.CT.mean, sym_chk = 34.883-SYM.CT.mean)%>%
  filter(!cam_chk<8 & !sym_chk<8)

bb11<-mTemp%>%
  filter(Filename == "acer_20230215_1.txt")%>%
  filter(!str_detect(Sample.Name, 'NEC'))%>%
  mutate(cam_chk = 31.42 -CAM.CT.mean, sym_chk = 34.883-SYM.CT.mean)%>%
  filter(cam_chk<8 | sym_chk<8)

fbatch<-batch10%>%
  full_join(batch9)%>%
  full_join(batch8)

remove(batch10,batch8, batch9,mlt)


#histogram for Cam NEC 9 
ggplot(fbatch, aes(x = CAM<.CT.mean, fill = str_detect(Sample.Name, pattern = "NEC9")))+
  geom_histogram()+
  scale_x_continuous(breaks = 17:40)

#histogram for Sym NEC9   
ggplot(fbatch, aes(x = SYM.CT.mean, fill = str_detect(Sample.Name, pattern = "NEC9")))+
  geom_histogram()

ggplot(fbatch, aes(x = CAM.CT.mean,fill = str_detect(Sample.Name, pattern = "NEC8")))+
  geom_histogram()
ggplot(fbatch, aes(x = SYM.CT.mean, fill = str_detect(Sample.Name, pattern = "NEC8")))+
  geom_histogram()


#filters out melt temperature to (sm)  be combined with timepoint and treatment in master2

sm<-distinct(master2, Sample.Name,Sample.Plate,CT,Target.Name)

master<-master%>%mutate_if(is.numeric,as.character)
sm<-sm%>%
  group_by(Target.Name, Sample.Name,Sample.Plate)%>%
  summarise(Ct.Mean = mean(CT))

final_frame<-full_join(sm,master, by = "Sample.Name")

final_frame<-final_frame%>%pivot_wider(names_from = Target.Name, values_from = Ct.Mean)%>%
  mutate(SYM.CAM = ((2^(CAM-SYM))*2)/.828)



q<-final_frame%>%
  group_by(TimePoint,SYM.CAM)%>%
  mutate(nlog = log10(SYM.CAM))


summary(q$SYM.CAM)


p<-ggplot(final_frame, aes(x = Date, y = log10(SYM.CAM), color = Treatment, group = interaction(Date,Treatment)))+
  geom_boxplot()
p
ggplotly(p)

str(master)

count3<-final_frame%>%
  filter(is.na(Sample.Plate) == TRUE)





