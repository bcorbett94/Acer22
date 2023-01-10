library(steponeR)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)
p2 <- list.files(path = "qPCR", pattern = ".txt", full.names = TRUE)
p2

master<-read.csv(header = TRUE, file ="acer_master.csv")

df <- steponeR(files = p2, 
              delim = "\t", 
              target.ratios = c("Sym.Cam"),
              fluor.norm = list(Cam = 0, Sym = 0), 
              ploidy = list(Cam = 2, Sym = 1),
              extract = list(Cam = .982, Sym = .813),
              copy.number = list(Cam = 1, Sym = 1) )
acer <- df$result

acer<-acer %>%
  filter(!Sample.Name== "positive",!Sample.Name == "negative") %>%
  filter(!Sample.Name == "neg", !Sample.Name == "pos", !Sample.Name == "NEC6")%>%
  filter(!Sample.Name == "277", !Sample.Name == "279")
nrow(acer)
#convets numeric timepoint in master sheet to character in order to combine them across the qpcr data
master<-master%>%mutate_if(is.numeric,as.character)%>%
    full_join(acer)
    #select(master,Sample.Name,TimePoint,Cam.CT.mean,Sym.CT.mean,Sym.Cam,Sym.CT.sd,Cam.CT.sd)
final<-select(master,Sample.Name,TimePoint,Cam.CT.mean,Sym.CT.mean,Sym.Cam,Sym.CT.sd,Cam.CT.sd )
final<-transform(final, TimePoint = as.numeric(TimePoint)) #converts Timpoint which was originally a "character" type to numeric type

master<-master %>%
  mutate(Date = as_date(Date,format = "%m/%e/%y"))



ggplot(master, aes(x = Date, y = log10(Sym.Cam), color = Treatment, group = interaction(Date,Treatment)))+
  geom_point()+
  geom_violin()

str(master)
