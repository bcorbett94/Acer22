#Host/Symbiont qPCR data Code comparing which primers between CAM & Pax 
library(steponeR)
library(tidyverse)
library(dplyr)
plates <- list.files(path = "qPCR", pattern = ".txt", full.names = TRUE)
plates

# Read in data and calculate target ratios
df <- steponeR(files = plates, 
               delim = "\t", 
               target.ratios = c("Sym.Cam", "Sym.Pax","Pax.Cam"), 
               #To determine which (cam or Pax) amplified the best and which one will be used moving forward
               fluor.norm = list(Cam = 0, Pax = 0, Sym = 0), 
               ploidy = list(Cam = 2, Pax = 2, Sym = 1),
               extract = list(Cam = .982, Pax = .982, Sym = .813),
               copy.number = list(Cam = 1, Pax = 1, Sym = 1) )
qpcr <- df$result

qpcr<-qpcr %>%
  filter(!Sample.Name== "positive",!Sample.Name == "negative") %>%
  filter(!Sample.Name == "277", !Sample.Name == "279")#Only two that didn't work

summary(qpcr)

CAM<-qpcr$Sym.Cam

hist(CAM)

ggplot(qpcr, aes(x= Sym.Pax)) + 
  geom_histogram()

ggplot(qpcr, aes(x= log10(Pax.Cam))) + 
  geom_histogram()
# tail to the left, Pax is amplifying earlier than Cam

ggplot(qpcr, aes(x= log10(Sym.Pax))) + 
  geom_histogram()

ggplot(qpcr, aes(x= log10(Sym.Cam))) + 
  geom_histogram()

exp(mean(log10(qpcr$Sym.Pax)))
exp(mean(log10(qpcr$Sym.Cam)))
#average o the log of the sym to cam, avg of log of sym to pax 

#CAM marker, palacio castro 2021 because it has been cited previously 
#histogram sym to pax ratio and sym to cam ratio, or box plot, 
#calculate average, symbiont to host ratio, cam to pax of box plot 
#averages, do a plot log function sym to cam 

#As a whole looking into


