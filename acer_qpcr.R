library(steponeR)
library(tidyverse)
library(dplyr)
p2 <- list.files(path = "qPCR", pattern = ".txt", full.names = TRUE)
p2

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
  filter(!Sample.Name == "277", !Sample.Name == "279")

