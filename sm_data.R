library(steponeR)
library(tidyverse)
library(dplyr)
pq <- list.files(path = ".", pattern = ".txt", full.names = TRUE)
p
#pl
# Read in data and calculate target ratios
df <- steponeR(files = pq, 
               delim = "\t", 
               target.ratios = c("C.D"), 
               fluor.norm = list(C = 2.234, D = 1),
               copy.number = list(C = 20, D = 3),
               ploidy = list(C = 1, D = 1),
               extract = list(C = 0.813, D = 0.813))
qpcr <- df$result

ampboth<-filter(qpcr,C.reps ==2 & D.reps ==2)

ampone<-filter(qpcr,C.reps ==2 | D.reps ==2)
