is.installed <- function(mypkg) is.element(mypkg, installed.packages()[,1]) 
is.installed('readxl') 
library(tidyverse)
library(readxl)

infile <- "../input/data/Meta_Reformatted.xlsx"
data_meta <- read_xlsx(infile, sheet = NULL, col_names = TRUE)

################### Rename colums to sample id ######	
data_meta <- data_meta %>%
  select(ID, Contraceptive, 'DAY 1':'DAY 42')

col.to <- str_split(c(names(data_meta[3:44])), " ", n=2) %>% map_chr(2)
col.from <- c(names(data_meta[3:44]))
data_meta <-data_meta %>% rename_at(vars(col.from), ~c(col.to))

################### Create Meta only tables with Sex and Bleed in Numeric/Charachter###### 
#With Charchters
data_meta <- data_meta %>%
	gather(Day, Info, -ID, -Contraceptive) %>%
	separate(Info, c("Bleeding", "Vag_swab", "Sexual_activity"), sep = ",") %>%
  mutate(ID = as.numeric(ID)) %>%
  mutate(Day = as.numeric(Day))

write.csv(data_meta, "../input/data/Meta_data_char.csv", row.names=FALSE, quote=FALSE)

#With numeric 
data_Meta <- data_meta %>%  
  mutate(dat =replace(Sexual_activity, Sexual_activity==0, NA),
         Sex =replace(dat, Sexual_activity=="S", 1)) %>%
  mutate(bat =replace(Bleeding, Bleeding==0, NA),
         Bleed =replace(bat, Bleeding=="B", 1)) %>%
  select(ID, Contraceptive, Day, Bleed, Sex) %>%
  mutate(Bleed = as.numeric(Bleed)) %>%
  mutate(Sex = as.numeric(Sex))

write.csv(data_Meta, "../input/data/Meta_data.csv", row.names=FALSE, quote=FALSE)
