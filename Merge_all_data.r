library(tidyverse)

#infile1 <- "../input/data/Meta_data.csv"
infile1 <- "../input/data/Meta_data_char.csv"
infile2 <- "../input/data/16S_data.csv"
infile3 <- "../input/data/sample_type_id.csv"
data_meta <- read_csv(infile1, col_names = TRUE)
data_16S <- read_csv(infile2, col_names = TRUE)
data_id <- read_csv(infile3, col_names = TRUE)

###########################Merge all data, 16S and Meta into a single tibble##########
data_all <- data_id %>%
  right_join(data_16S, by = c('Sample_id'='Sample_ID'))

data_all <- data_all %>%
  separate(Sample_type, c(NA, "Day"), sep = "cd") %>%
  mutate(Day = as.numeric(Day))

data_all <- data_meta %>% 
  right_join(data_all, by = c('ID'='Subject', 'Day'='Day'))
  
data_all <- data_all %>%
  separate(Event, c(NA, "Sample_type"), sep = "_") 

#write.csv(data_all, "../input/data/All_data.csv", row.names=FALSE, quote=FALSE)
write.csv(data_all, "../input/data/All_data_char.csv", row.names=FALSE, quote=FALSE)
