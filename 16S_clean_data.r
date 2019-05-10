library(tidyverse)

infile <- "../input/data/3x5_runs_cleaned.tsv"
data_16S <- read_tsv(infile, col_names = TRUE)

################### Filter out controlls #######
data_16S <- data_16S %>%
  select(Seq_ID, Taxonomy, contains("MiMens"), -contains("Pos"), -contains("Neg"), -contains("pos"), -contains("neg"))

################### Rename colums to sample id #######
new.names <- str_split(c(names(data_16S[3:490])), "__006", n=2) %>% map_chr(2)
old.names <- c(names(data_16S[3:490]))
data_16S <-data_16S %>% rename_at(vars(old.names), ~c(new.names))

################### Remove OTUs with reads > 42 #######
data_16S <- data_16S %>%
  mutate(Reads = rowSums(select(., -Seq_ID, -Taxonomy))) %>%
  filter(Reads > 42)

Hits <- data_16S %>% 
  pull(Reads) %>%
  summary(Hits)
Hits

################### Split Taxonomy column #######
data_16S <- data_16S %>%
  gather(Sample_ID, Values, -Seq_ID, -Taxonomy, -Reads) %>%
  separate(Taxonomy, c(NA, NA, NA, NA, "Family", "Genus", "Species"), sep = ";") %>%
  separate(Sample_ID, c("Plate", "Sample_ID"), sep = "__") %>%
  separate(Plate, c(NA, "Plate"), sep = "_") %>%
  select(-Reads)

write.csv(data_16S, "../input/data/16S_data.csv", row.names=FALSE, quote=FALSE)

################### Spread taxonomy (family) This was unnescceary (not tidy) #################
data_16S_spread <- data_16S %>%
  unite(temp, Seq_ID, Family) %>%
  spread(temp, Values)
################### Collapse all Family counts across all samples #####################################  
Taxa <- data_16S %>% 
  pull(Family) %>%
  unique()

data_16S_family <- data_16S %>%
  group_by(Family = str_extract(Family, paste(Taxa, collapse="|"))) %>% 
  summarise(Reads = sum(Values)) %>%
  na.omit()
