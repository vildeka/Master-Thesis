library(tidyverse)
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')
library(phyloseq)
install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")
install_github("gauravsk/ranacapa")
library(ranacapa)
library(RColorBrewer)
library(ggbiplot)
library(vegan)

infile <- "../input/data/3x5_runs_cleaned.tsv"
infile2 <- "../input/data/All_data_char.csv"
data_16S <- read_tsv(infile, col_names = TRUE)
data_all <- read_csv(infile2, col_names = TRUE)

###################### Normalizing the data ####### 
data_OTU_norm <- data_all %>%
  select(Sample_id, Seq_ID, Values) %>%
  group_by(Sample_id) %>%
  mutate(Percent = Values/sum(Values)) %>%
  ungroup() %>%
  select(-Values) %>% spread(Sample_id, Percent) %>%
  rename(rowname = Seq_ID) %>%
  column_to_rownames()
###################### Creating the phyloseq object ###########
data_reads <- data_all %>% #mark out the 15 samples that have less than 2000 reads
  select(Sample_id, Values) %>%
  group_by(Sample_id) %>% nest() %>%
  mutate(Counts = map(data, ~sum(.$Values))) %>%
  unnest(Counts) %>% select(Sample_id, Counts) %>% 
  filter(Counts < 2000) %>%
  filter(Counts > 0)

data_meta_ <- data_all %>%
  select(Sample_id, ID, Contraceptive, Day, Bleeding, Sexual_activity) %>%
  distinct(Sample_id, ID, Contraceptive, Day, Bleeding, Sexual_activity) %>%
  mutate(dat =replace(Sexual_activity, Sexual_activity=="0", "No"),
         Sex =replace(dat, dat=="S", "Yes")) %>%
  mutate(bat =replace(Bleeding, Bleeding=="0", "No"),
         Bleeding =replace(bat, Bleeding=="B", "Yes")) %>%
  select(Sample_id, ID, Contraceptive, Day, Bleeding, Sex)

data_meta <- data_meta_ %>%
  arrange(ID) %>%
  mutate(IDs = as.character(.$ID)) %>%
  full_join(., data_reads) %>%
  rename("Subject ID" = ID) %>%
  rename("Sexual activity" = Sex) %>%
  mutate(Counts = as.character(.$Counts)) %>%
  mutate(Counts =replace(Counts, Counts!="NA", "< 2000"),
         Counts =replace(Counts, Counts=="NA", "> 2000")) %>%
  rename(rowname = Sample_id) %>%
  column_to_rownames()
  
data_family <- data_all %>%
  select(Sample_id, Family, Values) %>%
  group_by(Sample_id, Family)  %>%
  summarise(Reads = sum(Values)) %>%
  ungroup() %>%
  spread(Sample_id, Reads) %>%
  rename(rowname = Family) %>%
  column_to_rownames()

data_OTU <- data_all %>%
  select(Sample_id, Seq_ID, Values) %>%
  spread(Sample_id, Values) %>%
  rename(rowname = Seq_ID) %>%
  column_to_rownames()

data_Taxa <- data_16S %>%
  select(Seq_ID, Taxonomy) %>%
  separate(Taxonomy, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  rename(rowname = Seq_ID) %>%
  column_to_rownames()
  
FAM = otu_table(data_family, taxa_are_rows = TRUE)
OTU = otu_table(data_OTU, taxa_are_rows = TRUE) 
OTU_norm = otu_table(data_OTU_norm, taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(data_Taxa))
MET = sample_data(data_meta)
physeq = phyloseq(OTU_norm, TAX, MET)
physeq_fam = phyloseq(OTU_norm, TAX, MET)

head(sample_data(data_meta))
###################### Rarefraction curve ####### 
rare <- ggrare(physeq, step = 500, label = NULL, color = "Counts",
       plot = TRUE, parallel = FALSE, se = FALSE)

rare <- rare + scale_x_continuous(name ="Sequence Sample Size", breaks =c(0,2000,10000,20000,30000))
rare <- rare + ylim(c(0,100))
###################### Creating the distance matrix ################
DistJSD <- distance(physeq, method="jsd") 
DistJSD_fam <- distance(physeq_fam, method="jsd") 
DistBC = distance(physeq, method = "bray")

#Root of Jensen-Shannon:
DistrJSD <- sqrt(DistJSD)
DistrJSD_fam <- sqrt(DistJSD_fam)

ordrJSD = ordinate(physeq, method = "PCoA", distance = DistrJSD)
ordBC = ordinate(physeq, method = "PCoA", distance = DistBC)

#looking at distribution
plot_O <- plot_ordination(physeq, ordrJSD, type = "scree", title="Scree Plot: Jensen-Shannon MDS", axes=c(3,1), color="Set")
plot_O <- plot_ordination(physeq, ordBC, type = "scree", title="Scree Plot: Bray-Curtis MDS", axes=c(3,1), color="Set")

###################### Colour pallet and custom theme ########

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sb <- scale_color_manual("Subject ID", values= c("#c9512e","#3fbcc1","#d94467","#5bb84a","#c34fb0","#7b64cc","#dd943d","#6587ca","#5c782d","#d189c4","#51a876","#a1486e","#adb46d","#d07869"))
sc <- scale_color_manual(MET$Contraceptive, name ="Contraceptive", values= c("#c9512e","#3fbcc1","#adb46d"))
sf <- scale_colour_gradientn(colours = myPalette(100), limits=c(1, 42))
scale_color_m

cleanup = theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(color="black"))

custom <- theme_light() + theme(legend.position="left", 
                                strip.text.x = element_text(colour='#999999'),
                                strip.background = element_rect(colour ='#e5e5e5', fill ='#e5e5e5')) 

###################### PCoA plotting ########
#Bray Curtis not used:
plot_ordination(physeq, ordBC,shape="Sexual activity") + 
  geom_point(mapping = aes(color= IDs, size=Day, shape=Sex)) + 
  ggtitle("PCoA: Bray-Curtis") + sb + custom

plot_ordination(physeq, ordBC, color="Day", shape="Sexual activity") + 
  geom_point(mapping = aes(shape=Sex)) + 
  ggtitle("PCoA: Bray-Curtis") + sc + facet_wrap( ~ID)

#sqrt of Jensen Shannon used:
PcOA_rJSD <- plot_ordination(physeq, ordrJSD, shape="Sexual.activity") + theme_light() + 
             geom_point(mapping = aes(color= IDs, size=Day, shape=Sexual.activity)) + 
             ggtitle("PCoA: sqrt of Jensen-Shannon") + sb + custom + 
             guides(col = guide_legend(ncol = 2))

PCOA_facet <- plot_ordination(physeq, ordrJSD, color="Day", shape="Sexual.activity") + 
              geom_point(mapping = aes(shape=Sexual.activity)) + custom +
              ggtitle("PCoA: sqrt of Jensen-Shannon") + sf + facet_wrap( ~Subject.ID)

PcOA_C <- plot_ordination(physeq, ordrJSD, shape="Sexual.activity") + theme_light() + 
             geom_point(mapping = aes(color= Contraceptive, size=Day, shape=Sexual.activity)) + 
             ggtitle("PCoA: sqrt of Jensen-Shannon") + sc + custom 

###################### Save all as pdf ########
pdf(file = "../input/plots/Beta_diversity.pdf",
    width = 11, # The width of the plot in inches
    height = 7 # The height of the plot in inches
)
PcOA_rJSD
PCOA_facet
dev.off()

###################### Luisa's mess ####
for(i in 2:(ncol(DistBC_matr)-1)){
  for(j in (i+1):ncol(DistBC_matr)){
    res[i]<-rownames(DistBC_matr)[i]
  }
}
#create a matrix of the distance matrix that can be saved 
DistBC_matr<-as.matrix(DistBC)

DistJSD_matr<-as.matrix(DistrJSD)
write.csv(DistJSD_matr, file = "../input/data/DistJSD.csv")
