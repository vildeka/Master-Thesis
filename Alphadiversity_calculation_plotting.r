library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(vegan)
library(cowplot)

infile <- "../input/data/All_data_char.csv"
data_ <- read_csv(infile, col_names = TRUE, guess_max = min(10000, 100000))

################################## Calculate diversity and Sum taxa per day #####
data_Diversity <- data_ %>%
  select(ID, Day, Genus, Family, Values) %>%
  group_by(ID, Day) %>%
  nest() %>%
  mutate(Tax_G = map(data, ~ aggregate(.$Values, by=list(.$Genus), sum))) %>%
  unnest(Tax_G) %>% spread(Group.1, x) %>%
  mutate(Diversity = diversity(.[-2], index = "invsimpson")) %>%
  select(ID, Day, Diversity)
  
data_Taxa <- data_ %>%
  select(ID, Day, Genus, Bleeding, Sexual_activity, Family, Values) %>%
  nest(-ID, -Day, -Bleeding, -Sexual_activity) %>%
  mutate(Tax_F = map(data, ~ aggregate(.$Values, by=list(.$Family), sum))) %>%
  unnest(Tax_F) %>%
  group_by(ID, Day) %>%
  mutate(Percent = x/sum(x)) %>%
  filter(Percent > 0.01) %>%
  arrange(ID)
  
data_Meta <- data_Taxa %>%  
  mutate(dat =replace(Sexual_activity, Sexual_activity=="0", as.numeric(NA)),
         "Sexual activity" =replace(dat, dat=="S", as.numeric(-0.1))) %>%
  mutate(bat =replace(Bleeding, Bleeding=="0", as.numeric(NA)),
         Bleeding =replace(bat, Bleeding=="B", as.numeric(-0.05))) %>%
  gather(Meta,value, Bleeding, "Sexual activity") %>%
  rename(Family = Group.1) %>%
  mutate(value = as.numeric(value)) %>%
  select(ID, Day, Meta, value, Family, Percent)

################################## Colour pallet ############
dd <- sort(unique(data_Meta$Family))
#dd.col <- c('#5B8A9A', '#77E8DC', '#254C53', '#8AD0F9', '#D27AE9', '#867AC4', '#E8CAEE', '#72329E', '#BEE091', '#48853A', '#8BEC6E', "#E77A9F", '#673B11', '#F3A852', '#E1502F', '#932846', '#A17A65', "#8B9662", '#EC3686', '#414B17', '#c85d3f', '#F0D447', '#4141F5', '#E743F4', '#F39680', '#F2A3BE', '#F2EC81')
dd.col <- c("#afb7ee","#ffc8d9","#c4a597","#7fe2e9","#c4ce96","#F8D0A4","#f6fff9","#e1caff","#9aacce","#A8EDFC","#ffd9b8","#f1a6b1","#c8ffd5","#E3E6AD","#d1b9ee","#88c29c","#6ececc","#83dafb","#7cb6b6","#fbdea9","#91c6f7","#f5f8bd","#8db1c5","#fab0aa","#96f3eb","#a1b37d","#c0a2c1")
names(dd.col)  <- dd

################################## Splitting of facets ######
custom <- theme_light() + theme(legend.position="left", 
                strip.text.x = element_text(colour='#999999'),
                strip.background = element_rect(colour ='#e5e5e5', fill ='#e5e5e5')) 

plot_t <- data_Meta %>%
  ungroup() %>%
  split(ceiling(group_indices(.,.$ID)/3)) %>%
  map(., ~ggplot(., aes(x=Day, y=Percent/2, fill=Family)) +
               geom_col(aes(fill=Family), width=1) + 
               guides(col = guide_legend(nrow = 2), fill = guide_legend(override.aes = list(shape = NA)), shape = guide_legend(override.aes = list(size = 1))) +
               facet_wrap(~ID) + scale_fill_manual(values = dd.col) +
               scale_y_continuous(labels = scales::percent) + xlim(c(0,42)) +
               custom + theme(axis.title.y = element_blank(), 
                              legend.position="bottom",
                              #legend.title = element_blank(), 
                              #legend.box.margin=unit(c(-3,0,0,0),"cm"),
                              strip.background = element_blank(),
                              strip.text.x=element_blank(),
                              plot.margin=unit(c(0,1,0.5,0), "cm"))) %>%
  map(., ~.x + geom_point(aes(x=Day, y=value, colour=Meta), size = 1))


plot_d <- data_Diversity %>%
  split(ceiling(group_indices(.,.$ID)/3)) %>%
  map(~ggplot(., aes(x=Day, y=Diversity, colour=Diversity)) + xlim(c(0,42))  +
        geom_line() + facet_wrap(~ID) +  
        labs(title = "Diversity and Taxonomy") + theme_light() + custom +
        theme(axis.text.x = element_blank(), 
              axis.title.y = element_blank(),
              axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(),
              legend.box.margin=unit(c(0,0,-10,0),"cm"),
              plot.margin=unit(c(0.5,1,0,0), "cm"),
              ))

################################## Put facet plots together ########
plot_Div_Tax <- map2(plot_d, plot_t, ~plot_grid(.x, .y, ncol = 1, align = 'v', axis = 'l', rel_heights = c(1, 3)))

################################## Save all as pdf ########
pdf(file = "../input/plots/Tax_Div.pdf",
    width = 8, # The width of the plot in inches
    height = 6 # The height of the plot in inches
    )
plot_Div_Tax
dev.off()

################################## Save as individual png ########
name <- map(plot_t, ~paste((unique(.x[["data"]][["ID"]])), collapse = ","))
plotnames = map2(plot_Div_Tax, name, ~paste0("Taxonomy_", .y, ".png")) %>%
  flatten()
plotnames

walk2(plotnames, plot_Div_Tax, ~ggsave(filename = .x, plot = .y, 
                                       path = "../input/plots/", 
                                       width = 6, height = 6,
                                             ))

################################## For all IDs in one facet plot #####

plot_T = ggplot(data_Meta, aes(x=Day, y=Percent/2, fill=Family)) + 
  geom_col(aes(fill=Family), width=1) + ggtitle("Taxonomy") + 
  custom + theme(legend.position="bottom", axis.title.y = element_blank()) + 
  guides(color=guide_legend(nrow=2), fill=guide_legend(nrow=4, override.aes = list(shape = NA))) +
  facet_wrap(~ID, scales = "free_y") + scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = dd.col) + 
  geom_point(aes(x=Day, y=value, colour=Meta), size = 1)

plot_D = ggplot(data_Diversity, aes(x=Day, y=Diversity, colour=Diversity)) + 
  geom_line() + custom + theme(legend.position="left") + 
  facet_wrap( ~ID)

################################# How to view the plots ############
data_plot$plot_I
print(data_plot$plot_M[3])
print(data_plot$plot_I) 
