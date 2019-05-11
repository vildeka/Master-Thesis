library(tidyverse)
library(ggplot2)

infile <- "../input/data/table_distances.csv"
data_Luisa <- read_csv(infile, col_names = TRUE)

######### Including only Subjects using Hormonon/non-Hormonal Contraceptive  #############
target <- c("11","15","29","51","71","80","128","139","150","155")

data_a <- data_Luisa %>%
  mutate(X1 = as.character(X1)) %>%
  filter(X1 %in% target) %>%
  rename("Subject ID"=X1) %>%
  gather(., "1":"42", key=Day, value=Values) %>%
  mutate(Day = as.numeric(as.character(.$Day))) %>% 
  ####### Median (not used in final plot) #######
  group_by(Day) %>%
  nest() %>%
  mutate(Median=map(data, ~log(median(.$Values, na.rm = TRUE)))) %>%
  unnest(Day, Median) %>%
  select(Day, Median) %>%

################### Defining time of Menses #############
rect <- data.frame(xmin=c(0,28), xmax=c(5,33), ymin=-Inf, ymax=Inf)
################### Plotting #######
data_a %>%
  ggplot(., aes(x=Day, y=log(Values))) + geom_point() + 
  geom_smooth(aes(colour = I("black"))) + geom_vline(xintercept=c(28), linetype="dotted") +
  theme(legend.position = "none") + scale_x_continuous(breaks =c(0,7,14,21,28,35,42)) + 
  labs(y= expression(paste("log(",Delta*D[J][S]," ", Delta*italic("t)"), sep="")),
       x= "Menatrual cycle (Days)") +
  geom_rect(data=rect, inherit.aes=FALSE, aes(xmin=xmin, xmax=xmax, ymin=ymin,
            ymax=ymax), color="transparent", fill="tomato", alpha=0.3) 
