library(tidyverse)
library(cowplot)

#reading files

read_csv('data/sample_oilcontent.csv')

oilcontent <- read_csv('data/sample_oilcontent.csv')


#some plots

ggplot(
  data = oilcontent,
  mapping = aes(x = stage, y = oil_contentDW, colour = genotype, size = weight_mg)
) + geom_point() 
