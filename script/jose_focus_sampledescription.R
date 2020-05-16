library(tidyverse)
library(cowplot)

#reading files

read_csv('data/sample_oilcontent.csv')

oilcontent <- read_csv('data/sample_oilcontent.csv')


#some plots

Oilcontent_plot <- ggplot(data = oilcontent,
  mapping = aes(x = stage, y = oil_contentDW, colour = genotype, size = weight_mg)
) + geom_point() 

Oilcontent_plot <- Oilcontent_plot +  
  labs(title = "Oil content in transgenic and wildtype cowpeas",
       #caption = 'ADD CAPTION',
       x = "Seed developmental stage",
       y = "Oil content (% DW)",
       size = "Seed weight (mg)",
       colour = "Genotype") + 
  #scale_x_continuous())+
  scale_y_continuous(limits=c(0, 10))+
  
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    strip.background = element_blank(),
    panel.grid.major = element_line(size = 1),
    axis.title = element_text(size = 10, face = "bold"),
    )

ggsave(filename = "results/oilcontent_plot.png", plot = Oilcontent_plot,
       width = 15, height = 12, dpi = 300, units = "cm")