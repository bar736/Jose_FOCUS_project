library(tidyverse)
library(cowplot)

#reading files

read_csv('data/embryo_all_normalised.csv')

embryo_expression <- read_csv('data/embryo_all_normalised.csv')

read_csv('data/cowpea_geneatlas.csv')

gene_atlas <- read_csv('data/cowpea_geneatlas.csv')


#merging files

gene_atlas <- gene_atlas <- rename(gene_atlas, Id=1)

colnames(gene_atlas)
  
gene_annotation <- select(gene_atlas, 1, 25)

view(full_join(embryo_expression, gene_annotation, by = c('Id')))

embryo_expression_annotation <- full_join(embryo_expression, gene_annotation, by = c('Id'))
embryo_expression_annotation <- rename(embryo_expression_annotation, annotation=10)

#tidying expression table

tidy_embryo_expression <- embryo_expression_annotation %>% 
  gather(key = Sample_name, value = ammount, -Id, -annotation)

tidy_embryo_expression <- separate(tidy_embryo_expression, Sample_name, into = c('sample_name', 'tissue', 'stage'), sep = "_")

tidy_embryo_expression <- separate(tidy_embryo_expression, sample_name, into = c('sample', 'genotype'), sep = "m.")

tidy_embryo_expression <- rename(tidy_embryo_expression, gene=1, expression=7)

tidy_embryo_expression <- tidy_embryo_expression %>% select(-sample)


#writing tidy combined data out
write_csv(tidy_embryo_expression, 'results/tidy_embryo_expression.csv')


