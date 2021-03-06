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
  
gene_annotation <- select(gene_atlas, 1,2, 25)

view(full_join(embryo_expression, gene_annotation, by = c('Id')))

embryo_expression_annotation <- full_join(embryo_expression, gene_annotation, by = c('Id'))
embryo_expression_annotation <- rename(embryo_expression_annotation, annotation=11, gene_name=10)

#tidying expression table

tidy_embryo_expression <- embryo_expression_annotation %>% 
  gather(key = Sample_name, value = ammount, -Id, -annotation, -gene_name)

tidy_embryo_expression <- tidy_embryo_expression %>% 
  mutate(sample = Sample_name)

tidy_embryo_expression <- separate(tidy_embryo_expression, Sample_name, into = c('sample_name', 'tissue', 'stage'), sep = "_")

tidy_embryo_expression <- separate(tidy_embryo_expression, sample_name, into = c('todelete', 'genotype'), sep = "m.")

tidy_embryo_expression <- rename(tidy_embryo_expression, gene_id=1, expression=8)


tidy_embryo_expression <- tidy_embryo_expression %>% select(-todelete)


#getting expression in log values
log10tidy_embryo_expression <- tidy_embryo_expression %>%
  filter(expression > 1)

log10tidy_embryo_expression <- log10tidy_embryo_expression %>% 
  mutate(expression = log10(expression))




#writing tidy combined data out
write_csv(tidy_embryo_expression, 'results/tidy_embryo_expression.csv')


#looking at the data

ggplot(
  data = tidy_embryo_expression,
  mapping = aes(x = gene_name, y = expression, colour = sample)
) + geom_point() 






#heat map

heatmap  <-  ggplot(data = log10tidy_embryo_expression, mapping = aes(x = sample,
                                                     y = gene_name,
                                                     fill = expression))+
  geom_tile()

  
#Heatmap with matrix


# Remove first column from seqdata and reduce gene numbers


#embryo_expression200 <- embryo_expression %>%
  #filter_all(all_vars(. > 200))

countmatrix <- embryo_expression[,-(1)]
countmatrix <- as.matrix(countmatrix)

#heatmap(countmatrix200)  #Error: cannot allocate vector of size 6.5 Gb


# Obtain CPMs
expressionCPM <- cpm(countmatrix)





# Store GeneID as rownames
#rownames(countmatrix) <- embryo_expression[,1]





#another way

embryo_expression_numeric <- embryo_expression %>% 
  mutate(norm.Wt_Embryo_1 = as.numeric(norm.WT_Embryo_1)) %>% 
  mutate(norm.Wt_Embryo_2 = as.numeric(norm.WT_Embryo_2)) %>% 
  mutate(norm.Wt_Embryo_3 = as.numeric(norm.WT_Embryo_3)) %>% 
  mutate(norm.Wt_Embryo_4 = as.numeric(norm.WT_Embryo_4)) %>% 
  mutate(norm.T_Embryo_1 = as.numeric(norm.T_Embryo_1)) %>% 
  mutate(norm.T_Embryo_2 = as.numeric(norm.T_Embryo_2)) %>% 
  mutate(norm.T_Embryo_3 = as.numeric(norm.T_Embryo_3)) %>% 
  mutate(norm.T_Embryo_4 = as.numeric(norm.T_Embryo_4)) 


install.packages("heatmap.plus")
library("heatmap.plus")


embryo_expression_matrix <- as.matrix(embryo_expression)

head(embryo_expression_matrix)

embryo_expression_matrix <- embryo_expression_matrix[,-1]

standard_expression <- scale(embryo_expression_matrix)



##more tets


row.names(embryo_expression) <- embryo_expression$Id
embryo_expression <- embryo_expression[1:41922, -1]
embryo_expression <- data.matrix(embryo_expression)

heatmap(embryo_expression)

