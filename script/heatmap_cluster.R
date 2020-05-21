library(tidyverse)

expn_data <- read_csv("results/tidy_embryo_expression.csv")

# Data processing
# Extremely basic, not ideal for a proper analysis
# 1: Remove genes with 0 expression
# 2: Convert to log scale
# 3: Normalise each gene to mean=0, sd=1

expn_scaled <- expn_data %>% 
	group_by(gene_name) %>% 
	filter(sum(expression) > 100) %>% 
  #filter(var(expression) > 50)
	mutate(
		log_expn = log2(expression + 1),
		norm_expn = (log_expn - mean(log_expn))/sd(log_expn)
	) %>% 
	ungroup()

# Check a plot
# Looks pretty good, but alphabetical order for genes not helpful

expn_scaled %>% 
	ggplot(aes(y = gene_name, x = sample, fill = norm_expn)) +
	geom_tile() +
	# Too many genes to show the names or the axis ticks so remove them
	scale_y_discrete(labels = NULL, breaks = NULL) +
	# Convert the Brewer RdBu palette into a continuous colour scale
	scale_fill_distiller(palette = "RdBu")

# Need expression in wide format for clustering
# The final step converts it from a tibble into a data.frame with
# the gene names as the 'rownames' of the data.frame
# This removes the gene names from the data (so we just have numeric data) 
# and makes them more metadata for each row (which will be important later)

expn_table <- expn_scaled %>%
	select(gene_name, norm_expn, sample) %>% 
	spread(key = sample, value = norm_expn) %>% 
	column_to_rownames("gene_name") 
	
# Can't use hierarchical clustering on this many genes because you need to
# make all-to-all pairwise comparisons (this is probably why the heatmap
# function didn't work). For this sort of preliminary data, I like to use k-means
# clustering with an arbitrary choice of the number of clusters (k) to get a quick
# look at what's in the data

expn_clustered <- cluster::clara(expn_table, k = 10)



# The cluster result (expn_clustered) has an element called clustering that 
# shows what cluster each gene is assigned to. We can access this with expn_clustered$clustering
# We sort this so that the genes from each cluster are together, then get 
# the gene names with names() -- the rownames from the expn_table data.frame
# are used to give names to the clustering output.

gene_order <- names(sort(expn_clustered$clustering))


#getting cluster numbers

cluster_info <- tibble(cluster = expn_clustered$clustering, gene_name = names(expn_clustered$clustering))

view(full_join(expn_scaled, cluster_info, by = c('gene_name')))

tidy_expression_cluster <- full_join(expn_scaled, cluster_info, by = c('gene_name'))
tidy_expression_cluster <- tidy_expression_cluster %>% select(-gene_id)

tidy_expression_cluster %>% 
  group_by(cluster) %>% 
  summarise(num_row = n())


# Now redo the plot with one change, the gene names are first converted to a
# factor, with the levels of the factor being our gene order.
# This makes the genes be plotted in the correct order by ggplot.

heatmap_plot <- tidy_expression_cluster %>% 
	mutate(gene_name = factor(gene_name, levels = gene_order)) %>% 
	ggplot(aes(y = gene_name, x = sample, fill = norm_expn)) +
	geom_tile() +
  facet_grid(~ genotype, switch = "x", scales = "free_x", space = "free_x") +
	scale_y_discrete(labels = NULL, breaks = NULL) +
	scale_fill_distiller(palette = "RdBu") +
  scale_x_discrete(labels=c("norm.T_Embryo_1" = "1", 
                            "norm.T_Embryo_2" = "2", 
                            "norm.T_Embryo_3" = "3", 
                            "norm.T_Embryo_4" = "4", 
                            "norm.WT_Embryo_1" = "1",
                            "norm.WT_Embryo_2" = "2",
                            "norm.WT_Embryo_3" = "3",
                            "norm.WT_Embryo_4" = "4"))


heatmap_plot <- heatmap_plot +  
  labs(title = "Gene expression heatmap during seed developmet",
       subtitle = "Comparison between the wildtype and the high-oil transgenic",
       caption = '(Raw RNAseq data analysed by J. Verdier)',
       x = "Sample developmental stage",
       y = "Genes",
       fill = "Normalized expression") +
  
  theme(
    plot.title = element_text(face = "bold"),
    #strip.background = element_blank(),
    axis.title = element_text(size = 10, face = "bold"),
    #axis.text.x = element_text(angle=45),
    legend.position="bottom")
  

ggsave(filename = "results/heatmap_plot.png", plot = heatmap_plot,
       width = 15, height = 20, dpi = 300, units = "cm")




####getting heat map og high expressed genes

expn_scaled_top <- expn_data %>% 
  group_by(gene_name) %>% 
  filter(sum(expression) > 100000) %>% 
  #filter(var(expression) > 50)
  mutate(
    log_expn = log2(expression + 1),
    norm_expn = (log_expn - mean(log_expn))/sd(log_expn)
  ) %>% 
  ungroup()

# Need expression in wide format for clustering
# The final step converts it from a tibble into a data.frame with
# the gene names as the 'rownames' of the data.frame
# This removes the gene names from the data (so we just have numeric data) 
# and makes them more metadata for each row (which will be important later)

expn_table_top <- expn_scaled_top %>%
  select(gene_name, norm_expn, sample) %>% 
  spread(key = sample, value = norm_expn) %>% 
  column_to_rownames("gene_name") 

# Can't use hierarchical clustering on this many genes because you need to
# make all-to-all pairwise comparisons (this is probably why the heatmap
# function didn't work). For this sort of preliminary data, I like to use k-means
# clustering with an arbitrary choice of the number of clusters (k) to get a quick
# look at what's in the data

expn_clustered_top <- cluster::clara(expn_table_top, k = 10)



# The cluster result (expn_clustered) has an element called clustering that 
# shows what cluster each gene is assigned to. We can access this with expn_clustered$clustering
# We sort this so that the genes from each cluster are together, then get 
# the gene names with names() -- the rownames from the expn_table data.frame
# are used to give names to the clustering output.

gene_order_top <- names(sort(expn_clustered_top$clustering))


#getting cluster numbers

cluster_info_top <- tibble(cluster = expn_clustered_top$clustering, gene_name = names(expn_clustered_top$clustering))

view(full_join(expn_scaled_top, cluster_info_top, by = c('gene_name')))

tidy_expression_cluster_top <- full_join(expn_scaled, cluster_info, by = c('gene_name'))
tidy_expression_cluster_top <- tidy_expression_cluster %>% select(-gene_id)

tidy_expression_cluster_top %>% 
  group_by(cluster) %>% 
  summarise(num_row = n())


# Now redo the plot with one change, the gene names are first converted to a
# factor, with the levels of the factor being our gene order.
# This makes the genes be plotted in the correct order by ggplot.

heatmap_plot_top <- tidy_expression_cluster_top %>% 
  mutate(gene_name = factor(gene_name, levels = gene_order)) %>% 
  ggplot(aes(y = gene_name, x = sample, fill = norm_expn)) +
  geom_tile() +
  facet_grid(~ genotype, switch = "x", scales = "free_x", space = "free_x") +
  scale_y_discrete(labels = NULL, breaks = NULL) +
  scale_fill_distiller(palette = "RdBu") +
  scale_x_discrete(labels=c("norm.T_Embryo_1" = "1", 
                            "norm.T_Embryo_2" = "2", 
                            "norm.T_Embryo_3" = "3", 
                            "norm.T_Embryo_4" = "4", 
                            "norm.WT_Embryo_1" = "1",
                            "norm.WT_Embryo_2" = "2",
                            "norm.WT_Embryo_3" = "3",
                            "norm.WT_Embryo_4" = "4"))


heatmap_plot_top <- heatmap_plot_top +  
  labs(title = "Gene expression heatmap during seed developmet",
       subtitle = "Comparison between the wildtype and the high-oil transgenic",
       x = "Sample developmental stage",
       y = "Genes",
       fill = "Normalized expression") +
  
  theme(
    plot.title = element_text(face = "bold"),
    #strip.background = element_blank(),
    axis.title = element_text(size = 10, face = "bold"),
    #axis.text.x = element_text(angle=45),
    legend.position="bottom")


ggsave(filename = "results/heatmap_plot_top.png", plot = heatmap_plot_top,
       width = 15, height = 20, dpi = 300, units = "cm")

