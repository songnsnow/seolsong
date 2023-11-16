# interpret results ------------------------------------------------------------------------------------------------
# load libraries
library(ggplot2)
library(numbat)
library(dplyr)
library(glue)
library(data.table)
library(ggtree)
library(stringr)
library(tidygraph)
library(patchwork)


# summarize output files to Numbat object
nb_test = Numbat$new(out_dir = '/mnt/mydata/output/rcc_1')

# plot
mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")
  
nb_test$plot_phylo_heatmap(
  clone_bar = TRUE, 
  p_min = 0.9,
  pal_clone = mypal
)
ggsave('/mnt/mydata/plot/heatmap_rcc.1.png',width=10, height=3,scale=1.5)

# pseudobulk
nb_test$bulk_clones %>% 
  filter(n_cells > 50) %>%
  plot_bulks(
    min_LLR = 10, # filtering CNVs by evidence
    legend = TRUE
  )
ggsave('/mnt/mydata/plot/bulk_rcc.1.png',width=10, height=3,scale=1.8)

# consensus
nb_test$plot_consensus()
ggsave('/mnt/mydata/plot/consensus_rcc.1.png',width=10, height=3,scale=1.5)




