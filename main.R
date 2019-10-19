library(tercen)
library(tidyverse)
library(flowCore)
library(SummarizedExperiment)
library(diffcyt)

ctx = tercenCtx()
data = ctx$as.matrix() 

# put row names on matrix
rownames(data) <- 1:nrow(data)
colnames(data) <- 1:ncol(data)

# transpose for flowCore
data = t(data)

marker_info <- as.data.frame(ctx$rselect())

experiment_info <- as.data.frame(ctx$cselect())
experiment_info <- experiment_info %>% select(c(1,2,3)) %>% distinct()

design <- createDesignMatrix(
  experiment_info, cols_design = c("group_id", "patient_id")
)

contrast <- createContrast(c(0, 1, rep(0, 7)))

# initialize 
seed = as.integer(ctx$op.value('seed'))
if (seed>0) set.seed(seed)

analysis_type =  as.character(ctx$op.value('analysis_type'))

d_input <- ctx$cselect() %>% 
  group_by(patient_id, group_id) %>% 
  group_rows %>% 
  sapply(function(x) data[x,])

if (analysis_type == "DA") {
  
  out_DA <- diffcyt(
    d_input = d_input, 
    experiment_info = experiment_info, 
    marker_info = marker_info, 
    design = design, 
    contrast = contrast, 
    analysis_type = "DA", 
    seed_clustering = seed
  )
  
  cluster_id    = rowData(out_DA$d_se)$cluster_id
  cluster_table = data.frame(cluster_id, .ci = seq_len(nrow(data))-1)
  stats_table   = rowData(out_DA$res) %>% as.tibble
  results_table = left_join(cluster_table, stats_table, by ="cluster_id")
  
  results_table = ctx$addNamespace(results_table)
  ctx$save(results_table)
}

if (analysis_type == "DS") {
  out_DS <- diffcyt(
    d_input = d_input, 
    experiment_info = experiment_info, 
    marker_info = marker_info, 
    design = design, 
    contrast = contrast, 
    analysis_type = "DS", 
    seed_clustering = seed, 
    plot = FALSE
  )
  
  cluster_id    = rowData(out_DS$d_se)$cluster_id
  cluster_table = data.frame(cluster_id, .ci = seq_len(nrow(data))-1)
  stats_table   = rowData(out_DS$res) %>% as.tibble
  results_table = left_join(cluster_table, stats_table, by ="cluster_id")
  marker_table  = marker_info %>% select(marker_name)  %>% mutate (.ri = seq_len(nrow(marker_info))-1)
  results_table = left_join(results_table, marker_table, by = c("marker_id" = "marker_name"))
  results_table = results_table %>% select(-marker_id, -ID)
  
  results_table = ctx$addNamespace(results_table)
  ctx$save(results_table)
}

