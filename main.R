library(tercen)
library(tercenApi)
library(dplyr, warn.conflicts = FALSE)
library(diffcyt)
library(SummarizedExperiment)

ctx = tercenCtx()

method <- ctx$op.value("method", as.character, "DA_edgeR")

counts <- ctx$as.matrix()
colData <- ctx$cselect() %>% 
  tidyr::unite("sample_id", sep = " - ") %>%
  mutate(.ci = seq_len(nrow(.)) - 1)
rowData <- ctx$rselect() %>%
  mutate(.ri = seq_len(nrow(.)) - 1) %>%
  mutate(cluster_id = .ri + 1) %>%
  as.data.frame()

first_color <- ctx$colors[[1]]
if(length(ctx$labels) > 0) {
  first_label <- ctx$labels[[1]]
  cols_design <- c("group_id", "patient_id")
} else {
  first_label <- NULL
  cols_design <- c("group_id")
}


### experimental design
experiment_info <- ctx$select(c(".ci", all_of(first_color), all_of(first_label))) %>%
  unique() %>%
  dplyr::rename(group_id = first_color, patient_id = first_label) %>%
  full_join(colData, by = ".ci") %>%
  as.data.frame()
rownames(experiment_info) <- experiment_info$sample_id
colnames(counts) <- rownames(experiment_info)
rownames(counts) <- rownames(rowData)

rowData$n_cells <- rowSums(counts)

se <- SummarizedExperiment(
  assays = list(counts = counts),
  rowData = rowData,
  colData = experiment_info[,"group_id"]
)
colnames(se) <- rownames(experiment_info)
rownames(se) <- rownames(rowData)

design <- createDesignMatrix(experiment_info, cols_design = cols_design)
contrast <- createContrast(as.numeric(grepl("group_id", colnames(design))))

if(method == "DA_edgeR") {
  res <- testDA_edgeR(se, design, contrast)
} else if(method == "DA_GLMM") {
  if(length(ctx$labels) > 0) {
    formula <- createFormula(experiment_info, cols_fixed = "group_id", cols_random = "patient_id")
  } else {
    formula <- createFormula(experiment_info, cols_fixed = "group_id")
  }
  res <- testDA_GLMM(se, formula, contrast)
}

rowData(res) %>% 
  as_tibble() %>%
  mutate(.ri = as.integer(cluster_id - 1)) %>%
  select(-cluster_id) %>%
  ctx$addNamespace() %>%
  ctx$save()
