suppressPackageStartupMessages({
  library(tercen)
  library(tercenApi)
  library(dplyr, warn.conflicts = FALSE)
  library(diffcyt)
  library(SummarizedExperiment)
})

ctx = tercenCtx()

method <- ctx$op.value("method", as.character, "DA_edgeR")

first_color <- ctx$colors[[1]][[1]]

if(length(ctx$colors) > 1) {
  second_color <- ctx$colors[[2]][[1]]
} else {
  second_color <- NULL
}

if(length(ctx$labels) > 0) {
  first_label <- ctx$labels[[1]]
  cols_design <- c("group_id", "patient_id")
} else {
  first_label <- NULL
  cols_design <- c("group_id")
}

df <- ctx$select(unique(c(
        ".ri", ".ci", ".y",
        ".axisIndex", ".colorLevels",
        all_of(first_color), all_of(second_color), all_of(first_label)
      )))

counts <- df %>%
  filter(.axisIndex == 0) %>%
  tidyr::pivot_wider(
    id_cols = ".ri",
    names_from = ".ci",
    values_from = ".y",
    values_fill = 0
  ) %>%
  select(-.ri) %>%
  as.data.frame()


colData <- ctx$cselect() %>% 
  tidyr::unite("sample_id", sep = " - ") %>%
  mutate(.ci = seq_len(nrow(.)) - 1)
rowData <- ctx$rselect() %>%
  mutate(.ri = seq_len(nrow(.)) - 1) %>%
  mutate(cluster_id = as.factor(.ri + 1)) %>%
  as.data.frame()

### experimental design
experiment_info <- df %>% select(c(".ci", all_of(first_color), all_of(first_label))) %>%
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
  colData = experiment_info[, "group_id"]
)
colnames(se) <- rownames(experiment_info)
rownames(se) <- rownames(rowData)


## Prepare design, contrast and model formula
design <- createDesignMatrix(experiment_info, cols_design = cols_design)
contrast <- createContrast(as.numeric(grepl("group_id", colnames(design))))

if(length(ctx$labels) > 0) {
  formula <- createFormula(
    experiment_info,
    cols_fixed = "group_id",
    cols_random = "patient_id"
  )
} else {
  formula <- createFormula(
    experiment_info,
    cols_fixed = "group_id"
  )
}

if(method == "DA_edgeR") {
  res <- testDA_edgeR(se, design, contrast)
} else if(method == "DA_GLMM") {
  res <- testDA_GLMM(se, formula, contrast)
} else {
  if(length(unique(df$.axisIndex)) < 2) stop("Two layers are needed to use a DS method.")
  
  list_medians <-  df %>% 
    filter(.axisIndex == 1) %>%
    mutate(marker = !!sym(second_color)) %>%
    group_by(.ci, .ri, marker) %>%
    summarise(median = median(.y, na.rm = TRUE)) %>%
    tidyr::pivot_wider(
      id_cols = c("marker", ".ri"),
      names_from = ".ci",
      values_from = "median",
      values_fill = NA
    ) %>%
    ungroup %>%
    split(.$marker)
  
  medians <- lapply(list_medians, function(x) {
    x %>% select(-marker, -.ri) %>% as.data.frame()
  })
  
  se_meds <- SummarizedExperiment(
    assays = medians,
    rowData = rowData,
    colData = experiment_info[, "group_id"]
  )
  colnames(se_meds) <- rownames(experiment_info)
  rownames(se_meds) <- rownames(rowData)
  
  # assumption: all the provided markers are state markers
  metadata(se_meds) <- list(
    id_type_markers = rep(FALSE, length(medians)),
    id_state_markers = rep(TRUE, length(medians))
  )
  
  if(method == "DS_limma") {
    res <- testDS_limma(se, se_meds, design, contrast)
  } else if(method == "DS_LMM") {
    res <- testDS_LMM(se, se_meds, formula, contrast)
  }
  
}

rowData(res) %>% 
  as_tibble() %>%
  mutate(cluster_id = as.integer(cluster_id)) %>%
  mutate(.ri = cluster_id - 1L) %>%
  select(-cluster_id) %>%
  ctx$addNamespace() %>%
  ctx$save()
