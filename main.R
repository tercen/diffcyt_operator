suppressPackageStartupMessages({
  library(tercen)
  library(tercenApi)
  library(dplyr, warn.conflicts = FALSE)
  library(diffcyt)
  library(SummarizedExperiment)
})

ctx = tercenCtx()

if(length(ctx$rnames) < 2) stop("Two row factors are required.")

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
  ".ri", ".ci", ".y", ".axisIndex", ".colorLevels",
  all_of(first_color), all_of(second_color), all_of(first_label)
)))

row_dat <- ctx$rselect() %>%
  mutate(.ri = seq_len(nrow(.)) - 1)

first_channel <- unique(row_dat[[2]])[1]
.ri_idx <- row_dat[row_dat[[2]] == first_channel, ] %>%
  select(.ri)

counts <- df %>%
  filter(.axisIndex == 0) %>%
  inner_join(.ri_idx, ".ri") %>%
  tidyr::pivot_wider(
    id_cols = ".ri",
    names_from = ".ci",
    values_from = ".y",
    values_fn = list(.y = length),
    values_fill = 0
  ) %>%
  select(-.ri) %>%
  as.data.frame()


colData <- ctx$cselect() %>% 
  tidyr::unite("sample_id", sep = " - ") %>%
  mutate(.ci = seq_len(nrow(.)) - 1)
rowData <- row_dat[row_dat[[2]] == first_channel, ] %>%
  mutate(cluster_id = as.factor(as.numeric(factor(.[[1]])))) %>%
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
  
  
  # colnames(row_dat)[1] <- "cluster_label"
  colnames(row_dat)[2] <- "marker"
  
  list_medians <- df %>%
    group_by(.ci, .ri) %>%
    summarise(median = median(.y, na.rm = TRUE)) %>%
    left_join(row_dat, ".ri") %>%
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
