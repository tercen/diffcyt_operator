suppressPackageStartupMessages({
  library(tercen)
  library(tercenApi)
  library(dplyr, warn.conflicts = FALSE)
  library(diffcyt)
  library(SummarizedExperiment)
})

ctx = tercenCtx()

method <- ctx$op.value("method", as.character, "DA_edgeR")
reference.index <- ctx$op.value("reference.index", as.double, 2)

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
  ".ri", ".ci", ".y", ".x", ".colorLevels",
  all_of(first_color), all_of(second_color), all_of(first_label)
)))

row_dat <- ctx$cselect() %>%
  mutate(.ci = seq_len(nrow(.)) - 1)

counts <- df %>%
  filter(.ri == 0) %>%
  tidyr::pivot_wider(
    id_cols = ".ci",
    names_from = ".x",
    values_from = ".y",
    values_fn = list(.y = length),
    values_fill = 0
  ) %>%
  select(-.ci) %>%
  as.data.frame()

rowData <- row_dat %>%
  mutate(cluster_id = as.factor(as.numeric(factor(.[[1]])))) %>%
  as.data.frame()

### experimental design
experiment_info <- df %>% select(c(".x", all_of(first_color), all_of(first_label))) %>%
  unique() %>%
  dplyr::rename(sample_id = .x, group_id = first_color, patient_id = first_label) %>%
  as.data.frame()
rownames(experiment_info) <- experiment_info$sample_id

# Reorder factors
lev <- unique(experiment_info$group_id)
lev_idx <- c(seq_along(lev)[reference.index], seq_along(lev)[-reference.index]) 
experiment_info$group_id <- factor(
  experiment_info$group_id, levels = lev[lev_idx]
)

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

group_idx <- as.numeric(grepl("group_id", colnames(design)))
contrast <- createContrast(group_idx)

cont_all <- contrast[, rep(1, sum(contrast))]
for(i in seq_len(ncol(cont_all))) {
  cont_all[, i] <- 0
  cont_all[i + 1, i] <- 1
}

##### iterate over possible contrasts
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
  list_medians <- df %>%
    group_by(.ci, .ri, .x) %>%
    summarise(median = median(.y, na.rm = TRUE)) %>%
    # left_join(row_dat, ".ri") %>%
    tidyr::pivot_wider(
      id_cols = c(".ri", ".ci"),
      names_from = ".x",
      values_from = "median",
      values_fill = NA
    ) %>%
    ungroup %>%
    split(.$.ri)
  
  medians <- lapply(list_medians, function(x) {
    x %>% select(-.ci, -.ri) %>% as.data.frame()
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
######
df_out <- rowData(res) %>% 
  as_tibble() %>%
  mutate(.ci = as.integer(cluster_id) - 1L) %>%
  mutate_at(vars(matches("marker_id")), function(x) return(as.integer(x) - 1L)) %>%
  rename_at(vars(matches("marker_id")), function(x) ".ri") %>%
  select(-cluster_id) %>%
  ctx$addNamespace()

ctx$save(df_out)
