# Script to generate expected test output for diffcyt operator
# This simulates what the operator does but with local test data

suppressPackageStartupMessages({
  library(dplyr, warn.conflicts = FALSE)
  library(diffcyt)
  library(SummarizedExperiment)
  library(tidyr)
})

# ============================================
# Test 2: DA_edgeR (simpler case without labels)
# ============================================
cat("Generating test 2 output (DA_edgeR)...\n")

df2 <- read.csv("test_2_in.csv", stringsAsFactors = FALSE)

# Simulate the Tercen ctx structure
# Map to crosstab view: row = cluster_id, col = sample_id, y = counts, color = group_id
df2_mapped <- df2 %>%
  mutate(
    .ci = as.integer(factor(cluster_id)) - 1L,  # row index
    .ri = 0L,  # single row for DA
    .x = sample_id,  # column
    .y = y_values  # y-axis value
  )

# Build row data (clusters)
row_dat <- df2_mapped %>%
  select(.ci, cluster_id) %>%
  unique() %>%
  arrange(.ci)

# Build counts matrix
counts <- df2_mapped %>%
  filter(.ri == 0) %>%
  pivot_wider(
    id_cols = ".ci",
    names_from = ".x",
    values_from = ".y",
    values_fn = list(.y = length),
    values_fill = 0
  ) %>%
  select(-.ci) %>%
  as.data.frame()

rowData <- row_dat %>%
  mutate(cluster_id = as.factor(as.numeric(factor(cluster_id)))) %>%
  as.data.frame()
rownames(rowData) <- 1:nrow(rowData)

# Experimental design
experiment_info <- df2 %>%
  select(sample_id, group_id) %>%
  unique() %>%
  as.data.frame()
rownames(experiment_info) <- experiment_info$sample_id

# Set reference level
lev <- unique(experiment_info$group_id)
experiment_info$group_id <- factor(experiment_info$group_id, levels = lev)

# Reorder columns
colnames(counts) <- rownames(experiment_info)
rownames(counts) <- rownames(rowData)

rowData$n_cells <- rowSums(counts)

# Create SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts)),
  rowData = rowData,
  colData = experiment_info[, "group_id", drop = FALSE]
)
colnames(se) <- rownames(experiment_info)
rownames(se) <- rownames(rowData)

# Design and contrast
design <- createDesignMatrix(experiment_info, cols_design = c("group_id"))
group_idx <- as.numeric(grepl("group_id", colnames(design)))
contrast <- createContrast(group_idx)

# Run DA_edgeR
res <- testDA_edgeR(se, design[colnames(se), ], contrast)

# Format output
df_out2 <- rowData(res) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cluster") %>%
  mutate(.ci = as.integer(cluster_id) - 1L) %>%
  mutate(group_1 = as.character(lev[1]), group_2 = as.character(lev[2])) %>%
  select(.ci, group_1, group_2, logFC, logCPM, LR, p_val, p_adj)

# Add namespace prefix
colnames(df_out2) <- c(".ci", "ds0.group_1", "ds0.group_2", "ds0.logFC", "ds0.logCPM", "ds0.LR", "ds0.p_val", "ds0.p_adj")

write.csv(df_out2, "test_2_out_1.csv", row.names = FALSE)
cat("Test 2 output written to test_2_out_1.csv\n")
print(df_out2)

# ============================================
# Test 1: DS_limma (with markers, no random effects)
# ============================================
cat("\nGenerating test 1 output (DS_limma)...\n")

df1 <- read.csv("test_1_in.csv", stringsAsFactors = FALSE)

# Simulate the Tercen ctx structure
# Map to crosstab view: row = (cluster_id, marker_id), col = sample_id, y = values, color = group_id
df1_mapped <- df1 %>%
  mutate(
    .ci = as.integer(factor(cluster_id)) - 1L,  # cluster row index
    .ri = as.integer(factor(marker_id)) - 1L,   # marker row index
    .x = sample_id,
    .y = y_values
  )

# Build row data (clusters)
row_dat1 <- df1_mapped %>%
  select(.ci, cluster_id) %>%
  unique() %>%
  arrange(.ci)

# Build counts matrix (from first row level only .ri == 0)
counts1 <- df1_mapped %>%
  filter(.ri == 0) %>%
  pivot_wider(
    id_cols = ".ci",
    names_from = ".x",
    values_from = ".y",
    values_fn = list(.y = length),
    values_fill = 0
  ) %>%
  select(-.ci) %>%
  as.data.frame()

rowData1 <- row_dat1 %>%
  mutate(cluster_id = as.factor(as.numeric(factor(cluster_id)))) %>%
  as.data.frame()
rownames(rowData1) <- 1:nrow(rowData1)

# Experimental design
experiment_info1 <- df1 %>%
  select(sample_id, group_id) %>%
  unique() %>%
  as.data.frame()
rownames(experiment_info1) <- experiment_info1$sample_id

# Set reference level
lev1 <- unique(experiment_info1$group_id)
experiment_info1$group_id <- factor(experiment_info1$group_id, levels = lev1)

# Reorder columns
colnames(counts1) <- rownames(experiment_info1)
rownames(counts1) <- rownames(rowData1)

rowData1$n_cells <- rowSums(counts1)

# Create SummarizedExperiment for counts
se1 <- SummarizedExperiment(
  assays = list(counts = as.matrix(counts1)),
  rowData = rowData1,
  colData = experiment_info1[, "group_id", drop = FALSE]
)
colnames(se1) <- rownames(experiment_info1)
rownames(se1) <- rownames(rowData1)

# Build medians list - similar to how main.R does it
list_medians <- df1_mapped %>%
  group_by(.ci, .ri, .x) %>%
  summarise(median = median(.y, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(".ri", ".ci"),
    names_from = ".x",
    values_from = "median",
    values_fill = NA
  ) %>%
  ungroup %>%
  split(.$.ri)

medians <- lapply(list_medians, function(x) {
  mat <- x %>% select(-.ci, -.ri) %>% as.data.frame()
  rownames(mat) <- rownames(rowData1)
  as.matrix(mat)
})

# Create medians SummarizedExperiment
se_meds <- SummarizedExperiment(
  assays = medians,
  rowData = rowData1,
  colData = experiment_info1[, "group_id", drop = FALSE]
)
rownames(se_meds) <- rownames(rowData1)

# Set marker metadata (all state markers for DS)
metadata(se_meds) <- list(
  id_type_markers = rep(FALSE, length(medians)),
  id_state_markers = rep(TRUE, length(medians))
)

# Design and contrast
design1 <- createDesignMatrix(experiment_info1, cols_design = c("group_id"))
group_idx1 <- as.numeric(grepl("group_id", colnames(design1)))
contrast1 <- createContrast(group_idx1)

# Run DS_limma with min_cells=1 for our small test data
res1 <- testDS_limma(se1, se_meds, design1[colnames(se_meds), ], contrast1, min_cells = 1)

# Format output
df_out1 <- rowData(res1) %>%
  as.data.frame() %>%
  mutate(.ci = as.integer(cluster_id) - 1L) %>%
  mutate(.ri = as.integer(marker_id) - 1L) %>%
  mutate(group_1 = as.character(lev1[1]), group_2 = as.character(lev1[2])) %>%
  select(.ci, .ri, group_1, group_2, logFC, p_val, p_adj)

# Add namespace prefix
colnames(df_out1) <- c(".ci", ".ri", "ds0.group_1", "ds0.group_2", "ds0.logFC", "ds0.p_val", "ds0.p_adj")

write.csv(df_out1, "test_1_out_1.csv", row.names = FALSE)
cat("Test 1 output written to test_1_out_1.csv\n")
print(df_out1)

cat("\nDone!\n")
