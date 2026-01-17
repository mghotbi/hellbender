# ============================================================
# OTU-level core + Chitinophagales/Chitinophagaceae OTU heatmap
# ============================================================
# Author: Mitra Ghotbi
# Purpose:
#  - Build phyloseq object from OTU counts + taxonomy + metadata
#  - Core OTUs per group (prevalence)
#  - OTUs drive Chitinophagales/Chitinophagaceae patterns?
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(phyloseq)
  library(pheatmap)
  library(RColorBrewer)
  library(stringr)
})

# ---------------------------
# 0) User settings (NO setwd)
# ---------------------------
otu_count_file <- "250508_16S_OTU_count_table_decontam_Lluvia_network.csv"
tax_table_file <- "250508_16S_taxonomy_table_decontam_Lluvia_network.csv"
metadata_file  <- "250508_16S_metadata_decontam_Lluvia_network.csv"

out_dir <- file.path(getwd(), "sons_reviewer_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Core definition
prev_core_threshold <- 0.30  # e.g. core = present in >=30% of samples within group

# Heatmap selection strategy (choose ONE)
heatmap_strategy <- "topN"   # "prevalence" or "topN"

# If strategy == "prevalence"
prev_heatmap_filter <- 0.05  # keep OTUs present in >=5% of samples overall (within subset)

# If strategy == "topN"
topN_otus <- 50              # show top N OTUs (by prevalence, tie-broken by variance)

# Transform settings for display
transform_type <- "log10p"   # "log10p" (log10(relab + pseudocount)) or "log10_no_pseudo" or "clr_like"
pseudocount    <- 1e-20       # only used for "log10p"
clip_min_log10 <- -6
clip_max_log10 <- -2

target_order  <- "Chitinophagales"
target_family <- "Chitinophagaceae"

focus_groups <- c("Chattanooga Zoo", "East TN Wild")   
focus_prev   <- 0.50

# ---------------------------
# 1) Read files
# ---------------------------
otu_count <- fread(otu_count_file)
tax_table <- fread(tax_table_file)
metadata  <- fread(metadata_file)

# ---------------------------
# 2) OTU matrix (OTUs x samples)
# ---------------------------
stopifnot("index" %in% names(otu_count))

otu_mat_samples_x_otus <- otu_count |>
  as.data.frame() |>
  column_to_rownames("index") |>
  as.matrix()

mode(otu_mat_samples_x_otus) <- "numeric"
otu_mat <- t(otu_mat_samples_x_otus)  # OTUs x samples

# ---------------------------
# 3) Taxonomy matrix (OTUs x ranks)
# ---------------------------
stopifnot("OTU" %in% names(tax_table))
tax_df <- as.data.frame(tax_table)

known_ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
rank_cols_present <- intersect(known_ranks, names(tax_df))
if (length(rank_cols_present) < 2) {
  stop("Taxonomy table missing expected rank columns. Found: ",
       paste(names(tax_df), collapse = ", "))
}

tax_mat <- tax_df |>
  select(OTU, all_of(rank_cols_present)) |>
  column_to_rownames("OTU") |>
  as.matrix()

colnames(tax_mat) <- str_to_title(colnames(tax_mat))

# ---------------------------
# 4) Metadata (samples x variables)
# ---------------------------
stopifnot("index" %in% names(metadata))
meta_df <- metadata |>
  as.data.frame() |>
  column_to_rownames("index")

# ---------------------------
# 5) Align samples + OTUs
# ---------------------------
common_samples <- intersect(colnames(otu_mat), rownames(meta_df))
if (length(common_samples) < 2) stop("Too few overlapping samples between OTU and metadata.")
otu_mat <- otu_mat[, common_samples, drop = FALSE]
meta_df <- meta_df[common_samples, , drop = FALSE]

common_otus <- intersect(rownames(otu_mat), rownames(tax_mat))
if (length(common_otus) < 2) stop("Too few overlapping OTUs between OTU and taxonomy.")
otu_mat <- otu_mat[common_otus, , drop = FALSE]
tax_mat <- tax_mat[common_otus, , drop = FALSE]

# ---------------------------
# 6) Build phyloseq
# ---------------------------
ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta_df)
)

cat("\nphyloseq object:\n"); print(ps)
cat("\nrank_names(ps):\n"); print(rank_names(ps))
cat("\nsample_variables(ps):\n"); print(sample_variables(ps))

# ---------------------------
# 7) Build group_broad if possible
# ---------------------------
needed_cols <- c("env_broad_scale", "site", "ecoregion_III", "bh_id_prefix")

if (all(needed_cols %in% sample_variables(ps))) {
  meta2 <- as(sample_data(ps), "data.frame") |>
    mutate(
      group_broad = case_when(
        env_broad_scale == "wild" &
          site == "Hiwassee River Picnic Area" ~ "East TN Wild",
        env_broad_scale == "environmental" &
          site == "Hiwassee River Picnic Area" ~ "East TN Environmental",
        env_broad_scale == "wild" &
          ecoregion_III == "Interior Plateau" &
          bh_id_prefix %in% c("BSC", "SB") ~ "Middle TN Wild",
        env_broad_scale == "environmental" &
          ecoregion_III == "Interior Plateau" ~ "Middle TN Environmental",
        site == "Chattanooga Zoo" ~ "Chattanooga Zoo",
        site == "Nashville Zoo" ~ "Nashville Zoo",
        TRUE ~ "Middle TN Recapture"
      )
    )
  sample_data(ps) <- sample_data(meta2)
  cat("\nAdded group_broad.\n")
}

# Choose grouping variable
if ("group_broad" %in% sample_variables(ps)) {
  group_var <- "group_broad"
} else if ("Category" %in% sample_variables(ps)) {
  group_var <- "Category"
} else if ("group" %in% sample_variables(ps)) {
  group_var <- "group"
} else {
  stop("No grouping variable found (group_broad / Category / group).")
}
cat("\nUsing grouping variable:", group_var, "\n")

# ---------------------------
# 8) Optional filtering
# ---------------------------
if ("Otu00001" %in% taxa_names(ps)) {
  ps <- prune_taxa(taxa_names(ps) != "Otu00001", ps)
  cat("Removed Otu00001\n")
}
if ("Genus" %in% rank_names(ps)) {
  ps <- subset_taxa(ps, is.na(Genus) | Genus != "Escherichia")
  cat("Removed Genus == Escherichia (if annotated)\n")
}

# ---------------------------
# 9) Core OTUs per group (prevalence)
# ---------------------------
otu_pa <- as(otu_table(ps) > 0, "matrix") # OTUs x samples
meta   <- as(sample_data(ps), "data.frame")

groups <- sort(unique(meta[[group_var]]))
groups <- groups[!is.na(groups)]

core_by_group <- function(group_name, prev = prev_core_threshold) {
  samp <- rownames(meta)[meta[[group_var]] == group_name]
  if (length(samp) == 0) return(character(0))
  mat <- otu_pa[, samp, drop = FALSE]
  prev_vec <- rowMeans(mat)
  names(prev_vec)[prev_vec >= prev]
}

core_list <- lapply(groups, core_by_group)
names(core_list) <- groups
core_shared <- Reduce(intersect, core_list)

cat("\nCore OTU counts per group (prev >=", prev_core_threshold, "):\n")
print(sapply(core_list, length))
cat("\nShared core OTUs across ALL groups:", length(core_shared), "\n")

# Save core shared taxonomy (nice supplement file)
tax <- as(tax_table(ps), "matrix")
core_shared_df <- data.frame(OTU = core_shared) |>
  left_join(as.data.frame(tax) |> rownames_to_column("OTU"), by = "OTU")
fwrite(core_shared_df, file.path(out_dir, "supp_table_core_shared_otu_taxonomy.tsv"), sep = "\t")

# ---------------------------
# 10) Target OTUs (Chitinophagales/Chitinophagaceae)
# ---------------------------
tax <- as(tax_table(ps), "matrix")
stopifnot(all(c("Order","Family") %in% colnames(tax)))

chitin_otus <- rownames(tax)[
  (!is.na(tax[, "Order"])  & tax[, "Order"]  == target_order) |
  (!is.na(tax[, "Family"]) & tax[, "Family"] == target_family)
]
chitin_otus <- intersect(chitin_otus, taxa_names(ps))
if (length(chitin_otus) == 0) stop("No OTUs matched target order/family.")

cat("\nTarget OTUs found:", length(chitin_otus), "\n")

# Prevalence by group for target OTUs
prev_by_group <- sapply(groups, function(g) {
  samp <- rownames(meta)[meta[[group_var]] == g]
  mat  <- otu_pa[chitin_otus, samp, drop = FALSE]
  rowMeans(mat)
})
prev_by_group <- as.matrix(prev_by_group)
colnames(prev_by_group) <- groups

prev_out <- as.data.frame(prev_by_group) |>
  rownames_to_column("OTU") |>
  left_join(as.data.frame(tax) |> rownames_to_column("OTU"), by = "OTU")

prev_file <- file.path(out_dir, "chitinophagales_chitinophagaceae_prevalence_by_group.tsv")
fwrite(prev_out, prev_file, sep = "\t")
cat("\nWrote prevalence table:\n", prev_file, "\n", sep="")

# Optional: OTUs shared between two focus groups
if (all(focus_groups %in% colnames(prev_by_group))) {
  keep <- prev_by_group[, focus_groups[1]] >= focus_prev &
    prev_by_group[, focus_groups[2]] >= focus_prev
  shared_focus <- rownames(prev_by_group)[keep]
  fwrite(
    data.frame(OTU = shared_focus),
    file.path(out_dir, "focus_shared_chitin_otus.tsv"),
    sep = "\t"
  )
  cat("\nOTUs with prevalence >=", focus_prev, " in BOTH ",
      focus_groups[1], " and ", focus_groups[2], ": ",
      length(shared_focus), "\n", sep = "")
}

# ---------------------------
# 11) Heatmap matrix (relative abundance -> transform)
# ---------------------------
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
relab  <- as(otu_table(ps_rel), "matrix")  # OTUs x samples

mat_sub <- relab[chitin_otus, , drop = FALSE]

# Choose OTUs to show (prevents your “one red stripe” weirdness)
prev_all <- rowMeans(mat_sub > 0)
var_all  <- apply(mat_sub, 1, var)

if (heatmap_strategy == "prevalence") {
  keep_otus <- names(prev_all)[prev_all >= prev_heatmap_filter]
  mat_sub  <- mat_sub[keep_otus, , drop = FALSE]
} else if (heatmap_strategy == "topN") {
  # rank by prevalence, then variance
  rank_df <- data.frame(OTU = rownames(mat_sub), prev = prev_all, var = var_all)
  rank_df <- rank_df[order(-rank_df$prev, -rank_df$var), ]
  keep_otus <- head(rank_df$OTU, min(topN_otus, nrow(rank_df)))
  mat_sub <- mat_sub[keep_otus, , drop = FALSE]
} else {
  stop("heatmap_strategy must be 'prevalence' or 'topN'")
}

cat("\nHeatmap OTUs selected:", nrow(mat_sub), "\n")
if (nrow(mat_sub) == 0) stop("No OTUs remain for heatmap after filtering/selection.")

# Transform for display
if (transform_type == "log10p") {
  mat_disp <- log10(mat_sub + pseudocount)
  mat_disp[mat_disp < clip_min_log10] <- clip_min_log10
  mat_disp[mat_disp > clip_max_log10] <- clip_max_log10
  legend_breaks <- c(-6, -5, -4, -3, -2)
  legend_labels <- c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2")
  subtitle <- paste0("Values = log10(relative abundance + ", pseudocount, ")")
} else if (transform_type == "log10_no_pseudo") {
  # Will set zeros to NA (otherwise -Inf)
  mat_disp <- log10(mat_sub)
  mat_disp[!is.finite(mat_disp)] <- NA
  subtitle <- "Values = log10(relative abundance); zeros shown as NA"
  legend_breaks <- NA
  legend_labels <- NA
} else if (transform_type == "clr_like") {
  # robust alternative (no pseudocount shown explicitly): log(x) - mean(log(x)) per sample
  # still needs a tiny offset internally to avoid log(0), but interpretation is “centered log”
  eps <- 1e-6
  tmp <- log(mat_sub + eps)
  mat_disp <- sweep(tmp, 2, colMeans(tmp), "-")
  subtitle <- "Values = centered log (CLR-like) of relab (internal eps=1e-6)"
  legend_breaks <- NA
  legend_labels <- NA
} else {
  stop("transform_type must be 'log10p', 'log10_no_pseudo', or 'clr_like'")
}

# ---------------------------
# 12) Column annotation
# ---------------------------
anno_col <- meta[colnames(mat_disp), , drop = FALSE] |> select(all_of(group_var))
colnames(anno_col) <- "Group"

group_levels <- sort(unique(anno_col$Group))
pal <- setNames(
  colorRampPalette(brewer.pal(min(8, length(group_levels)), "Set2"))(length(group_levels)),
  group_levels
)
anno_colors <- list(Group = pal)

# ---------------------------
# 13) Heatmap aesthetics
# ---------------------------
hm_cols <- colorRampPalette(
  c("#08306B", "#2171B5", "#6BAED6", "#FEE090", "#FC8D59", "#D73027")
)(200)

cluster_rows <- nrow(mat_disp) >= 2
cluster_cols <- ncol(mat_disp) >= 2

# ---------------------------
# 14) Save heatmap PDF reliably
# ---------------------------
pdf_file <- file.path(out_dir, "heatmap_chitinophagales_chitinophagaceae_otus.pdf")

pheatmap(
  mat = mat_disp,
  color = hm_cols,
  annotation_col = anno_col,
  annotation_colors = anno_colors,
  show_colnames = FALSE,
  fontsize_row = 7,
  border_color = NA,
  clustering_method = "average",
  cluster_rows = cluster_rows,
  cluster_cols = cluster_cols,
  main = paste0(target_order, " / ", target_family, " OTUs"),
  filename = pdf_file,
  width = 12,
  height = 8,
  legend_breaks = if (all(is.na(legend_breaks))) NULL else legend_breaks,
  legend_labels = if (all(is.na(legend_labels))) NULL else legend_labels
)

cat("\nSaved heatmap PDF to:\n", pdf_file, "\n", sep = "")
cat("\nSubtitle:\n", subtitle, "\n", sep = "")

# ---------------------------
# 15) Explanation  
# ---------------------------
writeLines(
  c(
    " note:",
    "- We repeated order/family-level patterns at OTU resolution by extracting all OTUs assigned to the target order/family.",
    "- We computed prevalence (fraction of samples with OTU > 0) within each group and exported a prevalence-by-group table.",
    "- The heatmap visualizes OTU-level relative abundances across samples, using a transform chosen for interpretability and sparsity.",
    "",
    paste0("Transform used: ", transform_type),
    subtitle,
    "",
    "If the plot looks 'blocky' or dominated by one OTU, use heatmap_strategy='topN' (recommended) to show multiple informative OTUs."
  ),
  file.path(out_dir, "heatmap_explanation_for_reviewer.txt")
)
