# ============================================================
# Core-family OTU heatmaps (ALL core families)
# ============================================================
# Author: Mitra Ghotbi
# What it does:
#  1) Builds a phyloseq object from your CSVs/rds
#  2) Computes OTU “core” per group (prevalence threshold)
#  3) Identifies “core families” = families that contain >=1 core OTU
#  4) For EACH core family, saves a PDF heatmap of OTU-level relative abundance
#     (log10(relab + pseudocount)) with dendrograms  
#
# Outputs (written into out_dir):
#  - core_family_summary.tsv
#  - heatmap_core_family_<FamilyName>_otus.pdf   (one per family)
#
# Notes:
#  - To run: put the three CSVs in your working directory 
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
# 0) User settings
# ---------------------------
otu_count_file <- "250508_16S_OTU_count_table_decontam_Lluvia_network.csv"
tax_table_file <- "250508_16S_taxonomy_table_decontam_Lluvia_network.csv"
metadata_file  <- "250508_16S_metadata_decontam_Lluvia_network.csv"

# Output directory
out_dir <- file.path(getwd(), "sons_reviewer_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Core definition (OTU-level)
prev_core_threshold <- 0.30   # OTU is "core" in a group if present in >=30% of that group's samples

# Heatmap display choices
pseudocount    <- 1e-20       # needed for log10; cannot log10(0)
clip_min_log10 <- -6
clip_max_log10 <- -2

# Which OTUs to include in each family heatmap
# - "all": all OTUs in the family (can be huge)
# - "core_only": only OTUs that are core (in at least one group)
otus_in_heatmap <- "core_only"   # recommended for reviewer response

# Further limit OTUs for readability (top by mean relative abundance across all samples)
top_n_otus_per_family <- 50       # set NA to disable and plot all selected OTUs

# Optional: exclude a problematic OTU / genus
drop_otu_id <- "Otu00001"         # set NULL if you don't want this filter
drop_genus  <- "Escherichia"      # set NULL if you don't want this filter

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
otu_mat <- t(otu_mat_samples_x_otus) # OTUs x samples

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

colnames(tax_mat) <- str_to_title(colnames(tax_mat))  # kingdom -> Kingdom etc.

# ---------------------------
# 4) Metadata matrix (samples x variables)
# ---------------------------
stopifnot("index" %in% names(metadata))
meta_df <- metadata |>
  as.data.frame() |>
  column_to_rownames("index")

# ---------------------------
# 5) Align samples + OTUs across tables
# ---------------------------
common_samples <- intersect(colnames(otu_mat), rownames(meta_df))
if (length(common_samples) < 2) stop("Too few overlapping samples between OTU table and metadata.")
otu_mat <- otu_mat[, common_samples, drop = FALSE]
meta_df <- meta_df[common_samples, , drop = FALSE]

common_otus <- intersect(rownames(otu_mat), rownames(tax_mat))
if (length(common_otus) < 2) stop("Too few overlapping OTUs between OTU table and taxonomy.")
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
# 7) Build group_broad if possible (optional)
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
# 8) Optional filtering (keep your behavior, but safe)
# ---------------------------
if (!is.null(drop_otu_id) && drop_otu_id %in% taxa_names(ps)) {
  ps <- prune_taxa(taxa_names(ps) != drop_otu_id, ps)
  cat("Removed OTU: ", drop_otu_id, "\n", sep = "")
}

if (!is.null(drop_genus) && "Genus" %in% rank_names(ps)) {
  ps <- subset_taxa(ps, is.na(Genus) | Genus != drop_genus)
  cat("Removed Genus == ", drop_genus, " (where annotated)\n", sep = "")
}

# ---------------------------
# 9) Core OTUs per group (prevalence) DspikeIn
# ---------------------------
otu_pa <- as(otu_table(ps) > 0, "matrix")  # OTUs x samples
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

core_any_group <- sort(unique(unlist(core_list)))   # OTUs core in >=1 group
core_all_groups <- Reduce(intersect, core_list)     # OTUs core in ALL groups

cat("\nCore OTU counts per group (prev >=", prev_core_threshold, "):\n")
print(sapply(core_list, length))
cat("\nCore OTUs in >=1 group:", length(core_any_group), "\n")
cat("Core OTUs in ALL groups:", length(core_all_groups), "\n")

# ---------------------------
# 10) Identify CORE FAMILIES (families containing >=1 core OTU)
# ---------------------------
tax <- as(tax_table(ps), "matrix")
if (!("Family" %in% colnames(tax))) stop("Tax table has no 'Family' column.")

tax_df2 <- as.data.frame(tax) |>
  rownames_to_column("OTU")

# Core OTUs -> families
core_tax <- tax_df2 |>
  filter(OTU %in% core_any_group) |>
  filter(!is.na(Family), Family != "", Family != "NA")

core_families <- sort(unique(core_tax$Family))
cat("\nNumber of core families (contain >=1 core OTU): ", length(core_families), "\n", sep = "")

if (length(core_families) == 0) {
  stop("No core families found (check Family annotations or core threshold).")
}

# Write a simple summary table you can cite in SI
core_family_summary <- core_tax |>
  count(Family, name = "n_core_otus_in_family") |>
  arrange(desc(n_core_otus_in_family))

summary_file <- file.path(out_dir, "core_family_summary.tsv")
fwrite(core_family_summary, summary_file, sep = "\t")
cat("\nWrote: ", summary_file, "\n", sep = "")

# ---------------------------
# 11) Relative abundance matrix (OTUs x samples)
# ---------------------------
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
relab  <- as(otu_table(ps_rel), "matrix")  # OTUs x samples

# ---------------------------
# 12) Sample annotation colors (stable + safe)
# ---------------------------
# IMPORTANT: anno_col rownames must match matrix column names EXACTLY
anno_col <- meta[colnames(relab), , drop = FALSE] |>
  select(all_of(group_var))

colnames(anno_col) <- "Group"
anno_col$Group <- as.character(anno_col$Group)

group_levels <- sort(unique(anno_col$Group))
pal <- setNames(
  colorRampPalette(brewer.pal(min(8, length(group_levels)), "Set2"))(length(group_levels)),
  group_levels
)
anno_colors <- list(Group = pal)

# Heatmap color palette
hm_cols <- colorRampPalette(c("#08306B", "#2171B5", "#6BAED6", "#FEE090", "#FC8D59", "#D73027"))(200)

# Legend (in relative abundance scale)
breaks_log <- c(-6, -5, -4, -3, -2)
labels_rel <- c("1e-6", "1e-5", "1e-4", "1e-3", "1e-2")

# ---------------------------
# 13) Helper: safe filename
# ---------------------------
safe_name <- function(x) {
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-]+", "", x)
  x
}

# ---------------------------
# 14) Loop through core families & save heatmaps
# ---------------------------
for (fam in core_families) {
  
  fam_safe <- safe_name(fam)
  
  # OTUs in this family
  fam_otus_all <- tax_df2$OTU[!is.na(tax_df2$Family) & tax_df2$Family == fam]
  fam_otus_all <- intersect(fam_otus_all, rownames(relab))
  
  if (length(fam_otus_all) == 0) {
    cat("Skipping family (no OTUs after alignment): ", fam, "\n", sep = "")
    next
  }
  
  # Decide which OTUs to show
  if (otus_in_heatmap == "core_only") {
    fam_otus <- intersect(fam_otus_all, core_any_group)
  } else {
    fam_otus <- fam_otus_all
  }
  
  if (length(fam_otus) == 0) {
    cat("Skipping family (no selected OTUs): ", fam, "\n", sep = "")
    next
  }
  
  # Subset matrix: OTUs x samples
  mat_sub <- relab[fam_otus, , drop = FALSE]
  
  # Optionally-> keep top N by mean relative abundance  
  if (!is.na(top_n_otus_per_family) && nrow(mat_sub) > top_n_otus_per_family) {
    means <- rowMeans(mat_sub, na.rm = TRUE)
    keep <- names(sort(means, decreasing = TRUE))[1:top_n_otus_per_family]
    mat_sub <- mat_sub[keep, , drop = FALSE]
  }
  
  # Log transform (cannot remove pseudocount if you keep log10)
  mat_log <- log10(mat_sub + pseudocount)
  
  # Clip for consistent color scaling
  mat_log[mat_log < clip_min_log10] <- clip_min_log10
  mat_log[mat_log > clip_max_log10] <- clip_max_log10
  
  # SAFE clustering (avoid hclust errors)
  cluster_rows <- nrow(mat_log) >= 2
  cluster_cols <- ncol(mat_log) >= 2
  
  # (Optional) order samples by group for readability
  ord <- order(anno_col[colnames(mat_log), "Group"])
  mat_log <- mat_log[, ord, drop = FALSE]
  
  # IMPORTANT: keep annotation aligned after reordering
  anno_plot <- anno_col[colnames(mat_log), , drop = FALSE]
  
  # Output PDF
  pdf_file <- file.path(out_dir, paste0("heatmap_core_family_", fam_safe, "_otus.pdf"))
  
   pheatmap(
    mat = mat_log,
    color = hm_cols,
    breaks = seq(clip_min_log10, clip_max_log10, length.out = 201),
    annotation_col = anno_plot,
    annotation_colors = anno_colors,
    show_colnames = FALSE,
    fontsize_row = 7,
    border_color = NA,
    clustering_method = "average",
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    main = paste0(
      "Core family: ", fam, "\n",
      "OTUs shown: ", nrow(mat_log),
      ifelse(!is.na(top_n_otus_per_family), paste0(" (top ", min(top_n_otus_per_family, nrow(mat_sub)), ")"), ""),
      "\nlog10(relative abundance + ", pseudocount, ")"
    ),
    legend_breaks = breaks_log,
    legend_labels = labels_rel,
    filename = pdf_file,
    width = 12,
    height = 8
  )
  
  cat("Saved family heatmap: ", pdf_file, "\n", sep = "")
}

cat("\nDONE.\nAll outputs are in:\n", out_dir, "\n", sep = "")

# ---------------------------
# 15) Quick explanation of Methods 
# ---------------------------
cat("\n--- Explanation (copy/paste) ---\n")
cat(paste0(
  "For each bacterial family containing at least one prevalence-defined core OTU, we generated OTU-level heatmaps.\n",
  "Relative abundances were computed per sample and visualized as log10(p + ", pseudocount, "), where p is OTU relative abundance.\n",
  "Because p < 1, log10(p) is negative (e.g., 1e-2 -> -2; 1e-6 -> -6). The pseudocount avoids log10(0).\n",
  "Row/column clustering is automatically disabled if fewer than 2 OTUs or 2 samples are available (hclust requires n >= 2).\n"
))
