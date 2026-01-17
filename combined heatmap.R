# ============================================================
# ONE heatmap for ALL core families  
# ============================================================

# --- Settings for the combined heatmap ---
combined_prev_min <- 0.10     # keep OTUs present in >=10% of all samples
combined_top_n    <- 300      # cap total OTUs 

pdf_width         <- 20       # BIG canvas
pdf_height        <- 12        
combined_pdf      <- file.path(out_dir, "heatmap_ALL_core_families_coreOTUs.pdf")

# --- Build the candidate OTU set (core OTUs only) ---
otu_candidates <- intersect(core_any_group, rownames(relab))

# Keep only OTUs that have a Family label
otu_to_family <- tax_df2 |>
  dplyr::select(OTU, Family) |>
  dplyr::filter(!is.na(Family), Family != "", Family != "NA")

otu_candidates <- intersect(otu_candidates, otu_to_family$OTU)

# --- Subset abundance matrix ---
mat0 <- relab[otu_candidates, , drop = FALSE]

# Prevalence filter across all samples
prev_all <- rowMeans(mat0 > 0, na.rm = TRUE)
mat1 <- mat0[names(prev_all)[prev_all >= combined_prev_min], , drop = FALSE]

if (nrow(mat1) < 2) {
  stop("Too few OTUs after combined prevalence filter. Lower combined_prev_min.")
}

# Cap to top N OTUs by prevalence (then mean abundance)
prev_all2 <- rowMeans(mat1 > 0, na.rm = TRUE)
mean_all2 <- rowMeans(mat1, na.rm = TRUE)

ord_keep <- order(prev_all2, mean_all2, decreasing = TRUE)
keep_otus <- rownames(mat1)[ord_keep][1:min(combined_top_n, nrow(mat1))]
mat_show <- mat1[keep_otus, , drop = FALSE]

# --- Log transform + clipping for display ---
mat_log <- log10(mat_show + pseudocount)
mat_log[mat_log < clip_min_log10] <- clip_min_log10
mat_log[mat_log > clip_max_log10] <- clip_max_log10

# --- Row annotation: Family ---
row_anno <- otu_to_family |>
  dplyr::filter(OTU %in% rownames(mat_log)) |>
  dplyr::distinct(OTU, .keep_all = TRUE) |>
  tibble::column_to_rownames("OTU")

row_anno <- row_anno[rownames(mat_log), , drop = FALSE]

fam_levels <- sort(unique(row_anno$Family))
fam_pal <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(fam_levels)),
  fam_levels
)

anno_colors2 <- c(anno_colors, list(Family = fam_pal))

# --- Order samples by group ---
ord_samp <- order(anno_col[colnames(mat_log), "Group"])
mat_log  <- mat_log[, ord_samp, drop = FALSE]
anno_plot <- anno_col[colnames(mat_log), , drop = FALSE]

# --- Order rows by Family ---
ord_rows <- order(row_anno$Family)
mat_log  <- mat_log[ord_rows, , drop = FALSE]
row_anno <- row_anno[ord_rows, , drop = FALSE]

# --- Safe clustering ---
cluster_rows <- nrow(mat_log) >= 2
cluster_cols <- ncol(mat_log) >= 2

# ============================================================
#   HEATMAP 
# ============================================================

pheatmap(
  mat = mat_log,
  color = hm_cols,
  breaks = seq(clip_min_log10, clip_max_log10, length.out = 201),
  
  annotation_col = anno_plot,
  annotation_row = row_anno,
  annotation_colors = anno_colors2,
  
  show_colnames = FALSE,
  
  fontsize        = 14,    
  fontsize_row    = 8,    
  fontsize_col    = 9,   
  
  border_color = NA,
  clustering_method = "average",
  cluster_rows = cluster_rows,
  cluster_cols = cluster_cols,
  filename = combined_pdf,
  width  = pdf_width,
  height = pdf_height
)

cat("\nSaved combined heatmap PDF:\n", combined_pdf, "\n", sep = "")
