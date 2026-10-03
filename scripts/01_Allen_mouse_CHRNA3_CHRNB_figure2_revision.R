# ============================================================
# Allen Institute Mouse Whole Cortex and Hippocampus 10x
# Chrna3 / Chrnb2 co-expression analysis
# Yao et al., 2021
#
# Project:
# Allen_Mouse_CHRNA3_CHRNB2
#
# Goal:
# Identify Chrna3 / Chrnb2 co-expression in adult mouse
# hippocampal GABAergic cells and examine inhibitory
# neuron subpopulations.
# ============================================================


# -----------------------------
# 1. Load required packages
# -----------------------------

library(data.table)
library(rhdf5)
library(ggplot2)
library(patchwork)
library(ggtext)


# -----------------------------
# 2. Define file locations
# -----------------------------

metadata_file <- "data/metadata.csv"
expression_file <- "data/expression_matrix.hdf5"


# -----------------------------
# 3. Inspect metadata structure
# -----------------------------

meta_header <- fread(
  metadata_file,
  nrows = 0
)

names(meta_header)


# -----------------------------
# 4. Load only metadata columns
#    needed for this analysis
# -----------------------------

cols_needed <- c(
  "sample_name",
  "external_donor_name_label",
  "donor_sex_label",
  "region_label",
  "class_label",
  "subclass_label",
  "cluster_label",
  "neighborhood_label",
  "cell_type_alias_label",
  "cell_type_designation_label",
  "cortical_layer_label"
)

meta <- fread(
  metadata_file,
  select = cols_needed
)

dim(meta)


# -----------------------------
# 5. Restrict to hippocampal
#    GABAergic cells
# -----------------------------

hip_gaba <- meta[
  region_label == "HIP" &
    class_label == "GABAergic"
]

dim(hip_gaba)


# -----------------------------
# 6. Summarize inhibitory
#    subpopulations and donors
# -----------------------------

sort(
  table(hip_gaba$subclass_label),
  decreasing = TRUE
)

length(
  unique(hip_gaba$external_donor_name_label)
)

table(
  hip_gaba$donor_sex_label
)

table(
  hip_gaba$external_donor_name_label,
  hip_gaba$donor_sex_label
)

table(
  hip_gaba$subclass_label,
  hip_gaba$external_donor_name_label
)


# -----------------------------
# 7. Inspect gene names in the
#    HDF5 expression matrix
# -----------------------------

genes <- h5read(
  expression_file,
  "/data/gene"
)

gene_idx <- match(
  c("Chrna3", "Chrnb2"),
  genes
)

gene_idx
genes[gene_idx]


# Safety check:
# Stop if either gene is not found

if (any(is.na(gene_idx))) {
  stop("Chrna3 and/or Chrnb2 could not be found in the HDF5 gene list.")
}


# -----------------------------
# 8. Read HDF5 cell identifiers
# -----------------------------

h5_samples <- h5read(
  expression_file,
  "/data/samples"
)

length(h5_samples)


# -----------------------------
# 9. Match hippocampal
#    GABAergic cells to HDF5
#    expression-matrix rows
# -----------------------------

sample_idx <- match(
  hip_gaba$sample_name,
  h5_samples
)

sum(is.na(sample_idx))
length(sample_idx)


# Safety check:
# Every selected metadata cell should match the HDF5 matrix

if (any(is.na(sample_idx))) {
  stop(
    paste(
      sum(is.na(sample_idx)),
      "hippocampal GABAergic cells did not match the HDF5 sample list."
    )
  )
}


# -----------------------------
# 10. Final checkpoint before
#     extracting expression data
# -----------------------------

cat("\nAnalysis checkpoint:\n")
cat("Hippocampal GABAergic cells:", nrow(hip_gaba), "\n")
cat("Matched HDF5 cells:", length(sample_idx), "\n")
cat("Unmatched cells:", sum(is.na(sample_idx)), "\n")
cat("Chrna3 gene index:", gene_idx[1], "\n")
cat("Chrnb2 gene index:", gene_idx[2], "\n")
# -----------------------------
# 11. Extract Chrna3 and Chrnb2
#     raw counts for selected cells
# -----------------------------

expr_counts <- h5read(
  expression_file,
  "/data/counts",
  index = list(
    sample_idx,
    gene_idx
  )
)

dim(expr_counts)

# Give the two expression columns informative names

colnames(expr_counts) <- c(
  "Chrna3_count",
  "Chrnb2_count"
)

# Inspect the first few cells

head(expr_counts)
# -----------------------------
# 12. Verify expression-row
#     alignment before merging
# -----------------------------

# Read the first five selected cells individually
# and compare them with the bulk extraction.

alignment_check <- t(
  sapply(
    1:5,
    function(i) {
      h5read(
        expression_file,
        "/data/counts",
        index = list(
          sample_idx[i],
          gene_idx
        )
      )
    }
  )
)

colnames(alignment_check) <- c(
  "Chrna3_count",
  "Chrnb2_count"
)

alignment_check

expr_counts[1:5, ]

identical(
  unname(alignment_check),
  unname(expr_counts[1:5, ])
)
# -----------------------------
# 13. Combine expression counts
#     with cell metadata
# -----------------------------

mouse <- copy(hip_gaba)

mouse[, Chrna3_count := expr_counts[, "Chrna3_count"]]
mouse[, Chrnb2_count := expr_counts[, "Chrnb2_count"]]

# Verify dimensions

dim(mouse)

# Inspect first few rows

mouse[
  1:10,
  .(
    sample_name,
    subclass_label,
    Chrna3_count,
    Chrnb2_count
  )
]


# -----------------------------
# 14. Classify cells by
#     Chrna3 / Chrnb2 detection
# -----------------------------

mouse[
  ,
  expression_group := fifelse(
    Chrna3_count > 0 & Chrnb2_count > 0,
    "Double-positive",
    fifelse(
      Chrna3_count > 0 & Chrnb2_count == 0,
      "Chrna3 only",
      fifelse(
        Chrna3_count == 0 & Chrnb2_count > 0,
        "Chrnb2 only",
        "Neither"
      )
    )
  )
]


# -----------------------------
# 15. Overall expression
#     categories
# -----------------------------

overall_counts <- table(mouse$expression_group)

overall_counts

round(
  100 * prop.table(overall_counts),
  3
)


# -----------------------------
# 16. Expression categories
#     by inhibitory subclass
# -----------------------------

subclass_summary <- mouse[
  ,
  .(
    total_cells = .N,
    Chrna3_only = sum(expression_group == "Chrna3 only"),
    Chrnb2_only = sum(expression_group == "Chrnb2 only"),
    double_positive = sum(expression_group == "Double-positive"),
    neither = sum(expression_group == "Neither")
  ),
  by = subclass_label
]

subclass_summary[
  ,
  pct_double_positive :=
    round(100 * double_positive / total_cells, 3)
]

setorder(
  subclass_summary,
  -pct_double_positive
)

subclass_summary
# -----------------------------
# 17. Double-positive frequency
#     by donor and subclass
# -----------------------------

donor_subclass_summary <- mouse[
  ,
  .(
    total_cells = .N,
    double_positive =
      sum(expression_group == "Double-positive")
  ),
  by = .(
    external_donor_name_label,
    donor_sex_label,
    subclass_label
  )
]

donor_subclass_summary[
  ,
  pct_double_positive :=
    round(
      100 * double_positive / total_cells,
      3
    )
]

# Show all donor/subclass combinations containing
# at least one double-positive cell

donor_subclass_summary[
  double_positive > 0
]

# Examine the three subclasses with the strongest signal

donor_subclass_summary[
  subclass_label == "Pvalb"
]

donor_subclass_summary[
  subclass_label == "Sncg"
]

donor_subclass_summary[
  subclass_label == "Lamp5"
]

# Save the current subclass summary

fwrite(
  subclass_summary,
  "results/mouse_HIP_GABA_subclass_summary.csv"
)
# -----------------------------
# 18. Fine transcriptomic clusters
#     containing double-positive cells
# -----------------------------

cluster_summary <- mouse[
  ,
  .(
    total_cells = .N,
    double_positive =
      sum(expression_group == "Double-positive")
  ),
  by = .(
    subclass_label,
    cluster_label
  )
]

cluster_summary[
  ,
  pct_double_positive :=
    round(
      100 * double_positive / total_cells,
      3
    )
]

# Show only clusters containing at least one
# Chrna3/Chrnb2 double-positive cell

dp_cluster_summary <- cluster_summary[
  double_positive > 0
]

setorder(
  dp_cluster_summary,
  subclass_label,
  -pct_double_positive
)

dp_cluster_summary
# -----------------------------
# 19. Donor distribution of
#     Pvalb-Vipr2 double-positive cells
# -----------------------------

mouse[
  subclass_label == "Pvalb" &
    expression_group == "Double-positive",
  .N,
  by = .(
    cluster_label,
    external_donor_name_label,
    donor_sex_label
  )
][
  order(
    cluster_label,
    external_donor_name_label
  )
]
# -----------------------------
# 20. Statistical comparison of
#     double-positive frequency
#     across inhibitory subclasses
# -----------------------------

subclass_2xk <- table(
  mouse$subclass_label,
  mouse$expression_group == "Double-positive"
)

subclass_2xk

set.seed(12345)

subclass_fisher <- fisher.test(
  subclass_2xk,
  simulate.p.value = TRUE,
  B = 100000
)

subclass_fisher


# -----------------------------
# 21. Focused enrichment tests
# -----------------------------

# Pvalb-Vipr2 versus all other hippocampal GABAergic cells

mouse[
  ,
  Pvalb_Vipr2 := subclass_label == "Pvalb" &
    cluster_label %in% c(
      "122_Pvalb Vipr2",
      "123_Pvalb Vipr2"
    )
]

pvalb_vipr2_table <- table(
  mouse$Pvalb_Vipr2,
  mouse$expression_group == "Double-positive"
)

pvalb_vipr2_table

fisher.test(pvalb_vipr2_table)


# Lamp5-Lhx6 versus other Lamp5 cells

lamp5_only <- mouse[
  subclass_label == "Lamp5"
]

lamp5_only[
  ,
  Lamp5_Lhx6 := grepl(
    "Lamp5 Lhx6",
    cluster_label
  )
]

lamp5_lhx6_table <- table(
  lamp5_only$Lamp5_Lhx6,
  lamp5_only$expression_group == "Double-positive"
)

lamp5_lhx6_table

fisher.test(lamp5_lhx6_table)
# -----------------------------
# 22. Create final mouse
#     subclass results table
# -----------------------------

mouse_subclass_results <- mouse[
  ,
  .(
    total_cells = .N,
    double_positive =
      sum(expression_group == "Double-positive"),
    pct_double_positive =
      round(
        100 * sum(expression_group == "Double-positive") / .N,
        3
      ),
    total_donors =
      uniqueN(external_donor_name_label),
    donors_with_double_positive =
      uniqueN(
        external_donor_name_label[
          expression_group == "Double-positive"
        ]
      )
  ),
  by = subclass_label
]

setorder(
  mouse_subclass_results,
  -pct_double_positive
)

mouse_subclass_results

fwrite(
  mouse_subclass_results,
  "results/mouse_HIP_GABA_final_subclass_results.csv"
)
# -----------------------------
# 23. Chrna3 : Chrnb2 ratios
#     in double-positive cells
# -----------------------------

double_positive_cells <- mouse[
  expression_group == "Double-positive"
]

# Safety check
dim(double_positive_cells)

# Calculate continuous log2 expression ratio
# Positive values = more Chrna3
# Negative values = more Chrnb2
# Zero = equal counts

double_positive_cells[
  ,
  log2_Chrna3_Chrnb2 :=
    log2(Chrna3_count / Chrnb2_count)
]

# Classify cells according to which subunit has the higher count

double_positive_cells[
  ,
  dominance := fifelse(
    Chrna3_count > Chrnb2_count,
    "Chrna3-dominant",
    fifelse(
      Chrnb2_count > Chrna3_count,
      "Chrnb2-dominant",
      "Equal"
    )
  )
]

# Overall dominance categories

table(
  double_positive_cells$dominance
)

# Actual Chrna3 / Chrnb2 count combinations

pair_table <- double_positive_cells[
  ,
  .N,
  by = .(
    Chrna3_count,
    Chrnb2_count
  )
]

setorder(
  pair_table,
  -N
)

pair_table

# Dominance categories by inhibitory subclass

table(
  double_positive_cells$subclass_label,
  double_positive_cells$dominance
)

# Continuous log2 ratio summary

summary(
  double_positive_cells$log2_Chrna3_Chrnb2
)

# Save the cell-level double-positive data

fwrite(
  double_positive_cells,
  "results/mouse_HIP_GABA_double_positive_cells.csv"
)

fwrite(
  pair_table,
  "results/mouse_HIP_GABA_double_positive_count_pairs.csv"
)
# ============================================================
# 24. Create revised manuscript Figure 2:
#     Chrna3 / Chrnb2 transcript co-detection in adult mouse
#     hippocampal GABAergic populations
# ============================================================


# -----------------------------
# Panel A data
# Overall expression categories
# -----------------------------

panel_a_data <- mouse[
  ,
  .N,
  by = expression_group
]

panel_a_data[
  ,
  pct := 100 * N / sum(N)
]

panel_a_data[
  ,
  expression_group := factor(
    expression_group,
    levels = c(
      "Neither",
      "Chrnb2 only",
      "Chrna3 only",
      "Double-positive"
    )
  )
]


# # -----------------------------
# Panel A plot
# -----------------------------

pA <- ggplot(
  panel_a_data,
  aes(
    x = expression_group,
    y = pct
  )
) +
  geom_col(
    fill = "grey70",
    color = "black",
    linewidth = 0.6,
    width = 0.72
  ) +
  geom_text(
    aes(
      y = pct + 2,
      label = sprintf(
        "%.2f%%\n(n=%d)",
        pct,
        N
      )
    ),
    vjust = 0,
    size = 3.3
  ) +
  scale_x_discrete(
    labels = c(
      "Neither",
      "<i>Chrnb2</i><br>only",
      "<i>Chrna3</i><br>only",
      "Double-<br>positive"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 80),
    breaks = c(
      0, 20, 40, 60, 80
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    x = NULL,
    y = "Cells (%)",
    title = "Expression categories"
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    axis.text.x = ggtext::element_markdown(
      angle = 0,
      hjust = 0.5
    ),
    plot.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )
# -----------------------------
# Panel B data
# Double-positive frequency by
# inhibitory subclass
# -----------------------------

panel_b_data <- copy(
  mouse_subclass_results
)

panel_b_data[
  ,
  subclass_label := factor(
    subclass_label,
    levels = panel_b_data[
      order(pct_double_positive),
      subclass_label
    ]
  )
]


# -----------------------------
# Panel B plot
# -----------------------------

pB <- ggplot(
  panel_b_data,
  aes(
    x = subclass_label,
    y = pct_double_positive
  )
) +
  geom_col(
    fill = "grey55",
    color = "black",
    linewidth = 0.6,
    width = 0.70
  ) +
  geom_text(
    aes(
      y = pct_double_positive + 0.12,
      label = sprintf(
        "%.2f%% (%d/%d)",
        pct_double_positive,
        double_positive,
        total_cells
      )
    ),
    hjust = 0,
    size = 3.1
  ) +
  coord_flip(
    clip = "off"
  ) +
  scale_y_continuous(
    limits = c(0, 6.2),
    breaks = c(
      0, 2, 4, 6
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    x = NULL,
    y = "Double-positive cells (%)",
    title = "GABAergic transcriptomic classes"
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      8, 20, 8, 8
    )
  )


# -----------------------------
# Panel C data
# Pvalb-Vipr2 enrichment
# -----------------------------

pvalb_focus <- mouse[
  ,
  .(
    total_cells = .N,
    double_positive =
      sum(
        expression_group ==
          "Double-positive"
      )
  ),
  by = .(
    population = ifelse(
      subclass_label == "Pvalb" &
        cluster_label %in%
        c(
          "122_Pvalb Vipr2",
          "123_Pvalb Vipr2"
        ),
      "Pvalb-Vipr2",
      "Other HIP GABAergic"
    )
  )
]

pvalb_focus[
  ,
  pct_double_positive :=
    100 * double_positive /
    total_cells
]

pvalb_focus[
  ,
  population := factor(
    population,
    levels = c(
      "Other HIP GABAergic",
      "Pvalb-Vipr2"
    )
  )
]


# -----------------------------
# Panel C plot
# -----------------------------

pC <- ggplot(
  pvalb_focus,
  aes(
    x = population,
    y = pct_double_positive,
    fill = population
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.6,
    width = 0.65
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f%%\n(%d/%d)",
        pct_double_positive,
        double_positive,
        total_cells
      )
    ),
    vjust = -0.4,
    size = 3.3
  ) +
  scale_fill_manual(
    values = c(
      "grey80",
      "grey35"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "Other hippocampal\nGABAergic",
      "Pvalb-Vipr2"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 14),
    breaks = c(
      0, 2.5, 5, 7.5, 10, 12.5
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    x = NULL,
    y = "Double-positive cells (%)",
    title = "Pvalb-Vipr2 enrichment",
    subtitle = "OR = 12.87; Fisher P = 3.85e-8"
  ) +
  guides(
    fill = "none"
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    plot.title = element_text(
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 8
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )


# -----------------------------
# Panel D data
# Lamp5-Lhx6 enrichment
# -----------------------------

lamp5_focus <- mouse[
  subclass_label == "Lamp5",
  .(
    total_cells = .N,
    double_positive =
      sum(
        expression_group ==
          "Double-positive"
      )
  ),
  by = .(
    population = ifelse(
      grepl(
        "Lamp5 Lhx6",
        cluster_label
      ),
      "Lamp5-Lhx6",
      "Other Lamp5"
    )
  )
]

lamp5_focus[
  ,
  pct_double_positive :=
    100 * double_positive /
    total_cells
]

lamp5_focus[
  ,
  population := factor(
    population,
    levels = c(
      "Other Lamp5",
      "Lamp5-Lhx6"
    )
  )
]


# -----------------------------
# Panel D plot
# -----------------------------

pD <- ggplot(
  lamp5_focus,
  aes(
    x = population,
    y = pct_double_positive,
    fill = population
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.6,
    width = 0.65
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f%%\n(%d/%d)",
        pct_double_positive,
        double_positive,
        total_cells
      )
    ),
    vjust = -0.4,
    size = 3.3
  ) +
  scale_fill_manual(
    values = c(
      "grey80",
      "grey45"
    )
  ) +
  scale_x_discrete(
    labels = c(
      "Other\nLamp5",
      "Lamp5-Lhx6"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1.35),
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1.00,
      1.25
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    x = NULL,
    y = "Double-positive cells (%)",
    title = "Lamp5-Lhx6 enrichment",
    subtitle = "OR = 10.38; Fisher P = 0.00213"
  ) +
  guides(
    fill = "none"
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    plot.title = element_text(
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 8
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )


# -----------------------------
# Combine Panels A-D
# -----------------------------

mouse_figure_revised <- (
  pA | pB
) / (
  pC | pD
) +
  plot_annotation(
    title = "Chrna3 and Chrnb2 transcript co-detection in adult mouse\nhippocampal GABAergic populations",
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 12
    ),
    plot.title = element_text(
      face = "bold",
      size = 12
    )
  )


# -----------------------------
# Render title with italic
# gene symbols using ggtext
# -----------------------------

mouse_figure_revised <- mouse_figure_revised &
  theme(
    plot.title = ggtext::element_markdown(
      face = "bold",
      size = 12
    )
  )

# Use markdown italics in the title
mouse_figure_revised <- (
  pA | pB
) / (
  pC | pD
) +
  plot_annotation(
    title = "*Chrna3* and *Chrnb2* transcript co-detection in adult mouse<br>hippocampal GABAergic populations",
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 12
    ),
    plot.title = ggtext::element_markdown(
      face = "bold",
      size = 12
    )
  )


# -----------------------------
# Display revised Figure 2
# -----------------------------

mouse_figure_revised


# -----------------------------
# Save revised Figure 2
# -----------------------------

ggsave(
  filename =
    "figures/Allen_mouse_CHRNA3_CHRNB2_Figure2_revised.pdf",
  plot = mouse_figure_revised,
  width = 7.09,
  height = 7.2,
  units = "in"
)

ggsave(
  filename =
    "figures/Allen_mouse_CHRNA3_CHRNB2_Figure2_revised.tiff",
  plot = mouse_figure_revised,
  width = 7.09,
  height = 7.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename =
    "figures/Allen_mouse_CHRNA3_CHRNB2_Figure2_revised.png",
  plot = mouse_figure_revised,
  width = 7.09,
  height = 7.2,
  units = "in",
  dpi = 300
)