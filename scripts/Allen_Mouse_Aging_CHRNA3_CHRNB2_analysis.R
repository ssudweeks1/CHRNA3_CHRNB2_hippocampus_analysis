# ============================================================
# Allen Institute Aging Mouse Brain
# Chrna3 / Chrnb2 co-expression analysis
#
# Young-adult (~2 month) versus aged (~18 month)
# Hippocampal GABAergic interneurons
#
# Project:
# Allen_Mouse_Aging_CHRNA3_CHRNB2
#
# Allen manifest used for downloaded data:
# releases/20260711/manifest.json
#
# Primary analysis:
#   HPF - HIP
#   Allen WMB taxonomy classes:
#       06 CTX-CGE GABA
#       07 CTX-MGE GABA
#   unbiased sampling only
#
# Gene detection:
#   raw UMI count > 0
# ============================================================


# ------------------------------------------------------------
# 1. Project setup
# ------------------------------------------------------------

# Make sure data.table is installed
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop(
    "The R package 'data.table' is required. ",
    "Install it with install.packages('data.table')."
  )
}

library(data.table)


# Confirm project location
cat("\nPROJECT DIRECTORY\n")
cat("-----------------\n")
print(getwd())


# Make sure required output folders exist
dir.create(
  "data/derived",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)



# ------------------------------------------------------------
# 2. Paths to downloaded Allen metadata
# ------------------------------------------------------------

cell_path <- file.path(
  "data",
  "abc_atlas",
  "metadata",
  "Zeng-Aging-Mouse-10Xv3",
  "20250131",
  "cell_metadata.csv"
)

cluster_path <- file.path(
  "data",
  "abc_atlas",
  "metadata",
  "Zeng-Aging-Mouse-10Xv3",
  "20250131",
  "cluster.csv"
)

donor_path <- file.path(
  "data",
  "abc_atlas",
  "metadata",
  "Zeng-Aging-Mouse-10Xv3",
  "20250131",
  "donor.csv"
)

mapping_path <- file.path(
  "data",
  "abc_atlas",
  "metadata",
  "Zeng-Aging-Mouse-WMB-taxonomy",
  "20241130",
  "cluster_mapping.csv"
)


# Confirm required metadata files exist
metadata_paths <- c(
  cell_path,
  cluster_path,
  donor_path,
  mapping_path
)

if (any(!file.exists(metadata_paths))) {
  
  missing_files <- metadata_paths[
    !file.exists(metadata_paths)
  ]
  
  stop(
    paste(
      "Required metadata file(s) not found:",
      paste(missing_files, collapse = "\n")
    )
  )
}


cat("\nMETADATA FILE CHECK\n")
cat("-------------------\n")
print(file.exists(metadata_paths))



# ------------------------------------------------------------
# Inspect metadata column names
# ------------------------------------------------------------

cell_header <- fread(
  cell_path,
  nrows = 0
)

cluster_header <- fread(
  cluster_path,
  nrows = 0
)

donor_header <- fread(
  donor_path,
  nrows = 0
)

mapping_header <- fread(
  mapping_path,
  nrows = 0
)



# ------------------------------------------------------------
# 3. Load metadata needed for hippocampal
#    GABAergic interneuron analysis
# ------------------------------------------------------------

cell_cols <- c(
  "cell_label",
  "cell_barcode",
  "cluster_alias",
  "wmb_cluster_alias",
  "library_label",
  "region_of_interest_label",
  "anatomical_division_label",
  "donor_label",
  "population_sampling",
  "donor_genotype",
  "donor_sex",
  "donor_age",
  "donor_age_category",
  "donor_in_wmb_study",
  "feature_matrix_label",
  "abc_sample_id"
)


# Verify every requested column exists
missing_cell_cols <- setdiff(
  cell_cols,
  names(cell_header)
)

if (length(missing_cell_cols) > 0) {
  
  stop(
    paste(
      "Required cell metadata columns are missing:",
      paste(missing_cell_cols, collapse = ", ")
    )
  )
}


# Load selected cell metadata
cell <- fread(
  cell_path,
  select = cell_cols
)


cat("\nCELL METADATA\n")
cat("-------------\n")
print(dim(cell))


# Load cluster metadata
cluster <- fread(
  cluster_path,
  select = c(
    "cluster_alias",
    "cluster_label",
    "cluster_name",
    "number_of_cells",
    "neurotransmitter_combined_label",
    "neurotransmitter_label"
  )
)


cat("\nCLUSTER METADATA\n")
cat("----------------\n")
print(dim(cluster))


# cluster_alias must be unique in cluster table
if (anyDuplicated(cluster$cluster_alias) != 0) {
  stop(
    "cluster_alias is not unique in cluster metadata."
  )
}


# ------------------------------------------------------------
# Basic metadata QC
# ------------------------------------------------------------

cat("\nAGE CATEGORY COUNTS — ALL CELLS\n")
print(
  table(
    cell$donor_age_category,
    useNA = "ifany"
  )
)


cat("\nDONOR AGES — ALL CELLS\n")
print(
  sort(
    unique(cell$donor_age)
  )
)


cat("\nSEX COUNTS — ALL CELLS\n")
print(
  table(
    cell$donor_sex,
    useNA = "ifany"
  )
)


cat("\nSAMPLING COUNTS — ALL CELLS\n")
print(
  table(
    cell$population_sampling,
    useNA = "ifany"
  )
)


cat("\nGENOTYPE COUNTS — ALL CELLS\n")
print(
  table(
    cell$donor_genotype,
    useNA = "ifany"
  )
)


cat("\nHPF / HIP REGION LABELS\n")
print(
  sort(
    unique(
      cell$region_of_interest_label[
        grepl(
          "HIP|HPF",
          cell$region_of_interest_label,
          ignore.case = TRUE
        )
      ]
    )
  )
)


cat("\nALLEN NEUROTRANSMITTER LABELS\n")
print(
  table(
    cluster$neurotransmitter_label,
    useNA = "ifany"
  )
)



# ------------------------------------------------------------
# 4. Add cluster-level neurotransmitter annotation
#    to each cell
# ------------------------------------------------------------

n_before <- nrow(cell)

cell <- merge(
  cell,
  cluster[
    ,
    .(
      cluster_alias,
      cluster_label,
      cluster_name,
      neurotransmitter_combined_label,
      neurotransmitter_label
    )
  ],
  by = "cluster_alias",
  all.x = TRUE,
  sort = FALSE
)


# Verify no cells were lost or duplicated
if (nrow(cell) != n_before) {
  stop(
    "Cell count changed after cluster metadata join."
  )
}


# Verify every cell received annotation
if (any(is.na(cell$neurotransmitter_label))) {
  stop(
    "Some cells lack neurotransmitter annotations."
  )
}


cat("\nCLUSTER ANNOTATION JOIN SUCCESSFUL\n")
cat("Cells after join:", nrow(cell), "\n")



# ------------------------------------------------------------
# 5. Load Allen WMB taxonomy mapping
# ------------------------------------------------------------

mapping <- fread(
  mapping_path
)


cat("\nWMB TAXONOMY MAPPING\n")
cat("--------------------\n")
print(dim(mapping))


cat("\nTAXONOMY LEVELS\n")
print(
  sort(
    unique(
      mapping$cluster_annotation_term_set_name
    )
  )
)



# ------------------------------------------------------------
# 6. Convert WMB taxonomy mapping to wide format
#    and attach class/subclass/supertype to cells
# ------------------------------------------------------------

mapping_keep <- mapping[
  cluster_annotation_term_set_name %in%
    c(
      "class",
      "subclass",
      "supertype"
    ),
  .(
    cluster_alias,
    cluster_annotation_term_set_name,
    cluster_annotation_term_name
  )
]


mapping_wide <- dcast(
  mapping_keep,
  cluster_alias ~ cluster_annotation_term_set_name,
  value.var = "cluster_annotation_term_name"
)


# Verify one row per cluster
if (anyDuplicated(mapping_wide$cluster_alias) != 0) {
  stop(
    "cluster_alias is not unique after taxonomy pivot."
  )
}


cat("\nWIDE TAXONOMY TABLE\n")
print(dim(mapping_wide))


# Add taxonomy to cells
n_before_taxonomy <- nrow(cell)

cell <- merge(
  cell,
  mapping_wide,
  by = "cluster_alias",
  all.x = TRUE,
  sort = FALSE
)


if (nrow(cell) != n_before_taxonomy) {
  stop(
    "Cell count changed after taxonomy join."
  )
}


if (
  any(is.na(cell$class)) ||
  any(is.na(cell$subclass)) ||
  any(is.na(cell$supertype))
) {
  stop(
    "Some cells lack class/subclass/supertype taxonomy."
  )
}



# ------------------------------------------------------------
# Create complete HPF-HIP population
# ------------------------------------------------------------

hip_all <- cell[
  region_of_interest_label == "HPF - HIP"
]


cat("\nHPF-HIP CELLS\n")
cat("-------------\n")
print(dim(hip_all))


cat("\nHPF-HIP CELLS BY WMB CLASS\n")
print(
  hip_all[
    ,
    .(
      n_cells = .N,
      n_clusters = uniqueN(cluster_alias),
      n_donors = uniqueN(donor_label)
    ),
    by = class
  ][
    order(-n_cells)
  ]
)



# ------------------------------------------------------------
# 7. Define hippocampal GABAergic interneuron population
# ------------------------------------------------------------

# Primary biological definition:
#
# HPF - HIP cells belonging to the Allen WMB
# cortical-type GABAergic classes:
#
#   06 CTX-CGE GABA
#   07 CTX-MGE GABA
#
# This retains the Sncg Gaba_4 cluster that Allen
# labels Glut-GABA at the neurotransmitter level,
# because it maps taxonomically to CTX-CGE GABA.

hip_gaba <- hip_all[
  class %in% c(
    "06 CTX-CGE GABA",
    "07 CTX-MGE GABA"
  )
]


cat("\nFINAL HIPPOCAMPAL GABAERGIC POPULATION\n")
cat("--------------------------------------\n")
print(dim(hip_gaba))


cat("\nNEUROTRANSMITTER COMPOSITION\n")
print(
  hip_gaba[
    ,
    .(
      n_cells = .N,
      n_clusters = uniqueN(cluster_alias),
      n_donors = uniqueN(donor_label)
    ),
    by = neurotransmitter_label
  ][
    order(-n_cells)
  ]
)


cat("\nSUBCLASS COMPOSITION\n")
print(
  hip_gaba[
    ,
    .(
      n_cells = .N,
      n_clusters = uniqueN(cluster_alias),
      n_donors = uniqueN(donor_label)
    ),
    by = subclass
  ][
    order(-n_cells)
  ]
)



# ------------------------------------------------------------
# 8. Characterize donor structure of complete
#    hippocampal GABAergic population
# ------------------------------------------------------------

cat("\nAGE-GROUP COMPOSITION — ALL SAMPLING\n")
print(
  hip_gaba[
    ,
    .(
      n_cells = .N,
      n_donors = uniqueN(donor_label)
    ),
    by = donor_age_category
  ]
)


cat("\nEXACT DONOR AGES\n")
print(
  hip_gaba[
    ,
    .(
      n_cells = .N,
      n_donors = uniqueN(donor_label)
    ),
    by = .(
      donor_age_category,
      donor_age
    )
  ][
    order(
      donor_age_category,
      donor_age
    )
  ]
)


cat("\nSEX BY AGE\n")
print(
  hip_gaba[
    ,
    .(
      n_cells = .N,
      n_donors = uniqueN(donor_label)
    ),
    by = .(
      donor_age_category,
      donor_sex
    )
  ]
)


cat("\nSAMPLING STRATEGY BY AGE\n")
print(
  hip_gaba[
    ,
    .(
      n_cells = .N,
      n_donors = uniqueN(donor_label)
    ),
    by = .(
      donor_age_category,
      population_sampling
    )
  ]
)


cat("\nGENOTYPE BY AGE\n")
print(
  hip_gaba[
    ,
    .(
      n_cells = .N,
      n_donors = uniqueN(donor_label)
    ),
    by = .(
      donor_age_category,
      donor_genotype
    )
  ]
)


cat("\nWMB STUDY OVERLAP BY AGE\n")
print(
  hip_gaba[
    ,
    .(
      n_cells = .N,
      n_donors = uniqueN(donor_label)
    ),
    by = .(
      donor_age_category,
      donor_in_wmb_study
    )
  ]
)



# ------------------------------------------------------------
# 9. Define primary age-comparison cohort
# ------------------------------------------------------------

# Primary analysis uses only unbiased sampling.
#
# Reason:
# All adult HPF-HIP GABAergic samples are unbiased,
# while Snap25+ neuron-enriched sampling occurs only
# among aged mice.

hip_gaba_unbiased <- hip_gaba[
  population_sampling == "unbiased"
]


cat("\nPRIMARY UNBIASED COHORT\n")
cat("-----------------------\n")
print(dim(hip_gaba_unbiased))


primary_age_summary <- hip_gaba_unbiased[
  ,
  .(
    n_cells = .N,
    n_donors = uniqueN(donor_label)
  ),
  by = donor_age_category
]


print(primary_age_summary)


# Save cohort summary
fwrite(
  primary_age_summary,
  "results/tables/primary_unbiased_age_summary.csv"
)



# ------------------------------------------------------------
# 10. Record primary donor cohort
# ------------------------------------------------------------

primary_donors <- unique(
  hip_gaba_unbiased[
    ,
    .(
      donor_label,
      donor_age_category,
      donor_age,
      donor_sex,
      donor_genotype
    )
  ]
)


setorder(
  primary_donors,
  donor_age_category,
  donor_label
)


cat("\nPRIMARY BIOLOGICAL REPLICATES\n")
cat("-----------------------------\n")
print(primary_donors)


cat(
  "\nNumber of primary donors:",
  nrow(primary_donors),
  "\n"
)


fwrite(
  primary_donors,
  "results/tables/primary_unbiased_donor_cohort.csv"
)



# ------------------------------------------------------------
# 11. Donor-by-subclass coverage in primary cohort
# ------------------------------------------------------------

# Count observed cells for each donor x subclass
donor_subclass_counts <- hip_gaba_unbiased[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    donor_label,
    donor_age_category,
    subclass
  )
]


# Create every possible donor x subclass combination
# so that zero-cell combinations are explicitly retained.
donor_subclass_grid <- CJ(
  donor_label = unique(
    hip_gaba_unbiased$donor_label
  ),
  subclass = unique(
    hip_gaba_unbiased$subclass
  ),
  unique = TRUE
)


# Donor-to-age lookup
donor_age_lookup <- unique(
  hip_gaba_unbiased[
    ,
    .(
      donor_label,
      donor_age_category
    )
  ]
)


donor_subclass_grid <- merge(
  donor_subclass_grid,
  donor_age_lookup,
  by = "donor_label",
  all.x = TRUE,
  sort = FALSE
)


# Add observed cell counts
donor_subclass_complete <- merge(
  donor_subclass_grid,
  donor_subclass_counts,
  by = c(
    "donor_label",
    "donor_age_category",
    "subclass"
  ),
  all.x = TRUE,
  sort = FALSE
)


# Missing combinations represent zero recovered cells
donor_subclass_complete[
  is.na(n_cells),
  n_cells := 0L
]


cat("\nDONOR x SUBCLASS GRID\n")
print(dim(donor_subclass_complete))


# Summarize coverage by age and subclass
subclass_qc <- donor_subclass_complete[
  ,
  .(
    n_donors_total = as.integer(.N),
    n_donors_with_cells =
      as.integer(sum(n_cells > 0)),
    min_cells_per_donor =
      as.integer(min(n_cells)),
    median_cells_per_donor =
      as.numeric(median(n_cells)),
    max_cells_per_donor =
      as.integer(max(n_cells)),
    total_cells =
      as.integer(sum(n_cells))
  ),
  by = .(
    donor_age_category,
    subclass
  )
][
  order(
    subclass,
    donor_age_category
  )
]


cat("\nPRIMARY SUBCLASS QC\n")
cat("-------------------\n")
print(subclass_qc)


fwrite(
  subclass_qc,
  "results/tables/primary_unbiased_subclass_QC.csv"
)



# ------------------------------------------------------------
# 12. Load Chrna3 / Chrnb2 raw-count expression data
#     and join to primary unbiased cohort
# ------------------------------------------------------------

expression_path <- file.path(
  "data",
  "derived",
  "Zeng_Aging_Mouse_Chrna3_Chrnb2_raw_counts.csv"
)


if (!file.exists(expression_path)) {
  
  stop(
    paste(
      "Two-gene raw-count file was not found:",
      expression_path
    )
  )
}


# Load two-gene expression table
expression_raw <- fread(
  expression_path
)


cat("\nTWO-GENE RAW EXPRESSION TABLE\n")
cat("-----------------------------\n")
print(dim(expression_raw))
print(head(expression_raw))


# Verify expected columns
required_expression_cols <- c(
  "cell_label",
  "Chrna3_raw",
  "Chrnb2_raw"
)

missing_expression_cols <- setdiff(
  required_expression_cols,
  names(expression_raw)
)

if (length(missing_expression_cols) > 0) {
  
  stop(
    paste(
      "Required expression columns are missing:",
      paste(
        missing_expression_cols,
        collapse = ", "
      )
    )
  )
}


# Cell labels must be unique
if (
  anyDuplicated(
    expression_raw$cell_label
  ) != 0
) {
  
  stop(
    "Duplicate cell_label values found in expression table."
  )
}


# Raw counts must be non-negative
if (
  any(expression_raw$Chrna3_raw < 0) ||
  any(expression_raw$Chrnb2_raw < 0)
) {
  
  stop(
    "Negative values found in raw expression counts."
  )
}


# ------------------------------------------------------------
# Join expression to primary cohort
# ------------------------------------------------------------

n_primary_before <- nrow(
  hip_gaba_unbiased
)


primary_expr <- merge(
  hip_gaba_unbiased,
  expression_raw,
  by = "cell_label",
  all.x = TRUE,
  sort = FALSE
)


if (
  nrow(primary_expr) !=
  n_primary_before
) {
  
  stop(
    "Primary cohort cell count changed during expression join."
  )
}


if (
  any(is.na(primary_expr$Chrna3_raw)) ||
  any(is.na(primary_expr$Chrnb2_raw))
) {
  
  stop(
    "Some primary cells did not receive Chrna3/Chrnb2 values."
  )
}


cat("\nPRIMARY COHORT + EXPRESSION\n")
cat("---------------------------\n")
print(dim(primary_expr))



# ------------------------------------------------------------
# Define expression / detection categories
#
# Gene detected = raw UMI count > 0
# ------------------------------------------------------------

primary_expr[
  ,
  Chrna3_detected :=
    Chrna3_raw > 0
]


primary_expr[
  ,
  Chrnb2_detected :=
    Chrnb2_raw > 0
]


primary_expr[
  ,
  expression_category :=
    fcase(
      
      Chrna3_detected &
        Chrnb2_detected,
      "Chrna3+Chrnb2+",
      
      Chrna3_detected &
        !Chrnb2_detected,
      "Chrna3 only",
      
      !Chrna3_detected &
        Chrnb2_detected,
      "Chrnb2 only",
      
      default = "Neither"
    )
]


primary_expr[
  ,
  expression_category :=
    factor(
      expression_category,
      levels = c(
        "Neither",
        "Chrna3 only",
        "Chrnb2 only",
        "Chrna3+Chrnb2+"
      )
    )
]


# ------------------------------------------------------------
# Basic expression QC
# ------------------------------------------------------------

cat("\nPRIMARY EXPRESSION CATEGORIES\n")
cat("-----------------------------\n")

print(
  table(
    primary_expr$expression_category,
    useNA = "ifany"
  )
)


primary_expression_age_counts <- primary_expr[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    donor_age_category,
    expression_category
  )
][
  order(
    donor_age_category,
    expression_category
  )
]


cat("\nEXPRESSION CATEGORIES BY AGE\n")
print(primary_expression_age_counts)


# Save age-level raw category counts
fwrite(
  primary_expression_age_counts,
  "results/tables/primary_expression_categories_by_age.csv"
)


# Save complete primary two-gene cell-level dataset
fwrite(
  primary_expr,
  file.path(
    "data",
    "derived",
    "primary_unbiased_hip_gaba_Chrna3_Chrnb2.csv"
  )
)



# ------------------------------------------------------------
# Final checkpoint
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT COMPLETED THROUGH SECTION 12\n")
cat("============================================================\n")

cat(
  "Primary cells:",
  nrow(primary_expr),
  "\n"
)

cat(
  "Primary donors:",
  uniqueN(primary_expr$donor_label),
  "\n"
)

cat(
  "Adult donors:",
  uniqueN(
    primary_expr[
      donor_age_category == "adult",
      donor_label
    ]
  ),
  "\n"
)

cat(
  "Aged donors:",
  uniqueN(
    primary_expr[
      donor_age_category == "aged",
      donor_label
    ]
  ),
  "\n"
)

cat("\nReady for donor-aware Chrna3/Chrnb2 analysis.\n")
# ------------------------------------------------------------
# 13. Donor-level Chrna3 / Chrnb2 expression summary
# ------------------------------------------------------------

# The mouse/donor is the biological replicate.
#
# For each donor calculate:
#   - total analyzed cells
#   - Chrna3-positive cells
#   - Chrnb2-positive cells
#   - double-positive cells
#   - four expression categories
#   - corresponding percentages


donor_expression_summary <- primary_expr[
  ,
  .(
    n_cells = .N,
    
    n_chrna3_pos =
      sum(Chrna3_detected),
    
    n_chrnb2_pos =
      sum(Chrnb2_detected),
    
    n_double_pos =
      sum(
        Chrna3_detected &
          Chrnb2_detected
      ),
    
    n_chrna3_only =
      sum(
        Chrna3_detected &
          !Chrnb2_detected
      ),
    
    n_chrnb2_only =
      sum(
        !Chrna3_detected &
          Chrnb2_detected
      ),
    
    n_neither =
      sum(
        !Chrna3_detected &
          !Chrnb2_detected
      ),
    
    pct_chrna3 =
      100 * mean(
        Chrna3_detected
      ),
    
    pct_chrnb2 =
      100 * mean(
        Chrnb2_detected
      ),
    
    pct_double =
      100 * mean(
        Chrna3_detected &
          Chrnb2_detected
      )
  ),
  
  by = .(
    donor_label,
    donor_age_category,
    donor_age,
    donor_sex,
    donor_genotype
  )
]


# ------------------------------------------------------------
# Conditional co-expression measures
#
# These answer:
#
# Among Chrna3-positive cells,
# what percentage also express Chrnb2?
#
# Among Chrnb2-positive cells,
# what percentage also express Chrna3?
# ------------------------------------------------------------

donor_expression_summary[
  ,
  pct_chrnb2_among_chrna3 :=
    fifelse(
      n_chrna3_pos > 0,
      100 * n_double_pos /
        n_chrna3_pos,
      NA_real_
    )
]


donor_expression_summary[
  ,
  pct_chrna3_among_chrnb2 :=
    fifelse(
      n_chrnb2_pos > 0,
      100 * n_double_pos /
        n_chrnb2_pos,
      NA_real_
    )
]


# ------------------------------------------------------------
# Sort donors
# ------------------------------------------------------------

setorder(
  donor_expression_summary,
  donor_age_category,
  donor_label
)


cat("\nDONOR-LEVEL EXPRESSION SUMMARY\n")
cat("------------------------------\n")

print(
  donor_expression_summary
)


# ------------------------------------------------------------
# Verify exactly nine biological replicates
# ------------------------------------------------------------

cat(
  "\nNumber of donor-level observations:",
  nrow(donor_expression_summary),
  "\n"
)


cat(
  "Adult donors:",
  donor_expression_summary[
    donor_age_category == "adult",
    .N
  ],
  "\n"
)


cat(
  "Aged donors:",
  donor_expression_summary[
    donor_age_category == "aged",
    .N
  ],
  "\n"
)


# ------------------------------------------------------------
# Cross-check donor totals against cell-level results
# ------------------------------------------------------------

cat("\nDONOR-TOTAL CROSS-CHECK\n")
cat("-----------------------\n")

cat(
  "Total cells:",
  sum(
    donor_expression_summary$n_cells
  ),
  "\n"
)

cat(
  "Total Chrna3-positive:",
  sum(
    donor_expression_summary$n_chrna3_pos
  ),
  "\n"
)

cat(
  "Total Chrnb2-positive:",
  sum(
    donor_expression_summary$n_chrnb2_pos
  ),
  "\n"
)

cat(
  "Total double-positive:",
  sum(
    donor_expression_summary$n_double_pos
  ),
  "\n"
)


# ------------------------------------------------------------
# Summary of donor-level percentages by age group
#
# These summaries are descriptive only.
# Formal donor-aware statistics will be performed next.
# ------------------------------------------------------------

age_donor_summary <- donor_expression_summary[
  ,
  .(
    n_donors = .N,
    
    mean_pct_chrna3 =
      mean(pct_chrna3),
    
    median_pct_chrna3 =
      median(pct_chrna3),
    
    mean_pct_chrnb2 =
      mean(pct_chrnb2),
    
    median_pct_chrnb2 =
      median(pct_chrnb2),
    
    mean_pct_double =
      mean(pct_double),
    
    median_pct_double =
      median(pct_double),
    
    mean_pct_chrnb2_among_chrna3 =
      mean(
        pct_chrnb2_among_chrna3,
        na.rm = TRUE
      ),
    
    median_pct_chrnb2_among_chrna3 =
      median(
        pct_chrnb2_among_chrna3,
        na.rm = TRUE
      )
  ),
  by = donor_age_category
]


cat("\nDESCRIPTIVE DONOR-LEVEL AGE SUMMARY\n")
cat("-----------------------------------\n")

print(
  age_donor_summary
)


# ------------------------------------------------------------
# Save donor-level results
# ------------------------------------------------------------

fwrite(
  donor_expression_summary,
  "results/tables/primary_donor_expression_summary.csv"
)


fwrite(
  age_donor_summary,
  "results/tables/primary_donor_expression_age_summary.csv"
)


cat(
  "\nSection 13 donor-level summaries saved successfully.\n"
)
# ------------------------------------------------------------
# 14. Donor-aware adult vs aged statistical analysis
#
# Exact enumeration permutation tests
# Biological replicate = donor / mouse
# ------------------------------------------------------------


# ------------------------------------------------------------
# Function for exact two-group permutation test
#
# Test statistic:
#   mean(aged donors) - mean(adult donors)
#
# All possible group assignments are enumerated exactly.
# The number of permutations depends on the number of
# donors with a defined value for the outcome being tested.
# ------------------------------------------------------------

exact_donor_permutation_test <- function(
    data,
    outcome,
    group_col = "donor_age_category",
    adult_label = "adult",
    aged_label = "aged"
) {
  
  # Keep only required columns
  d <- data[
    !is.na(get(outcome)) &
      get(group_col) %in% c(
        adult_label,
        aged_label
      ),
    .(
      group = get(group_col),
      value = get(outcome)
    )
  ]
  
  
  # Number of donors in each group
  n_total <- nrow(d)
  
  n_adult <- sum(
    d$group == adult_label
  )
  
  n_aged <- sum(
    d$group == aged_label
  )
  
  
  # Observed donor-level values
  adult_values <- d[
    group == adult_label,
    value
  ]
  
  aged_values <- d[
    group == aged_label,
    value
  ]
  
  
  # Observed difference:
  # aged minus adult
  observed_difference <- (
    mean(aged_values) -
      mean(adult_values)
  )
  
  
  # ----------------------------------------------------------
  # Enumerate every possible assignment of n_adult donors
  # ----------------------------------------------------------
  
  adult_combinations <- combn(
    n_total,
    n_adult,
    simplify = FALSE
  )
  
  
  permutation_differences <- vapply(
    adult_combinations,
    FUN.VALUE = numeric(1),
    
    FUN = function(adult_index) {
      
      aged_index <- setdiff(
        seq_len(n_total),
        adult_index
      )
      
      mean(
        d$value[aged_index]
      ) -
        mean(
          d$value[adult_index]
        )
    }
  )
  
  
  # ----------------------------------------------------------
  # Exact two-sided permutation P value
  # ----------------------------------------------------------
  
  tolerance <- 1e-12
  
  p_value <- mean(
    abs(permutation_differences) >=
      (
        abs(observed_difference) -
          tolerance
      )
  )
  
  
  # ----------------------------------------------------------
  # Return results
  # ----------------------------------------------------------
  
  data.table(
    outcome = outcome,
    
    n_adult = n_adult,
    n_aged = n_aged,
    
    adult_mean = mean(
      adult_values
    ),
    
    adult_median = median(
      adult_values
    ),
    
    aged_mean = mean(
      aged_values
    ),
    
    aged_median = median(
      aged_values
    ),
    
    difference_aged_minus_adult =
      observed_difference,
    
    exact_permutation_p =
      p_value,
    
    n_permutations =
      length(
        permutation_differences
      )
  )
}



# ------------------------------------------------------------
# Outcomes to test
# ------------------------------------------------------------

outcomes_to_test <- c(
  
  # Individual gene detection
  "pct_chrna3",
  "pct_chrnb2",
  
  # Co-expression among all analyzed cells
  "pct_double",
  
  # Conditional co-expression
  "pct_chrnb2_among_chrna3",
  "pct_chrna3_among_chrnb2"
)



# ------------------------------------------------------------
# Run exact permutation test for each outcome
# ------------------------------------------------------------

donor_age_tests <- rbindlist(
  lapply(
    outcomes_to_test,
    function(x) {
      
      exact_donor_permutation_test(
        donor_expression_summary,
        x
      )
      
    }
  ),
  fill = TRUE
)



# ------------------------------------------------------------
# Add human-readable outcome labels
# ------------------------------------------------------------

donor_age_tests[
  ,
  outcome_label := fcase(
    
    outcome == "pct_chrna3",
    "% Chrna3-positive",
    
    outcome == "pct_chrnb2",
    "% Chrnb2-positive",
    
    outcome == "pct_double",
    "% Chrna3+Chrnb2+",
    
    outcome ==
      "pct_chrnb2_among_chrna3",
    "% Chrnb2+ among Chrna3+",
    
    outcome ==
      "pct_chrna3_among_chrnb2",
    "% Chrna3+ among Chrnb2+",
    
    default = outcome
  )
]



# ------------------------------------------------------------
# Multiple-comparison correction
#
# Holm adjustment is reported across the five
# donor-level outcomes.
# ------------------------------------------------------------

donor_age_tests[
  ,
  holm_adjusted_p :=
    p.adjust(
      exact_permutation_p,
      method = "holm"
    )
]


# Put columns in convenient order
setcolorder(
  donor_age_tests,
  c(
    "outcome",
    "outcome_label",
    "n_adult",
    "n_aged",
    "adult_mean",
    "adult_median",
    "aged_mean",
    "aged_median",
    "difference_aged_minus_adult",
    "exact_permutation_p",
    "holm_adjusted_p",
    "n_permutations"
  )
)



# ------------------------------------------------------------
# Display results
# ------------------------------------------------------------

cat("\n")
cat("============================================================\n")
cat("DONOR-AWARE EXACT AGE COMPARISONS\n")
cat("============================================================\n")

print(
  donor_age_tests
)



# ------------------------------------------------------------
# Confirm complete enumeration for each outcome
# ------------------------------------------------------------

cat(
  "\nPERMUTATION ENUMERATION BY OUTCOME\n"
)

print(
  donor_age_tests[
    ,
    .(
      outcome_label,
      n_adult,
      n_aged,
      n_permutations
    )
  ]
)

# ------------------------------------------------------------
# Save statistical results
# ------------------------------------------------------------

fwrite(
  donor_age_tests,
  "results/tables/primary_donor_exact_permutation_tests.csv"
)


cat(
  "\nSection 14 donor-aware statistical tests ",
  "saved successfully.\n"
)
# ------------------------------------------------------------
# 15. Sensitivity analyses
#
# Sensitivity 1:
#   All sampling methods, all 11 donors
#
# Sensitivity 2:
#   Unbiased sampling only, restricted to genotypes
#   represented in both adult and aged groups
#
# The primary analysis remains the unbiased-sampling
# 4-adult vs 5-aged donor comparison from Section 14.
# ------------------------------------------------------------


# ============================================================
# Helper function:
# Create donor-level Chrna3 / Chrnb2 summary for any
# cell-level expression dataset
# ============================================================

make_donor_expression_summary <- function(expr_data) {
  
  donor_summary <- expr_data[
    ,
    .(
      n_cells = .N,
      
      n_chrna3_pos =
        sum(Chrna3_detected),
      
      n_chrnb2_pos =
        sum(Chrnb2_detected),
      
      n_double_pos =
        sum(
          Chrna3_detected &
            Chrnb2_detected
        ),
      
      n_chrna3_only =
        sum(
          Chrna3_detected &
            !Chrnb2_detected
        ),
      
      n_chrnb2_only =
        sum(
          !Chrna3_detected &
            Chrnb2_detected
        ),
      
      n_neither =
        sum(
          !Chrna3_detected &
            !Chrnb2_detected
        ),
      
      pct_chrna3 =
        100 * mean(
          Chrna3_detected
        ),
      
      pct_chrnb2 =
        100 * mean(
          Chrnb2_detected
        ),
      
      pct_double =
        100 * mean(
          Chrna3_detected &
            Chrnb2_detected
        )
    ),
    
    by = .(
      donor_label,
      donor_age_category,
      donor_age,
      donor_sex,
      donor_genotype
    )
  ]
  
  
  # Conditional co-expression:
  # undefined if a donor has zero cells positive for
  # the conditioning gene.
  
  donor_summary[
    ,
    pct_chrnb2_among_chrna3 :=
      fifelse(
        n_chrna3_pos > 0,
        100 * n_double_pos /
          n_chrna3_pos,
        NA_real_
      )
  ]
  
  
  donor_summary[
    ,
    pct_chrna3_among_chrnb2 :=
      fifelse(
        n_chrnb2_pos > 0,
        100 * n_double_pos /
          n_chrnb2_pos,
        NA_real_
      )
  ]
  
  
  setorder(
    donor_summary,
    donor_age_category,
    donor_label
  )
  
  
  return(
    donor_summary
  )
}



# ============================================================
# Helper function:
# Run the five donor-aware exact permutation tests
# ============================================================

run_sensitivity_tests <- function(
    donor_summary,
    analysis_name
) {
  
  results <- rbindlist(
    lapply(
      outcomes_to_test,
      
      function(x) {
        
        exact_donor_permutation_test(
          donor_summary,
          x
        )
        
      }
    ),
    fill = TRUE
  )
  
  
  # Add readable outcome labels
  
  results[
    ,
    outcome_label := fcase(
      
      outcome == "pct_chrna3",
      "% Chrna3-positive",
      
      outcome == "pct_chrnb2",
      "% Chrnb2-positive",
      
      outcome == "pct_double",
      "% Chrna3+Chrnb2+",
      
      outcome ==
        "pct_chrnb2_among_chrna3",
      "% Chrnb2+ among Chrna3+",
      
      outcome ==
        "pct_chrna3_among_chrnb2",
      "% Chrna3+ among Chrnb2+",
      
      default = outcome
    )
  ]
  
  
  # Holm correction within this sensitivity analysis
  
  results[
    ,
    holm_adjusted_p :=
      p.adjust(
        exact_permutation_p,
        method = "holm"
      )
  ]
  
  
  # Add analysis identifier
  
  results[
    ,
    analysis := analysis_name
  ]
  
  
  setcolorder(
    results,
    c(
      "analysis",
      "outcome",
      "outcome_label",
      "n_adult",
      "n_aged",
      "adult_mean",
      "adult_median",
      "aged_mean",
      "aged_median",
      "difference_aged_minus_adult",
      "exact_permutation_p",
      "holm_adjusted_p",
      "n_permutations"
    )
  )
  
  
  return(
    results
  )
}



# ============================================================
# 15A. Sensitivity analysis 1
#     All sampling methods / all donors
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SENSITIVITY 1: ALL SAMPLING METHODS\n")
cat("============================================================\n")


# Join the raw two-gene counts to all taxonomy-defined
# hippocampal GABAergic cells.

n_all_before <- nrow(
  hip_gaba
)


all_expr <- merge(
  hip_gaba,
  expression_raw,
  by = "cell_label",
  all.x = TRUE,
  sort = FALSE
)


# Verify no cells were lost or duplicated

if (
  nrow(all_expr) !=
  n_all_before
) {
  
  stop(
    "All-sampling cohort changed size during expression join."
  )
}


if (
  any(is.na(all_expr$Chrna3_raw)) ||
  any(is.na(all_expr$Chrnb2_raw))
) {
  
  stop(
    "Missing expression values in all-sampling cohort."
  )
}


# Define detection from raw UMI counts

all_expr[
  ,
  Chrna3_detected :=
    Chrna3_raw > 0
]


all_expr[
  ,
  Chrnb2_detected :=
    Chrnb2_raw > 0
]


cat(
  "Cells:",
  nrow(all_expr),
  "\n"
)

cat(
  "Donors:",
  uniqueN(all_expr$donor_label),
  "\n"
)


cat("\nALL-SAMPLING AGE STRUCTURE\n")

print(
  all_expr[
    ,
    .(
      n_cells = .N,
      n_donors = uniqueN(donor_label)
    ),
    by = donor_age_category
  ]
)


# Create one row per donor

all_sampling_donor_summary <-
  make_donor_expression_summary(
    all_expr
  )


cat("\nALL-SAMPLING DONOR EXPRESSION SUMMARY\n")

print(
  all_sampling_donor_summary
)


# Run exact donor-level tests

all_sampling_tests <-
  run_sensitivity_tests(
    all_sampling_donor_summary,
    "All sampling methods"
  )


cat("\nALL-SAMPLING EXACT AGE TESTS\n")

print(
  all_sampling_tests
)


cat("\nALL-SAMPLING PERMUTATION ENUMERATION\n")

print(
  all_sampling_tests[
    ,
    .(
      outcome_label,
      n_adult,
      n_aged,
      n_permutations
    )
  ]
)


# Save sensitivity-1 outputs

fwrite(
  all_sampling_donor_summary,
  "results/tables/sensitivity_all_sampling_donor_summary.csv"
)


fwrite(
  all_sampling_tests,
  "results/tables/sensitivity_all_sampling_exact_tests.csv"
)



# ============================================================
# 15B. Sensitivity analysis 2
#     Genotype-matched unbiased cohort
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SENSITIVITY 2: GENOTYPE-MATCHED UNBIASED COHORT\n")
cat("============================================================\n")


# ------------------------------------------------------------
# Determine which genotypes occur in BOTH age groups
# rather than hard-coding the genotype names.
# ------------------------------------------------------------

adult_genotypes <- unique(
  hip_gaba_unbiased[
    donor_age_category == "adult",
    donor_genotype
  ]
)


aged_genotypes <- unique(
  hip_gaba_unbiased[
    donor_age_category == "aged",
    donor_genotype
  ]
)


common_genotypes <- intersect(
  adult_genotypes,
  aged_genotypes
)


cat("\nGENOTYPES REPRESENTED IN BOTH AGE GROUPS\n")

print(
  common_genotypes
)


# Restrict the unbiased cohort to common genotypes

hip_gaba_genotype_matched <-
  hip_gaba_unbiased[
    donor_genotype %in%
      common_genotypes
  ]


cat("\nGENOTYPE-MATCHED COHORT SIZE\n")

print(
  dim(
    hip_gaba_genotype_matched
  )
)


cat("\nGENOTYPE-MATCHED AGE STRUCTURE\n")

print(
  hip_gaba_genotype_matched[
    ,
    .(
      n_cells = .N,
      n_donors = uniqueN(donor_label)
    ),
    by = donor_age_category
  ]
)


cat("\nGENOTYPE-MATCHED DONORS\n")

print(
  unique(
    hip_gaba_genotype_matched[
      ,
      .(
        donor_label,
        donor_age_category,
        donor_age,
        donor_sex,
        donor_genotype
      )
    ]
  )[
    order(
      donor_age_category,
      donor_label
    )
  ]
)


# ------------------------------------------------------------
# Expression table
#
# primary_expr already contains expression values for
# all unbiased cells, so simply subset it by genotype.
# ------------------------------------------------------------

genotype_matched_expr <- primary_expr[
  donor_genotype %in%
    common_genotypes
]


# Verify dimensions correspond to the filtered metadata

if (
  nrow(genotype_matched_expr) !=
  nrow(hip_gaba_genotype_matched)
) {
  
  stop(
    "Genotype-matched expression and metadata sizes differ."
  )
}


cat(
  "\nGenotype-matched expression cells:",
  nrow(genotype_matched_expr),
  "\n"
)


cat(
  "Genotype-matched donors:",
  uniqueN(
    genotype_matched_expr$donor_label
  ),
  "\n"
)


# Create donor-level expression summary

genotype_matched_donor_summary <-
  make_donor_expression_summary(
    genotype_matched_expr
  )


cat("\nGENOTYPE-MATCHED DONOR EXPRESSION SUMMARY\n")

print(
  genotype_matched_donor_summary
)


# Run exact donor-level tests

genotype_matched_tests <-
  run_sensitivity_tests(
    genotype_matched_donor_summary,
    "Genotype-matched unbiased"
  )


cat("\nGENOTYPE-MATCHED EXACT AGE TESTS\n")

print(
  genotype_matched_tests
)


cat("\nGENOTYPE-MATCHED PERMUTATION ENUMERATION\n")

print(
  genotype_matched_tests[
    ,
    .(
      outcome_label,
      n_adult,
      n_aged,
      n_permutations
    )
  ]
)


# Save sensitivity-2 outputs

fwrite(
  genotype_matched_donor_summary,
  "results/tables/sensitivity_genotype_matched_donor_summary.csv"
)


fwrite(
  genotype_matched_tests,
  "results/tables/sensitivity_genotype_matched_exact_tests.csv"
)



# ============================================================
# 15C. Combine primary and sensitivity analyses
#     for easy comparison
# ============================================================

primary_tests_for_comparison <- copy(
  donor_age_tests
)


primary_tests_for_comparison[
  ,
  analysis :=
    "Primary unbiased"
]


setcolorder(
  primary_tests_for_comparison,
  c(
    "analysis",
    setdiff(
      names(
        primary_tests_for_comparison
      ),
      "analysis"
    )
  )
)


combined_age_tests <- rbindlist(
  list(
    primary_tests_for_comparison,
    all_sampling_tests,
    genotype_matched_tests
  ),
  fill = TRUE
)


cat("\n")
cat("============================================================\n")
cat("PRIMARY + SENSITIVITY ANALYSES\n")
cat("============================================================\n")


# Display the central co-expression result first

cat("\nChrna3 + Chrnb2 DOUBLE-POSITIVE RESULT\n")

print(
  combined_age_tests[
    outcome == "pct_double",
    .(
      analysis,
      n_adult,
      n_aged,
      adult_mean,
      aged_mean,
      difference_aged_minus_adult,
      exact_permutation_p,
      holm_adjusted_p,
      n_permutations
    )
  ]
)


cat("\nALL FIVE OUTCOMES\n")

print(
  combined_age_tests[
    ,
    .(
      analysis,
      outcome_label,
      n_adult,
      n_aged,
      adult_mean,
      aged_mean,
      difference_aged_minus_adult,
      exact_permutation_p,
      holm_adjusted_p
    )
  ]
)


# Save combined comparison table

fwrite(
  combined_age_tests,
  "results/tables/primary_and_sensitivity_age_tests.csv"
)


cat("\n")
cat("============================================================\n")
cat("SECTION 15 SENSITIVITY ANALYSES COMPLETE\n")
cat("============================================================\n")
# ------------------------------------------------------------
# 16. Subclass-level Chrna3 / Chrnb2 co-expression analysis
#
# Primary cohort only:
#   HPF - HIP
#   CTX-CGE + CTX-MGE GABA
#   unbiased sampling
#
# Primary subclass endpoint:
#   percentage of cells that are Chrna3+Chrnb2+
#
# Biological replicate:
#   donor / mouse
# ------------------------------------------------------------


cat("\n")
cat("============================================================\n")
cat("SECTION 16: SUBCLASS-LEVEL CO-EXPRESSION ANALYSIS\n")
cat("============================================================\n")


# ============================================================
# 16A. Count cells and double-positive cells
#      for each donor x subclass combination
# ============================================================

donor_subclass_expression_observed <- primary_expr[
  ,
  .(
    n_cells = .N,
    
    n_chrna3_pos =
      sum(
        Chrna3_detected
      ),
    
    n_chrnb2_pos =
      sum(
        Chrnb2_detected
      ),
    
    n_double_pos =
      sum(
        Chrna3_detected &
          Chrnb2_detected
      )
  ),
  by = .(
    donor_label,
    donor_age_category,
    subclass
  )
]


# ------------------------------------------------------------
# Create the complete donor x subclass grid.
#
# This preserves donors with zero recovered cells in a
# particular subclass.
# ------------------------------------------------------------

all_primary_subclasses <- sort(
  unique(
    primary_expr$subclass
  )
)


subclass_expression_grid <- CJ(
  donor_label =
    primary_donors$donor_label,
  
  subclass =
    all_primary_subclasses,
  
  unique = TRUE
)


# Add donor age group

subclass_expression_grid <- merge(
  subclass_expression_grid,
  
  primary_donors[
    ,
    .(
      donor_label,
      donor_age_category
    )
  ],
  
  by = "donor_label",
  all.x = TRUE,
  sort = FALSE
)


# Add observed expression counts

donor_subclass_expression <- merge(
  subclass_expression_grid,
  
  donor_subclass_expression_observed,
  
  by = c(
    "donor_label",
    "donor_age_category",
    "subclass"
  ),
  
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# Explicitly represent zero recovered cells
# ------------------------------------------------------------

donor_subclass_expression[
  is.na(n_cells),
  `:=`(
    n_cells = 0L,
    n_chrna3_pos = 0L,
    n_chrnb2_pos = 0L,
    n_double_pos = 0L
  )
]


# ------------------------------------------------------------
# Calculate donor-level percentages.
#
# IMPORTANT:
#
# If n_cells == 0, pct_double is undefined (NA),
# not zero.
#
# If n_cells > 0 but n_double_pos == 0,
# pct_double is legitimately 0%.
# ------------------------------------------------------------

donor_subclass_expression[
  ,
  pct_double :=
    fifelse(
      n_cells > 0,
      100 * n_double_pos /
        n_cells,
      NA_real_
    )
]


donor_subclass_expression[
  ,
  pct_chrna3 :=
    fifelse(
      n_cells > 0,
      100 * n_chrna3_pos /
        n_cells,
      NA_real_
    )
]


donor_subclass_expression[
  ,
  pct_chrnb2 :=
    fifelse(
      n_cells > 0,
      100 * n_chrnb2_pos /
        n_cells,
      NA_real_
    )
]


# Sort for inspection

setorder(
  donor_subclass_expression,
  subclass,
  donor_age_category,
  donor_label
)


cat("\nDONOR x SUBCLASS EXPRESSION TABLE\n")
cat("---------------------------------\n")

print(
  donor_subclass_expression
)


# Save complete donor-level subclass table

fwrite(
  donor_subclass_expression,
  "results/tables/primary_donor_subclass_expression.csv"
)



# ============================================================
# 16B. Descriptive subclass summary by age
# ============================================================

subclass_age_summary <- donor_subclass_expression[
  !is.na(pct_double),
  .(
    n_donors_with_cells = .N,
    
    total_cells =
      sum(n_cells),
    
    total_double_pos =
      sum(n_double_pos),
    
    mean_pct_double =
      mean(pct_double),
    
    median_pct_double =
      median(pct_double),
    
    min_pct_double =
      min(pct_double),
    
    max_pct_double =
      max(pct_double)
  ),
  by = .(
    subclass,
    donor_age_category
  )
][
  order(
    subclass,
    donor_age_category
  )
]


cat("\nSUBCLASS CO-EXPRESSION SUMMARY BY AGE\n")
cat("-------------------------------------\n")

print(
  subclass_age_summary
)


fwrite(
  subclass_age_summary,
  "results/tables/primary_subclass_double_positive_age_summary.csv"
)



# ============================================================
# 16C. Exact donor-level permutation test
#      for double-positive percentage within each subclass
# ============================================================


# ------------------------------------------------------------
# Wrapper function for one subclass
# ------------------------------------------------------------

test_one_subclass <- function(
    subclass_name
) {
  
  d <- donor_subclass_expression[
    subclass == subclass_name &
      !is.na(pct_double)
  ]
  
  
  # Count evaluable donors
  n_adult_available <- d[
    donor_age_category == "adult",
    .N
  ]
  
  n_aged_available <- d[
    donor_age_category == "aged",
    .N
  ]
  
  
  # Require at least two evaluable donors in each group.
  #
  # All current subclasses are expected to satisfy this,
  # but the safeguard makes the script robust.
  if (
    n_adult_available < 2 ||
    n_aged_available < 2
  ) {
    
    return(
      data.table(
        subclass = subclass_name,
        
        n_adult = n_adult_available,
        n_aged = n_aged_available,
        
        adult_mean = NA_real_,
        adult_median = NA_real_,
        
        aged_mean = NA_real_,
        aged_median = NA_real_,
        
        difference_aged_minus_adult =
          NA_real_,
        
        exact_permutation_p =
          NA_real_,
        
        n_permutations =
          NA_integer_
      )
    )
  }
  
  
  # Run the same exact permutation procedure used
  # for the overall analysis.
  
  result <- exact_donor_permutation_test(
    data = d,
    outcome = "pct_double"
  )
  
  
  result[
    ,
    subclass := subclass_name
  ]
  
  
  # Remove redundant outcome column
  result[
    ,
    outcome := NULL
  ]
  
  
  # Put subclass first
  setcolorder(
    result,
    c(
      "subclass",
      setdiff(
        names(result),
        "subclass"
      )
    )
  )
  
  
  return(
    result
  )
}



# ------------------------------------------------------------
# Run the test across all eight subclasses
# ------------------------------------------------------------

subclass_double_tests <- rbindlist(
  lapply(
    all_primary_subclasses,
    test_one_subclass
  ),
  fill = TRUE
)



# ============================================================
# 16D. Add cell-count and donor-coverage information
# ============================================================

subclass_test_qc <- donor_subclass_expression[
  ,
  .(
    total_cells =
      sum(n_cells),
    
    donors_total =
      .N,
    
    donors_with_cells =
      sum(n_cells > 0),
    
    min_cells_per_donor =
      min(n_cells),
    
    median_cells_per_donor =
      as.numeric(
        median(n_cells)
      ),
    
    max_cells_per_donor =
      max(n_cells)
  ),
  by = .(
    subclass,
    donor_age_category
  )
]


# Adult QC columns

adult_subclass_qc <- subclass_test_qc[
  donor_age_category == "adult",
  .(
    subclass,
    
    adult_total_cells =
      total_cells,
    
    adult_donors_with_cells =
      donors_with_cells,
    
    adult_min_cells_per_donor =
      min_cells_per_donor,
    
    adult_median_cells_per_donor =
      median_cells_per_donor,
    
    adult_max_cells_per_donor =
      max_cells_per_donor
  )
]


# Aged QC columns

aged_subclass_qc <- subclass_test_qc[
  donor_age_category == "aged",
  .(
    subclass,
    
    aged_total_cells =
      total_cells,
    
    aged_donors_with_cells =
      donors_with_cells,
    
    aged_min_cells_per_donor =
      min_cells_per_donor,
    
    aged_median_cells_per_donor =
      median_cells_per_donor,
    
    aged_max_cells_per_donor =
      max_cells_per_donor
  )
]


# Add QC information to statistical results

subclass_double_tests <- merge(
  subclass_double_tests,
  adult_subclass_qc,
  by = "subclass",
  all.x = TRUE,
  sort = FALSE
)


subclass_double_tests <- merge(
  subclass_double_tests,
  aged_subclass_qc,
  by = "subclass",
  all.x = TRUE,
  sort = FALSE
)



# ============================================================
# 16E. Holm correction across the eight subclass tests
# ============================================================

subclass_double_tests[
  ,
  holm_adjusted_p :=
    p.adjust(
      exact_permutation_p,
      method = "holm"
    )
]



# ------------------------------------------------------------
# Put key result columns first
# ------------------------------------------------------------

setcolorder(
  subclass_double_tests,
  c(
    "subclass",
    
    "n_adult",
    "n_aged",
    
    "adult_mean",
    "adult_median",
    
    "aged_mean",
    "aged_median",
    
    "difference_aged_minus_adult",
    
    "exact_permutation_p",
    "holm_adjusted_p",
    
    "n_permutations",
    
    "adult_total_cells",
    "aged_total_cells",
    
    "adult_donors_with_cells",
    "aged_donors_with_cells",
    
    "adult_min_cells_per_donor",
    "aged_min_cells_per_donor",
    
    "adult_median_cells_per_donor",
    "aged_median_cells_per_donor",
    
    "adult_max_cells_per_donor",
    "aged_max_cells_per_donor"
  )
)



# Sort by Allen subclass number / name

setorder(
  subclass_double_tests,
  subclass
)



cat("\n")
cat("============================================================\n")
cat("SUBCLASS-LEVEL Chrna3+Chrnb2+ AGE TESTS\n")
cat("============================================================\n")

print(
  subclass_double_tests
)



# ============================================================
# 16F. Compact publication-oriented results table
# ============================================================

subclass_double_compact <- subclass_double_tests[
  ,
  .(
    subclass,
    
    adult_donors =
      n_adult,
    
    aged_donors =
      n_aged,
    
    adult_cells =
      adult_total_cells,
    
    aged_cells =
      aged_total_cells,
    
    adult_mean_pct_double =
      adult_mean,
    
    aged_mean_pct_double =
      aged_mean,
    
    difference_pp =
      difference_aged_minus_adult,
    
    exact_p =
      exact_permutation_p,
    
    holm_p =
      holm_adjusted_p,
    
    n_permutations
  )
]


cat("\nCOMPACT SUBCLASS RESULTS\n")
cat("------------------------\n")

print(
  subclass_double_compact
)



# ============================================================
# 16G. Identify sparse subclasses for cautious interpretation
# ============================================================

subclass_double_compact[
  ,
  interpretation_flag :=
    fifelse(
      
      adult_cells < 50 |
        aged_cells < 50,
      
      "Sparse cell representation — interpret cautiously",
      
      "Adequate for donor-level exploratory comparison"
    )
]


cat("\nSUBCLASS INTERPRETATION FLAGS\n")
cat("-----------------------------\n")

print(
  subclass_double_compact[
    ,
    .(
      subclass,
      adult_cells,
      aged_cells,
      adult_donors,
      aged_donors,
      interpretation_flag
    )
  ]
)



# ============================================================
# 16H. Save results
# ============================================================

fwrite(
  subclass_double_tests,
  "results/tables/primary_subclass_double_positive_exact_tests.csv"
)


fwrite(
  subclass_double_compact,
  "results/tables/primary_subclass_double_positive_compact.csv"
)


cat("\n")
cat("============================================================\n")
cat("SECTION 16 SUBCLASS ANALYSIS COMPLETE\n")
cat("============================================================\n")
# ------------------------------------------------------------
# 17. Publication-quality figures
#
# Primary cohort:
#   HPF - HIP
#   CTX-CGE + CTX-MGE GABA
#   unbiased sampling
#
# Biological replicate:
#   donor / mouse
#
# Figures:
#   17A. Overall donor-level expression by age
#   17B. Double-positive co-expression by subclass
#   17C. Subclass expression profiles
# ------------------------------------------------------------


# ============================================================
# 17A. Plotting setup
# ============================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "The ggplot2 package is required for Section 17. ",
    "Install it with: install.packages('ggplot2')"
  )
}

library(ggplot2)


# Make sure output folder exists
dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Consistent age-group labels and order
# ------------------------------------------------------------

age_levels <- c(
  "adult",
  "aged"
)


age_labels <- c(
  adult = "Young adult\n(~2 months)",
  aged = "Aged\n(~18 months)"
)


# ------------------------------------------------------------
# Publication theme
# ------------------------------------------------------------

theme_allen <- theme_classic(
  base_size = 12
) +
  theme(
    
    axis.title = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 10
    ),
    
    plot.title = element_text(
      size = 13,
      face = "bold"
    ),
    
    plot.subtitle = element_text(
      size = 10
    ),
    
    plot.margin = margin(
      8,
      10,
      8,
      8
    )
  )



# ============================================================
# 17B. FIGURE 17A
# Overall donor-level expression by age
# ============================================================


# ------------------------------------------------------------
# Convert donor summary to long format
# ------------------------------------------------------------

overall_plot_data <- melt(
  copy(
    donor_expression_summary
  ),
  
  id.vars = c(
    "donor_label",
    "donor_age_category"
  ),
  
  measure.vars = c(
    "pct_chrna3",
    "pct_chrnb2",
    "pct_double"
  ),
  
  variable.name = "measure",
  value.name = "percent"
)


# ------------------------------------------------------------
# Friendly labels
# ------------------------------------------------------------

overall_plot_data[
  ,
  measure_label := fcase(
    
    measure == "pct_chrna3",
    "Chrna3+",
    
    measure == "pct_chrnb2",
    "Chrnb2+",
    
    measure == "pct_double",
    "Chrna3+Chrnb2+",
    
    default = measure
  )
]


overall_plot_data[
  ,
  donor_age_category := factor(
    donor_age_category,
    levels = age_levels
  )
]


overall_plot_data[
  ,
  measure_label := factor(
    measure_label,
    levels = c(
      "Chrna3+",
      "Chrnb2+",
      "Chrna3+Chrnb2+"
    )
  )
]



# ------------------------------------------------------------
# Mean and SD for display
# ------------------------------------------------------------

overall_summary_plot <- overall_plot_data[
  ,
  .(
    mean_percent =
      mean(
        percent
      ),
    
    sd_percent =
      sd(
        percent
      )
  ),
  by = .(
    donor_age_category,
    measure_label
  )
]



# ------------------------------------------------------------
# Exact permutation P-value labels from Section 14
# ------------------------------------------------------------

overall_pvalues <- donor_age_tests[
  outcome %in% c(
    "pct_chrna3",
    "pct_chrnb2",
    "pct_double"
  ),
  
  .(
    measure_label = fcase(
      
      outcome == "pct_chrna3",
      "Chrna3+",
      
      outcome == "pct_chrnb2",
      "Chrnb2+",
      
      outcome == "pct_double",
      "Chrna3+Chrnb2+"
    ),
    
    p_label =
      sprintf(
        "Exact P = %.3f",
        exact_permutation_p
      )
  )
]


overall_pvalues[
  ,
  measure_label := factor(
    measure_label,
    levels = c(
      "Chrna3+",
      "Chrnb2+",
      "Chrna3+Chrnb2+"
    )
  )
]


# ------------------------------------------------------------
# Set P-value annotation height separately for each facet
# ------------------------------------------------------------

overall_annotation_y <- overall_plot_data[
  ,
  .(
    annotation_y =
      max(
        percent,
        na.rm = TRUE
      ) * 1.18 +
      0.5
  ),
  by = measure_label
]


overall_pvalues <- merge(
  overall_pvalues,
  overall_annotation_y,
  by = "measure_label",
  all.x = TRUE
)



# ------------------------------------------------------------
# Build Figure 17A
# ------------------------------------------------------------

fig17A <- ggplot(
  overall_plot_data,
  
  aes(
    x = donor_age_category,
    y = percent
  )
) +
  
  # Mean ± SD
  geom_errorbar(
    data = overall_summary_plot,
    
    mapping = aes(
      x = donor_age_category,
      
      ymin =
        pmax(
          mean_percent -
            sd_percent,
          0
        ),
      
      ymax =
        mean_percent +
        sd_percent
    ),
    
    width = 0.12,
    linewidth = 0.6,
    inherit.aes = FALSE
  ) +
  
  # Horizontal mean marker
  geom_point(
    data = overall_summary_plot,
    
    aes(
      x = donor_age_category,
      y = mean_percent
    ),
    
    shape = 95,
    size = 10,
    inherit.aes = FALSE
  ) +
  
  # Individual donor mice
  geom_point(
    aes(
      shape = donor_age_category,
      fill = donor_age_category
    ),
    
    size = 3.1,
    stroke = 0.9,
    
    position = position_jitter(
      width = 0.06,
      height = 0,
      seed = 123
    )
  ) +
  
  # Exact permutation P values
  geom_text(
    data = overall_pvalues,
    
    aes(
      x = 1.5,
      y = annotation_y,
      label = p_label
    ),
    
    inherit.aes = FALSE,
    size = 3.2
  ) +
  
  facet_wrap(
    ~ measure_label,
    scales = "free_y",
    nrow = 1
  ) +
  
  # Force all panels to include zero
  expand_limits(
    y = 0
  ) +
  
  scale_x_discrete(
    labels = age_labels
  ) +
  
  scale_shape_manual(
    values = c(
      adult = 21,
      aged = 24
    )
  ) +
  
  scale_fill_manual(
    values = c(
      adult = "white",
      aged = "grey40"
    )
  ) +
  
  labs(
    x = NULL,
    
    y = "Cells positive (%)",
    
    title =
      "Chrna3 and Chrnb2 expression in hippocampal GABAergic interneurons",
    
    subtitle =
      paste0(
        "Points represent individual mice ",
        "(young adult n = 4; aged n = 5); ",
        "horizontal bars indicate means and whiskers indicate ± SD"
      )
  ) +
  
  theme_allen +
  
  theme(
    legend.position = "none"
  )


print(
  fig17A
)



# ------------------------------------------------------------
# Save Figure 17A
# ------------------------------------------------------------

ggsave(
  filename =
    "figures/Figure17A_overall_donor_age_expression.png",
  
  plot = fig17A,
  
  width = 8.5,
  height = 4.2,
  units = "in",
  dpi = 300
)


ggsave(
  filename =
    "figures/Figure17A_overall_donor_age_expression.tiff",
  
  plot = fig17A,
  
  width = 8.5,
  height = 4.2,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename =
    "figures/Figure17A_overall_donor_age_expression.pdf",
  
  plot = fig17A,
  
  width = 8.5,
  height = 4.2,
  units = "in"
)



# ============================================================
# 17C. FIGURE 17B
# Chrna3+Chrnb2+ co-expression by subclass
# ============================================================


# ------------------------------------------------------------
# Cleaner subclass labels
# ------------------------------------------------------------

subclass_label_lookup <- data.table(
  
  subclass = c(
    "046 Vip Gaba",
    "047 Sncg Gaba",
    "048 RHP-COA Ndnf Gaba",
    "049 Lamp5 Gaba",
    "050 Lamp5 Lhx6 Gaba",
    "051 Pvalb chandelier Gaba",
    "052 Pvalb Gaba",
    "053 Sst Gaba"
  ),
  
  subclass_short = c(
    "Vip",
    "Sncg",
    "RHP-COA Ndnf",
    "Lamp5",
    "Lamp5 Lhx6",
    "Pvalb chandelier",
    "Pvalb",
    "Sst"
  )
)



# ------------------------------------------------------------
# Add short subclass labels to donor-level data
# ------------------------------------------------------------

subclass_plot_data <- merge(
  
  copy(
    donor_subclass_expression
  ),
  
  subclass_label_lookup,
  
  by = "subclass",
  all.x = TRUE,
  sort = FALSE
)


subclass_plot_data[
  ,
  donor_age_category := factor(
    donor_age_category,
    levels = age_levels
  )
]


subclass_plot_data[
  ,
  subclass_short := factor(
    subclass_short,
    levels =
      subclass_label_lookup$subclass_short
  )
]


# Only donors with recovered cells in a subclass
# have a defined percentage.

subclass_plot_data_defined <-
  subclass_plot_data[
    !is.na(pct_double)
  ]



# ------------------------------------------------------------
# Mean and SD by subclass and age
# ------------------------------------------------------------

subclass_double_summary_plot <-
  subclass_plot_data_defined[
    ,
    .(
      mean_pct_double =
        mean(
          pct_double
        ),
      
      sd_pct_double =
        sd(
          pct_double
        )
    ),
    
    by = .(
      subclass_short,
      donor_age_category
    )
  ]



# ------------------------------------------------------------
# Exact subclass P values and donor sample sizes
# ------------------------------------------------------------

subclass_pvalues <- merge(
  
  subclass_double_tests[
    ,
    .(
      subclass,
      exact_permutation_p,
      n_adult,
      n_aged
    )
  ],
  
  subclass_label_lookup,
  
  by = "subclass",
  all.x = TRUE
)



# ------------------------------------------------------------
# Position P-value labels above each subclass panel
# ------------------------------------------------------------

subclass_annotation_y <-
  subclass_plot_data_defined[
    ,
    .(
      annotation_y =
        max(
          pct_double,
          na.rm = TRUE
        ) * 1.14 +
        1
    ),
    
    by = subclass_short
  ]


subclass_pvalues <- merge(
  subclass_pvalues,
  subclass_annotation_y,
  by = "subclass_short",
  all.x = TRUE
)


subclass_pvalues[
  ,
  p_label :=
    sprintf(
      "P = %.3f   (n = %d/%d)",
      exact_permutation_p,
      n_adult,
      n_aged
    )
]


subclass_pvalues[
  ,
  subclass_short := factor(
    subclass_short,
    levels =
      subclass_label_lookup$subclass_short
  )
]



# ------------------------------------------------------------
# Build Figure 17B
# ------------------------------------------------------------

fig17B <- ggplot(
  subclass_plot_data_defined,
  
  aes(
    x = donor_age_category,
    y = pct_double
  )
) +
  
  # Mean ± SD
  geom_errorbar(
    data =
      subclass_double_summary_plot,
    
    mapping = aes(
      x = donor_age_category,
      
      ymin =
        pmax(
          mean_pct_double -
            sd_pct_double,
          0
        ),
      
      ymax =
        mean_pct_double +
        sd_pct_double
    ),
    
    width = 0.12,
    linewidth = 0.55,
    inherit.aes = FALSE
  ) +
  
  # Mean marker
  geom_point(
    data =
      subclass_double_summary_plot,
    
    mapping = aes(
      x = donor_age_category,
      y = mean_pct_double
    ),
    
    shape = 95,
    size = 8,
    inherit.aes = FALSE
  ) +
  
  # Individual donor values
  geom_point(
    aes(
      shape = donor_age_category,
      fill = donor_age_category
    ),
    
    size = 2.7,
    stroke = 0.8,
    
    position =
      position_jitter(
        width = 0.06,
        height = 0,
        seed = 456
      )
  ) +
  
  # Exact permutation P values + evaluable donor counts
  geom_text(
    data = subclass_pvalues,
    
    aes(
      x = 1.5,
      y = annotation_y,
      label = p_label
    ),
    
    inherit.aes = FALSE,
    size = 2.6
  ) +
  
  facet_wrap(
    ~ subclass_short,
    scales = "free_y",
    ncol = 4
  ) +
  
  # All panels must include zero
  expand_limits(
    y = 0
  ) +
  
  scale_x_discrete(
    labels = age_labels
  ) +
  
  scale_shape_manual(
    values = c(
      adult = 21,
      aged = 24
    )
  ) +
  
  scale_fill_manual(
    values = c(
      adult = "white",
      aged = "grey40"
    )
  ) +
  
  labs(
    x = NULL,
    
    y = "Chrna3+Chrnb2+ cells (%)",
    
    title =
      "Chrna3/Chrnb2 co-expression by hippocampal GABAergic subclass",
    
    subtitle =
      paste0(
        "Points represent individual mice; ",
        "horizontal bars indicate donor means and whiskers indicate ± SD; ",
        "P values are exact permutation tests"
      )
  ) +
  
  theme_allen +
  
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      size = 7
    ),
    
    strip.text = element_text(
      size = 9,
      face = "bold"
    )
  )


print(
  fig17B
)



# ------------------------------------------------------------
# Save Figure 17B
# ------------------------------------------------------------

ggsave(
  filename =
    "figures/Figure17B_subclass_double_positive_age.png",
  
  plot = fig17B,
  
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 300
)


ggsave(
  filename =
    "figures/Figure17B_subclass_double_positive_age.tiff",
  
  plot = fig17B,
  
  width = 10,
  height = 6.5,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename =
    "figures/Figure17B_subclass_double_positive_age.pdf",
  
  plot = fig17B,
  
  width = 10,
  height = 6.5,
  units = "in"
)



# ============================================================
# 17D. FIGURE 17C
# Subclass expression profile:
# Chrna3+, Chrnb2+, and Chrna3+Chrnb2+
#
# NOTE:
# Subclasses are categorical populations.
# Therefore, points are NOT connected by lines.
# ============================================================


# ------------------------------------------------------------
# Convert donor x subclass table to long format
# ------------------------------------------------------------

subclass_profile_data <- melt(
  
  subclass_plot_data[
    n_cells > 0
  ],
  
  id.vars = c(
    "donor_label",
    "donor_age_category",
    "subclass",
    "subclass_short"
  ),
  
  measure.vars = c(
    "pct_chrna3",
    "pct_chrnb2",
    "pct_double"
  ),
  
  variable.name =
    "expression_measure",
  
  value.name =
    "percent"
)



# ------------------------------------------------------------
# Friendly expression labels
# ------------------------------------------------------------

subclass_profile_data[
  ,
  expression_label := fcase(
    
    expression_measure ==
      "pct_chrna3",
    "Chrna3+",
    
    expression_measure ==
      "pct_chrnb2",
    "Chrnb2+",
    
    expression_measure ==
      "pct_double",
    "Chrna3+Chrnb2+",
    
    default =
      expression_measure
  )
]


subclass_profile_data[
  ,
  expression_label := factor(
    expression_label,
    
    levels = c(
      "Chrna3+",
      "Chrnb2+",
      "Chrna3+Chrnb2+"
    )
  )
]


subclass_profile_data[
  ,
  donor_age_category := factor(
    donor_age_category,
    levels = age_levels
  )
]


subclass_profile_data[
  ,
  subclass_short := factor(
    subclass_short,
    levels =
      subclass_label_lookup$subclass_short
  )
]



# ------------------------------------------------------------
# Calculate donor means and SD
# ------------------------------------------------------------

subclass_profile_summary <-
  subclass_profile_data[
    ,
    .(
      mean_percent =
        mean(
          percent,
          na.rm = TRUE
        ),
      
      sd_percent =
        sd(
          percent,
          na.rm = TRUE
        ),
      
      n_donors =
        sum(
          !is.na(percent)
        )
    ),
    
    by = .(
      subclass_short,
      donor_age_category,
      expression_label
    )
  ]



# ------------------------------------------------------------
# Build Figure 17C
#
# No connecting lines:
# each subclass is a categorical population.
# ------------------------------------------------------------

fig17C <- ggplot(
  
  subclass_profile_summary,
  
  aes(
    x = subclass_short,
    y = mean_percent,
    group = donor_age_category,
    shape = donor_age_category
  )
) +
  
  # Mean ± SD
  geom_errorbar(
    aes(
      ymin =
        pmax(
          mean_percent -
            sd_percent,
          0
        ),
      
      ymax =
        mean_percent +
        sd_percent
    ),
    
    width = 0.12,
    
    position = position_dodge(
      width = 0.35
    ),
    
    linewidth = 0.45
  ) +
  
  # Mean points
  geom_point(
    aes(
      fill = donor_age_category
    ),
    
    position = position_dodge(
      width = 0.35
    ),
    
    size = 2.8,
    stroke = 0.8
  ) +
  
  facet_wrap(
    ~ expression_label,
    scales = "free_y",
    ncol = 1
  ) +
  
  # Make zero visible in every expression panel
  expand_limits(
    y = 0
  ) +
  
  scale_shape_manual(
    values = c(
      adult = 21,
      aged = 24
    ),
    
    labels = c(
      adult = "Young adult (~2 months)",
      aged = "Aged (~18 months)"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      adult = "white",
      aged = "grey40"
    ),
    
    labels = c(
      adult = "Young adult (~2 months)",
      aged = "Aged (~18 months)"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      
      "Vip" =
        "Vip",
      
      "Sncg" =
        "Sncg",
      
      "RHP-COA Ndnf" =
        "RHP-COA\nNdnf",
      
      "Lamp5" =
        "Lamp5",
      
      "Lamp5 Lhx6" =
        "Lamp5\nLhx6",
      
      "Pvalb chandelier" =
        "Pvalb\nchandelier",
      
      "Pvalb" =
        "Pvalb",
      
      "Sst" =
        "Sst"
    )
  ) +
  
  labs(
    x =
      "Allen WMB GABAergic subclass",
    
    y =
      "Cells positive (%)",
    
    title =
      "Chrna3 and Chrnb2 expression profiles across hippocampal GABAergic subclasses",
    
    subtitle =
      "Points show donor means; error bars indicate ± SD",
    
    shape = NULL,
    fill = NULL
  ) +
  
  theme_allen +
  
  theme(
    legend.position = "top",
    
    axis.text.x = element_text(
      size = 8
    ),
    
    strip.text = element_text(
      size = 9,
      face = "bold"
    )
  )


print(
  fig17C
)



# ------------------------------------------------------------
# Save Figure 17C
# ------------------------------------------------------------

ggsave(
  filename =
    "figures/Figure17C_subclass_expression_profiles.png",
  
  plot = fig17C,
  
  width = 9.5,
  height = 9,
  units = "in",
  dpi = 300
)


ggsave(
  filename =
    "figures/Figure17C_subclass_expression_profiles.tiff",
  
  plot = fig17C,
  
  width = 9.5,
  height = 9,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


ggsave(
  filename =
    "figures/Figure17C_subclass_expression_profiles.pdf",
  
  plot = fig17C,
  
  width = 9.5,
  height = 9,
  units = "in"
)



# ============================================================
# 17E. Save plotting data
#
# Saving the exact values used to construct each figure
# improves reproducibility and manuscript review.
# ============================================================

fwrite(
  overall_plot_data,
  "results/tables/Figure17A_plot_data.csv"
)


fwrite(
  overall_summary_plot,
  "results/tables/Figure17A_plot_summary.csv"
)


fwrite(
  subclass_plot_data_defined,
  "results/tables/Figure17B_plot_data.csv"
)


fwrite(
  subclass_double_summary_plot,
  "results/tables/Figure17B_plot_summary.csv"
)


fwrite(
  subclass_profile_summary,
  "results/tables/Figure17C_plot_summary.csv"
)



# ============================================================
# 17F. Completion message
# ============================================================

cat("\n")

cat(
  "============================================================\n"
)

cat(
  "SECTION 17 FIGURE GENERATION COMPLETE\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nFigures saved in the 'figures' folder:\n"
)


cat(
  "  Figure17A_overall_donor_age_expression\n"
)


cat(
  "  Figure17B_subclass_double_positive_age\n"
)


cat(
  "  Figure17C_subclass_expression_profiles\n"
)


cat(
  "\nEach figure was saved as PNG, TIFF, and PDF.\n"
)
# ------------------------------------------------------------
# ------------------------------------------------------------
# 18. Donor-aware transcript-count intensity analysis
#
# Question:
# Among cells in which Chrna3 and/or Chrnb2 are detected,
# does transcript-count intensity differ between
# young-adult and aged mice?
#
# Primary cohort:
#   HPF - HIP
#   CTX-CGE + CTX-MGE GABA
#   unbiased sampling
#
# Biological replicate:
#   donor / mouse
#
# IMPORTANT:
# Raw UMI counts are used here.
# These values are not library-size normalized.
# Therefore, this analysis is complementary / exploratory.
#
# Primary statistical summaries:
#   donor mean log1p(raw UMI count)
#
# log1p(x) = log(1 + x)
# ------------------------------------------------------------


cat("\n")
cat("============================================================\n")
cat("SECTION 18: POSITIVE-CELL TRANSCRIPT-COUNT INTENSITY\n")
cat("============================================================\n")



# ============================================================
# 18A. Robust helper functions
#
# Explicit as.numeric() prevents data.table from receiving
# integer values for some donors and double values for others.
# ============================================================

safe_mean <- function(x) {
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  as.numeric(
    mean(x)
  )
}


safe_median <- function(x) {
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  as.numeric(
    median(x)
  )
}


safe_mean_log1p <- function(x) {
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  as.numeric(
    mean(
      log1p(x)
    )
  )
}



# ============================================================
# 18B. Create donor-level positive-cell summaries
# ============================================================

donor_intensity_summary <- primary_expr[
  ,
  {
    
    # --------------------------------------------------------
    # Chrna3-positive cells
    # --------------------------------------------------------
    
    chrna3_values <-
      Chrna3_raw[
        Chrna3_raw > 0
      ]
    
    
    # --------------------------------------------------------
    # Chrnb2-positive cells
    # --------------------------------------------------------
    
    chrnb2_values <-
      Chrnb2_raw[
        Chrnb2_raw > 0
      ]
    
    
    # --------------------------------------------------------
    # Double-positive cells
    # --------------------------------------------------------
    
    double_index <-
      Chrna3_raw > 0 &
      Chrnb2_raw > 0
    
    
    chrna3_double_values <-
      Chrna3_raw[
        double_index
      ]
    
    
    chrnb2_double_values <-
      Chrnb2_raw[
        double_index
      ]
    
    
    # --------------------------------------------------------
    # Return donor-level summaries
    # --------------------------------------------------------
    
    list(
      
      # Cell counts
      n_cells =
        as.integer(.N),
      
      n_chrna3_positive =
        as.integer(
          length(
            chrna3_values
          )
        ),
      
      n_chrnb2_positive =
        as.integer(
          length(
            chrnb2_values
          )
        ),
      
      n_double_positive =
        as.integer(
          sum(
            double_index
          )
        ),
      
      
      # ------------------------------------------------------
      # Chrna3 among Chrna3-positive cells
      # ------------------------------------------------------
      
      chrna3_pos_mean_raw =
        safe_mean(
          chrna3_values
        ),
      
      chrna3_pos_median_raw =
        safe_median(
          chrna3_values
        ),
      
      chrna3_pos_mean_log1p =
        safe_mean_log1p(
          chrna3_values
        ),
      
      
      # ------------------------------------------------------
      # Chrnb2 among Chrnb2-positive cells
      # ------------------------------------------------------
      
      chrnb2_pos_mean_raw =
        safe_mean(
          chrnb2_values
        ),
      
      chrnb2_pos_median_raw =
        safe_median(
          chrnb2_values
        ),
      
      chrnb2_pos_mean_log1p =
        safe_mean_log1p(
          chrnb2_values
        ),
      
      
      # ------------------------------------------------------
      # Chrna3 among double-positive cells
      # ------------------------------------------------------
      
      chrna3_double_mean_raw =
        safe_mean(
          chrna3_double_values
        ),
      
      chrna3_double_median_raw =
        safe_median(
          chrna3_double_values
        ),
      
      chrna3_double_mean_log1p =
        safe_mean_log1p(
          chrna3_double_values
        ),
      
      
      # ------------------------------------------------------
      # Chrnb2 among double-positive cells
      # ------------------------------------------------------
      
      chrnb2_double_mean_raw =
        safe_mean(
          chrnb2_double_values
        ),
      
      chrnb2_double_median_raw =
        safe_median(
          chrnb2_double_values
        ),
      
      chrnb2_double_mean_log1p =
        safe_mean_log1p(
          chrnb2_double_values
        )
    )
  },
  
  by = .(
    donor_label,
    donor_age_category,
    donor_age,
    donor_sex,
    donor_genotype
  )
]



# ------------------------------------------------------------
# Sort donors
# ------------------------------------------------------------

setorder(
  donor_intensity_summary,
  donor_age_category,
  donor_label
)


cat("\nDONOR-LEVEL TRANSCRIPT-COUNT SUMMARY\n")
cat("------------------------------------\n")

print(
  donor_intensity_summary
)



# ============================================================
# 18C. Verify positive-cell totals
# ============================================================

cat("\nPOSITIVE-CELL COUNT CHECK\n")
cat("-------------------------\n")


print(
  donor_intensity_summary[
    ,
    .(
      donor_label,
      donor_age_category,
      n_cells,
      n_chrna3_positive,
      n_chrnb2_positive,
      n_double_positive
    )
  ]
)


cat(
  "\nTotal analyzed cells:",
  sum(
    donor_intensity_summary$n_cells
  ),
  "\n"
)


cat(
  "Total Chrna3-positive cells:",
  sum(
    donor_intensity_summary$n_chrna3_positive
  ),
  "\n"
)


cat(
  "Total Chrnb2-positive cells:",
  sum(
    donor_intensity_summary$n_chrnb2_positive
  ),
  "\n"
)


cat(
  "Total double-positive cells:",
  sum(
    donor_intensity_summary$n_double_positive
  ),
  "\n"
)



# ============================================================
# 18D. Define the four planned intensity outcomes
#
# Primary statistic:
# donor mean log1p(raw UMI count)
# ============================================================

intensity_outcomes <- c(
  
  "chrna3_pos_mean_log1p",
  
  "chrnb2_pos_mean_log1p",
  
  "chrna3_double_mean_log1p",
  
  "chrnb2_double_mean_log1p"
)



# ============================================================
# 18E. Exact donor-aware permutation tests
#
# Reuses the exact permutation function created in Section 14.
# ============================================================

intensity_tests <- rbindlist(
  
  lapply(
    
    intensity_outcomes,
    
    function(x) {
      
      exact_donor_permutation_test(
        
        data =
          donor_intensity_summary,
        
        outcome =
          x
      )
    }
  ),
  
  fill = TRUE
)



# ------------------------------------------------------------
# Add readable outcome labels
# ------------------------------------------------------------

intensity_tests[
  ,
  outcome_label := fcase(
    
    outcome ==
      "chrna3_pos_mean_log1p",
    "Chrna3 count among Chrna3+ cells",
    
    outcome ==
      "chrnb2_pos_mean_log1p",
    "Chrnb2 count among Chrnb2+ cells",
    
    outcome ==
      "chrna3_double_mean_log1p",
    "Chrna3 count among double-positive cells",
    
    outcome ==
      "chrnb2_double_mean_log1p",
    "Chrnb2 count among double-positive cells",
    
    default =
      outcome
  )
]



# ------------------------------------------------------------
# Holm correction across the four planned tests
# ------------------------------------------------------------

intensity_tests[
  ,
  holm_adjusted_p :=
    p.adjust(
      exact_permutation_p,
      method = "holm"
    )
]



# ------------------------------------------------------------
# Put columns in convenient order
# ------------------------------------------------------------

setcolorder(
  intensity_tests,
  
  c(
    "outcome",
    "outcome_label",
    
    "n_adult",
    "n_aged",
    
    "adult_mean",
    "adult_median",
    
    "aged_mean",
    "aged_median",
    
    "difference_aged_minus_adult",
    
    "exact_permutation_p",
    "holm_adjusted_p",
    
    "n_permutations"
  )
)



cat("\n")
cat("============================================================\n")
cat("DONOR-AWARE TRANSCRIPT-COUNT INTENSITY TESTS\n")
cat("============================================================\n")


print(
  intensity_tests
)



# ============================================================
# 18F. Check permutation enumeration for each outcome
# ============================================================

cat(
  "\nPERMUTATION ENUMERATION BY INTENSITY OUTCOME\n"
)


print(
  intensity_tests[
    ,
    .(
      outcome_label,
      n_adult,
      n_aged,
      n_permutations
    )
  ]
)



# ============================================================
# 18G. Descriptive raw-count summaries by age
#
# These are descriptive only.
# Statistical inference above uses donor mean log1p counts.
# ============================================================

raw_count_descriptive <- donor_intensity_summary[
  ,
  .(
    
    n_donors =
      .N,
    
    
    # --------------------------------------------------------
    # Chrna3-positive cells
    # --------------------------------------------------------
    
    mean_donor_chrna3_pos_raw =
      mean(
        chrna3_pos_mean_raw,
        na.rm = TRUE
      ),
    
    median_donor_chrna3_pos_raw =
      median(
        chrna3_pos_median_raw,
        na.rm = TRUE
      ),
    
    
    # --------------------------------------------------------
    # Chrnb2-positive cells
    # --------------------------------------------------------
    
    mean_donor_chrnb2_pos_raw =
      mean(
        chrnb2_pos_mean_raw,
        na.rm = TRUE
      ),
    
    median_donor_chrnb2_pos_raw =
      median(
        chrnb2_pos_median_raw,
        na.rm = TRUE
      ),
    
    
    # --------------------------------------------------------
    # Chrna3 within double-positive cells
    # --------------------------------------------------------
    
    mean_donor_chrna3_double_raw =
      mean(
        chrna3_double_mean_raw,
        na.rm = TRUE
      ),
    
    median_donor_chrna3_double_raw =
      median(
        chrna3_double_median_raw,
        na.rm = TRUE
      ),
    
    
    # --------------------------------------------------------
    # Chrnb2 within double-positive cells
    # --------------------------------------------------------
    
    mean_donor_chrnb2_double_raw =
      mean(
        chrnb2_double_mean_raw,
        na.rm = TRUE
      ),
    
    median_donor_chrnb2_double_raw =
      median(
        chrnb2_double_median_raw,
        na.rm = TRUE
      )
  ),
  
  by =
    donor_age_category
]



cat("\nDESCRIPTIVE RAW-COUNT SUMMARY BY AGE\n")
cat("------------------------------------\n")


print(
  raw_count_descriptive
)



# ============================================================
# 18H. Prepare plotting data
# ============================================================

intensity_plot_data <- melt(
  
  copy(
    donor_intensity_summary
  ),
  
  id.vars = c(
    "donor_label",
    "donor_age_category"
  ),
  
  measure.vars =
    intensity_outcomes,
  
  variable.name =
    "outcome",
  
  value.name =
    "mean_log1p_umi"
)



# ------------------------------------------------------------
# Add readable facet labels
# ------------------------------------------------------------

intensity_plot_data[
  ,
  outcome_label := fcase(
    
    outcome ==
      "chrna3_pos_mean_log1p",
    "Chrna3 within\nChrna3+ cells",
    
    outcome ==
      "chrnb2_pos_mean_log1p",
    "Chrnb2 within\nChrnb2+ cells",
    
    outcome ==
      "chrna3_double_mean_log1p",
    "Chrna3 within\ndouble-positive cells",
    
    outcome ==
      "chrnb2_double_mean_log1p",
    "Chrnb2 within\ndouble-positive cells",
    
    default =
      outcome
  )
]



intensity_plot_data[
  ,
  donor_age_category := factor(
    donor_age_category,
    levels = age_levels
  )
]



intensity_plot_data[
  ,
  outcome_label := factor(
    
    outcome_label,
    
    levels = c(
      "Chrna3 within\nChrna3+ cells",
      "Chrnb2 within\nChrnb2+ cells",
      "Chrna3 within\ndouble-positive cells",
      "Chrnb2 within\ndouble-positive cells"
    )
  )
]



# ============================================================
# 18I. Calculate plotting means and SD
# ============================================================

intensity_plot_summary <- intensity_plot_data[
  !is.na(mean_log1p_umi),
  
  .(
    
    mean_value =
      mean(
        mean_log1p_umi
      ),
    
    sd_value =
      sd(
        mean_log1p_umi
      ),
    
    n_donors =
      .N
  ),
  
  by = .(
    donor_age_category,
    outcome_label
  )
]



# ============================================================
# 18J. Prepare exact P-value annotations
# ============================================================

intensity_pvalues <- intensity_tests[
  ,
  .(
    
    outcome_label = fcase(
      
      outcome ==
        "chrna3_pos_mean_log1p",
      "Chrna3 within\nChrna3+ cells",
      
      outcome ==
        "chrnb2_pos_mean_log1p",
      "Chrnb2 within\nChrnb2+ cells",
      
      outcome ==
        "chrna3_double_mean_log1p",
      "Chrna3 within\ndouble-positive cells",
      
      outcome ==
        "chrnb2_double_mean_log1p",
      "Chrnb2 within\ndouble-positive cells"
    ),
    
    p_label =
      sprintf(
        "Exact P = %.3f",
        exact_permutation_p
      ),
    
    n_adult,
    
    n_aged
  )
]



intensity_pvalues[
  ,
  outcome_label := factor(
    
    outcome_label,
    
    levels = c(
      "Chrna3 within\nChrna3+ cells",
      "Chrnb2 within\nChrnb2+ cells",
      "Chrna3 within\ndouble-positive cells",
      "Chrnb2 within\ndouble-positive cells"
    )
  )
]



# ------------------------------------------------------------
# Annotation height independently for each facet
# ------------------------------------------------------------

intensity_annotation_y <- intensity_plot_data[
  !is.na(mean_log1p_umi),
  
  .(
    annotation_y =
      max(
        mean_log1p_umi,
        na.rm = TRUE
      ) * 1.12 +
      0.05
  ),
  
  by =
    outcome_label
]



intensity_pvalues <- merge(
  
  intensity_pvalues,
  
  intensity_annotation_y,
  
  by =
    "outcome_label",
  
  all.x = TRUE
)



# ============================================================
# 18K. Build Figure 18
# ============================================================

fig18 <- ggplot(
  
  intensity_plot_data[
    !is.na(mean_log1p_umi)
  ],
  
  aes(
    x =
      donor_age_category,
    
    y =
      mean_log1p_umi
  )
) +
  
  # ----------------------------------------------------------
# Mean ± SD
# ----------------------------------------------------------

geom_errorbar(
  
  data =
    intensity_plot_summary,
  
  mapping = aes(
    
    x =
      donor_age_category,
    
    ymin =
      pmax(
        mean_value -
          sd_value,
        0
      ),
    
    ymax =
      mean_value +
      sd_value
  ),
  
  width = 0.12,
  linewidth = 0.55,
  
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Horizontal mean marker
# ----------------------------------------------------------

geom_point(
  
  data =
    intensity_plot_summary,
  
  mapping = aes(
    
    x =
      donor_age_category,
    
    y =
      mean_value
  ),
  
  shape = 95,
  size = 8,
  
  inherit.aes = FALSE
) +
  
  # ----------------------------------------------------------
# Individual donor mice
# ----------------------------------------------------------

geom_point(
  
  aes(
    shape =
      donor_age_category,
    
    fill =
      donor_age_category
  ),
  
  size = 2.9,
  stroke = 0.8,
  
  position =
    position_jitter(
      width = 0.06,
      height = 0,
      seed = 789
    )
) +
  
  # ----------------------------------------------------------
# Exact permutation P values
# ----------------------------------------------------------

geom_text(
  
  data =
    intensity_pvalues,
  
  mapping = aes(
    
    x = 1.5,
    
    y =
      annotation_y,
    
    label =
      p_label
  ),
  
  inherit.aes = FALSE,
  
  size = 2.9
) +
  
  facet_wrap(
    
    ~ outcome_label,
    
    scales =
      "free_y",
    
    ncol = 2
  ) +
  
  # Force zero into every panel
  expand_limits(
    y = 0
  ) +
  
  scale_x_discrete(
    labels =
      age_labels
  ) +
  
  scale_shape_manual(
    
    values = c(
      adult = 21,
      aged = 24
    )
  ) +
  
  scale_fill_manual(
    
    values = c(
      adult = "white",
      aged = "grey40"
    )
  ) +
  
  labs(
    
    x = NULL,
    
    y =
      "Donor mean log(1 + raw UMI count)",
    
    title =
      paste0(
        "Transcript-count intensity in Chrna3- and Chrnb2-positive ",
        "hippocampal interneurons"
      ),
    
    subtitle =
      paste0(
        "Points represent individual mice; ",
        "horizontal bars indicate donor means and ",
        "whiskers indicate ± SD"
      )
  ) +
  
  theme_allen +
  
  theme(
    
    legend.position =
      "none",
    
    axis.text.x =
      element_text(
        size = 8
      ),
    
    strip.text =
      element_text(
        size = 9,
        face = "bold"
      )
  )



print(
  fig18
)



# ============================================================
# 18L. Save Figure 18
# ============================================================

ggsave(
  
  filename =
    "figures/Figure18_positive_cell_transcript_intensity.png",
  
  plot =
    fig18,
  
  width = 8.5,
  height = 7,
  
  units = "in",
  
  dpi = 300
)



ggsave(
  
  filename =
    "figures/Figure18_positive_cell_transcript_intensity.tiff",
  
  plot =
    fig18,
  
  width = 8.5,
  height = 7,
  
  units = "in",
  
  dpi = 600,
  
  compression =
    "lzw"
)



ggsave(
  
  filename =
    "figures/Figure18_positive_cell_transcript_intensity.pdf",
  
  plot =
    fig18,
  
  width = 8.5,
  height = 7,
  
  units = "in"
)



# ============================================================
# 18M. Save Section 18 data and statistical results
# ============================================================

fwrite(
  
  donor_intensity_summary,
  
  "results/tables/primary_donor_transcript_intensity_summary.csv"
)



fwrite(
  
  intensity_tests,
  
  "results/tables/primary_donor_transcript_intensity_exact_tests.csv"
)



fwrite(
  
  raw_count_descriptive,
  
  "results/tables/primary_transcript_intensity_raw_descriptive.csv"
)



fwrite(
  
  intensity_plot_data,
  
  "results/tables/Figure18_plot_data.csv"
)



fwrite(
  
  intensity_plot_summary,
  
  "results/tables/Figure18_plot_summary.csv"
)



# ============================================================
# 18N. Final checkpoint
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SECTION 18 TRANSCRIPT-COUNT ANALYSIS COMPLETE\n")
cat("============================================================\n")


cat(
  "\nPrimary Section 18 outcomes use ",
  "donor mean log1p(raw UMI count).\n"
)


cat(
  "Raw-count summaries are descriptive and ",
  "should be interpreted cautiously.\n"
)


cat(
  "\nFigure saved as:\n"
)


cat(
  "  Figure18_positive_cell_transcript_intensity.png\n"
)


cat(
  "  Figure18_positive_cell_transcript_intensity.tiff\n"
)


cat(
  "  Figure18_positive_cell_transcript_intensity.pdf\n"
)
# ============================================================
# 19. REVISED MAIN FIGURE 3 + SUPPLEMENTARY FIGURE S2
#
# Main Figure 3:
#   A = overall donor-level age comparison  (current fig17A)
#   B = subclass-level double-positive age comparison (current fig17B)
#
# Supplementary Figure S2:
#   positive-cell transcript-count intensity (current fig18)
# ============================================================

# ------------------------------------------------------------
# 19A. Required package for combining plots
# ------------------------------------------------------------

if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop(
    "The patchwork package is required for Section 19. ",
    "Install it with: install.packages('patchwork')"
  )
}

library(patchwork)

# Make sure output folder exists
dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 19B. Short age labels for publication figure versions
# ------------------------------------------------------------

short_age_labels <- c(
  adult = "Young",
  aged  = "Aged"
)

# ------------------------------------------------------------
# 19C. Clean publication versions of current fig17A and fig17B
#      for assembly into revised Figure 3
# ------------------------------------------------------------

fig17A_main <- suppressWarnings(
  fig17A +
    scale_x_discrete(labels = short_age_labels) +
    labs(
      title = NULL,
      subtitle = NULL
    ) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9),
      strip.text = element_text(size = 10, face = "bold"),
      plot.margin = margin(4, 6, 2, 6)
    )
)

fig17B_main <- suppressWarnings(
  fig17B +
    scale_x_discrete(labels = short_age_labels) +
    labs(
      title = NULL,
      subtitle = NULL
    ) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      strip.text = element_text(size = 9, face = "bold"),
      plot.margin = margin(2, 6, 4, 6)
    )
)

# ------------------------------------------------------------
# 19D. Assemble revised main Figure 3
# ------------------------------------------------------------

Figure3_revised <- (
  fig17A_main /
    fig17B_main
) +
  plot_layout(
    ncol = 1,
    heights = c(0.78, 1.22)
  ) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      size = 16,
      face = "bold"
    )
  )

print(Figure3_revised)

# ------------------------------------------------------------
# 19E. Save revised main Figure 3
#      Final width ~180 mm (7.09 in)
# ------------------------------------------------------------

ggsave(
  filename = "figures/Figure3_Aging_Mouse_Donor_Aware_REVISED.png",
  plot = Figure3_revised,
  width = 7.09,
  height = 8.70,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "figures/Figure3_Aging_Mouse_Donor_Aware_REVISED.tiff",
  plot = Figure3_revised,
  width = 7.09,
  height = 8.70,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "figures/Figure3_Aging_Mouse_Donor_Aware_REVISED.pdf",
  plot = Figure3_revised,
  width = 7.09,
  height = 8.70,
  units = "in"
)

# ------------------------------------------------------------
# 19F. Clean version of current fig18 for Supplementary Figure S2
# ------------------------------------------------------------

Supplementary_Figure_S2 <- suppressWarnings(
  fig18 +
    scale_x_discrete(labels = short_age_labels) +
    labs(
      title = NULL,
      subtitle = NULL
    ) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 8.5),
      axis.text.x = element_text(size = 8.5),
      strip.text = element_text(size = 9, face = "bold"),
      plot.margin = margin(4, 6, 4, 6)
    )
)

print(Supplementary_Figure_S2)

# ------------------------------------------------------------
# 19G. Save Supplementary Figure S2
#      (former Figure 3C / current fig18)
# ------------------------------------------------------------

ggsave(
  filename = "figures/Supplementary_Figure_S2_Transcript_Intensity.png",
  plot = Supplementary_Figure_S2,
  width = 7.09,
  height = 5.80,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "figures/Supplementary_Figure_S2_Transcript_Intensity.tiff",
  plot = Supplementary_Figure_S2,
  width = 7.09,
  height = 5.80,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "figures/Supplementary_Figure_S2_Transcript_Intensity.pdf",
  plot = Supplementary_Figure_S2,
  width = 7.09,
  height = 5.80,
  units = "in"
)

cat("\nSection 19 complete.\n")
cat("Files saved as:\n")
cat("  figures/Figure3_Aging_Mouse_Donor_Aware_REVISED.png\n")
cat("  figures/Figure3_Aging_Mouse_Donor_Aware_REVISED.tiff\n")
cat("  figures/Figure3_Aging_Mouse_Donor_Aware_REVISED.pdf\n")
cat("  figures/Supplementary_Figure_S2_Transcript_Intensity.png\n")
cat("  figures/Supplementary_Figure_S2_Transcript_Intensity.tiff\n")
cat("  figures/Supplementary_Figure_S2_Transcript_Intensity.pdf\n")