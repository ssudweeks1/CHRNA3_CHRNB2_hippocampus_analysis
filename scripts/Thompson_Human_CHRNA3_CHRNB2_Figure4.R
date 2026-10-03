# ============================================================
# Thompson Human Hippocampus CHRNA3 / CHRNB2 Analysis
#
# FINAL REPRODUCIBLE SCRIPT
# ============================================================


# ============================================================
# 1. PROJECT SETUP
# ============================================================

project_dir <-
  "C:/Users/sns34/Documents/Thompson_Human_CHRNA3_CHRNB2"


if (!dir.exists(project_dir)) {
  
  stop(
    paste0(
      "Project directory not found:\n",
      project_dir
    )
  )
}


setwd(project_dir)


cat("\n")
cat("============================================================\n")
cat("THOMPSON HUMAN CHRNA3 / CHRNB2 ANALYSIS\n")
cat("============================================================\n")


cat(
  "Project directory: ",
  getwd(),
  "\n",
  sep = ""
)



# ============================================================
# 2. LOAD SAVED THOMPSON WORKSPACE
# ============================================================

workspace_file <-
  file.path(
    project_dir,
    "Thompson_human_CHRNA3_CHRNB2_workspace.RData"
  )


if (!file.exists(workspace_file)) {
  
  stop(
    paste0(
      "Saved Thompson workspace not found:\n",
      workspace_file
    )
  )
}


cat("\nLoading saved Thompson workspace...\n")


loaded_objects <-
  load(
    workspace_file,
    envir = .GlobalEnv
  )


cat(
  "Workspace loaded successfully.\n"
)


cat(
  "Number of objects loaded: ",
  length(loaded_objects),
  "\n",
  sep = ""
)



# ------------------------------------------------------------
# Confirm that the key objects were actually contained in the
# saved workspace.
# ------------------------------------------------------------

required_objects <- c(
  "meta",
  "ratio_df",
  "donor_results",
  "paired_results"
)


missing_objects <-
  required_objects[
    !vapply(
      required_objects,
      function(x) {
        
        exists(
          x,
          envir = .GlobalEnv,
          inherits = FALSE
        )
      },
      logical(1)
    )
  ]


if (length(missing_objects) > 0) {
  
  cat("\nObjects actually loaded from the .RData file:\n")
  
  print(
    loaded_objects
  )
  
  
  stop(
    paste0(
      "\nRequired workspace object(s) missing: ",
      paste(
        missing_objects,
        collapse = ", "
      ),
      "\n\n",
      "The list above shows exactly what was stored in the workspace."
    )
  )
}


cat("\nRequired workspace objects found:\n")

print(
  required_objects
)



# ============================================================
# 3. REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "ggplot2",
  "patchwork"
)


missing_packages <-
  required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]


if (length(missing_packages) > 0) {
  
  cat(
    "\nInstalling missing package(s): ",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
  
  
  install.packages(
    missing_packages
  )
}


library(ggplot2)
library(patchwork)



# ============================================================
# 4. OUTPUT DIRECTORIES
# ============================================================

dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  "results",
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  file.path(
    "results",
    "tables"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)



# ============================================================
# 5. HELPER FUNCTIONS
# ============================================================

fmt_int <- function(x) {
  
  format(
    x,
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}



safe_ggsave <- function(
    filename,
    plot,
    ...
) {
  
  tryCatch(
    
    {
      
      ggsave(
        filename = filename,
        plot = plot,
        ...
      )
      
      cat(
        "Saved: ",
        filename,
        "\n",
        sep = ""
      )
    },
    
    error = function(e) {
      
      extension <-
        tools::file_ext(
          filename
        )
      
      
      base_name <-
        substr(
          filename,
          1,
          nchar(filename) -
            nchar(extension) -
            1
        )
      
      
      timestamp <-
        format(
          Sys.time(),
          "%Y%m%d_%H%M%S"
        )
      
      
      alternate_filename <-
        paste0(
          base_name,
          "_NEW_",
          timestamp,
          ".",
          extension
        )
      
      
      cat(
        "\nCould not overwrite:\n",
        filename,
        "\n",
        sep = ""
      )
      
      
      cat(
        "Saving instead as:\n",
        alternate_filename,
        "\n",
        sep = ""
      )
      
      
      ggsave(
        filename = alternate_filename,
        plot = plot,
        ...
      )
    }
  )
}



safe_write_csv <- function(
    x,
    filename
) {
  
  tryCatch(
    
    {
      
      write.csv(
        x,
        filename,
        row.names = FALSE
      )
    },
    
    error = function(e) {
      
      timestamp <-
        format(
          Sys.time(),
          "%Y%m%d_%H%M%S"
        )
      
      
      alternate_filename <-
        sub(
          "\\.csv$",
          paste0(
            "_NEW_",
            timestamp,
            ".csv"
          ),
          filename
        )
      
      
      write.csv(
        x,
        alternate_filename,
        row.names = FALSE
      )
    }
  )
}



# ============================================================
# 6. PUBLICATION THEME
# ============================================================

theme_human <- theme_classic(
  base_size = 12
) +
  
  theme(
    
    axis.title =
      element_text(
        size = 11
      ),
    
    axis.text =
      element_text(
        size = 10
      ),
    
    plot.title =
      element_text(
        size = 12,
        face = "bold"
      ),
    
    plot.subtitle =
      element_text(
        size = 9
      ),
    
    legend.title =
      element_blank(),
    
    legend.text =
      element_text(
        size = 9
      ),
    
    plot.margin =
      margin(
        8,
        8,
        8,
        8
      )
  )



# ============================================================
# 7. VERIFY REQUIRED COLUMN STRUCTURE
# ============================================================

required_donor_columns <- c(
  "donor",
  "sex",
  "age",
  "GABA_MGE_n",
  "double_positive_n",
  "double_positive_percent"
)


missing_donor_columns <-
  setdiff(
    required_donor_columns,
    names(
      donor_results
    )
  )


if (length(missing_donor_columns) > 0) {
  
  stop(
    paste0(
      "donor_results is missing required column(s): ",
      paste(
        missing_donor_columns,
        collapse = ", "
      )
    )
  )
}



required_ratio_columns <- c(
  "cell",
  "donor",
  "fine_class",
  "superfine_class",
  "CHRNA3_count",
  "CHRNB2_count"
)


missing_ratio_columns <-
  setdiff(
    required_ratio_columns,
    names(
      ratio_df
    )
  )


if (length(missing_ratio_columns) > 0) {
  
  stop(
    paste0(
      "ratio_df is missing required column(s): ",
      paste(
        missing_ratio_columns,
        collapse = ", "
      )
    )
  )
}



required_paired_columns <- c(
  "donor",
  "sex",
  "age",
  "CORT_n",
  "CORT_double_positive_n",
  "CORT_double_positive_percent",
  "nonCORT_n",
  "nonCORT_double_positive_n",
  "nonCORT_double_positive_percent"
)


missing_paired_columns <-
  setdiff(
    required_paired_columns,
    names(
      paired_results
    )
  )


if (length(missing_paired_columns) > 0) {
  
  stop(
    paste0(
      "paired_results is missing required column(s): ",
      paste(
        missing_paired_columns,
        collapse = ", "
      )
    )
  )
}


cat(
  "\nRequired column structure verified successfully.\n"
)
# ============================================================
# 5. BASIC INHIBITORY DATASET
# ============================================================

inhibitory_classes <- c("GABA.MGE", "GABA.PENK", "GABA.CGE", "GABA.LAMP5")

inh_meta <- meta[
  as.character(meta$fine.cell.class) %in% inhibitory_classes,
]

inh_subtype_meta <- inh_meta[!is.na(inh_meta$superfine.cell.class), ]

n_total_human_donors <- length(unique(as.character(donor_results$donor)))

overall_total <- nrow(inh_meta)
overall_double <- nrow(ratio_df)
overall_double_percent <- 100 * overall_double / overall_total

cat("\nBASIC DATA CHECKS\n")
cat("-----------------\n")
cat("Inhibitory nuclei: ", overall_total, "\n", sep = "")
cat("Double-positive nuclei: ", overall_double, "\n", sep = "")
cat("Human donors: ", n_total_human_donors, "\n", sep = "")

if (overall_total != 10747L) warning("Expected 10,747 inhibitory nuclei.")
if (overall_double != 66L) warning("Expected 66 double-positive inhibitory nuclei.")
if (n_total_human_donors != 10L) warning("Expected 10 human donors.")


# ============================================================
# 6. FIGURE 4 DATA AND STATISTICS
# ============================================================

figure4A_data <- data.frame(
  category = "CHRNA3+\nCHRNB2+",
  n_double = overall_double,
  n_total = overall_total,
  percent = overall_double_percent
)

class_data <- data.frame(
  inhibitory_class = inhibitory_classes,
  n_total = vapply(
    inhibitory_classes,
    function(x) sum(as.character(inh_meta$fine.cell.class) == x, na.rm = TRUE),
    numeric(1)
  ),
  n_double = vapply(
    inhibitory_classes,
    function(x) sum(as.character(ratio_df$fine_class) == x, na.rm = TRUE),
    numeric(1)
  )
)
class_data$n_total <- as.integer(class_data$n_total)
class_data$n_double <- as.integer(class_data$n_double)
class_data$percent <- 100 * class_data$n_double / class_data$n_total
class_data$inhibitory_class <- factor(
  class_data$inhibitory_class,
  levels = inhibitory_classes
)

class_fisher_matrix <- cbind(
  Double_positive = class_data$n_double,
  Not_double_positive = class_data$n_total - class_data$n_double
)
rownames(class_fisher_matrix) <- as.character(class_data$inhibitory_class)
class_fisher <- fisher.test(class_fisher_matrix)

# Previously verified GABA.MGE raw-count expression categories.
mge_category_data <- data.frame(
  expression_category = c(
    "Neither", "CHRNB2\nonly", "CHRNA3\nonly", "Double-\npositive"
  ),
  n = c(2576L, 1089L, 64L, 47L)
)
mge_category_total <- sum(mge_category_data$n)
if (mge_category_total != 3776L) stop("GABA.MGE counts do not sum to 3,776.")
mge_category_data$percent <- 100 * mge_category_data$n / mge_category_total
mge_category_data$expression_category <- factor(
  mge_category_data$expression_category,
  levels = c("Neither", "CHRNB2\nonly", "CHRNA3\nonly", "Double-\npositive")
)

donor_plot_data <- as.data.frame(donor_results)
donor_plot_data$double_positive_percent <- as.numeric(donor_plot_data$double_positive_percent)
donor_plot_data$double_positive_n <- as.integer(donor_plot_data$double_positive_n)
donor_plot_data <- donor_plot_data[order(donor_plot_data$double_positive_percent), ]
donor_plot_data$donor <- factor(donor_plot_data$donor, levels = donor_plot_data$donor)
human_donor_median <- median(donor_plot_data$double_positive_percent)
human_donor_min <- min(donor_plot_data$double_positive_percent)
human_donor_max <- max(donor_plot_data$double_positive_percent)


# ============================================================
# 7. FIGURE 4
# ============================================================

fig4A <- ggplot(figure4A_data, aes(x = category, y = percent)) +
  geom_col(width = 0.58, fill = "grey45", color = "black", linewidth = 0.6) +
  geom_text(
    aes(label = sprintf("%.3f%%\n(%s/%s)", percent, fmt_int(n_double), fmt_int(n_total))),
    vjust = -0.45,
    size = 3.8
  ) +
  scale_y_continuous(limits = c(0, 0.85), expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = NULL,
    y = "Double-positive nuclei (%)",
    title = "Overall inhibitory population",
    subtitle = "10,747 GABAergic inhibitory nuclei"
  ) +
  theme_human +
  theme(axis.text.x = element_text(face = "bold"))

fig4B <- ggplot(class_data, aes(x = inhibitory_class, y = percent)) +
  geom_col(width = 0.63, fill = "grey55", color = "black", linewidth = 0.6) +
  geom_text(
    aes(label = sprintf("%.3f%%\n(%s/%s)", percent, fmt_int(n_double), fmt_int(n_total))),
    vjust = -0.35,
    size = 3.0
  ) +
  scale_y_continuous(limits = c(0, 1.55), expand = expansion(mult = c(0, 0.03))) +
  labs(
    x = "Broad inhibitory class",
    y = "CHRNA3+CHRNB2+ nuclei (%)",
    title = "Broad inhibitory classes",
    subtitle = paste0(
      "Fisher exact P = ",
      format(class_fisher$p.value, scientific = TRUE, digits = 4)
    )
  ) +
  theme_human +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

fig4C <- ggplot(mge_category_data, aes(x = expression_category, y = percent)) +
  geom_col(width = 0.65, fill = "grey65", color = "black", linewidth = 0.6) +
  geom_text(
    aes(label = sprintf("%.2f%%\n(n = %s)", percent, fmt_int(n))),
    vjust = -0.35,
    size = 3.0
  ) +
  scale_y_continuous(limits = c(0, 76), expand = expansion(mult = c(0, 0.02))) +
  labs(
    x = "Expression category",
    y = "GABA.MGE nuclei (%)",
    title = "Expression categories within GABA.MGE",
    subtitle = "n = 3,776 nuclei"
  ) +
  theme_human

fig4D <- ggplot(donor_plot_data, aes(x = donor, y = double_positive_percent)) +
  geom_hline(yintercept = human_donor_median, linetype = "dashed", linewidth = 0.6) +
  geom_point(shape = 21, size = 3.2, stroke = 0.9, fill = "white") +
  geom_text(aes(label = double_positive_n), nudge_y = 0.18, size = 2.7) +
  scale_y_continuous(
    limits = c(0, max(3.5, human_donor_max + 0.5)),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Human donor",
    y = "GABA.MGE double-positive nuclei (%)",
    title = "Donor-level recurrence",
    subtitle = sprintf(
      "Median = %.2f%%; range = %.2f-%.2f%%\nNumbers above points indicate double-positive nuclei",
      human_donor_median, human_donor_min, human_donor_max
    )
  ) +
  theme_human +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figure4_human <- (fig4A | fig4B) / (fig4C | fig4D) +
  plot_annotation(
    title = "CHRNA3 and CHRNB2 transcript co-detection in human hippocampal inhibitory neurons",
    subtitle = "Thompson et al. anterior human hippocampus single-nucleus RNA-sequencing dataset",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.tag = element_text(size = 14, face = "bold")
    )
  )

print(figure4_human)

safe_ggsave(
  file.path("figures", "Figure4_Thompson_human_CHRNA3_CHRNB2.png"),
  figure4_human, width = 10, height = 8, units = "in", dpi = 300
)
safe_ggsave(
  file.path("figures", "Figure4_Thompson_human_CHRNA3_CHRNB2.tiff"),
  figure4_human, width = 10, height = 8, units = "in", dpi = 600,
  compression = "lzw"
)
safe_ggsave(
  file.path("figures", "Figure4_Thompson_human_CHRNA3_CHRNB2.pdf"),
  figure4_human, width = 10, height = 8, units = "in"
)

safe_write_csv(figure4A_data, file.path("results", "tables", "Figure4A_overall_human_double_positive.csv"))
safe_write_csv(class_data, file.path("results", "tables", "Figure4B_human_inhibitory_classes.csv"))
safe_write_csv(mge_category_data, file.path("results", "tables", "Figure4C_GABA_MGE_expression_categories.csv"))
safe_write_csv(donor_plot_data, file.path("results", "tables", "Figure4D_GABA_MGE_donor_results.csv"))


# ============================================================
# 8. SUPERFINE INHIBITORY-SUBTYPE ANALYSIS
# ============================================================

subtype_denominators <- aggregate(
  list(total_n = rep(1L, nrow(inh_subtype_meta))),
  by = list(
    fine_class = as.character(inh_subtype_meta$fine.cell.class),
    superfine_class = as.character(inh_subtype_meta$superfine.cell.class)
  ),
  FUN = sum
)

double_subtypes <- aggregate(
  list(double_positive_n = rep(1L, nrow(ratio_df))),
  by = list(
    fine_class = as.character(ratio_df$fine_class),
    superfine_class = as.character(ratio_df$superfine_class)
  ),
  FUN = sum
)

human_subtype_results <- merge(
  subtype_denominators,
  double_subtypes,
  by = c("fine_class", "superfine_class"),
  all.x = TRUE
)
human_subtype_results$double_positive_n[is.na(human_subtype_results$double_positive_n)] <- 0L
human_subtype_results$total_n <- as.integer(human_subtype_results$total_n)
human_subtype_results$double_positive_n <- as.integer(human_subtype_results$double_positive_n)
human_subtype_results$double_positive_percent <-
  100 * human_subtype_results$double_positive_n / human_subtype_results$total_n

double_donor_counts <- aggregate(
  donor ~ fine_class + superfine_class,
  data = ratio_df,
  FUN = function(x) length(unique(x))
)
names(double_donor_counts)[3] <- "donors_with_double_positive"

human_subtype_results <- merge(
  human_subtype_results,
  double_donor_counts,
  by = c("fine_class", "superfine_class"),
  all.x = TRUE
)
human_subtype_results$donors_with_double_positive[
  is.na(human_subtype_results$donors_with_double_positive)
] <- 0L
human_subtype_results$donors_with_double_positive <-
  as.integer(human_subtype_results$donors_with_double_positive)

human_subtype_results <- human_subtype_results[
  order(
    -human_subtype_results$double_positive_percent,
    -human_subtype_results$double_positive_n
  ),
]
rownames(human_subtype_results) <- NULL

cat("\nALL HUMAN INHIBITORY SUPERFINE SUBTYPES\n")
cat("---------------------------------------\n")
print(human_subtype_results, row.names = FALSE)


# ============================================================
# 9. GABA.MGE HETEROGENEITY AND CORT ENRICHMENT
# ============================================================

mge_fine <- human_subtype_results[
  human_subtype_results$fine_class == "GABA.MGE",
]

mge_matrix <- cbind(
  Double_positive = mge_fine$double_positive_n,
  Not_double_positive = mge_fine$total_n - mge_fine$double_positive_n
)
rownames(mge_matrix) <- mge_fine$superfine_class
mge_fisher <- fisher.test(mge_matrix)

cort_row <- human_subtype_results[
  human_subtype_results$fine_class == "GABA.MGE" &
    human_subtype_results$superfine_class == "CORT",
]
if (nrow(cort_row) != 1L) stop("Expected exactly one GABA.MGE/CORT row.")

pooled_cort_double <- cort_row$double_positive_n
pooled_cort_total <- cort_row$total_n
pooled_noncort_double <- sum(mge_fine$double_positive_n) - pooled_cort_double
pooled_noncort_total <- sum(mge_fine$total_n) - pooled_cort_total

cort_matrix <- matrix(
  c(
    pooled_cort_double,
    pooled_cort_total - pooled_cort_double,
    pooled_noncort_double,
    pooled_noncort_total - pooled_noncort_double
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    c("CORT", "Other_GABA_MGE"),
    c("Double_positive", "Not_double_positive")
  )
)
cort_fisher <- fisher.test(cort_matrix)

# Exploratory LAMP5 and CGE heterogeneity tests.
lamp5_results <- human_subtype_results[human_subtype_results$fine_class == "GABA.LAMP5", ]
lamp5_matrix <- cbind(
  Double_positive = lamp5_results$double_positive_n,
  Not_double_positive = lamp5_results$total_n - lamp5_results$double_positive_n
)
rownames(lamp5_matrix) <- lamp5_results$superfine_class
lamp5_fisher <- fisher.test(lamp5_matrix)

cge_results <- human_subtype_results[human_subtype_results$fine_class == "GABA.CGE", ]
cge_matrix <- cbind(
  Double_positive = cge_results$double_positive_n,
  Not_double_positive = cge_results$total_n - cge_results$double_positive_n
)
rownames(cge_matrix) <- cge_results$superfine_class
cge_fisher <- fisher.test(cge_matrix)

cat("\nGABA.MGE superfine heterogeneity Fisher P = ",
    format(mge_fisher$p.value, scientific = TRUE, digits = 6), "\n", sep = "")
cat("CORT pooled frequency = ", sprintf("%.3f%%", 100 * pooled_cort_double / pooled_cort_total), "\n", sep = "")
cat("Other GABA.MGE pooled frequency = ", sprintf("%.3f%%", 100 * pooled_noncort_double / pooled_noncort_total), "\n", sep = "")
cat("CORT pooled OR = ", sprintf("%.3f", unname(cort_fisher$estimate)), "\n", sep = "")
cat("CORT pooled 95% CI = ", sprintf("%.3f to %.3f", cort_fisher$conf.int[1], cort_fisher$conf.int[2]), "\n", sep = "")
cat("CORT pooled Fisher P = ", format(cort_fisher$p.value, scientific = TRUE, digits = 6), "\n", sep = "")
cat("LAMP5 subtype Fisher P = ", format(lamp5_fisher$p.value, scientific = TRUE, digits = 6), "\n", sep = "")
cat("CGE subtype Fisher P = ", format(cge_fisher$p.value, scientific = TRUE, digits = 6), "\n", sep = "")


# ============================================================
# 10. DONOR-AWARE CORT ANALYSIS
# ============================================================

cort_paired <- as.data.frame(paired_results)
cort_paired$CORT_double_positive_percent <- as.numeric(cort_paired$CORT_double_positive_percent)
cort_paired$nonCORT_double_positive_percent <- as.numeric(cort_paired$nonCORT_double_positive_percent)
cort_paired$paired_difference_pp <-
  cort_paired$CORT_double_positive_percent -
  cort_paired$nonCORT_double_positive_percent
cort_paired$fold_difference <- ifelse(
  cort_paired$nonCORT_double_positive_percent > 0,
  cort_paired$CORT_double_positive_percent /
    cort_paired$nonCORT_double_positive_percent,
  NA_real_
)
cort_paired <- cort_paired[order(cort_paired$donor), ]
rownames(cort_paired) <- NULL

n_donors <- nrow(cort_paired)
n_donors_cort_positive <- sum(cort_paired$CORT_double_positive_n > 0)
n_cort_higher <- sum(cort_paired$paired_difference_pp > 0)
n_noncort_higher <- sum(cort_paired$paired_difference_pp < 0)
n_ties <- sum(cort_paired$paired_difference_pp == 0)
mean_cort_percent <- mean(cort_paired$CORT_double_positive_percent)
median_cort_percent <- median(cort_paired$CORT_double_positive_percent)
mean_noncort_percent <- mean(cort_paired$nonCORT_double_positive_percent)
median_noncort_percent <- median(cort_paired$nonCORT_double_positive_percent)
paired_differences <- cort_paired$paired_difference_pp
observed_mean_difference <- mean(paired_differences)
median_difference_pp <- median(paired_differences)

# Exact paired sign-flip permutation test.
sign_grid <- expand.grid(rep(list(c(-1, 1)), length(paired_differences)))
sign_matrix <- as.matrix(sign_grid)
permuted_mean_differences <- apply(
  sign_matrix,
  1,
  function(signs) mean(paired_differences * signs)
)
paired_permutation_p <- mean(
  abs(permuted_mean_differences) >=
    abs(observed_mean_difference) - 1e-12
)

# Direction-only exact sign test.
non_tied_differences <- paired_differences[paired_differences != 0]
sign_test_successes <- sum(non_tied_differences > 0)
sign_test <- binom.test(
  x = sign_test_successes,
  n = length(non_tied_differences),
  p = 0.5,
  alternative = "two.sided"
)

# Donor-stratified CMH sensitivity analysis.
cmh_array <- array(
  0,
  dim = c(2, 2, n_donors),
  dimnames = list(
    Subtype = c("CORT", "nonCORT"),
    Expression = c("Double_positive", "Not_double_positive"),
    Donor = as.character(cort_paired$donor)
  )
)

for (i in seq_len(n_donors)) {
  cmh_array["CORT", "Double_positive", i] <-
    cort_paired$CORT_double_positive_n[i]
  cmh_array["CORT", "Not_double_positive", i] <-
    cort_paired$CORT_n[i] - cort_paired$CORT_double_positive_n[i]
  cmh_array["nonCORT", "Double_positive", i] <-
    cort_paired$nonCORT_double_positive_n[i]
  cmh_array["nonCORT", "Not_double_positive", i] <-
    cort_paired$nonCORT_n[i] - cort_paired$nonCORT_double_positive_n[i]
}

cmh_result <- mantelhaen.test(cmh_array, correct = FALSE)

cort_donor_aware_summary <- data.frame(
  analysis = c(
    "Exact paired sign-flip permutation",
    "Exact sign test",
    "Donor-stratified CMH",
    "Pooled Fisher exact"
  ),
  statistic_description = c(
    "Mean paired CORT minus non-CORT percentage",
    "Donors with CORT frequency greater than non-CORT",
    "Common odds ratio stratified by donor",
    "Pooled cell-level odds ratio"
  ),
  statistic_value = c(
    observed_mean_difference,
    sign_test_successes,
    unname(cmh_result$estimate),
    unname(cort_fisher$estimate)
  ),
  p_value = c(
    paired_permutation_p,
    sign_test$p.value,
    cmh_result$p.value,
    cort_fisher$p.value
  )
)

cat("\nDONOR-AWARE CORT RESULTS\n")
cat("------------------------\n")
cat("CORT-positive donors: ", n_donors_cort_positive, "/", n_donors, "\n", sep = "")
cat("CORT > non-CORT in ", n_cort_higher, "/", n_donors, " donors\n", sep = "")
cat("Mean donor CORT frequency: ", sprintf("%.3f%%", mean_cort_percent), "\n", sep = "")
cat("Mean donor non-CORT frequency: ", sprintf("%.3f%%", mean_noncort_percent), "\n", sep = "")
cat("Mean paired difference: ", sprintf("%.3f percentage points", observed_mean_difference), "\n", sep = "")
cat("Exact paired permutation P = ", sprintf("%.6f", paired_permutation_p), "\n", sep = "")
cat("Exact sign-test P = ", sprintf("%.6f", sign_test$p.value), "\n", sep = "")
cat("CMH common OR = ", sprintf("%.3f", unname(cmh_result$estimate)), "\n", sep = "")
cat("CMH P = ", format(cmh_result$p.value, scientific = TRUE, digits = 6), "\n", sep = "")


# ============================================================
# 11. SAVE STATISTICAL TABLES
# ============================================================

safe_write_csv(
  human_subtype_results,
  file.path("results", "tables", "Human_superfine_CHRNA3_CHRNB2_results.csv")
)
safe_write_csv(
  lamp5_results,
  file.path("results", "tables", "Human_LAMP5_superfine_results.csv")
)
safe_write_csv(
  cge_results,
  file.path("results", "tables", "Human_CGE_superfine_results.csv")
)
safe_write_csv(
  cort_paired,
  file.path("results", "tables", "CORT_donor_paired_results.csv")
)
safe_write_csv(
  cort_donor_aware_summary,
  file.path("results", "tables", "CORT_donor_aware_statistical_summary.csv")
)
safe_write_csv(
  data.frame(permuted_mean_difference = permuted_mean_differences),
  file.path("results", "tables", "CORT_exact_signflip_permutation_distribution.csv")
)


# ============================================================
# 12. FIGURE 5 DATA
# ============================================================

fig5A_data <- human_subtype_results[
  order(
    -human_subtype_results$double_positive_percent,
    -human_subtype_results$double_positive_n
  ),
]

# Reader-friendly labels retaining parent family where useful.
fig5A_data$display_subtype <- fig5A_data$superfine_class

idx_mge <- fig5A_data$fine_class == "GABA.MGE"
fig5A_data$display_subtype[idx_mge] <-
  paste0(fig5A_data$superfine_class[idx_mge], " (MGE)")

idx_cge <- fig5A_data$fine_class == "GABA.CGE"
fig5A_data$display_subtype[idx_cge] <-
  paste0(fig5A_data$superfine_class[idx_cge], " (CGE)")

idx_cxcl14 <-
  fig5A_data$fine_class == "GABA.LAMP5" &
  fig5A_data$superfine_class == "CXCL14"
fig5A_data$display_subtype[idx_cxcl14] <- "CXCL14 (LAMP5)"

fig5A_data$display_subtype <- factor(
  fig5A_data$display_subtype,
  levels = rev(fig5A_data$display_subtype)
)

fig5A_data$highlight <- ifelse(
  fig5A_data$superfine_class == "CORT",
  "CORT",
  "Other"
)

# IMPORTANT FINAL LABEL CHANGE:
# report only the number of donors containing >=1 double-positive
# nucleus for the subtype. Do not display /10, because not every
# subtype necessarily has recovered nuclei in every donor.
fig5A_data$plot_label <- sprintf(
  "%.2f%% (%s/%s; %d donor%s)",
  fig5A_data$double_positive_percent,
  fmt_int(fig5A_data$double_positive_n),
  fmt_int(fig5A_data$total_n),
  fig5A_data$donors_with_double_positive,
  ifelse(fig5A_data$donors_with_double_positive == 1, "", "s")
)

fig5A_xmax <- max(
  10,
  ceiling(max(fig5A_data$double_positive_percent) * 1.05)
)

fig5B_data <- rbind(
  data.frame(
    donor = as.character(cort_paired$donor),
    group = "non-CORT GABA.MGE",
    percent = cort_paired$nonCORT_double_positive_percent
  ),
  data.frame(
    donor = as.character(cort_paired$donor),
    group = "CORT",
    percent = cort_paired$CORT_double_positive_percent
  )
)
fig5B_data$group <- factor(
  fig5B_data$group,
  levels = c("non-CORT GABA.MGE", "CORT")
)
fig5B_summary <- aggregate(percent ~ group, data = fig5B_data, FUN = mean)
names(fig5B_summary)[2] <- "mean_percent"
fig5B_summary$group <- factor(
  fig5B_summary$group,
  levels = c("non-CORT GABA.MGE", "CORT")
)

fig5C_data <- cort_paired[
  order(cort_paired$CORT_double_positive_percent),
]
fig5C_data$donor <- factor(fig5C_data$donor, levels = fig5C_data$donor)
cort_donor_median <- median(fig5C_data$CORT_double_positive_percent)
cort_donor_min <- min(fig5C_data$CORT_double_positive_percent)
cort_donor_max <- max(fig5C_data$CORT_double_positive_percent)

bottom_panel_max <- max(
  fig5B_data$percent,
  fig5C_data$CORT_double_positive_percent,
  na.rm = TRUE
)
bottom_ymax <- ceiling(bottom_panel_max * 1.16 / 5) * 5
bottom_ymax <- max(bottom_ymax, 35)


# ============================================================
# 13. FIGURE 5
# ============================================================

fig5A <- ggplot(
  fig5A_data,
  aes(
    x = double_positive_percent,
    y = display_subtype,
    fill = highlight
  )
) +
  geom_col(width = 0.68, color = "black", linewidth = 0.5) +
  geom_text(aes(label = plot_label), hjust = -0.04, size = 2.85) +
  scale_fill_manual(values = c("CORT" = "grey25", "Other" = "grey70")) +
  scale_x_continuous(
    limits = c(0, fig5A_xmax),
    breaks = pretty(c(0, fig5A_xmax), n = 5),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "CHRNA3+CHRNB2+ nuclei (%)",
    y = "Human inhibitory superfine subtype",
    title = "Fine inhibitory subtypes",
    subtitle = paste0(
      "Labels: double-positive frequency (positive/total nuclei; ",
      "donors with >=1 double-positive nucleus)"
    )
  ) +
  theme_human +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "bold", size = 9),
    plot.margin = margin(5, 185, 6, 8)
  )

fig5B <- ggplot(
  fig5B_data,
  aes(x = group, y = percent, group = donor)
) +
  geom_line(color = "grey60", linewidth = 0.55) +
  geom_point(aes(fill = group), shape = 21, size = 3.0, stroke = 0.8) +
  geom_point(
    data = fig5B_summary,
    aes(x = group, y = mean_percent),
    inherit.aes = FALSE,
    shape = 95,
    size = 11
  ) +
  annotate(
    "text",
    x = 1.5,
    y = bottom_ymax - 2,
    label = sprintf("Exact paired permutation P = %.4f", paired_permutation_p),
    size = 3.2
  ) +
  scale_fill_manual(
    values = c("non-CORT GABA.MGE" = "white", "CORT" = "grey35")
  ) +
  scale_y_continuous(
    limits = c(0, bottom_ymax),
    breaks = seq(0, bottom_ymax, by = 10),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Double-positive nuclei (%)",
    title = "Matched donor CORT enrichment"
  ) +
  theme_human +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 9),
    plot.margin = margin(8, 20, 8, 8)
  )

fig5C <- ggplot(
  fig5C_data,
  aes(x = donor, y = CORT_double_positive_percent)
) +
  geom_hline(
    yintercept = cort_donor_median,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_point(shape = 21, size = 3.1, stroke = 0.85, fill = "grey35") +
  geom_text(aes(label = CORT_double_positive_n), nudge_y = 0.9, size = 2.8) +
  scale_y_continuous(
    limits = c(0, bottom_ymax),
    breaks = seq(0, bottom_ymax, by = 10),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Human donor",
    y = "CORT double-positive nuclei (%)",
    title = "CORT recurrence across donors",
    subtitle = sprintf(
      "Double-positive CORT nuclei detected in %d/%d donors; dashed line = median %.2f%%",
      n_donors_cort_positive,
      n_donors,
      cort_donor_median
    )
  ) +
  theme_human +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(8, 8, 8, 20)
  )

figure5_human <- fig5A / (fig5B | fig5C) +
  plot_layout(heights = c(1.15, 1)) +
  plot_annotation(
    title = "CHRNA3 and CHRNB2 transcript co-detection across human hippocampal inhibitory subtypes",
    subtitle = "Thompson et al. anterior human hippocampus single-nucleus RNA-sequencing dataset",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 15, face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.tag = element_text(size = 14, face = "bold"),
      plot.margin = margin(6, 8, 6, 8)
    )
  )

print(figure5_human)

safe_ggsave(
  file.path("figures", "Figure5_Thompson_human_fine_subtype_CORT.png"),
  figure5_human, width = 12, height = 9.5, units = "in", dpi = 300
)
safe_ggsave(
  file.path("figures", "Figure5_Thompson_human_fine_subtype_CORT.tiff"),
  figure5_human, width = 12, height = 9.5, units = "in", dpi = 600,
  compression = "lzw"
)
safe_ggsave(
  file.path("figures", "Figure5_Thompson_human_fine_subtype_CORT.pdf"),
  figure5_human, width = 12, height = 9.5, units = "in"
)

safe_write_csv(
  fig5A_data,
  file.path("results", "tables", "Figure5A_human_superfine_subtypes.csv")
)
safe_write_csv(
  fig5B_data,
  file.path("results", "tables", "Figure5B_CORT_paired_long_format.csv")
)
safe_write_csv(
  fig5B_summary,
  file.path("results", "tables", "Figure5B_CORT_group_means.csv")
)
safe_write_csv(
  fig5C_data,
  file.path("results", "tables", "Figure5C_CORT_donor_recurrence.csv")
)


# ============================================================
# 14. SAVE ANALYSIS OBJECTS AND SESSION INFO
# ============================================================

save(
  human_subtype_results,
  class_fisher,
  mge_matrix,
  mge_fisher,
  cort_matrix,
  cort_fisher,
  lamp5_results,
  lamp5_fisher,
  cge_results,
  cge_fisher,
  cort_paired,
  paired_permutation_p,
  permuted_mean_differences,
  sign_test,
  cmh_result,
  cort_donor_aware_summary,
  figure4_human,
  figure5_human,
  file = file.path(
    "results",
    "Thompson_CHRNA3_CHRNB2_fine_subtype_analysis.RData"
  )
)

capture.output(
  sessionInfo(),
  file = file.path("results", "Thompson_CHRNA3_CHRNB2_sessionInfo.txt")
)


# ============================================================
# 15. FINAL REPORT
# ============================================================

cat("\n============================================================\n")
cat("THOMPSON HUMAN ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Overall double-positive: ", overall_double, "/", overall_total,
    " (", sprintf("%.3f%%", overall_double_percent), ")\n", sep = "")
cat("Broad-class Fisher P = ",
    format(class_fisher$p.value, scientific = TRUE, digits = 6), "\n", sep = "")
cat("GABA.MGE superfine Fisher P = ",
    format(mge_fisher$p.value, scientific = TRUE, digits = 6), "\n", sep = "")
cat("CORT: ", pooled_cort_double, "/", pooled_cort_total,
    " (", sprintf("%.3f%%", 100 * pooled_cort_double / pooled_cort_total), ")\n", sep = "")
cat("Other GABA.MGE: ", pooled_noncort_double, "/", pooled_noncort_total,
    " (", sprintf("%.3f%%", 100 * pooled_noncort_double / pooled_noncort_total), ")\n", sep = "")
cat("CORT pooled OR = ", sprintf("%.3f", unname(cort_fisher$estimate)),
    " (95% CI ", sprintf("%.3f-%.3f", cort_fisher$conf.int[1], cort_fisher$conf.int[2]), ")\n", sep = "")
cat("CORT pooled Fisher P = ",
    format(cort_fisher$p.value, scientific = TRUE, digits = 6), "\n", sep = "")
cat("CORT double-positive nuclei in ", n_donors_cort_positive, "/", n_donors,
    " donors\n", sep = "")
cat("Mean donor CORT frequency = ", sprintf("%.3f%%", mean_cort_percent), "\n", sep = "")
cat("Mean donor non-CORT frequency = ", sprintf("%.3f%%", mean_noncort_percent), "\n", sep = "")
cat("Mean paired difference = ", sprintf("%.3f percentage points", observed_mean_difference), "\n", sep = "")
cat("Exact paired permutation P = ", sprintf("%.6f", paired_permutation_p), "\n", sep = "")
cat("Exact sign-test P = ", sprintf("%.6f", sign_test$p.value), "\n", sep = "")
cat("CMH common OR = ", sprintf("%.3f", unname(cmh_result$estimate)), "\n", sep = "")
cat("CMH P = ", format(cmh_result$p.value, scientific = TRUE, digits = 6), "\n", sep = "")
cat("LAMP5 subtype Fisher P = ",
    format(lamp5_fisher$p.value, scientific = TRUE, digits = 6), "\n", sep = "")
cat("CGE subtype Fisher P = ",
    format(cge_fisher$p.value, scientific = TRUE, digits = 6), "\n", sep = "")
cat("\nSCRIPT FINISHED SUCCESSFULLY\n")
cat("============================================================\n")