#!/usr/bin/env Rscript
#################### Improved Admixture Plot R script ####################
### Usage: Rscript plot_admix.R [admixture_dir] [output_prefix] [cv_file] [fam_file]
### Example: Rscript plot_admix.R admixture_results/ my_study cross_validation.txt hapmap3.fam

options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Load libraries
library(tidyverse)
if (!require("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr", lib = personal_lib, dependencies = TRUE)
  library(ggpubr, lib.loc = personal_lib)
}
if (!require("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer", lib = personal_lib)
  library(RColorBrewer, lib.loc = personal_lib)
}


# Parse command line arguments
cmd_args <- commandArgs(trailingOnly = TRUE)
if (length(cmd_args) < 2) {
  stop("Usage: Rscript plot_admix.R [admixture_dir] [output_prefix] [cv_file] [fam_file]
  
Arguments:
  admixture_dir:  Directory containing .Q files (default: current directory)
  output_prefix:  Prefix for output files (required)
  cv_file:       Cross-validation file (default: cross_validation.txt) 
  fam_file:      .fam file for sample info (optional, for ordering)")
}

# Set parameters with defaults
admixture_dir <- ifelse(length(cmd_args) >= 1, cmd_args[1], ".")
out.name <- cmd_args[2]
cv_file <- ifelse(length(cmd_args) >= 3, cmd_args[3], "cross_validation.txt")
fam_file <- ifelse(length(cmd_args) >= 4, cmd_args[4], NULL)

# Ensure directory ends with /
if (!endsWith(admixture_dir, "/")) {
  admixture_dir <- paste0(admixture_dir, "/")
}

# Set up logging
log_file <- paste0(out.name, "_plot_admix.log")
sink(log_file, append = FALSE, split = TRUE)

cat("=== ADMIXTURE Plotting Script ===\n")
cat("Admixture directory:", admixture_dir, "\n")
cat("Output prefix:", out.name, "\n")
cat("CV file:", cv_file, "\n")
if (!is.null(fam_file)) cat("FAM file:", fam_file, "\n")

#######################
# 1. CROSS VALIDATION #
#######################

if (file.exists(cv_file)) {
  cat("\n--- Processing Cross Validation ---\n")
  
  cv <- read.delim(cv_file, header = FALSE)
  cv <- cv %>% 
    separate(V1, into = c("stuff", "stuff2", "K", 'CV'), sep = " ") %>%  
    mutate(
      K = as.numeric(str_extract(K, "\\d+")),  # Extract just the number
      CV = as.numeric(CV)
    ) %>%
    filter(K >= 2)  # Filter out K=1
  
  bestK <- cv %>% slice_min(CV) %>% pull(K)
  
  cv.p <- ggplot(cv, aes(x = K, y = CV)) +
    geom_line(col = "gray20", size = 1) +
    geom_point(size = 4, col = "gray20") +
    geom_point(data = cv %>% filter(K == bestK), 
               size = 4, col = "red") +
    labs(title = paste("Cross validation: best K =", bestK), 
         x = "K (number of clusters)",
         y = "Cross validation error") +
    theme_bw(base_size = 14) +
    scale_x_continuous(breaks = cv$K)
  
  ggsave(cv.p, filename = paste0(out.name, "_cross_validation.png"), 
         width = 10, height = 6, bg = "white", dpi = 300)
  
  cat("Best number of clusters (K) based on cross validation:", bestK, "\n")
  cat("CV plot saved as:", paste0(out.name, "_cross_validation.png"), "\n")
} else {
  cat("Warning: CV file not found:", cv_file, "\n")
  bestK <- NULL
}

###################
# 2. ADMIXTURE Q  #
###################

cat("\n--- Processing ADMIXTURE Q files ---\n")

# Find all Q files in the directory
q_files <- list.files(path = admixture_dir, pattern = '\\.Q$', full.names = TRUE)

if (length(q_files) == 0) {
  stop("No .Q files found in directory: ", admixture_dir)
}

cat("Found", length(q_files), "Q files:\n")
for (f in q_files) cat("  ", f, "\n")

# Extract K values from filenames dynamically
get_k_from_filename <- function(filename) {
  # Try different patterns to extract K value
  basename_file <- basename(filename)
  
  # Pattern 1: filename.K.Q (e.g., hapmap3.3.Q)
  k_match <- str_extract(basename_file, "\\.(\\d+)\\.Q$")
  if (!is.na(k_match)) {
    return(as.numeric(str_extract(k_match, "\\d+")))
  }
  
  # Pattern 2: filename_K.Q (e.g., hapmap3_3.Q)  
  k_match <- str_extract(basename_file, "_(\\d+)\\.Q$")
  if (!is.na(k_match)) {
    return(as.numeric(str_extract(k_match, "\\d+")))
  }
  
  # Pattern 3: filenameK.Q (e.g., hapmap3K3.Q)
  k_match <- str_extract(basename_file, "K(\\d+)\\.Q$")
  if (!is.na(k_match)) {
    return(as.numeric(str_extract(k_match, "\\d+")))
  }
  
  # If no pattern matches, return NA
  return(NA)
}

# Process each Q file
q_data_list <- list()
plot_list <- list()

# Load sample information if fam file provided
sample_info <- NULL
if (!is.null(fam_file) && file.exists(fam_file)) {
  sample_info <- read.table(fam_file, col.names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))
  cat("Loaded sample info for", nrow(sample_info), "individuals\n")
}

for (q_file in q_files) {
  k_val <- get_k_from_filename(q_file)
  
  if (is.na(k_val)) {
    cat("Warning: Could not extract K value from filename:", q_file, "\n")
    next
  }
  
  # Skip K=1
  if (k_val == 1) {
    cat("Skipping K=1 from", q_file, "\n")
    next
  }
  
  cat("Processing K =", k_val, "from", q_file, "\n")
  
  # Read Q matrix
  q_matrix <- read.table(q_file, header = FALSE)
  
  # Add sample IDs
  q_matrix$SampleID <- 1:nrow(q_matrix)
  if (!is.null(sample_info)) {
    q_matrix$IID <- sample_info$IID[1:nrow(q_matrix)]
  } else {
    q_matrix$IID <- paste0("Sample_", 1:nrow(q_matrix))
  }
  
  # Convert to long format
  value_cols <- paste0("V", 1:k_val)
  q_long <- q_matrix %>%
    pivot_longer(cols = all_of(value_cols), 
                 names_to = "Ancestry", 
                 values_to = "Proportion") %>%
    mutate(
      K = k_val,
      Ancestry = factor(Ancestry, levels = value_cols)
    )
  
  # Calculate dominant ancestry for ordering
  q_matrix_dominant <- q_matrix %>%
    select(all_of(value_cols)) %>%
    mutate(
      SampleID = 1:nrow(.),
      DominantAncestry = max.col(.[value_cols]),
      MaxProportion = apply(.[value_cols], 1, max)
    ) %>%
    arrange(DominantAncestry, desc(MaxProportion))
  
  # Reorder samples based on clustering
  q_long$SampleID <- factor(q_long$SampleID, levels = q_matrix_dominant$SampleID)
  
  # Store data
  q_data_list[[paste0("K", k_val)]] <- q_long
  
  # Create color palette
  if (k_val <= 8) {
    colors <- brewer.pal(max(3, k_val), "Set2")[1:k_val]
  } else {
    colors <- rainbow(k_val)
  }
  
  # Determine label frequency based on sample size
  n_samples <- nrow(q_matrix)
  if (n_samples <= 50) {
    label_freq <- 1  # Show every sample
  } else if (n_samples <= 200) {
    label_freq <- 5  # Show every 5th sample
  } else if (n_samples <= 500) {
    label_freq <- 10  # Show every 10th sample
  } else {
    label_freq <- 20  # Show every 20th sample
  }
  
  # Get unique samples in their clustered order
  sample_order <- unique(q_long[c("SampleID", "IID")])
  
  # Select samples to show labels for
  samples_to_show <- seq(1, nrow(sample_order), by = label_freq)
  
  # Create breaks and labels for the x-axis
  x_breaks <- sample_order$SampleID[samples_to_show]
  x_labels <- sample_order$IID[samples_to_show]
  
  # Create labeled individual plot
  p_labeled <- ggplot(q_long, aes(x = SampleID, y = Proportion, fill = Ancestry)) +
    geom_col(width = 1, size = 0) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_x_discrete(breaks = x_breaks, labels = x_labels) +
    labs(
      title = paste("K =", k_val),
      x = "Individual Sample ID",
      y = "Ancestry proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.spacing = unit(0, "lines"),
      plot.margin = margin(5, 5, 20, 5)  # Extra bottom margin for rotated labels
    )
  
  # Create clean plot for combined display (no labels)
  p_clean <- ggplot(q_long, aes(x = SampleID, y = Proportion, fill = Ancestry)) +
    geom_col(width = 1, size = 0) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    labs(
      title = paste("K =", k_val),
      x = "Individuals",
      y = "Ancestry proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      panel.spacing = unit(0, "lines"),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # Store clean plot for combined plot
  plot_list[[paste0("K", k_val)]] <- p_clean
  
  # Adjust width based on number of samples for readability
  plot_width <- max(12, min(24, n_samples * 0.1))
  
  # Save labeled individual plot
  ggsave(p_labeled, filename = paste0(out.name, "_K", k_val, "_admixture.png"), 
         width = plot_width, height = 6, bg = "white", dpi = 300)
}

# Sort by K value for combined plot (excluding K=1)
k_values <- sort(as.numeric(str_extract(names(plot_list), "\\d+")))
k_values <- k_values[k_values >= 2]  # Ensure K=1 is excluded
sorted_plots <- plot_list[paste0("K", k_values)]

cat("\nProcessed K values:", paste(k_values, collapse = ", "), "\n")

# Create combined plot
if (length(sorted_plots) > 0) {
  cat("\n--- Creating combined plot ---\n")
  
  # Determine optimal arrangement
  n_plots <- length(sorted_plots)
  if (n_plots <= 4) {
    ncol_val <- 1
    height_val <- 2.5 * n_plots
  } else {
    ncol_val <- 2
    height_val <- 2.5 * ceiling(n_plots / 2)
  }
  
  # Clean up plots for combined display
  for (i in 1:length(sorted_plots)) {
    # Remove both x-axis and y-axis titles from ALL plots (will add globally)
    sorted_plots[[i]] <- sorted_plots[[i]] + 
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank())
  }
  
  # Combine plots - Fix the error by ensuring proper plot list
  combined_plot <- do.call(ggarrange, c(sorted_plots, 
                                       list(ncol = ncol_val, 
                                            align = "v",
                                            heights = rep(1, length(sorted_plots)))))
  
  # Add overall labels - only one "Ancestry proportion" and one "Individuals"
  final_plot <- annotate_figure(combined_plot,
                               left = text_grob("Ancestry proportion", rot = 90, size = 14),
                               bottom = text_grob("Individuals", size = 14))
  
  # Save combined plot
  ggsave(final_plot, 
         filename = paste0(out.name, "_admixture_combined.png"), 
         width = 14, height = height_val, bg = "white", dpi = 300)
  
  cat("Combined plot saved as:", paste0(out.name, "_admixture_combined.png"), "\n")
  
  # Highlight best K if available
  if (!is.null(bestK) && paste0("K", bestK) %in% names(plot_list)) {
    best_k_plot <- plot_list[[paste0("K", bestK)]] +
      labs(title = paste("Best K =", bestK, "(lowest CV error)"),
           x = "Individuals", y = "Ancestry proportion") +
      theme(legend.position = "bottom")
    
    ggsave(best_k_plot, 
           filename = paste0(out.name, "_best_K", bestK, "_admixture.png"), 
           width = 14, height = 5, bg = "white", dpi = 300)
    
    cat("Best K plot saved as:", paste0(out.name, "_best_K", bestK, "_admixture.png"), "\n")
    
    # Create a version with legible individual labels using the preserved order
    best_k_data_for_labels <- q_data_list[[paste0("K", bestK)]]
    n_samples <- length(unique(best_k_data_for_labels$SampleID))
    
    # Determine label frequency based on sample size
    if (n_samples <= 50) {
      label_freq <- 1  # Show every sample
    } else if (n_samples <= 200) {
      label_freq <- 5  # Show every 5th sample
    } else if (n_samples <= 500) {
      label_freq <- 10  # Show every 10th sample
    } else {
      label_freq <- 20  # Show every 20th sample
    }
    
    # Get unique samples in their clustered order
    sample_order <- unique(best_k_data_for_labels[c("SampleID", "IID")])
    # The SampleID is already a factor with the correct clustered order
    sample_order$display_position <- as.numeric(sample_order$SampleID)
    
    # Select samples to show labels for
    samples_to_show <- seq(1, nrow(sample_order), by = label_freq)
    show_labels <- rep("", nrow(sample_order))
    show_labels[samples_to_show] <- sample_order$IID[samples_to_show]
    
    # Add the label information back to the plot data
    label_mapping <- data.frame(
      SampleID = sample_order$SampleID,
      show_label = show_labels
    )
    
    plot_data_with_labels <- best_k_data_for_labels %>%
      left_join(label_mapping, by = "SampleID")
    
    # Create breaks and labels for the x-axis
    x_breaks <- sample_order$SampleID[samples_to_show]
    x_labels <- sample_order$IID[samples_to_show]
    
    best_k_plot_labeled <- ggplot(plot_data_with_labels, aes(x = SampleID, y = Proportion, fill = Ancestry)) +
      geom_col(width = 1, size = 0) +
      scale_fill_manual(values = if (bestK <= 8) brewer.pal(max(3, bestK), "Set2")[1:bestK] else rainbow(bestK)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
      scale_x_discrete(breaks = x_breaks, labels = x_labels) +
      labs(
        title = paste("Best K =", bestK, "(lowest CV error) - Individual Labels"),
        x = "Individual Sample ID",
        y = "Ancestry proportion"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(5, 5, 20, 5)  # Extra bottom margin for rotated labels
      )
    
    # Adjust width based on number of samples for readability
    plot_width <- max(12, min(24, n_samples * 0.1))
    
    ggsave(best_k_plot_labeled, 
           filename = paste0(out.name, "_best_K", bestK, "_labeled.png"), 
           width = plot_width, height = 6, bg = "white", dpi = 300)
    
    cat("Best K plot with individual labels saved as:", paste0(out.name, "_best_K", bestK, "_labeled.png"), "\n")
  }
}

############################
# 3. SAMPLE ASSIGNMENTS    #
############################

# Create sample assignment file for best K
if (!is.null(bestK) && paste0("K", bestK) %in% names(q_data_list)) {
  cat("\n--- Creating sample assignments for best K ---\n")
  
  # Get the Q data for best K
  best_k_data <- q_data_list[[paste0("K", bestK)]]
  
  # Create assignment table
  assignments <- best_k_data %>%
    group_by(SampleID, IID) %>%
    slice_max(Proportion, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      PopulationID = paste0("Pop", as.numeric(str_extract(Ancestry, "\\d+"))),
      MaxProportion = round(Proportion, 3)
    ) %>%
    select(SampleID, IID, PopulationID, MaxProportion) %>%
    arrange(as.numeric(SampleID))
  
  # Save assignment file
  assignment_file <- paste0(out.name, "_K", bestK, "_assignments.txt")
  write.table(assignments, file = assignment_file, 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  cat("Sample assignments for K =", bestK, "saved as:", assignment_file, "\n")
  cat("Assignment summary:\n")
  assignment_summary <- assignments %>% 
    count(PopulationID, name = "n_samples") %>%
    mutate(percentage = round(100 * n_samples / sum(n_samples), 1))
  print(assignment_summary)
  
  # Also save a detailed version with all ancestry proportions
  detailed_assignments <- best_k_data %>%
    select(SampleID, IID, Ancestry, Proportion) %>%
    pivot_wider(names_from = Ancestry, values_from = Proportion, names_prefix = "Pop") %>%
    mutate(across(starts_with("Pop"), ~ round(.x, 3))) %>%
    arrange(as.numeric(SampleID))
  
  detailed_file <- paste0(out.name, "_K", bestK, "_detailed_assignments.txt")
  write.table(detailed_assignments, file = detailed_file, 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  cat("Detailed assignments saved as:", detailed_file, "\n")
}

# Summary
cat("\n=== SUMMARY ===\n")
cat("Processed", length(q_files), "Q files (excluding K=1)\n")
cat("K values found:", paste(k_values, collapse = ", "), "\n")
if (!is.null(bestK)) {
  cat("Best K from CV:", bestK, "\n")
}
cat("Output files created:\n")
if (file.exists(cv_file)) {
  cat("  -", paste0(out.name, "_cross_validation.png"), "\n")
}
cat("  -", paste0(out.name, "_admixture_combined.png"), "\n")
for (k in k_values) {
  cat("  -", paste0(out.name, "_K", k, "_admixture.png"), "\n")
}
if (!is.null(bestK)) {
  cat("  -", paste0(out.name, "_best_K", bestK, "_admixture.png"), "\n")
  cat("  -", paste0(out.name, "_K", bestK, "_assignments.txt"), "\n")
  cat("  -", paste0(out.name, "_K", bestK, "_detailed_assignments.txt"), "\n")
}
cat("  -", log_file, "\n")

# Close logging
sink()

cat("Analysis complete! Check", log_file, "for details.\n")