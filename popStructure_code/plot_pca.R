#!/usr/bin/env Rscript
#################### Dynamic PCA Plot R script ####################
### Usage: Rscript plot_pca.R <eigenvalues> <eigenvectors> <metadata> <output_name> [sample_col] [color_cols...]
### Example: Rscript plot_pca.R pca.eigenval pca.eigenvec metadata.txt output ID pop seq_run location
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Create a personal library directory if it doesn't exist
personal_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(personal_lib)) {
  dir.create(personal_lib, recursive = TRUE)
}

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
if (length(cmd_args) < 4) {
  stop("Usage: Rscript plot_pca.R <eigenvalues> <eigenvectors> <metadata> <output_name> [sample_col] [color_cols...]")
}

val.file <- cmd_args[1]
vec.file <- cmd_args[2]
meta.file <- cmd_args[3]
out.name <- cmd_args[4]
sample_col <- if(length(cmd_args) >= 5) cmd_args[5] else NULL
color_cols <- if(length(cmd_args) > 5) cmd_args[6:length(cmd_args)] else NULL

# Set up logging and output directory
output_dir <- "plink_pca_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

log_file <- file.path(output_dir, paste0(out.name, "_plot_pca.log"))
sink(log_file, append = FALSE, split = TRUE)

cat("=== Dynamic PCA Plot Started ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input files:", paste(cmd_args[1:4], collapse = ", "), "\n")
if (!is.null(sample_col)) cat("Sample column:", sample_col, "\n")
if (!is.null(color_cols)) cat("Color columns:", paste(color_cols, collapse = ", "), "\n\n")

# Read PCA data
cat("Reading PCA files...\n")
pca.val <- read.table(val.file)
pca.vec <- read.table(vec.file)

# Calculate variance explained
pca.val.var <- pca.val %>% 
  mutate(varex = V1/sum(V1), 
         percVarex = round(V1/sum(V1)*100, 2))

# Check eigenvector file structure and assign column names
cat("Checking eigenvector file structure...\n")
cat("Eigenvector file dimensions:", nrow(pca.vec), "x", ncol(pca.vec), "\n")

# Check if first column(s) contain sample IDs vs PCs
# Sample ID columns are typically character or have limited unique values
# PC columns are typically numeric with many decimal values
col1_is_sample_id <- is.character(pca.vec[,1]) || is.factor(pca.vec[,1]) || 
                     (!is.numeric(pca.vec[,1])) ||
                     (is.numeric(pca.vec[,1]) && all(pca.vec[,1] == as.integer(pca.vec[,1])))

col2_is_sample_id <- if(ncol(pca.vec) > 1) {
  is.character(pca.vec[,2]) || is.factor(pca.vec[,2]) || 
  (!is.numeric(pca.vec[,2])) ||
  (is.numeric(pca.vec[,2]) && all(pca.vec[,2] == as.integer(pca.vec[,2])))
} else {
  FALSE
}

# Determine file structure
if (col1_is_sample_id && col2_is_sample_id && ncol(pca.vec) > 2) {
  # Format: FID IID PC1 PC2 PC3 ...
  cat("Detected format: FID IID PC1 PC2 ...\n")
  max_pcs <- min(ncol(pca.vec) - 2, 20)
  pc_names <- paste0("PC", 1:max_pcs)
  colnames(pca.vec) <- c("FID", "IID", pc_names)
} else if (col1_is_sample_id && ncol(pca.vec) > 1) {
  # Format: IID PC1 PC2 PC3 ...
  cat("Detected format: IID PC1 PC2 ...\n")
  max_pcs <- min(ncol(pca.vec) - 1, 20)
  pc_names <- paste0("PC", 1:max_pcs)
  colnames(pca.vec) <- c("IID", pc_names)
  # Add FID column for consistency (set to same as IID)
  pca.vec$FID <- pca.vec$IID
  pca.vec <- pca.vec[, c("FID", "IID", pc_names)]
} else {
  # Fallback - assume all columns are PCs, use row numbers as IID
  cat("Warning: Could not detect sample ID columns. Using row numbers as IID\n")
  max_pcs <- min(ncol(pca.vec), 20)
  pc_names <- paste0("PC", 1:max_pcs)
  pca.vec <- pca.vec[, 1:max_pcs]
  colnames(pca.vec) <- pc_names
  pca.vec$IID <- as.character(1:nrow(pca.vec))
  pca.vec$FID <- pca.vec$IID
  pca.vec <- pca.vec[, c("FID", "IID", pc_names)]
}

cat("Final eigenvector structure: FID, IID, and", length(pc_names), "PCs\n")

# Read and process metadata
cat("Reading metadata...\n")
# Detect file format and read accordingly
if (grepl("\\.csv$", meta.file, ignore.case = TRUE)) {
  cat("Detected CSV format\n")
  meta <- read.csv(meta.file, stringsAsFactors = FALSE)
} else {
  cat("Detected tab-delimited format\n")
  meta <- read.delim(meta.file, stringsAsFactors = FALSE)
}
cat("Metadata dimensions:", nrow(meta), "x", ncol(meta), "\n")
cat("Metadata columns:", paste(colnames(meta), collapse = ", "), "\n")

# Auto-detect sample column if not provided
if (is.null(sample_col)) {
  # Try common sample column names
  possible_cols <- c("IID", "ID", "sample", "Sample", "sample_id", "Sample_ID", 
                    "individual", "Individual", "SampleID", "sample.id")
  sample_col <- intersect(possible_cols, colnames(meta))[1]
  
  if (is.null(sample_col)) {
    cat("Available columns in metadata:\n")
    print(colnames(meta))
    stop("Could not auto-detect sample column. Please specify it as the 5th argument.")
  }
}

cat("Using sample column:", sample_col, "\n")

# Rename sample column to IID for joining
if (sample_col != "IID") {
  meta <- meta %>% rename(IID = !!sym(sample_col))
}

# Check and harmonize data types before joining
cat("Checking data types for joining...\n")
cat("PCA IID type:", class(pca.vec$IID), "\n")
cat("Metadata IID type:", class(meta$IID), "\n")

# Function to harmonize IID columns
harmonize_iid_types <- function(pca_data, meta_data) {
  pca_iid_class <- class(pca_data$IID)
  meta_iid_class <- class(meta_data$IID)
  
  if (pca_iid_class != meta_iid_class) {
    cat("Data type mismatch detected. Harmonizing...\n")
    
    # Try to convert both to character first (safest approach)
    pca_data$IID <- as.character(pca_data$IID)
    meta_data$IID <- as.character(meta_data$IID)
    
    cat("Converted both IID columns to character\n")
    
    # Check if we can convert to numeric (in case original data was meant to be numeric)
    pca_numeric_test <- suppressWarnings(as.numeric(pca_data$IID))
    meta_numeric_test <- suppressWarnings(as.numeric(meta_data$IID))
    
    # If both can be converted to numeric without introducing NAs, use numeric
    if (!any(is.na(pca_numeric_test)) && !any(is.na(meta_numeric_test))) {
      pca_data$IID <- pca_numeric_test
      meta_data$IID <- meta_numeric_test
      cat("Successfully converted both IID columns to numeric\n")
    } else {
      cat("Keeping both IID columns as character\n")
    }
  } else {
    cat("Data types match - no conversion needed\n")
  }
  
  return(list(pca = pca_data, meta = meta_data))
}

# Harmonize the data types
harmonized <- harmonize_iid_types(pca.vec, meta)
pca.vec <- harmonized$pca
meta <- harmonized$meta

# Verify the types match now
cat("After harmonization:\n")
cat("PCA IID type:", class(pca.vec$IID), "\n")
cat("Metadata IID type:", class(meta$IID), "\n")

# Auto-detect color columns if not provided
if (is.null(color_cols)) {
  # Exclude the sample column and any obvious non-grouping columns
  exclude_cols <- c("IID", "FID", "depth", "missing", "coverage", "reads", 
                   "date", "Date", "notes", "Notes", "file", "File")
  
  potential_cols <- setdiff(colnames(meta), exclude_cols)
  
  # Filter for columns with reasonable number of categories (2-20) and check data types
  categorical_cols <- c()
  for (col in potential_cols) {
    col_data <- meta[[col]][!is.na(meta[[col]])]
    unique_vals <- length(unique(col_data))
    
    # Check if column is suitable for grouping
    is_character_or_factor <- is.character(col_data) || is.factor(col_data)
    is_reasonable_categories <- unique_vals >= 2 && unique_vals <= 20
    is_mostly_non_numeric <- !is.numeric(col_data) || 
                            (is.numeric(col_data) && unique_vals <= 10)
    
    if (is_reasonable_categories && (is_character_or_factor || is_mostly_non_numeric)) {
      categorical_cols <- c(categorical_cols, col)
    }
  }
  
  color_cols <- categorical_cols
  cat("Auto-detected grouping columns:", paste(color_cols, collapse = ", "), "\n")
}

# Join PCA data with metadata
cat("Attempting to join PCA data with metadata...\n")
pca.vec.meta <- left_join(pca.vec, meta, by = "IID")
cat("Join successful!\n")

missing_samples <- sum(is.na(pca.vec.meta[[color_cols[1]]]))
if (missing_samples > 0) {
  cat("Warning:", missing_samples, "samples missing metadata\n")
}

# Function to create color palette
get_colors <- function(n) {
  if (n <= 8) {
    brewer.pal(max(3, n), "Set2")
  } else if (n <= 12) {
    brewer.pal(n, "Set3")
  } else {
    rainbow(n)
  }
}

# Function to create a PCA plot
create_pca_plot <- function(data, color_var = NULL, title = "PCA") {
  p <- ggplot(data, aes(PC1, PC2))
  
  if (!is.null(color_var) && color_var %in% colnames(data)) {
    # Remove rows with NA values for the color variable
    data_clean <- data[!is.na(data[[color_var]]), ]
    
    if (nrow(data_clean) == 0) {
      cat("Warning: All values are NA for", color_var, "- creating basic plot\n")
      p <- p + geom_point(size = 3, alpha = 0.8, color = "steelblue")
    } else {
      # Check if variable is numeric and continuous
      col_data <- data_clean[[color_var]]
      unique_vals <- length(unique(col_data))
      
      if (is.numeric(col_data) && unique_vals > 10) {
        # Continuous variable - use color gradient
        p <- ggplot(data_clean, aes(PC1, PC2)) +
             geom_point(aes(color = !!sym(color_var)), size = 3, alpha = 0.8) +
             scale_color_gradient(low = "blue", high = "red")
        cat("Using continuous color scale for", color_var, "\n")
      } else {
        # Categorical variable - convert to factor and use discrete colors
        data_clean[[color_var]] <- as.factor(data_clean[[color_var]])
        n_groups <- length(unique(data_clean[[color_var]]))
        colors <- get_colors(n_groups)
        
        p <- ggplot(data_clean, aes(PC1, PC2)) +
             geom_point(aes(color = !!sym(color_var)), size = 3, alpha = 0.8) +
             scale_color_manual(values = colors)
        cat("Using discrete color scale for", color_var, "with", n_groups, "groups\n")
      }
    }
  } else {
    p <- p + geom_point(size = 3, alpha = 0.8, color = "steelblue")
  }
  
  p <- p +
    xlab(paste0("PC1 (", pca.val.var$percVarex[1], "%)")) +
    ylab(paste0("PC2 (", pca.val.var$percVarex[2], "%)")) +
    ggtitle(title) +
    theme_minimal() +
    theme(text = element_text(size = 12),
          legend.position = "bottom",
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 9))
  
  return(p)
}

# Create plots
plots <- list()

# Basic PCA without coloring
plots[["basic"]] <- create_pca_plot(pca.vec.meta, title = "PCA - All Samples")

# Create plots for each specified grouping variable
for (col in color_cols) {
  if (col %in% colnames(pca.vec.meta)) {
    clean_col_name <- gsub("[^A-Za-z0-9]", "_", col)
    plots[[clean_col_name]] <- create_pca_plot(pca.vec.meta, col, 
                                               paste("PCA colored by", col))
    cat("Created plot for:", col, "\n")
  } else {
    cat("Warning: Column", col, "not found in metadata\n")
  }
}

# Save individual plots
timestamp <- format(Sys.time(), "%Y-%m-%d")
for (plot_name in names(plots)) {
  filename <- file.path(output_dir, paste0(out.name, "_pca_", plot_name, "_", timestamp, ".png"))
  ggsave(plots[[plot_name]], filename = filename, width = 10, height = 8, 
         units = "in", dpi = 300, bg="white")
  cat("Saved:", filename, "\n")
}

# Create combined plot if multiple grouping variables
if (length(plots) > 1 && length(plots) <= 6) {
  combined_plot <- ggarrange(plotlist = plots, ncol = 2, nrow = ceiling(length(plots)/2))
  combined_filename <- file.path(output_dir, paste0(out.name, "_pca_combined_", timestamp, ".png"))
  ggsave(combined_plot, filename = combined_filename, 
         width = 16, height = 6 * ceiling(length(plots)/2), 
         units = "in", dpi = 300)
  cat("Saved combined plot:", combined_filename, "\n")
}

# Save PCA coordinates with metadata for further analysis
output_data <- pca.vec.meta %>% 
  select(FID, IID, all_of(pc_names[1:min(10, length(pc_names))]), 
         any_of(color_cols))
csv_filename <- file.path(output_dir, paste0(out.name, "_pca_coords_with_metadata.csv"))
write.csv(output_data, csv_filename, row.names = FALSE)
cat("Saved PCA coordinates:", csv_filename, "\n")

# Print summary
cat("\n=== Summary ===\n")
cat("Total samples:", nrow(pca.vec.meta), "\n")
cat("PC1 variance explained:", pca.val.var$percVarex[1], "%\n")
cat("PC2 variance explained:", pca.val.var$percVarex[2], "%\n")
cat("Plots created:", length(plots), "\n")

# Print session info
cat("\n=== Session Info ===\n")
sessionInfo()

cat("\n=== Dynamic PCA Plot Completed ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()