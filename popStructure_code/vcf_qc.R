#################### ipyrad VCF QC R script ####################
### Usage: module load Rtidyverse; Rscript vcf_qc.R indiv_miss.imiss site_miss.lmiss indiv_depth.idepth

options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Load in libraries
library(tidyverse)
if (!require("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr")
  library(ggpubr)
}


# Load input files: individual missingness, individual mean depth, site missingness from VCF
cmd_args <- commandArgs(trailingOnly = TRUE)

imiss.file <- cmd_args[1]
smiss.file <- cmd_args[2]
idp.file <- cmd_args[3]

# Set up logging file
log_file <- "vcf_qc_r.log"
sink(log_file, append = FALSE, split = TRUE)  # split = TRUE shows output on console too

# Add timestamp to log
cat("=== VCF QC Analysis Started ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input files:", paste(cmd_args, collapse = ", "), "\n\n")


# Read in data
imiss <- read.delim(imiss.file)
smiss <- read.delim(smiss.file)
idp <- read.delim(idp.file)

# Calculate summary statistics
imiss.summary <- imiss %>% summarise(n = n(),
                                     samples = n(),
                                     mean_missing = mean(F_MISS, na.rm = T),
                                     median_missing = median(F_MISS),
                                     sd_missing = sd(F_MISS),
                                     min_missing = min(F_MISS),
                                     max_missing = max(F_MISS),
                                     q25_missing = quantile(F_MISS, 0.25, na.rm = T),
                                     q75_missing = quantile(F_MISS, 0.75, na.rm = T))

smiss.summary <- smiss %>% summarise(n = n(),
                                     samples = n(),
                                     mean_missing = mean(F_MISS, na.rm = T),
                                     median_missing = median(F_MISS),
                                     sd_missing = sd(F_MISS),
                                     min_missing = min(F_MISS),
                                     max_missing = max(F_MISS),
                                     q25_missing = quantile(F_MISS, 0.25, na.rm = T),
                                     q75_missing = quantile(F_MISS, 0.75, na.rm = T))

dp.summary <- idp %>% summarise(n = n(),
                                     samples = n(),
                                     mean_missing = mean(MEAN_DEPTH, na.rm = T),
                                     median_missing = median(MEAN_DEPTH),
                                     sd_missing = sd(MEAN_DEPTH),
                                     min_missing = min(MEAN_DEPTH),
                                     max_missing = max(MEAN_DEPTH),
                                     q25_missing = quantile(MEAN_DEPTH, 0.25, na.rm = T),
                                     q75_missing = quantile(MEAN_DEPTH, 0.75, na.rm = T))


# Create summary distribution plots
# Individual Missingness histogram 
hist_data <- hist(imiss$F_MISS, breaks=15, plot=FALSE)
x_min <- min(imiss$F_MISS)
x_max <- max(imiss$F_MISS)
x_range <- x_max - x_min
y_max <- max(hist_data$counts)

imiss.p <- ggplot(imiss, aes(x=F_MISS)) +
  geom_histogram(bins=15, fill="lightpink", color="pink4", alpha=0.8) +
  geom_vline(xintercept = imiss.summary$mean, linetype="dashed", color="red") +
  annotate("text",  x = imiss.summary$mean, y = Inf, hjust = -0.1, vjust = 1.5, color = "red",
           label = paste("Mean:", round(imiss.summary$mean, 3))) +
  annotate("text", 
           x = x_min + x_range * 0.815,
           y = y_max * 0.775,
           label = paste("Min: ", round(imiss.summary$min, 3), "\n",
                         "Max: ", round(imiss.summary$max, 3), "\n", 
                         "Median: ", round(imiss.summary$median, 3), "\n",
                         "Mean: ", round(imiss.summary$mean, 3)),
           hjust = 0.5, vjust = 0.5, size = 3.2) +
  labs(title = "Individual Missingness Distribution", x = "Freq. of missing", y = "Count") +
  theme_minimal()


# SNP Missingness historgram
hist_data_smiss <- hist(smiss$F_MISS, breaks=15, plot=FALSE)
x_min_smiss <- min(smiss$F_MISS)
x_max_smiss <- max(smiss$F_MISS)
x_range_smiss <- x_max_smiss - x_min_smiss
y_max_smiss <- max(hist_data_smiss$counts)

smiss.p <- ggplot(smiss, aes(x=F_MISS)) +
  geom_histogram(bins=15, fill="lightblue", color="darkblue", alpha=0.8) +
  geom_vline(xintercept = smiss.summary$mean_missing, linetype="dashed", color="red") +
  annotate("text", x = smiss.summary$mean_missing, y = Inf, hjust = -0.1, vjust = 1.5, color = "red",
           label = paste("Mean:", round(smiss.summary$mean_missing, 3))) +
  annotate("text", 
           x = x_min_smiss + x_range_smiss * 0.185,
           y = y_max_smiss * 0.775,
           label = paste("Min: ", round(smiss.summary$min_missing, 3), "\n",
                         "Max: ", round(smiss.summary$max_missing, 3), "\n", 
                         "Median: ", round(smiss.summary$median_missing, 3), "\n",
                         "Mean: ", round(smiss.summary$mean_missing, 3)),
           hjust = 0.5, vjust = 0.5, size = 3.2) +
  labs(title = "SNP Missingness Distribution", x = "Freq. of missing", y = "Count") +
  theme_minimal()


# Mean individual depth histogram
hist_data_dp <- hist(idp$MEAN_DEPTH, breaks=15, plot=FALSE)
x_min_dp <- min(idp$MEAN_DEPTH)
x_max_dp <- max(idp$MEAN_DEPTH)
x_range_dp <- x_max_dp - x_min_dp
y_max_dp <- max(hist_data_dp$counts)

dp.p <- ggplot(idp, aes(x=MEAN_DEPTH)) +
  geom_histogram(bins=15, fill="lightgreen", color="darkgreen", alpha=0.8) +
  geom_vline(xintercept = dp.summary$mean_missing, linetype="dashed", color="red") +
  annotate("text", x = dp.summary$mean_missing, y = Inf, hjust = -0.1, vjust = 1.5, color = "red",
           label = paste("Mean:", round(dp.summary$mean_missing, 3))) +
  annotate("text", 
           x = x_min_dp + x_range_dp * 0.815,
           y = y_max_dp * 0.775,
           label = paste("Min: ", round(dp.summary$min_missing, 3), "\n",
                         "Max: ", round(dp.summary$max_missing, 3), "\n", 
                         "Median: ", round(dp.summary$median_missing, 3), "\n",
                         "Mean: ", round(dp.summary$mean_missing, 3)),
           hjust = 0.5, vjust = 0.5, size = 3.2) +
  labs(title = "Individual Depth Distribution", x = "Mean depth", y = "Count") +
  theme_minimal()


# Create additional summary plots
# Missingness per individual
imiss.p2 <- ggplot(imiss, aes(x=reorder(INDV, F_MISS), y=F_MISS)) +
  geom_point(color="pink4", alpha=0.8) +
  geom_hline(yintercept = 0.5, linetype="dashed", color="darkred") +
  annotate("text", x = 0, y = 0.5, hjust = -0.1, vjust = 1.5, color = "darkred",
           label = paste("50%")) +
  labs(title = "Individual Missingness", x = "Sample", y = "Freq. of missingness") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), axis.text.x = element_blank())

# SNP missingness per chromosome/scaffold
# Find top 10 chromosomes with highest missingness
top_missing_chroms <- smiss %>%
  group_by(CHR) %>%
  summarise(
    n_snps = n(),
    mean_missingness = mean(F_MISS),
    median_missingness = median(F_MISS),
    max_missingness = max(F_MISS),
    .groups = 'drop'
  ) %>%
  arrange(desc(mean_missingness)) %>%
  slice_head(n = 10)

smiss.p2 <- ggplot(smiss %>% filter(CHR %in% top_missing_chroms$CHR), aes(x=reorder(CHR, F_MISS), y=F_MISS)) +
  geom_boxplot(fill="lightblue", color="darkblue") +
  labs(title = "Top 10 chroms with highest SNP missingness", x = "Chromosome", y = "Freq. of missingness") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 45))

# Depth per individual
dp.p2 <- ggplot(idp, aes(x=reorder(INDV, MEAN_DEPTH), y=MEAN_DEPTH)) +
  geom_point(color="darkgreen", alpha=0.8) +
  #geom_hline(yintercept = 10, linetype="dashed", color="red") +
  labs(title = "Individual Depth", x = "Sample", y = "Mean depth") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), axis.text.x = element_blank())


# Combine summary plots
toprow <- ggarrange(imiss.p, smiss.p, dp.p, nrow = 1)
botrow <- ggarrange(imiss.p2, smiss.p2, dp.p2, nrow = 1)
plots <- ggarrange(toprow, botrow, nrow = 2)


# Save the combined plots as PDF
pdf("vcf_qc_summary_plots.pdf", width = 15, height = 10)
print(plots)
dev.off()

cat("Summary plots saved to: vcf_qc_summary_plots.pdf\n")


# Identify samples with missing data at multiple cutoffs
missing_cutoffs <- c(0.8, 0.75, 0.5)
missing_cutoff_names <- c("80", "75", "50")

for(i in 1:length(missing_cutoffs)) {
  cutoff <- missing_cutoffs[i]
  cutoff_name <- missing_cutoff_names[i]
  
  # Filter samples above cutoff
  missing_samples <- imiss %>%
    filter(F_MISS > cutoff) %>%
    select(INDV, F_MISS) %>%
    arrange(desc(F_MISS))
  
  # Save to CSV
  filename <- paste0("samples_missing_", cutoff_name, "percent.csv")
  write.csv(missing_samples, filename, row.names = FALSE, quote = F)
  
  cat("Samples with >", cutoff_name, "% missing data saved to:", filename, "\n")
  cat("Number of samples with >", cutoff_name, "% missing data:", nrow(missing_samples), "\n")
}


# Identify samples with less than 5 mean depth
low_depth_samples <- idp %>%
  filter(MEAN_DEPTH < 5) %>%
  select(INDV, MEAN_DEPTH) %>%
  arrange(MEAN_DEPTH)

# Save low depth samples to CSV
write.csv(low_depth_samples, "samples_low_depth.csv", row.names = FALSE, quote = F)

cat("Samples with <5 mean depth saved to: samples_low_depth.csv\n")
cat("Number of samples with <5 mean depth:", nrow(low_depth_samples), "\n")


# Print summary of flagged samples at different cutoffs
cat("\n=== SUMMARY OF FLAGGED SAMPLES ===\n")
for(i in 1:length(missing_cutoffs)) {
  cutoff <- missing_cutoffs[i]
  cutoff_name <- missing_cutoff_names[i]
  n_samples <- nrow(imiss %>% filter(F_MISS > cutoff))
  cat("High missingness (>", cutoff_name, "%):", n_samples, "samples\n")
}
cat("Low depth (<5x):", nrow(low_depth_samples), "samples\n")

# Check for overlap between high missingness (80%) and low depth samples
high_missing_80 <- imiss %>% filter(F_MISS > 0.8)
overlap_samples <- intersect(high_missing_80$INDV, low_depth_samples$INDV)
cat("Samples with both high missingness (>80%) AND low depth:", length(overlap_samples), "samples\n")

if(length(overlap_samples) > 0) {
  cat("Overlapping samples:", paste(overlap_samples, collapse = ", "), "\n")
}

# Print session info
sessionInfo()

# End log file
cat("\n=== VCF QC Analysis Completed ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
sink()  # Close the log file