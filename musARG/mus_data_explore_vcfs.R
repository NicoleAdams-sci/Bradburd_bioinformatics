# Explore Mus Freeze2 VCFs from Beth Dumont
# 10/15/2025
# Nicole Adams

library(tidyverse)
library(viridisLite)
library(ggpubr)
# Load mapping libraries
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(leaflet)


setwd("~/Documents/NicoleAdams/bioinfo_help/musARG/")

# saved excel file as csv to load easily in R
meta <- read.csv("data/mus_WholeGenomeSequencing_20251002.csv") #full data

###### Clean data ######
# make columns numeric
meta <- meta %>% mutate(across(c(Insert.Size, X..bases.sequenced, Approx.Depth, Latitude, Longitude), as.numeric))
meta$Subspecies <- gsub("Mus cypriacus", "cypriacus", meta$Subspecies)
meta$Subspecies <- gsub("Mus spretus", "spretus", meta$Subspecies)

# load in individual missingness and avg depth calc'd from vcftools --missing-indv and --depth
imiss <- read.delim("output/chr19_indiv_miss.imiss")
idp <- read.delim("output/chr19_indiv_depth.idepth")

# merge missingness and depth
imiss.idp <- full_join(imiss, idp)

# check to see the overlap in sample naming bt vcf and metadata
length(intersect(meta$Sample.ID, imiss.idp$INDV))

# add column to metadata so match vcf
meta$INDV <- meta$Sample.ID

# merge vcf samples with missingness and depth with metadata
mus <- left_join(imiss.idp, meta)

# Summary table and bar plot for Sample Type (wild-caught, inbred etc)
sampletype_summary <- mus %>%
  group_by(Sample.Type) %>%
  summarise(count = n())

samptype.p <- ggplot(sampletype_summary, aes(x = reorder(Sample.Type, -count), y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count of Samples by Sample type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


################# Identify & Remove unwanted samples ################# 
# save list of non-wild-caught samples to remove from VCF
inbreds <- mus %>% filter(Sample.Type != "wild-caught") %>% select(INDV)
write.csv(inbreds, file="data/inbred_mus_samples.csv", quote = FALSE, row.names = FALSE)

# save list of hybrid samples to remove from VCF
hybrids <- mus %>% filter(Subspecies == "hybrid") %>% select(INDV)
write.csv(hybrids, file="data/hybrids_mus_samples.csv", quote = FALSE, row.names = FALSE)

# save list of low depth samples - intersection of samples with low Approx.Depth and low vcf MEAN_DEPTH
adp10 <- mus %>% filter(Approx.Depth < 10) %>% select(INDV)
dp10 <- mus %>% filter(MEAN_DEPTH < 10) %>% select(INDV)
lowcov <- mus %>% filter(Notes == "low coverage") %>% select(INDV)
low.dp <- intersect(adp10, dp10)
write.csv(low.dp, file="data/low_depth_mus_samples.csv", quote = FALSE, row.names = FALSE)

# Remove non-wild-caught mus samples (N=80) and hybrids (N=19)
mus <- mus %>% filter(Sample.Type == "wild-caught") %>% filter(Subspecies != "hybrid") %>% filter(!INDV %in% low.dp$INDV)
#################################################################### 


# Summary table and bar plot for 'Subspecies'
subspecies_summary <- mus %>%
  group_by(Subspecies) %>%
  summarise(count = n())

subsp.p <- ggplot(subspecies_summary, aes(x = reorder(Subspecies, -count), y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count of Samples by Subspecies") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Summary table and bar plot for 'Sex'
sex_summary <- mus %>%
  group_by(Sex) %>%
  summarise(count = n())

sex.p <- ggplot(sex_summary, aes(x = reorder(Sex, -count), y=count)) +
  geom_col() +
  labs(title = "Count of Samples by Sex") +
  theme_minimal() 

# Summary table and plot for 'Approx Depth' (from metadata)
dp.p <- ggplot(mus %>% filter(!is.na(Approx.Depth)), aes(x = Approx.Depth)) +
  geom_histogram(binwidth = 2, fill = "lightblue", color = "navy") +
  labs(title = "Distribution of Sequencing Depth", x = "Approximate Depth", y = "Count") +
  theme_minimal()

# Summary table and bar plot for 'Library Layout'
library_layout_summary <- mus %>%
  group_by(Library.Layout) %>%
  summarise(count = n())

lib.p <- ggplot(library_layout_summary, aes(x = Library.Layout, y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count of Samples by Library Layout") +
  theme_minimal() 

# Summary table and bar plot for 'Sequencing Platform'
sequencing_platform_summary <- mus %>%
  group_by(Sequencing.Platform) %>%
  summarise(count = n())

seq.platform.p <- ggplot(sequencing_platform_summary %>% filter(count > 10), aes(x = reorder(Sequencing.Platform, -count), y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "Count of Samples by Sequencing Platform") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Summary table and plot of subspecies and sex
subspecies_sex_summary <- mus %>%
  group_by(Subspecies, Sex) %>%
  summarise(count = n()) %>%
  ungroup()

subsp.sex.p <- ggplot(subspecies_sex_summary, aes(x = Subspecies, y = count, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d() +
  labs(title = "Sample Count by Subspecies and Sex", fill = "Sex") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Summary plot for 'MEAN_DEPTH' (from vcftools)
dp.p2 <- ggplot(mus, aes(x = MEAN_DEPTH)) +
  geom_histogram(binwidth = 2, fill = "lightblue", color = "navy") +
  labs(title = "Average Individual Depth", x = "Mean Depth (vcf)", y = "Count") +
  theme_minimal()

# Summary plot for individual missningness (from vcftools)
imiss.p <- ggplot(mus, aes(x = F_MISS)) +
  geom_histogram( fill = "cornflowerblue", color = "navy") +
  labs(title = "Individual Missingness in filtered VCF", x = "Frequency of missingness", y = "Count") +
  theme_minimal()


## Map

# Clean your data by removing rows with missing lat/long
map_data <- mus %>%
  drop_na(Latitude, Longitude)

# Simple map using ggplot2
world_map <- ne_countries(scale = "medium", returnclass = "sf")

map.p <- ggplot(data = world_map) +
  geom_sf() +
  geom_point(data = map_data, aes(x = Longitude, y = Latitude), color = "red", size = 2, alpha = 0.5) +
  labs(title = "Sample Locations") +
  theme_void()

# Interactive map using leaflet
# This is great for exploring individual points
leaflet(map_data) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~Longitude,
    lat = ~Latitude,
    popup = ~paste("Subspecies:", `Subspecies`, "<br>", "Site:", `Site`)
  )

# Summarize data by Latitude and Longitude to get the sample size at each location
map_data_summary <- mus %>%
  drop_na(Latitude, Longitude) %>%
  group_by(Latitude, Longitude) %>%
  summarise(sample_size = n()) %>%
  ungroup()

# Plotting by point size
map.sampSz.p <- ggplot(data = world_map) +
  geom_sf() +
  geom_point(data = map_data_summary, aes(x = Longitude, y = Latitude, size = sample_size), color = "red", alpha = 0.5) +
  scale_size_continuous(name = "Sample Size") +
  labs(title = paste("Sample Locations by Size ( N =",dim(mus)[1], ")")) +
  theme_minimal() +
  theme(legend.position = "right")

# Plotting by subspecies
map_data$Subspecies <- map_data$Subspecies %>%
  na_if("") %>%
  replace_na("not_labeled")

map_data_summary2 <- map_data %>%
  group_by(Latitude, Longitude, Subspecies) %>%
  summarise(sample_size = n()) %>%
  ungroup()

map.sampSz.p2 <- ggplot(data = world_map) +
  geom_sf() +
  geom_point(data = map_data_summary2, aes(x = Longitude, y = Latitude, 
                                           shape = Subspecies, size = sample_size, color = Subspecies),
             alpha = 0.5) +
  scale_size_continuous(name = "Sample Size") +
  scale_color_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 18)) +
  labs(title = paste("Sample Locations by Subspecies ( N =",dim(mus)[1], ")")) +
  theme_minimal() +
  theme(legend.position = "right")


# ---  Combine the Plots  ---
bottom_plots <- ggarrange(subsp.p, dp.p2, imiss.p, ncol = 3)

# Stack the map on top of the combined summary plots
final_plot <- ggarrange(map.sampSz.p2, bottom_plots, nrow = 2,
                        heights = c(2, 1), # Make the top plot (map) taller
                        common.legend = FALSE)

# Add a single title to the entire figure
final_plot <- annotate_figure(final_plot,
                              top = text_grob("Summary of Mus Freeze2 filtered VCF", face = "bold", size = 18))


# Save summary figure
today <- format(Sys.Date(), format = "%Y%m%d")
ggsave(
  filename = paste0("mus_data_summary_", today, ".png"),
  plot = final_plot,
  width = 12, # Set the width of the saved image in inches
  height = 8, # Set the height of the saved image in inches
  bg = "white",
  dpi = 300 # Set the resolution for high quality
)
