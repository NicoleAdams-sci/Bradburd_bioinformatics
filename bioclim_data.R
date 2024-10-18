# Downloading bioclim data, extracting sample lat/longs from it, looking for correlations, making PCA on climate data
# 10/18/2024
# Nicole E. Adams

library(geodata) # Used to download worldclim/bioclim data
library(raster)
library(psych)    # Used to investigate correlations among variables
library(corrplot) # Used to investigate correlations among variables
library(ggbiplot) # plot pca biplot


# Download bioclim data - all variables, only the USA, at a 10 minutes of a degree resolution
bioclim_data <- worldclim_global(var = "bio", res = 10, path = "~/Documents/NicoleAdams/worldClimData/")


# Once downloaded, you don't need to run that first line, instead you can call the data locally
bioclim_dir <- "~/Documents/NicoleAdams/worldClimData/climate/wc2.1_10m/"
bioclim_files <- list.files(bioclim_dir, pattern = "\\.tif$", full.names = T)


# load all the files into a raster stack
bioclim_stack <- stack(bioclim_files)

# plot just for funsies
plot(bioclim_stack) # plots the whole world for each variable


# pull out specific lat longs
mice <- read.csv("~/Documents/NicoleAdams/pman/bradburd_pman_metdata.csv")

mice_clim_values <- extract(bioclim_stack, mice [, c("LON", "LAT")])

mice_clim <- cbind(mice, mice_clim_values)

# rename columns to just the bioclim variable (for variable descriptions see https://www.worldclim.org/data/bioclim.html#google_vignette)
colnames(mice_clim) <- gsub("wc2.1_10m_", "", colnames(mice_clim))

# Check for correlation between bioclim variables - generally |0.7| is considered correlated
pairs.panels(mice_clim[,c(8:26)], scale=T)

mice_clim.corr <- cor(mice_clim[,c(8:26)])
corrplot(mice_clim.corr,method = "ellipse", tl.srt = 0)

# PCA on clim variables. "You should set the scale. argument to TRUE, to scale the data matrix to unit variances."
pca <- prcomp(mice_clim[,c(8:26)], scale. = TRUE)
pca_scores <- data.frame(pca$x) # object stores the transformation and the transformed values, which are called scores
pairs(pca_scores) # the transformed variables (now called PC1, PC2, PC3, PC4...) should be now uncorrelated
summary(pca) 

ggbiplot(pca,
         circle = TRUE,
         varname.size = 2,
         ellipse = TRUE, ellipse.level = 0.5, ellipse.alpha = 0.1) +
  theme_minimal() 


# plot bioclim variables for the region of interest
buffer <- 1
xmin <- min(mice$LON, na.rm = T)
xmax <- max(mice$LON, na.rm = T)
ymin <- min(mice$LAT, na.rm = T)
ymax <- max(mice$LAT, na.rm = T)

crop_extent <- extent(xmin-buffer, xmax+buffer, ymin-buffer, ymax+buffer)

bioclim_crop <- crop(bioclim_stack, crop_extent)

plot(bioclim_crop)


# ~*~ Choose which bioclim variables to use in the GEA based on the correlations, biplot, and what you know about the system ~*~

# save the clim data
# make sure the pca variables are in the same order as the bamfile samples!
order_index <- match(bamlist$sample, mice_clim$GSB_ID) # the order w/in match matters!!
mice_clim_ordered <- mice_clim[order_index, ]
#write.table(mice_clim_ordered, "~/Documents/NicoleAdams/pman/gea/gsb_bioclim_ordered.txt", sep = "\t", col.names = T, row.names = F, quote = FALSE)


# If you want to keep bioclim PCs you can do the following
# save PCs (cumulative variance explained for PC1-4=98.2%)
pca.covar <- pca$x[,1:4]
# make sure the pca variables are in the same order as the bamfile samples!
rownames(pca.covar) <- bamlist$sample
order_index2 <- match(bamlist$sample, row.names(pca.covar)) # the order w/in match matters!!
pca.covar_ordered <- pca.covar[order_index2, ]
#write.table(pca.covar_ordered, "~/Documents/NicoleAdams/pman/gea/bioclim_pc.covar.txt", sep = "\t", col.names = F, row.names = F, quote = FALSE)




