# Create qq and manhattan plots and list of significant SNPs for each env variable using the output from angsd -doAssoc
# 2024
# Nicole E. Adams

## The doAssoc result files can be large so you can remove rows with missing data via command line before putting them into R
#for file in `ls *.lrt*`
#do
#zcat < $file | grep -vE "999.000000"  | gzip > sub.$file
#done


# Load in libraries
library(tidyverse)
library(fastman)

# Load in angsd -doAsso results for all env tests 
lrt.files <- list.files(path="~/Documents/NicoleAdams/pman/gea", pattern='*.lrt*', full.names = TRUE)

# Loop to produce qq and manhattan plots  
lrt.list <- list()
for (FILE in lrt.files){
 lrt <- read.table(gzfile(FILE), header=T, sep="\t")
 
 # get names depending on which environmental database I used 
 if(grepl("PRISM", FILE)==FALSE) {
   envA <- unlist(strsplit(FILE, "[/|,.,_]+"))[c(8,12)]
   env <- paste(envA[1], envA[2], sep = ".")
 }
 
 if(grepl("PRISM", FILE)==TRUE) {
   envA <- unlist(strsplit(FILE, "[/|,.,_]+"))[c(8:9,11)] # for PRISM 
   env <- paste(envA[1], envA[2], envA[3], sep = ".") # for PRISM
 }
 
 # assign names
 dfName <- paste( env, 'df', sep = '.' )
 qqName <- paste( env, 'qq', sep = '.' )
 manName <- paste( env, 'man', sep = '.' )
 
 # filter out sites where the p-value is infinity
 lrt.sub <- lrt[!(-log10(lrt$P) == Inf),]
 lrt.sub <- lrt.sub[!(-log10(lrt.sub$P) == Inf),]
 lrt.sub <- lrt.sub[!(lrt.sub$P == Inf),]
 
 lrt.sub$SNP<-paste("r",1:length(lrt.sub$Chromosome), sep="")
 
 # assign dataframe to list
 lrt.list[[dfName]] <- lrt.sub
 
	# save qq plot to PDF
 pdf(paste0("~/Documents/NicoleAdams/pman/gea/gsb_plots_", env, ".pdf"), width=8, height=10)
    par(mfrow=c(2, 1))
  qqman::qq(lrt.sub$P, main = paste0(qqName))
  dev.off()

	# set up chr values for manhattan plot
 par(mar=c(5.1, 4.1, 4.1, 2.1))
 lrt2plot <- lrt.sub %>% filter(P < 0.05)
 mylabs <- unique(lrt2plot$Chromosome) # preserve the order, otherwise it will order alphabetically
 lrt2plot$CHR <- as.numeric(factor(lrt2plot$Chromosome, levels = mylabs))
 lrt2plot$CHR2 <- paste0("chr", lrt2plot$CHR)
 
	# save manhattan plot as png
  png(paste0("~/Documents/NicoleAdams/pman/gea/gsb_mani_", env, ".png", sep = ""), width=600, height=500, res=120)
  fastman(lrt2plot, chr="CHR", bp="Position", p="P", cex=0.5, main=paste0(manName))
  dev.off()
 
	# identify and save list of significant SNPs at two thresholds
 lrt.sig <- lrt2plot %>% filter(P < 0.00001) # default manhattan() suggestiveline -log10(1e-5)
 lrt.sig.high <- lrt2plot %>% filter(P < 5e-8) # default manhattan() genomewideline -log10(5e-8)

	write.table(lrt.sig, file=paste0("~/Documents/NicoleAdams/pman/gea/gsb_", env, "_sigSNPs.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
	write.table(lrt.sig.high, file=paste0("~/Documents/NicoleAdams/pman/gea/gsb_", env, "_HiSigSNPs.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
 
}