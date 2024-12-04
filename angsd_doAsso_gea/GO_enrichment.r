# Perform GO enrichment on significant SNPs from doAssoc (or outlier tests) using the annotation file created from LD-annot
# 2024
# Nicole E. Adams modified from Rachael Bay's GO.R

# Load in libraries
library(biomaRt)
library(topGO)
library(tidyverse)
library(stringr)
library(vroom)

options(timeout = 30000)

st=format(Sys.time(), "%Y-%m-%d")

# Load in LD-annot output and reformat
ld <- vroom("~/Documents/NicoleAdams/pman/gea/gsb_doAssoc_impute/annot/pman.noLow.dwnSamp_0.9")
ld <- separate(ld,col="annotation",into=c("ID","Dbxref","Name","gbkey","gene","gene_biotype"),sep=";")
ld <- separate(ld,col="gene",into=c("g","gene"),sep="=")
ld$chr <- paste(gsub("chr","NC-",ld$chromosome),".1",sep="")
ld$chr <- gsub("scaff","NW-",ld$chr)
ld$snp <- paste(ld$chr,ld$region_start,sep="_")

#Grab pman annotations from biomaRt
 ensembl = useEnsembl(biomart="ensembl")
 
 pm = useEnsembl(biomart="ensembl", dataset="pmbairdii_gene_ensembl")
 listAttributes(pm)
 pm_genes <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "description",
                                     "chromosome_name", "start_position", "end_position",
                                     "external_gene_name", "gene_biotype", "go_id", "name_1006", "definition_1006",
                                     "namespace_1003","ensembl_transcript_id"),
                      values = T,
                      mart = pm)
 saveRDS(pm_genes,"~/Documents/NicoleAdams/pman/gea/pmangenes.rds")
 
#Grab mus musculus annotations from biomaRt
 mm = useEnsembl(biomart="ensembl",dataset="mmusculus_gene_ensembl")
 listAttributes(mm)
 mm_genes <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "description",
                                     "chromosome_name", "start_position", "end_position",
                                     "external_gene_name", "gene_biotype", "go_id", "name_1006", "definition_1006",
                                     "namespace_1003","ensembl_transcript_id"),
                      values = T,
                      mart = mm)
 saveRDS(mm_genes,"~/Documents/NicoleAdams/pman/gea/musmusgenes.rds")

# Once downloaded, you don't need those biomaRt lines, instead you can call the data locally 
pman_genes <- readRDS("~/Documents/NicoleAdams/pman/gea/pmangenes.rds")
mus_genes <- readRDS("~/Documents/NicoleAdams/pman/gea/musmusgenes.rds")



###Add GO terms
genes <- na.omit(unique(ld$gene))
gene_data <-data.frame(gene=genes,
                       zdesc=pman_genes$description[match(genes,pman_genes$external_gene_name)],
                       cdesc=mus_genes$description[match(genes,mus_genes$external_gene_name)])
gene_data$cdesc[gene_data$cdesc==""] <- NA
length(grep("LOC",gene_data$gene)) # Not really annotated
length(which(!is.na(gene_data$cdesc) | !is.na(gene_data$zdesc))) # Annotated with some gene name

gos <- data.frame(gene_id=genes,desc=NA,GO=NA)
for (g in genes) {
  if(!is.na(gene_data$zdesc[gene_data$gene==g])) {
    sub <- subset(pman_genes,external_gene_name==g)
    gostring <- paste(sub$go_id,collapse=",")
    gos$desc[gos$gene_id==g] <- sub$description[1]
  }
  else if (!is.na(gene_data$cdesc[gene_data$gene==g])) {
    sub <- subset(mus_genes,external_gene_name==g)
    gostring <- paste(sub$go_id,collapse=",")   
    gos$desc[gos$gene_id==g] <- sub$description[1]
  }
  else {gostring <- NA}
  gos$GO[gos$gene_id==g] <- gostring
}

write.table(gos[,c("gene_id","GO")],paste0("~/Documents/NicoleAdams/pman/gea/gsb_doAssoc_impute/annot/pman.noLow.dwnSamp_GOmap_", st, ".txt"),quote=F,row.names=F,col.names=F,sep="\t")





################################### 
# Load in significant SNP files
sig.files <- list.files(path="~/Documents/NicoleAdams/pman/gea/gsb_doAssoc_impute", pattern='gsb.imput.*_sigSNPs.txt', full.names = TRUE)
sig.files <- sig.files[1:3] #keep only CRU files

goList <- c()
for (file in sig.files){
  enamA <- unlist(strsplit(file, "[/|,.,_]+"))[c(13,14)]
  parname <- paste(enamA[1], enamA[2], sep = ".")
  
  sites <- read.delim(file,header=T)
  sites$marker <- paste(sites$Chromosome, sites$Position, sep = "_")
  
  ####Pull out genes linked to significant snps
  candSNP <- sites
  candGenes <- unique(ld$gene[ld$snp%in%candSNP[,15]])
  candGeneDesc <- sapply(candGenes,function(x) gos$desc[gos$gene_id==x])
  candGeneDescShort <- gsub(" \\[.*\\]","",candGeneDesc)
  geneFrame <- data.frame(GeneID=candGenes,Description=candGeneDescShort)
  write.csv(geneFrame,paste("~/Documents/NicoleAdams/pman/gea/gsb_doAssoc_impute/annot/",parname, "_", st,".genelist.csv",sep=""),row.names=F)

  ###For a given column, do GO enrichment for significant SNPs
  geneID2GO <- readMappings(file=paste0("~/Documents/NicoleAdams/pman/gea/gsb_doAssoc_impute/annot/pman.noLow.dwnSamp_GOmap_", st, ".txt"),sep="\t",IDsep=",")
  allgenes <- unique(ld$gene) #This is the 'background'
  myIG <- factor(as.numeric(allgenes%in%candGenes))
  names(myIG) <- allgenes
  
  ##BP
  GOdata <- new("topGOdata",ontology="BP",allGenes=myIG, nodeSize=10,
                annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
  test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
  resultFisher <- getSigGroups(GOdata,test.stat)
  res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
  filt <- res[res$classic<0.05 & res$Significant>=5,]
  write.csv(filt,paste("~/Documents/NicoleAdams/pman/gea/gsb_doAssoc_impute/annot/",parname, "_",st, ".GO_BP.csv",sep=""))
  
  ##MF
  GOdata <- new("topGOdata",ontology="MF",allGenes=myIG, nodeSize=10,
                annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
  test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
  resultFisher <- getSigGroups(GOdata,test.stat)
  res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
  filt <- res[res$classic<0.05 & res$Significant>=5,]
  write.csv(filt,paste("~/Documents/NicoleAdams/pman/gea/gsb_doAssoc_impute/annot/",parname, "_", st, ".GO_MF.csv",sep=""))
  
  ##CC
  GOdata <- new("topGOdata",ontology="CC",allGenes=myIG, nodeSize=10,
                annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
  test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
  resultFisher <- getSigGroups(GOdata,test.stat)
  res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
  filt <- res[res$classic<0.05 & res$Significant>=5,]
  write.csv(filt,paste("~/Documents/NicoleAdams/pman/gea/gsb_doAssoc_impute/annot/",parname, "_", st, ".GO_CC.csv",sep=""))
    
}
