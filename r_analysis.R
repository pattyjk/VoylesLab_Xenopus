##Voyles lab X. tropicalis data
setwd('Documents/GitHub/Abby_Xenopus/')

#read in OTU tables normalized by 16S operon count
otu1<-read.delim('asv_table1.txt', header=T, row.names=1)
otu2<-read.delim('asv_table2.txt', header=T, row.names = 1)

# Get all unique sample names
all_samples <- union(colnames(otu1), colnames(otu2))

# Get all unique OTU IDs
all_otus <- union(rownames(otu1), rownames(otu2))

# Create full tables with all OTUs and samples (fill missing with 0)
otu1_full <- matrix(0, nrow = length(all_otus), ncol = length(all_samples),
                    dimnames = list(all_otus, all_samples))
otu2_full <- otu1_full

# Fill in existing values from each OTU table
otu1_full[rownames(otu1), colnames(otu1)] <- as.matrix(otu1)
otu2_full[rownames(otu2), colnames(otu2)] <- as.matrix(otu2)

# Sum the two matrices to merge overlapping samples
otu_merged <- otu1_full + otu2_full

# Convert to data.frame and write to file if needed
otu_merged_df <- as.data.frame(otu_merged)

#write to file
write.table(otu_merged_df, "merged_otu_table.txt", sep = "\t", quote = FALSE)

#read in meta data
meta<-read.delim('abby_map.txt', header=T)

#################
#decontam table
library(phyloseq)
library(decontam)

# Keep only samples that exist in both files
shared.samples <- intersect(colnames(otu_merged_df), meta$SampleID)
otu <- otu[, shared.samples]
meta <- meta[match(shared.samples, meta$SampleID), ]

# Convert to phyloseq
OTU <- otu_table(as.matrix(otu), taxa_are_rows = TRUE)
SAM <- sample_data(meta)
sample_names(SAM) <- meta$SampleID  # ensure they match
sample_names(OTU) <- colnames(otu)
ps <- phyloseq(OTU, SAM)

# Identify negative controls
# Adjust this line if your controls are named differently
sample_data(ps)$is.neg <- sample_data(ps)$Type %in% c("Control", "NegExtraction", "NegPCR")

# Run decontam
contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg")

# Output contaminants
contaminants <- as.data.frame(rownames(contamdf.prev[contamdf.prev$contaminant, ]))
names(contaminants)<-'OTUS'
write.table(contaminants, "contaminants.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#remove contaminants
dim(otu_merged_df)
#[1] 6529  169
otu_no_contam<-otu_merged_df[-which(row.names(otu_merged_df) %in% contaminants$OTUS), ]
dim(otu_no_contam)
#[1] 6165  169

#write to file
write.table(otu_no_contam, 'no_contam_otu_table.txt', row.names = T, sep='\t', quote=F)

#read in plasmid table
library(readr)
plasmid_table_combined <- read_delim("plasmid_table_combined.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
plasmid<-as.data.frame(t(plasmid_table_combined))
plasmid$SampleID<-row.names(plasmid)
names(plasmid)<-c('PlasmidCount', 'SampleID')

#bind plasmid data to meta data
meta2<-merge(meta, plasmid, by.x='SampleID', by.y='SampleID')

#add total reads to meta data
tots_reads<-as.data.frame(colSums(otu_no_contam))
tots_reads$SampleID<-row.names(tots_reads)
names(tots_reads)<-c('TotalReads', 'SampleID')
meta2<-merge(meta2, tots_reads, by='SampleID')
write.table(meta2, 'metadata_withabun.txt',sep='\t', quote=F, row.names=F)

# Total 16S copies of community = 
#16S copies of plasmid spike-in X (1/ % of reads of plasmid in sequenced community [read count / total reads] - 1)
# spike in was 10,000

#calculate estimated total 16S abundance (10,000 spikein)
meta2$TotalAbun<-abs(10000 * (1/( (meta2$PlasmidCount/meta2$TotalReads)-1)))

#calculate true sOTU abundance (sOTU abun/total abun*total read count)
# Keep only samples present in both OTU and metadata
common_samples <- intersect(colnames(otu_no_contam), meta2$SampleID)
otu <- otu_no_contam[, common_samples]
meta2 <- meta2 %>% filter(SampleID %in% common_samples)

# Order metadata to match OTU columns
meta2 <- meta2[match(colnames(otu), meta2$SampleID), ]

# Convert OTU counts to relative abundance (column-wise)
otu_rel <- sweep(otu, 2, colSums(otu), FUN = "/")

# Multiply each sample's relative abundance by its TotalAbun value
otu_scaled <- sweep(otu_rel, 2, meta2$TotalAbun, FUN = "*")

#remove the controls
otu_scaled <- otu[, !(colnames(otu_scaled) %in% c("NegPCR", "NegExtraction2", "NegExtraction"))]

#write to file for future use
write.table(otu_scaled, 'otu_table_norm_abun_no_contam.txt', sep='\t', quote=F, row.names=T)

####################
#Diversity analysis
#read in table for analysis
otu_scaled<-read.delim('otu_table_norm_abun_no_contam.txt', header=T, row.names=1)

#depth
min(colSums(otu_scaled))
#[1] 5420

#rarefy data 
set.seed(515)
xt_rare<-rrarefy(t(round(otu_scaled, 0)), sample=5420)

#calculate PCoA based on BC similarity
ko_pcoa<-capscale(xt_rare  ~ 1, distance='bray')

#pull out x/y coordinates
ko.scores<-scores(ko_pcoa)

#grab only sample coordinates, write to data frame
ko.coords<-as.data.frame(ko.scores$sites)

#create sample names as a column
ko.coords$SampleID<-row.names(ko.coords)

#map back meta data
ko.coords<-merge(ko.coords, meta2, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained for first two axis
100*round(ko_pcoa$CA$eig[1]/sum(ko_pcoa$CA$eig), 3)
#23.1
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#10.3

#plot PCoA
library(ggplot2)
ggplot(ko.coords, aes(MDS1, MDS2, color=Treatment))+
  geom_point(size=2)+
  theme_bw()+
  xlab("PC1- 23.1%")+
  ylab("PC2- 10.3%")

###Calculate alpha diversity
#Richness
xt.alph<-as.data.frame(specnumber(rrarefy(t(otu_scaled), sample = 5420)))
xt.alph$Richness<-as.numeric(xt.alph$`specnumber(rrarefy(t(otu_scaled), sample = 5420))`)
xt.alph$SampleID<-row.names(xt.alph)
xt.alph<-merge(xt.alph, meta, by='SampleID')

#Shannon
xt.alpha2<-as.data.frame(vegan::diversity(rrarefy(t(otu_scaled), sample=5420), index = 'shannon'))
names(xt.alpha2)<-"Shannon"
cort.alph<-cbind(xt.alph, xt.alpha2)

#plot richness
ggplot(xt.alph, aes(AMX_Treatment, Richness, fill=Jliv_Treatment))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("sOTU Richness")
