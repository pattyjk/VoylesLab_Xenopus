##Voyles lab X. tropicalis data
setwd('Documents/GitHub/Abby_Xenopus/')

#read in OTU tables
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

min(colSums(otu_merged_df))
#5476

#rarefy data 
set.seed(515)
xt_rare<-rrarefy(t(otu_merged_df), sample=5476)

#calculate PCoA based on BC similarity
ko_pcoa<-capscale(xt_rare  ~ 1, distance='bray')

#pull out x/y coordinates
ko.scores<-scores(ko_pcoa)

#grab only sample coordinates, write to data frame
ko.coords<-as.data.frame(ko.scores$sites)

#create sample names as a column
ko.coords$SampleID<-row.names(ko.coords)

#map back meta data
ko.coords<-merge(ko.coords, meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained for first two axis
100*round(ko_pcoa$CA$eig[1]/sum(ko_pcoa$CA$eig), 3)
#22.9
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#10.2

#plot PCoA
library(ggplot2)
ggplot(ko.coords, aes(MDS1, MDS2, color=AMX_Treatment, Shape=Jliv_Treatment))+
  geom_point(size=2)+
  theme_bw()+
  xlab("PC1- 22.9%")+
  ylab("PC2- 10.2%")

###Calculate alpha diversity
#Richness
xt.alph<-as.data.frame(specnumber(rrarefy(t(otu_merged_df), sample=5476)))
xt.alph$SampleID<-row.names(xt.alph)
xt.alph<-merge(xt.alph, meta, by='SampleID')
xt.alph$Richness<-as.numeric(xt.alph$`specnumber(rrarefy(t(otu_merged_df), sample = 5476))`)

#Shannon
xt.alpha2<-as.data.frame(vegan::diversity(rrarefy(t(otu_merged_df), sample=5476), index = 'shannon'))
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



