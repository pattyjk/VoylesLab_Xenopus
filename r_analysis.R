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

# Optional: save to file
write.table(otu_merged_df, "merged_otu_table.txt", sep = "\t", quote = FALSE)

#read in meta data
meta<-read.delim('abby_map.txt', header=T)

min(colSums(otu_merged_df))
#5476



