# ============================================================
'R code for Microbiome analysis

Christina Pavloudi
christina.pavloudi@embrc.eu
https://cpavloud.github.io/mysite/

	Copyright (C) 2024 Christina Pavloudi
  
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
  
    This script is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.'

# =============================================================

############################LOAD LIBRARIES #######################################
library(vegan); packageVersion("vegan")
library(ecodist); packageVersion("ecodist")
library(GGally); packageVersion("GGally")
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
library(tidyverse); packageVersion("tidyverse")
library(RColorBrewer); packageVersion("RColorBrewer")
library(DECIPHER); packageVersion("DECIPHER")
library(microViz); packageVersion("microViz")
library(ape); packageVersion("ape")
library(decontam); packageVersion("decontam")
library(DT); packageVersion("DT")
#For the LefSe
library(microbiomeMarker); packageVersion("microbiomeMarker")
library(MicrobiotaProcess); packageVersion("MicrobiotaProcess")
library(patchwork)
#For the upset plot
library(UpSetR); packageVersion("UpSetR")
library(ComplexUpset); packageVersion("ComplexUpset")
# for the kruskal_test and wilcox_test
library(coin)
library(microbial); packageVersion("microbial")

# Define a default theme for ggplot graphics
theme_set(theme_bw()) 

########################################################################################
############################# 16S rRNA amplicons #######################################
########################################################################################

############################ Preparation of Data #######################################
#import the OTU table (or else biotic data)
bac <- read.csv("DADA2_outputs/seq_table.csv", sep = ",", header=TRUE, row.names = 1)
bac <- select(bac, -SC72C1, -SC72C2, -SH72C2)
new_bac <- read.csv("DADA2_new_outputs/seq_table.csv", sep = ",", header=TRUE, row.names = 1)

#convert the rownames into an ASV column
new_bac <- tibble::rownames_to_column(new_bac, "ASV")
#replace the ASV string 
new_bac$ASV <- str_replace(new_bac$ASV, "ASV_", "ASV0_")
#convert the ASV column into the rowname of the data frame
new_bac  <- new_bac  %>% remove_rownames %>% column_to_rownames(var="ASV")

#import the taxonomy table
taxonomybac <- read.csv("DADA2_outputs/seq_Taxonomy_silva_species.csv", sep = ",", header=FALSE, row.names = 1, na.strings=c("","NA"))
new_taxonomybac <- read.csv("DADA2_new_outputs/seq_Taxonomy_silva_species.csv", sep = ",", header=FALSE, row.names = 1, na.strings=c("","NA"))

#convert the rownames into an ASV column
new_taxonomybac <- tibble::rownames_to_column(new_taxonomybac, "ASV")
#replace the ASV string 
new_taxonomybac$ASV <- str_replace(new_taxonomybac$ASV, "ASV_", "ASV0_")
#convert the ASV column into the rowname of the data frame
new_taxonomybac  <- new_taxonomybac  %>% remove_rownames %>% column_to_rownames(var="ASV")

colnames(taxonomybac) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(new_taxonomybac) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#check where there are NA values in the taxonomy data frame
colSums(is.na(taxonomybac))
colSums(is.na(new_taxonomybac))

# You can either remove rows with only NAs
taxonomybac <- taxonomybac[rowSums(is.na(taxonomybac)) != ncol(taxonomybac), ]
new_taxonomybac <- new_taxonomybac[rowSums(is.na(new_taxonomybac)) != ncol(new_taxonomybac), ]
#ORRRRR
# Replace the NA Kingdom values with a specific word/phrase (and maybe remove them in a later stage)
# taxonomybac$Kingdom[is.na(taxonomybac$Kingdom)] <- "Uknown"

taxonomy <- taxonomybac
taxonomy$correct_species <- NA

#fill in correct_species column
  for (i in 1:nrow(taxonomy))
    if (!is.na(taxonomy$Species[i])){
      taxonomy$correct_species[i]=paste(taxonomy$Genus[i], taxonomy$Species[i], sep = " ")
    }

taxonomy <- select(taxonomy, -Species)
colnames(taxonomy)[colnames(taxonomy) == "correct_species"] ="Species"


#fill in Phylum column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Phylum[i])){
      taxonomy$Phylum[i]=paste("Unknown", taxonomy$Kingdom[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Phylum))==0) {
    break
  }
}

#fill in Class column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Class[i])){
      taxonomy$Class[i]=paste("Unknown", taxonomy$Phylum[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Class))==0) {
    break
  }
}

taxonomy$Class <- str_replace(taxonomy$Class, "Unknown Unknown ", "Unknown ")


#fill in Order column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Order[i])){
      taxonomy$Order[i]=paste("Unknown", taxonomy$Class[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Order))==0) {
    break
  }
}

taxonomy$Order <- str_replace(taxonomy$Order, "Unknown Unknown ", "Unknown ")


#fill in Family column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Family[i])){
      taxonomy$Family[i]=paste("Unknown", taxonomy$Order[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Family))==0) {
    break
  }
}

taxonomy$Family <- str_replace(taxonomy$Family, "Unknown Unknown ", "Unknown ")


#fill in Genus column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Genus[i])){
      taxonomy$Genus[i]=paste("Unknown", taxonomy$Family[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Genus))==0) {
    break
  }
}

taxonomy$Genus <- str_replace(taxonomy$Genus, "Unknown Unknown ", "Unknown ")

#fill in Species column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Species[i])){
      taxonomy$Species[i]=paste("Unknown", taxonomy$Genus[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Species))==0) {
    break
  }
}

taxonomy$Species <- str_replace(taxonomy$Species, "Unknown Unknown ", "Unknown ")

taxonomybac <- taxonomy

taxonomybac <- taxonomybac[!(taxonomybac$Order %in% "Chloroplast"),]
taxonomybac <- taxonomybac[!(taxonomybac$Family %in% "Mitochondria"),]

taxonomy <- new_taxonomybac
taxonomy$correct_species <- NA

#fill in correct_species column
for (i in 1:nrow(taxonomy))
  if (!is.na(taxonomy$Species[i])){
    taxonomy$correct_species[i]=paste(taxonomy$Genus[i], taxonomy$Species[i], sep = " ")
  }

taxonomy <- select(taxonomy, -Species)
colnames(taxonomy)[colnames(taxonomy) == "correct_species"] ="Species"


#fill in Phylum column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Phylum[i])){
      taxonomy$Phylum[i]=paste("Unknown", taxonomy$Kingdom[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Phylum))==0) {
    break
  }
}

#fill in Class column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Class[i])){
      taxonomy$Class[i]=paste("Unknown", taxonomy$Phylum[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Class))==0) {
    break
  }
}

taxonomy$Class <- str_replace(taxonomy$Class, "Unknown Unknown ", "Unknown ")


#fill in Order column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Order[i])){
      taxonomy$Order[i]=paste("Unknown", taxonomy$Class[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Order))==0) {
    break
  }
}

taxonomy$Order <- str_replace(taxonomy$Order, "Unknown Unknown ", "Unknown ")


#fill in Family column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Family[i])){
      taxonomy$Family[i]=paste("Unknown", taxonomy$Order[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Family))==0) {
    break
  }
}

taxonomy$Family <- str_replace(taxonomy$Family, "Unknown Unknown ", "Unknown ")


#fill in Genus column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Genus[i])){
      taxonomy$Genus[i]=paste("Unknown", taxonomy$Family[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Genus))==0) {
    break
  }
}

taxonomy$Genus <- str_replace(taxonomy$Genus, "Unknown Unknown ", "Unknown ")

#fill in Species column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Species[i])){
      taxonomy$Species[i]=paste("Unknown", taxonomy$Genus[i], sep = " ")
    }
  if (sum(is.na(taxonomy$Species))==0) {
    break
  }
}

taxonomy$Species <- str_replace(taxonomy$Species, "Unknown Unknown ", "Unknown ")

new_taxonomybac <- taxonomy

new_taxonomybac <- new_taxonomybac[!(new_taxonomybac$Order %in% "Chloroplast"),]
new_taxonomybac <- new_taxonomybac[!(new_taxonomybac$Family %in% "Mitochondria"),]


#check where there are NA values in the taxonomy data frame and in the biotic data
colSums(is.na(taxonomybac))
colSums(is.na(bac))
colSums(is.na(new_taxonomybac))
colSums(is.na(new_bac))

#convert the taxonomy data from data frame to matrix
taxonomy_matrix_bac <- as.matrix(taxonomybac)
taxonomy_matrix_new_bac <- as.matrix(new_taxonomybac)

# prepare the object for the phyloseq object
TAX_BAC = phyloseq::tax_table(taxonomy_matrix_bac)
TAX_NEW_BAC = phyloseq::tax_table(taxonomy_matrix_new_bac)

#convert the biotic data from data frame to matrix
biotic_matrix_bac <- as.matrix(bac)
biotic_matrix_new_bac <- as.matrix(new_bac)

#tranpose biotic data for the calculation of diversity indices
biotic_presabs <- decostand(bac, method = "pa")
biotic_presabs_matrix <- as.matrix(biotic_presabs)
biotic_presabs_new <- decostand(new_bac, method = "pa")
biotic_presabs_new_matrix <- as.matrix(biotic_presabs_new)

#prepare the object for the phyloseq object
OTU_PR = phyloseq::otu_table(biotic_presabs_matrix, taxa_are_rows = TRUE)
OTU_NEW_PR = phyloseq::otu_table(biotic_presabs_new_matrix, taxa_are_rows = TRUE)

OTU_BAC = phyloseq::otu_table(biotic_matrix_bac, taxa_are_rows = TRUE)
OTU_NEW_BAC = phyloseq::otu_table(biotic_matrix_new_bac, taxa_are_rows = TRUE)

#import the metadata of the samples
metadata_physeq <- read.csv("DADA2_outputs/metadata.csv", header=TRUE, row.names = 1) 
# prepare the objects for the phyloseq object
META = phyloseq::sample_data(metadata_physeq)
head(META)

#load the tree file
#load("fitGTR.RData")
#TREE <- phy_tree(fitGTR$tree)
#TREE <- read_tree("seqs_aligned.fa.treefile")
#rename tree tip labels
#TREE$tip.label <- taxa_names(TAX_BAC)

#check if the names are correct and the same
#setequal(taxa_names(TAX_BAC), taxa_names(TREE))

######################## PHYLOSEQ analysis #######################################

# combine them all to create the phyloseq object
physeq_bac = phyloseq(OTU_BAC, TAX_BAC, META)
physeq_bac_PR = phyloseq::phyloseq(OTU_PR, TAX_BAC, META)
physeq_new_bac = phyloseq(OTU_NEW_BAC, TAX_NEW_BAC, META)
physeq_new_bac_PR = phyloseq::phyloseq(OTU_NEW_PR, TAX_NEW_BAC, META)

physeq = merge_phyloseq(physeq_bac, physeq_new_bac)
physeq_PR = merge_phyloseq(physeq_bac_PR, physeq_new_bac_PR)

#Check if there are ASVs with no counts and how many there are
any(taxa_sums(physeq) == 0)
sum(taxa_sums(physeq) == 0)
#Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
physeq <- phyloseq::prune_taxa(taxa_sums(physeq) > 0, physeq)
# Remove samples that are now empty, ie. that have no counts
physeq <- phyloseq::prune_samples(sample_sums(physeq) > 0, physeq)

#Check if there are ASVs with no counts and how many there are
any(taxa_sums(physeq_PR) == 0)
sum(taxa_sums(physeq_PR) == 0)
#Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
physeq_PR <- phyloseq::prune_taxa(taxa_sums(physeq_PR) > 0, physeq_PR)
# Remove samples that are now empty, ie. that have no counts
physeq_PR <- phyloseq::prune_samples(sample_sums(physeq_PR) > 0, physeq_PR)

#physeq_bac_PR = phyloseq::phyloseq(OTU_PR, TAX_BAC, META, TREE)
#physeq_bac = phyloseq::phyloseq(OTU_BAC, TAX_BAC, META, TREE)

#physeq = physeq_bac
#physeq_PR = physeq_bac_PR

#Inspect Library Sizes
df <- as.data.frame(sample_data(physeq)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(physeq)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Treatment)) + geom_point()
ggsave("LibrarySize.png", width = 8, height = 6, dpi = 600)

#merge the OTUs at the Phylum level
physeq_merged_Phylum <- phyloseq::tax_glom(physeq, "Phylum")
ps0 <- transform_sample_counts(physeq_merged_Phylum, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#Remove the OTU column, it makes no sense
OTU_TAX_merged <- select(OTU_TAX_merged, -OTU)
OTU_TAX_merged <- select(OTU_TAX_merged, -Class, -Order, -Family, -Genus, -Species)
#save the final phyto otu table with the taxonomies
write.table(OTU_TAX_merged, "16S_OTU_TAX_Merged_Phylum.tsv", row.names = FALSE)


#merge the OTUs at the Class level
physeq_merged_Class <- phyloseq::tax_glom(physeq, "Class")
ps0 <- transform_sample_counts(physeq_merged_Class, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#Remove the OTU column, it makes no sense
OTU_TAX_merged <- select(OTU_TAX_merged, -OTU)
OTU_TAX_merged <- select(OTU_TAX_merged, -Order, -Family, -Genus, -Species)
#save the final phyto otu table with the taxonomies
write.table(OTU_TAX_merged, "16S_OTU_TAX_Merged_Class.tsv", row.names = FALSE)

#merge the OTUs at the Species level
physeq_merged_Species <- phyloseq::tax_glom(physeq, "Species")
ps0 <- transform_sample_counts(physeq_merged_Species, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#Remove the OTU column, it makes no sense
OTU_TAX_merged <- select(OTU_TAX_merged, -OTU)
#save the final phyto otu table with the taxonomies
write.table(OTU_TAX_merged, "16S_OTU_TAX_Merged_Species.tsv", row.names = FALSE)

#Check if there are ASVs with no counts and how many there are
any(taxa_sums(physeq_merged_Species) == 0)
sum(taxa_sums(physeq_merged_Species) == 0)
#Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
physeq_merged_Species <- phyloseq::prune_taxa(taxa_sums(physeq_merged_Species) > 0, physeq_merged_Species)
# Remove samples that are now empty, ie. that have no counts
physeq_merged_Species <- phyloseq::prune_samples(sample_sums(physeq_merged_Species) > 0, physeq_merged_Species)


#get the data frame from the phyloseq object
pd_species <- psmelt(physeq_merged_Species)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd_species[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)

# Changing order of levels
# You can use this according to your metadata
# Now I am just playing around
pd_species$Treatment <- factor(pd_species$Treatment,      # Reordering group factor levels
                       levels = c("Control", "Infection"))
sample_data(physeq_merged_Species)$Treatment <- factor(sample_data(physeq_merged_Species)$Treatment,
                                     levels =  c("Control", "Infection"), ordered = TRUE)
levels(sample_data(physeq_merged_Species)$Treatment)

pd_species$Nematode <- factor(pd_species$Nematode,      # Reordering group factor levels
                       levels = c("Steinernema carpocapsae", "Steinernema hermaphroditum"))
sample_data(physeq_merged_Species)$Nematode <- factor(sample_data(physeq_merged_Species)$Nematode,
                                        levels =  c("Steinernema carpocapsae", "Steinernema hermaphroditum"), ordered = TRUE)
levels(sample_data(physeq_merged_Species)$Nematode)

pd_species$Timepoint <- factor(pd_species$Timepoint,      # Reordering group factor levels
                      levels = c("24", "48", "72"))
sample_data(physeq_merged_Species)$Timepoint <- factor(sample_data(physeq_merged_Species)$Timepoint,
                                       levels =  c("24", "48", "72"), ordered = TRUE)
levels(sample_data(physeq_merged_Species)$Timepoint)


#and do the actual plotting
barchart_palette <- ggplot(pd_species, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18))
ggsave("phylum_Barchart_Silva.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_species, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Colony) +
  facet_wrap(~Treatment, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Treatment.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_species, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Colony) +
  facet_wrap(~Nematode, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Nematode.png", width = 22, height = 16, dpi = 600)

#plot_bar(physeq, x="Treatment", fill="Phylum") 

plot_bar(physeq_merged_Species, x="Treatment", fill="Phylum", title = "") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=17), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=20), 
        strip.text.x = element_text(size = 17))
ggsave("phylum_Barchart_Treatment_combined.png", width = 10, height = 8, dpi = 600)

plot_bar(physeq_merged_Species, x="Treatment", fill="Phylum", title = "") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  facet_grid(~Nematode, scales="free_x") +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=17), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=20), 
        strip.text.x = element_text(size = 17))
ggsave("phylum_Barchart_Treatment_Nematode_combined.png", width = 12, height = 8, dpi = 600)

# create the nmds plot colour coded by e.g. Type (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq_merged_Species, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq_merged_Species, ord.nmds.bray, color="Treatment", shape = "Nematode") +
  geom_text(aes(label = Timepoint), size = 5, vjust = 2) + 
  geom_point(size=7)
ggsave("16S_P1_Silva_Treatment.png", width = 8, height = 8, dpi = 600)

#PERMANOVA in phyloseq
#This library interferes with phyloseq, so I am loading it now
metadata_permanova <- as(sample_data(physeq_merged_Species), "data.frame")
#permanova.Nematode <- adonis2(distance(physeq, method="bray") ~ Nematode, data = metadata_permanova)
#permanova.Nematode
permanova.Nematode <- microbial::betatest(physeq_merged_Species, group ="Nematode", distance = "bray")
permanova.Nematode
permanova.Treatment <- microbial::betatest(physeq_merged_Species, group ="Treatment", distance = "bray")
permanova.Treatment

dist.bc <- phyloseq::distance(physeq_merged_Species, method = "bray")
permanova <- adonis2(dist.bc ~ Treatment, data = metadata_permanova, perm=999)
permanova

# plot the diversity indices with colour coding by e.g. dpw (info included in the metadata) 
richness <- plot_richness(physeq_merged_Species, measures=c("Observed", "Chao1", "ACE"), color="Treatment")
ggsave("richness.png", width = 10, height = 6, dpi = 600)

alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
p_alpha_meas <- plot_richness(physeq_merged_Species, "Nematode", "Treatment", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas + geom_boxplot(data=p_alpha_meas$data, aes(x=Nematode, y=value, color=NULL), alpha=0.1)
ggsave("richness_boxplot_all.png", width = 10, height = 8, dpi = 600)

# create a heatmap, for the Phylum rank 
#heatmap <- phyloseq::plot_heatmap(physeq, taxa.label="Phylum")
#ggsave("heatmap_phylum.png", width = 8, height = 8, dpi = 600)

# create a better heatmap, for the Phylum rank 
#p <- phyloseq::plot_heatmap(physeq, "NMDS", "bray", "Treatment", "Phylum")
#p$scales$scales[[1]]$name <- "Phylum"
#p$scales$scales[[2]]$name <- "Treatment"
#ggsave("heatmap_phylum_better.png", width = 8, height = 8, dpi = 600)

#reordering labels
#p <- phyloseq::plot_heatmap(physeq, "NMDS", "bray", "Phylum", sample.label="Treatment", sample.order	= "Treatment")
#p$data$Status <- factor(p$data$Status, levels = c("Healthy", "Diseased", "Recovered"), ordered=TRUE)
#levels(p$data$Status) <- c("Healthy", "Diseased", "Recovered")
#p$scales$scales[[1]]$name <- "Phylum"
#p$scales$scales[[2]]$name <- "Treatment"
#ggsave("heatmap_phylum_better_reïrder.png", width = 8, height = 8, dpi = 600)

#Workaround to make a new phyloseq object on the correct order 
#and then re-create the heatmap
#subset the phyloseq object and create a separate object for each level of the factor you want
physeq_carp <- subset_samples(physeq_merged_Species, Nematode=="Steinernema carpocapsae")
physeq_herm <- subset_samples(physeq_merged_Species, Nematode=="Steinernema hermaphroditum")

#Check if there are ASVs with no counts and how many there are
any(taxa_sums(physeq_carp) == 0)
sum(taxa_sums(physeq_carp) == 0)
#Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
physeq_carp <- phyloseq::prune_taxa(taxa_sums(physeq_carp) > 0, physeq_carp)
# Remove samples that are now empty, ie. that have no counts
physeq_carp <- phyloseq::prune_samples(sample_sums(physeq_carp) > 0, physeq_carp)

#Check if there are ASVs with no counts and how many there are
any(taxa_sums(physeq_herm) == 0)
sum(taxa_sums(physeq_herm) == 0)
#Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
physeq_herm <- phyloseq::prune_taxa(taxa_sums(physeq_herm) > 0, physeq_herm)
# Remove samples that are now empty, ie. that have no counts
physeq_herm <- phyloseq::prune_samples(sample_sums(physeq_herm) > 0, physeq_herm)



#merge the phyloseq objects in the order you want
#physeq_reordered <- merge_phyloseq(physeq_carp, physeq_herm)
#p <- phyloseq::plot_heatmap(physeq_reordered,  "NMDS", "bray", "Timepoint", "Phylum")
#p$scales$scales[[1]]$name <- "Phylum"
#p$scales$scales[[2]]$name <- "Timepoint"
#ggsave("heatmap_phylum_better_oncemore.png", width = 8, height = 8, dpi = 600)

#plot_net(physeq, distance = "bray", maxdist=0.4, 
#         type = "samples", color="Treatment", shape="Nematode")
#ggsave("network.png", width = 8, height = 8, dpi = 600)

#network <- make_network(physeq, type="samples", distance="bray", max.dist = 0.4, 
#             keep.isolates=TRUE)
#plot_network(network, physeq, type="samples", 
#             color="Treatment", shape="Timepoint")

#NMDS on Weighted Unifrac distances
#this needs a phylogenetic tree
#ordu_uni = ordinate(physeq, "NMDS", "unifrac", weighted=TRUE)
#p_uni <- plot_ordination(physeq, ordu_uni, color="Treatment", shape = "Nematode") +
#  geom_point(size=5, alpha=0.75) +
#  scale_colour_brewer(type="qual", palette="Set1")+ 
#  ggtitle("MDS on weighted-UniFrac distance")
#ggsave("16S_UNI_WEIGHTED_Silva_Status.png", width = 8, height = 8, dpi = 600)

#add ellipses
#p_uni + 
#  stat_ellipse(type = "norm", linetype = 2) +
#  stat_ellipse(type = "t") +
#  theme_bw()
#ggsave("16S_UNI_WEIGHTED_Silva_Status_ellipses.png", width = 8, height = 8, dpi = 600)

#################### Proteobacteria #############################################

physeq_Proteobacteria <- subset_taxa(physeq_merged_Species, Phylum =="Proteobacteria")

#get the data frame from the phyloseq object
pd_Proteobacteria <- psmelt(physeq_Proteobacteria)

#Count how many Phyla are there in your samples
HowManyClasses <- length(unique(unlist(pd_Proteobacteria[,c("Class")])))
HowManySpecies <- length(unique(unlist(pd_Proteobacteria[,c("Species")])))


# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
ClassesPalette = getPalette(HowManyClasses)
SpeciesPalette = getPalette(HowManySpecies)

#and do the actual plotting
barchart_palette <- ggplot(pd_Proteobacteria, aes(x = Sample, y = Abundance, factor(Class), fill = factor(Class))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ClassesPalette) +
  labs(fill = "Class") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18))
ggsave("phylum_Barchart_Silva_Proteo.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_Proteobacteria, aes(x = Sample, y = Abundance, factor(Class), fill = factor(Class))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ClassesPalette) +
  labs(fill = "Class") +
  #facet_grid(~Colony) +
  facet_wrap(~Treatment, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Treatment_Proteo.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_Proteobacteria, aes(x = Sample, y = Abundance, factor(Species), fill = factor(Species))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = SpeciesPalette) +
  labs(fill = "Species") +
  #facet_grid(~Colony) +
  facet_wrap(~Treatment, scales="free_x") +
  guides(fill=guide_legend(ncol=2)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Treatment_Proteo_Species.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_Proteobacteria, aes(x = Sample, y = Abundance, factor(Class), fill = factor(Class))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ClassesPalette) +
  labs(fill = "Class") +
  #facet_grid(~Colony) +
  facet_wrap(~Nematode, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Nematode_Proteo.png", width = 22, height = 16, dpi = 600)

#merge the OTUs at the Species level
physeq_merged_proteo_Species <- phyloseq::tax_glom(physeq_Proteobacteria, "Species")
ps0 <- transform_sample_counts(physeq_merged_proteo_Species, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#Remove the OTU column, it makes no sense
OTU_TAX_merged <- select(OTU_TAX_merged, -OTU)
#save the final phyto otu table with the taxonomies
write.table(OTU_TAX_merged, "16S_OTU_TAX_Merged_Proteo_Species.tsv", row.names = FALSE)


#################### Bacteroidota #############################################

physeq_Bacteroidota <- subset_taxa(physeq_merged_Species, Phylum =="Bacteroidota")

#get the data frame from the phyloseq object
pd_Bacteroidota <- psmelt(physeq_Bacteroidota)

#Count how many Phyla are there in your samples
HowManyClasses <- length(unique(unlist(pd_Bacteroidota[,c("Class")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
ClassesPalette = getPalette(HowManyClasses)

#and do the actual plotting
barchart_palette <- ggplot(pd_Bacteroidota, aes(x = Sample, y = Abundance, factor(Class), fill = factor(Class))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ClassesPalette) +
  labs(fill = "Class") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18))
ggsave("phylum_Barchart_Silva_Bacteroidota.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_Bacteroidota, aes(x = Sample, y = Abundance, factor(Class), fill = factor(Class))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ClassesPalette) +
  labs(fill = "Class") +
  #facet_grid(~Colony) +
  facet_wrap(~Treatment, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Treatment_Bacteroidota.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_Bacteroidota, aes(x = Sample, y = Abundance, factor(Class), fill = factor(Class))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ClassesPalette) +
  labs(fill = "Class") +
  #facet_grid(~Colony) +
  facet_wrap(~Nematode, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Nematode_Bacteroidota.png", width = 22, height = 16, dpi = 600)

#################### Firmicutes #############################################

physeq_Firmicutes <- subset_taxa(physeq_merged_Species, Phylum =="Firmicutes")

#get the data frame from the phyloseq object
pd_Firmicutes <- psmelt(physeq_Firmicutes)

#Count how many Phyla are there in your samples
HowManyClasses <- length(unique(unlist(pd_Firmicutes[,c("Class")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
ClassesPalette = getPalette(HowManyClasses)

#and do the actual plotting
barchart_palette <- ggplot(pd_Firmicutes, aes(x = Sample, y = Abundance, factor(Class), fill = factor(Class))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ClassesPalette) +
  labs(fill = "Class") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18))
ggsave("phylum_Barchart_Silva_Firmicutes.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_Firmicutes, aes(x = Sample, y = Abundance, factor(Class), fill = factor(Class))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ClassesPalette) +
  labs(fill = "Class") +
  #facet_grid(~Colony) +
  facet_wrap(~Treatment, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Treatment_Firmicutes.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_Firmicutes, aes(x = Sample, y = Abundance, factor(Class), fill = factor(Class))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = ClassesPalette) +
  labs(fill = "Class") +
  #facet_grid(~Colony) +
  facet_wrap(~Nematode, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Nematode_Firmicutes.png", width = 22, height = 16, dpi = 600)

###############################################################################


#get the data frame from the phyloseq object
pd_carp <- psmelt(physeq_carp)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd_carp[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)

# Changing order of levels
# You can use this according to your metadata
# Now I am just playing around
pd_carp$Treatment <- factor(pd_carp$Treatment,      # Reordering group factor levels
                            levels = c("Control", "Infection"))
sample_data(physeq_carp)$Treatment <- factor(sample_data(physeq_carp)$Treatment,
                                             levels =  c("Control", "Infection"), ordered = TRUE)
levels(sample_data(physeq_carp)$Treatment)

pd_carp$Timepoint <- factor(pd_carp$Timepoint,      # Reordering group factor levels
                            levels = c("24", "48", "72"))
sample_data(physeq_carp)$Timepoint <- factor(sample_data(physeq_carp)$Timepoint,
                                             levels =  c("24", "48", "72"), ordered = TRUE)
levels(sample_data(physeq_carp)$Timepoint)


#and do the actual plotting
barchart_palette <- ggplot(pd_carp, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18))
ggsave("phylum_Barchart_Silva_carp.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_carp, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Colony) +
  facet_wrap(~Treatment, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Treatment_carp.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_carp, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Colony) +
  facet_wrap(~Timepoint, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Timepoint_carp.png", width = 22, height = 16, dpi = 600)

plot_bar(physeq_carp, x="Treatment", fill="Phylum", title = "") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=17), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=20), 
        strip.text.x = element_text(size = 17))
ggsave("phylum_Barchart_Treatment_combined_carp.png", width = 10, height = 8, dpi = 600)

plot_bar(physeq_carp, x="Treatment", fill="Phylum", title = "") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  facet_grid(~Timepoint, scales="free_x") +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=17), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=20), 
        strip.text.x = element_text(size = 17))
ggsave("phylum_Barchart_Treatment_Timepoint_combined_carp.png", width = 10, height = 8, dpi = 600)

plot_bar(physeq_carp, x="Timepoint", fill="Phylum", title = "") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  facet_grid(~Treatment, scales="free_x") +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=17), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=20), 
        strip.text.x = element_text(size = 17))
ggsave("phylum_Barchart_Treatment_Timepoint_combined_carp_2.png", width = 10, height = 8, dpi = 600)

# create the nmds plot colour coded by e.g. Type (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq_carp, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq_carp, ord.nmds.bray, color="Treatment", shape = "Timepoint") +
  #geom_text(aes(label = Timepoint), size = 5, vjust = 2) + 
  geom_point(size=7)
ggsave("16S_P1_Silva_Treatment_carp.png", width = 8, height = 8, dpi = 600)

#PERMANOVA in phyloseq
#This library interferes with phyloseq, so I am loading it now
metadata_permanova <- as(sample_data(physeq_carp), "data.frame")
#permanova.Nematode <- adonis2(distance(physeq_carp, method="bray") ~ Nematode, data = metadata_permanova)
#permanova.Nematode
permanova.Timepoint <- microbial::betatest(physeq_carp, group ="Timepoint", distance = "bray")
permanova.Timepoint
permanova.Treatment <- microbial::betatest(physeq_carp, group ="Treatment", distance = "bray")
permanova.Treatment


# plot the diversity indices with colour coding by e.g. dpw (info included in the metadata) 
richness <- plot_richness(physeq_carp, measures=c("Observed", "Chao1", "ACE"), color="Treatment")
ggsave("richness_carp.png", width = 10, height = 6, dpi = 600)

alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
p_alpha_meas <- plot_richness(physeq_carp, "Timepoint" , "Treatment", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas + geom_boxplot(data=p_alpha_meas$data, aes(x=Timepoint, y=value, color=NULL), alpha=0.1)
ggsave("richness_boxplot_all_carp.png", width = 10, height = 8, dpi = 600)

##############################################################################


#get the data frame from the phyloseq object
pd_herm <- psmelt(physeq_herm)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd_herm[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)

# Changing order of levels
# You can use this according to your metadata
# Now I am just playing around
pd_herm$Treatment <- factor(pd_herm$Treatment,      # Reordering group factor levels
                            levels = c("Control", "Infection"))
sample_data(physeq_herm)$Treatment <- factor(sample_data(physeq_herm)$Treatment,
                                             levels =  c("Control", "Infection"), ordered = TRUE)
levels(sample_data(physeq_herm)$Treatment)

pd_herm$Timepoint <- factor(pd_herm$Timepoint,      # Reordering group factor levels
                            levels = c("24", "48", "72"))
sample_data(physeq_herm)$Timepoint <- factor(sample_data(physeq_herm)$Timepoint,
                                             levels =  c("24", "48", "72"), ordered = TRUE)
levels(sample_data(physeq_herm)$Timepoint)


#and do the actual plotting
barchart_palette <- ggplot(pd_herm, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Year, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18))
ggsave("phylum_Barchart_Silva_herm.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_herm, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Colony) +
  facet_wrap(~Treatment, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Treatment_herm.png", width = 22, height = 16, dpi = 600)

#and do the actual plotting
barchart_palette <- ggplot(pd_herm, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  #facet_grid(~Colony) +
  facet_wrap(~Timepoint, scales="free_x") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=18), 
        strip.text.x = element_text(size = 18))
ggsave("phylum_Barchart_Silva_Timepoint_herm.png", width = 22, height = 16, dpi = 600)

plot_bar(physeq_herm, x="Treatment", fill="Phylum", title = "") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=17), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=20), 
        strip.text.x = element_text(size = 17))
ggsave("phylum_Barchart_Treatment_combined_herm.png", width = 10, height = 8, dpi = 600)

plot_bar(physeq_herm, x="Treatment", fill="Phylum", title = "") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  facet_grid(~Timepoint, scales="free_x") +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=17), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=20), 
        strip.text.x = element_text(size = 17))
ggsave("phylum_Barchart_Treatment_Timepoint_combined_herm.png", width = 10, height = 8, dpi = 600)

plot_bar(physeq_herm, x="Timepoint", fill="Phylum", title = "") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  facet_grid(~Treatment, scales="free_x") +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.text=element_text(size=17), axis.title=element_text(size=20,face="bold"),
        legend.title=element_text(size=20), legend.text=element_text(size=20), 
        strip.text.x = element_text(size = 17))
ggsave("phylum_Barchart_Treatment_Timepoint_combined_herm_2.png", width = 10, height = 8, dpi = 600)

# create the nmds plot colour coded by e.g. Type (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq_herm, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq_herm, ord.nmds.bray, color="Treatment", shape = "Timepoint") +
  #geom_text(aes(label = Timepoint), size = 5, vjust = 2) + 
  geom_point(size=7)
ggsave("16S_P1_Silva_Treatment_herm.png", width = 8, height = 8, dpi = 600)

#PERMANOVA in phyloseq
#This library interferes with phyloseq, so I am loading it now
metadata_permanova <- as(sample_data(physeq_herm), "data.frame")
#permanova.Nematode <- adonis2(distance(physeq_herm, method="bray") ~ Nematode, data = metadata_permanova)
#permanova.Nematode
permanova.Timepoint <- microbial::betatest(physeq_herm, group ="Timepoint", distance = "bray")
permanova.Timepoint
permanova.Treatment <- microbial::betatest(physeq_herm, group ="Treatment", distance = "bray")
permanova.Treatment


# plot the diversity indices with colour coding by e.g. dpw (info included in the metadata) 
richness <- plot_richness(physeq_herm, measures=c("Observed", "Chao1", "ACE"), color="Treatment")
ggsave("richness_herm.png", width = 10, height = 6, dpi = 600)

alpha_meas = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson")
p_alpha_meas <- plot_richness(physeq_herm, "Timepoint" , "Treatment", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas + geom_boxplot(data=p_alpha_meas$data, aes(x=Timepoint, y=value, color=NULL), alpha=0.1)
ggsave("richness_boxplot_all_herm.png", width = 10, height = 8, dpi = 600)

##############################################################################

########### Baloon plot
#create a column in the taxonomy table and in the bac table with the ASVs
#taxonomy_asvs <- tibble::rownames_to_column(taxonomybac, "ASVs")
#bac_asvs <- tibble::rownames_to_column(bac, "ASVs")

#select from the bac_asvs table the ASVs that are in the taxonomy_bac table
#bac_asvs <- semi_join(bac_asvs %>% mutate, taxonomy_asvs, by = 'ASVs')

#convert the ASVs column into the rownames of the data frames
#taxonomy_asvs <- taxonomy_asvs %>% remove_rownames %>% column_to_rownames(var="ASVs")
#bac_asvs <- bac_asvs %>% remove_rownames %>% column_to_rownames(var="ASVs")

# merge data frames 
#extended_otu_table <- cbind(bac_asvs, taxonomy_asvs$Phylum, taxonomy_asvs$Class) 

#delete the phyla with empty rows
#extended_otu_table_non_na <- extended_otu_table %>% drop_na()

# rename the column name of the phylum 
#names(extended_otu_table_non_na)[names(extended_otu_table_non_na) == 'taxonomy_asvs$Phylum'] <- "Phylum"
#names(extended_otu_table_non_na)[names(extended_otu_table_non_na) == 'taxonomy_asvs$Class'] <- "Class"

# sum entries of each phylum
#only_phyla <- select(extended_otu_table_non_na, -Class)
#phylum_sum <- only_phyla %>% group_by(Phylum) %>% summarise_all(sum)

# make balloon plot 
#phylum_sum %>% gather(Variable, Value, -Phylum) %>%
#  ggplot(aes(x = Variable, y = Phylum, size = Value)) +
#  geom_point(colour = "sky blue") +
#  scale_size_continuous(range = c(-1, 10)) +
#  labs(x = "",
#       y = "") +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
#        axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
#        legend.title=element_text(size=20), legend.text=element_text(size=18))
#ggsave("ballot_plot.png", width = 12, height = 10, dpi = 600)


#phylum_sum %>% gather(Variable, Value, -Phylum) %>%
#  ggplot(aes(x = Variable, y = Phylum, size = Value)) +
#  geom_point(colour = "sky blue") +
#  scale_size_continuous(range = c(-1, 10)) +
#  labs(title = "Balloon Plot for x by y",
#       subtitle = "Area is proportional to Freq.",
#       x = "",
#       y = "") +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
#        legend.title=element_text(size=20), legend.text=element_text(size=18))
#ggsave("ballot_plot_2.png", width = 12, height = 10, dpi = 600)

##### Now, you can select specific Phyla of interest and create subsets of the 
##### physeq object and repeat the analyses

######################### LefSe etc. ###########################################

#physeq_carp_decont <- subset_samples(physeq_carp, Code != "SC72C1" | Code != "SC72C2")
#physeq_herm_decont<- subset_samples(physeq_herm, Code != "SH72C2")


#Check if there are ASVs with no counts and how many there are
#any(taxa_sums(physeq_carp_decont) == 0)
#sum(taxa_sums(physeq_carp_decont) == 0)
#Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
#physeq_carp_decont <- phyloseq::prune_taxa(taxa_sums(physeq_carp_decont) > 0, physeq_carp_decont)
# Remove samples that are now empty, ie. that have no counts
#physeq_carp_decont <- phyloseq::prune_samples(sample_sums(physeq_carp_decont) > 0, physeq_carp_decont)


#Check if there are ASVs with no counts and how many there are
#any(taxa_sums(physeq_herm_decont) == 0)
#sum(taxa_sums(physeq_herm_decont) == 0)
#Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
#physeq_herm_decont <- phyloseq::prune_taxa(taxa_sums(physeq_herm_decont) > 0, physeq_herm_decont)
# Remove samples that are now empty, ie. that have no counts
#physeq_herm_decont <- phyloseq::prune_samples(sample_sums(physeq_herm_decont) > 0, physeq_herm_decont)

###############################################################################

######### You will need to make modifications below, and choose the correct factor for your data

# #LefSe on ASV level
# OTU_lefse <- run_lefse(
#   physeq_carp_decont,
#   group = "Treatment",
#   subgroup = "Timepoint",
#   taxa_rank = "none", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
#   transform = c("identity"), #return the original data without any transformation (default).
#   norm = "CPM", # do not normalize.
#   kw_cutoff = 0.05,
#   lda_cutoff = 2,
#   bootstrap_n = 30,
#   bootstrap_fraction = 2/3,
#   wilcoxon_cutoff = 0.05,
#   multigrp_strat = FALSE,
#   strict = "0",
#   sample_min = 7,
#   only_same_subgrp = FALSE,
#   curv = FALSE
# )
# 
# markertable_status_type <- marker_table(OTU_lefse)
# markertable_status_type <- as.data.frame(markertable_status_type)
# 
# OTU_lefse
# head(marker_table(OTU_lefse))
# write.csv(marker_table(OTU_lefse), "lefse_status_carp_treatment_timepoint.csv")
# 
# p_abd <- plot_abundance(OTU_lefse, group = "Treatment")
# # customize the plot with ggplot2, modify the fill color manually
# p_abd + scale_fill_manual(values = c("Control" = "green",
#                                      "Infection" = "pink"))
# ggsave("lefse_abundance_status_carp_treatment_timepoint.png", width = 8, height = 8, dpi = 600)
# 
# 
# heatmap <- plot_heatmap(OTU_lefse, transform = "log10p", group = "Treatment", 
#                         sample_label = TRUE, annotation_col = c("Control" = "green",
#                                                                 "Infection" = "pink"))
# heatmap
# #ggsave("lefse_heatmap_status_type.png", width = 8, height = 8, dpi = 600)
# 
# 
# #better heatmap
# heatmap <- plot_heatmap(OTU_lefse, transform = "log10p", group = "Treatment", 
#                         cluster_marker = TRUE, row_dend_side = "right", 
#                         #cluster_sample = TRUE, column_dend_side = "bottom",
#                         column_order = c("D1", "D2", "D3", "D4", "R1", "R2", "R3", "R4", "H1", "H2", "H3", "H4", 
#                                          "SW_D1", "SW_D2", "SW_R1", "SW_R2", "SW_H1", "SW_H2"), 
#                         sample_label = TRUE, annotation_col = c("Control" = "green",
#                                                                 "Infection" = "pink"))
# heatmap
# #ggsave("lefse_heatmap_status_type_better.png", width = 8, height = 8, dpi = 600)
# 
# # bar plot
# plot_ef_bar(OTU_lefse) + scale_fill_manual(values = c("Control" = "green",
#                                                       "Infection" = "pink"))
# ggsave("lefse_bar_status_type.png", width = 8, height = 8, dpi = 600)
# 
# # dot plot
# p <- plot_ef_dot(OTU_lefse) 
# p + scale_fill_manual(values = c("Control" = "green",
#                                  "Infection" = "pink"))
# ggsave("lefse_dot_status_type.png", width = 8, height = 8, dpi = 600)
# 
# ##############
# 
# #LefSe on ASV level only group
# OTU_lefse_group <- run_lefse(
#   physeq_carp_decont,
#   group = "Treatment",
#   #subgroup = "Sample_type",
#   taxa_rank = "none", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
#   transform = c("identity"), #return the original data without any transformation (default).
#   norm = "CPM", # do not normalize.
#   kw_cutoff = 0.05,
#   lda_cutoff = 2,
#   bootstrap_n = 30,
#   bootstrap_fraction = 2/3,
#   wilcoxon_cutoff = 0.05,
#   multigrp_strat = FALSE,
#   strict = "0",
#   sample_min = 7,
#   only_same_subgrp = FALSE,
#   curv = FALSE
# )
# 
# markertable_status <- marker_table(OTU_lefse_group)
# markertable_status <- as.data.frame(markertable_status)
# 
# OTU_lefse_group
# head(marker_table(OTU_lefse_group))
# write.csv(marker_table(OTU_lefse_group), "lefse_status.csv")
# 
# p_abd <- plot_abundance(OTU_lefse_group, group = "Treatment")
# # customize the plot with ggplot2, modify the fill color manually
# p_abd + scale_fill_manual(values = c("Control" = "green",
#                                      "Infection" = "pink"))
# ggsave("lefse_abundance_status.png", width = 8, height = 8, dpi = 600)
# 
# 
# heatmap <- plot_heatmap(OTU_lefse_group, transform = "log10p", group = "Treatment", 
#                         sample_label = TRUE, annotation_col = c("Control" = "green",
#                                                                 "Infection" = "pink"))
# ggsave("lefse_heatmap_status.png", width = 8, height = 8, dpi = 600)
# 
# # bar plot
# plot_ef_bar(OTU_lefse_group) + scale_fill_manual(values = c("Control" = "green",
#                                                             "Infection" = "pink"))
# ggsave("lefse_bar_status.png", width = 8, height = 8, dpi = 600)
# 
# # dot plot
# p <- plot_ef_dot(OTU_lefse_group) 
# p + scale_fill_manual(values = c("Control" = "green",
#                                  "Infection" = "pink"))
# ggsave("lefse_dot_status.png", width = 8, height = 8, dpi = 600)

###############################################################################

#LefSe on Species level
OTU_lefse_species <- run_lefse(
  physeq_carp,
  group = "Treatment",
  subgroup = "Timepoint",
  taxa_rank = "Species", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
  transform = c("identity"), #return the original data without any transformation (default).
  norm = "CPM", # do not normalize.
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = "0",
  sample_min = 7,
  only_same_subgrp = FALSE,
  curv = FALSE
)

markertable_species <- marker_table(OTU_lefse_species)
markertable_species <- as.data.frame(markertable_species)

OTU_lefse_species
head(marker_table(OTU_lefse_species))
write.csv(marker_table(OTU_lefse_species), "lefse_species_carp_Treatment_Timepoint.csv")

p_abd <- plot_abundance(OTU_lefse_species, group = "Treatment")
# customize the plot with ggplot2, modify the fill color manually
p_abd + scale_fill_manual(values = c("Control" = "green",
                                     "Infection" = "pink"))
ggsave("lefse_abundance_species_carp_Treatment_Timepoint.png", width = 8, height = 6, dpi = 600)


heatmap <- plot_heatmap(OTU_lefse_species, transform = "log10p", group = "Treatment", 
                        sample_label = TRUE, annotation_col = c("Control" = "green",
                                                                "Infection" = "pink"))
heatmap
ggsave("lefse_heatmap_species_carp_Treatment_Timepoint.png", width = 8, height = 8, dpi = 600)

# bar plot
plot_ef_bar(OTU_lefse_species) + scale_fill_manual(values = c("Control" = "green",
                                                              "Infection" = "pink"))
ggsave("lefse_bar_species_carp_Treatment_Timepoint.png", width = 8, height = 6, dpi = 600)

# dot plot
p <- plot_ef_dot(OTU_lefse_species) 
p + scale_fill_manual(values = c("Control" = "green",
                                 "Infection" = "pink"))
ggsave("lefse_dot_species_carp_Treatment_Timepoint.png", width = 8, height = 6, dpi = 600)

##########

#LefSe on Species level only group
OTU_lefse_group_species <- run_lefse(
  physeq_carp,
  group = "Treatment",
  #subgroup = "Sample_type",
  taxa_rank = "Species", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
  transform = c("identity"), #return the original data without any transformation (default).
  norm = "CPM", # do not normalize.
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = "0",
  sample_min = 7,
  only_same_subgrp = FALSE,
  curv = FALSE
)

markertable_status_species <- marker_table(OTU_lefse_group_species)
markertable_status_species <- as.data.frame(markertable_status_species)

OTU_lefse_group_species
head(marker_table(OTU_lefse_group_species))
write.csv(marker_table(OTU_lefse_group_species), "lefse_species_carp_Treatment.csv")

p_abd <- plot_abundance(OTU_lefse_group_species, group = "Treatment")
# customize the plot with ggplot2, modify the fill color manually
p_abd + scale_fill_manual(values = c("Control" = "green",
                                     "Infection" = "pink"))
ggsave("lefse_abundance_species_carp_Treatment.png", width = 8, height = 8, dpi = 600)


heatmap <- plot_heatmap(OTU_lefse_group_species, transform = "log10p", group = "Treatment", 
                        sample_label = TRUE, annotation_col = c("Control" = "green",
                                                                "Infection" = "pink"))
heatmap
ggsave("lefse_heatmap_species_carp_Treatment.png", width = 8, height = 12, dpi = 600)

# bar plot
plot_ef_bar(OTU_lefse_group_species) + scale_fill_manual(values = c("Control" = "green",
                                                                    "Infection" = "pink"))
ggsave("lefse_bar_species_carp_Treatment.png", width = 8, height = 8, dpi = 600)

# dot plot
p <- plot_ef_dot(OTU_lefse_group_species) 
p + scale_fill_manual(values = c("Control" = "green",
                                 "Infection" = "pink"))
ggsave("lefse_dot_species_carp_Treatment.png", width = 8, height = 8, dpi = 600)

##############


#LefSe on Species level
OTU_lefse_species <- run_lefse(
  physeq_herm,
  group = "Treatment",
  subgroup = "Timepoint",
  taxa_rank = "Species", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
  transform = c("identity"), #return the original data without any transformation (default).
  norm = "CPM", # do not normalize.
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = "0",
  sample_min = 7,
  only_same_subgrp = FALSE,
  curv = FALSE
)

markertable_species <- marker_table(OTU_lefse_species)
markertable_species <- as.data.frame(markertable_species)

OTU_lefse_species
head(marker_table(OTU_lefse_species))
write.csv(marker_table(OTU_lefse_species), "lefse_species_herm_Treatment_Timepoint.csv")

p_abd <- plot_abundance(OTU_lefse_species, group = "Treatment")
# customize the plot with ggplot2, modify the fill color manually
p_abd + scale_fill_manual(values = c("Control" = "green",
                                     "Infection" = "pink"))
ggsave("lefse_abundance_species_herm_Treatment_Timepoint.png", width = 8, height = 8, dpi = 600)


heatmap <- plot_heatmap(OTU_lefse_species, transform = "log10p", group = "Treatment", 
                        sample_label = TRUE, annotation_col = c("Control" = "green",
                                                                "Infection" = "pink"))
heatmap
ggsave("lefse_heatmap_species_herm_Treatment_Timepoint.png", width = 8, height = 8, dpi = 600)

# bar plot
plot_ef_bar(OTU_lefse_species) + scale_fill_manual(values = c("Control" = "green",
                                                              "Infection" = "pink"))
ggsave("lefse_bar_species_herm_Treatment_Timepoint.png", width = 8, height = 8, dpi = 600)

# dot plot
p <- plot_ef_dot(OTU_lefse_species) 
p + scale_fill_manual(values = c("Control" = "green",
                                 "Infection" = "pink"))
ggsave("lefse_dot_species_herm_Treatment_Timepoint.png", width = 8, height = 8, dpi = 600)

##########

#LefSe on Species level only group
OTU_lefse_group_species <- run_lefse(
  physeq_herm,
  group = "Treatment",
  #subgroup = "Sample_type",
  taxa_rank = "Species", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
  transform = c("identity"), #return the original data without any transformation (default).
  norm = "CPM", # do not normalize.
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = "0",
  sample_min = 7,
  only_same_subgrp = FALSE,
  curv = FALSE
)

markertable_status_species <- marker_table(OTU_lefse_group_species)
markertable_status_species <- as.data.frame(markertable_status_species)

OTU_lefse_group_species
head(marker_table(OTU_lefse_group_species))
write.csv(marker_table(OTU_lefse_group_species), "lefse_species_herm_Treatment.csv")

p_abd <- plot_abundance(OTU_lefse_group_species, group = "Treatment")
# customize the plot with ggplot2, modify the fill color manually
p_abd + scale_fill_manual(values = c("Control" = "green",
                                     "Infection" = "pink"))
ggsave("lefse_abundance_species_herm_Treatment.png", width = 8, height = 8, dpi = 600)


heatmap <- plot_heatmap(OTU_lefse_group_species, transform = "log10p", group = "Treatment", 
                        sample_label = TRUE, annotation_col = c("Control" = "green",
                                                                "Infection" = "pink"))
heatmap
ggsave("lefse_heatmap_species_herm_Treatment.png", width = 8, height = 12, dpi = 600)

# bar plot
plot_ef_bar(OTU_lefse_group_species) + scale_fill_manual(values = c("Control" = "green",
                                                                    "Infection" = "pink"))
ggsave("lefse_bar_species_herm_Treatment.png", width = 8, height = 8, dpi = 600)

# dot plot
p <- plot_ef_dot(OTU_lefse_group_species) 
p + scale_fill_manual(values = c("Control" = "green",
                                 "Infection" = "pink"))
ggsave("lefse_dot_species_herm_Treatment.png", width = 8, height = 8, dpi = 600)

#############################################################

# #Lefse another way
# deres <- diff_analysis(obj = physeq_decont, classgroup = "Colony", 
#                        #subclass = "Sample_type",
#                        mlfun = "lda",
#                        filtermod = "pvalue",
#                        firstcomfun = "kruskal_test",
#                        firstalpha = 0.05,
#                        strictmod = TRUE,
#                        secondcomfun = "wilcox_test",
#                        subclmin = 3,
#                        subclwilc = TRUE,
#                        secondalpha = 0.01,
#                        lda=2)
# 
# cols = c("Healthy" = "green", "Recovered" = "blue", "Diseased" = "pink")
# 
# #visualization of differential results (with action = “get”) by ggdiffbox
# diffbox <- ggdiffbox(obj=deres, box_notch=FALSE,
#                      colorlist=cols, l_xlabtext="relative abundance")
# diffbox
# ggsave("diffbox.png", width = 8, height = 12, dpi = 600)
# 
# 
# #visualization of different results by ggdiffclade
# diffclade_p <- ggdiffclade(
#   obj=deres,
#   alpha=0.3,
#   linewd=0.15,
#   skpointsize=0.6,
#   layout="radial",
#   taxlevel=3,
#   removeUnknown = TRUE,
#   reduce = FALSE # This argument is to remove the branch of unknown taxonomy.
# ) +
#   scale_fill_manual(
#     values = cols
#   ) +
#   guides(color = guide_legend(
#     keywidth = 0.1,
#     keyheight = 0.2,
#     order = 3,
#     ncol=1)
#   ) +
#   theme(
#     panel.background = element_rect(fill=NA),
#     legend.position = "right",
#     plot.margin = margin(0,0,0,0),
#     legend.key.width = unit(0.2, "cm"),
#     legend.key.height = unit(0.2, "cm"),
#     legend.spacing.y = unit(0.02, "cm"),
#     legend.title = element_text(size=7),
#     legend.text = element_text(size=6),
#     legend.box.spacing = unit(0.02,"cm")
#   )
# diffclade_p
# ggsave("diffclade.png", width = 12, height = 12, dpi = 600)

# #convert the phyloseq object to MPSE object
# MPSE <- as.MPSE(physeq_decont)
# new_lefse <- mp_diff_analysis(
#   .data = MPSE,
#   .abundance=normalize,
#   .group= Colony,
#   #.sec.group = Condition,
#   ml.method = "lda",
#   ldascore = 2,
#   relative = TRUE)
# 
# p <- MPSE %>% new_lefse
# p <- p + 
#   scale_fill_manual(
#     aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
#     values = c("Healthy" = "green", "Recovered" = "blue", "Diseased" = "pink")
#   )  + 
#   scale_fill_manual(
#     values = c("Healthy" = "green", "Recovered" = "blue", "Diseased" = "pink") # The LDA barplot layer
#   )
# ### and the fill aes for hight light layer of tree was renamed to 'fill_new_new'
# p <- p + 
#   scale_fill_manual(
#     aesthetics = "fill_new_new",
#     values = c("Healthy" = "green", "Recovered" = "blue", "Diseased" = "pink")
#   )
# p
# 
# ### visualizing the differential taxa with cladogram
# f <- MPSE %>% 
#   mp_plot_diff_cladogram(
#     label.size = 2.5, 
#     hilight.alpha = .3, 
#     bg.tree.size = .5, 
#     bg.point.size = 2, 
#     bg.point.stroke = .25
#   ) + 
#   scale_fill_diff_cladogram(values = c("Healthy" = "green", 
#                                        "Recovered" = "blue", "Diseased" = "pink"
#   )) +
#     scale_size_continuous(range = c(1, 4))
# f
# 
# 

####################### UPSET PLOT #############################################


### Create upset diagram
upsetda <- get_upset(obj=physeq_carp, factorNames="Treatment")
#upset(upsetda, sets=unique(as.vector(sample_data(physeq)$Status)), 
#      sets.bar.color = "#56B4E9",
#      order.by = "freq", 
#      empty.intersections = "on",mainbar.y.label	= "Number of ASVs", 
#      point.size = 8, line.size	= 2)
#ggsave("upset_status.png", width = 12, height = 8, dpi = 600)
write.csv(upsetda, "upset_carp_Treatment.csv")


Treatment = c("Control", "Infection")

upset(
  upsetda,
  Treatment,
  name='Treatment',
  queries=list(
    upset_query(
      intersect=c('Control', 'Infection'),
      color='red',
      fill='red',
      only_components=c('intersections_matrix', 'Intersection size')
    )),
  base_annotations=list(
    'Number of Species'=intersection_size(
      text=list(
        vjust=-0.1,
        hjust=-0.1,
        angle=45, size=5
      )
    )
    + annotate(
      geom='text', size=5, x=Inf, y=Inf,
      label=paste('Total Species:', nrow(upsetda)),
      vjust=1, hjust=1
    )
  ),
  min_size=5,
  width_ratio=0.1,  stripes='white',
  set_sizes=FALSE,
  themes=upset_default_themes(text=element_text(size=15))
)
ggsave("upset_carp_Treatment.png", width = 6, height = 6, dpi = 600)


### Create upset diagram
upsetda <- get_upset(obj=physeq_herm, factorNames="Treatment")
#upset(upsetda, sets=unique(as.vector(sample_data(physeq)$Status)), 
#      sets.bar.color = "#56B4E9",
#      order.by = "freq", 
#      empty.intersections = "on",mainbar.y.label	= "Number of ASVs", 
#      point.size = 8, line.size	= 2)
#ggsave("upset_status.png", width = 12, height = 8, dpi = 600)
write.csv(upsetda, "upset_herm_Treatment.csv")


Treatment = c("Control", "Infection")

upset(
  upsetda,
  Treatment,
  name='Treatment',
  queries=list(
    upset_query(
      intersect=c('Control', 'Infection'),
      color='red',
      fill='red',
      only_components=c('intersections_matrix', 'Intersection size')
    )),
  base_annotations=list(
    'Number of Species'=intersection_size(
      text=list(
        vjust=-0.1,
        hjust=-0.1,
        angle=45, size=5
      )
    )
    + annotate(
      geom='text', size=5, x=Inf, y=Inf,
      label=paste('Total Species:', nrow(upsetda)),
      vjust=1, hjust=1
    )
  ),
  min_size=5,
  width_ratio=0.1,  stripes='white',
  set_sizes=FALSE,
  themes=upset_default_themes(text=element_text(size=15))
)
ggsave("upset_herm_Treatment.png", width = 6, height = 6, dpi = 600)

##################### Rarefaction #############################################

# for reproducibly random number
set.seed(1024)
#rareres <- get_rarecurve(obj=physeq, chunks=400)

#p_rare <- ggrarecurve(obj=rareres,
#                      indexNames=c("Observe","Chao1","ACE"),
#) +
#  theme(legend.spacing.y=unit(0.01,"cm"),
#        legend.text=element_text(size=4))

#prare1 <- ggrarecurve(obj=rareres, factorNames="Status",
#                      indexNames=c("Observe", "Chao1", "ACE")
#) +
#  scale_fill_manual(values = c("Healthy" = "green", "Recovered" = "blue",
#                               "Diseased" = "pink"))+
#  scale_color_manual(values = c("Healthy" = "green", "Recovered" = "blue",
#                                "Diseased" = "pink"))+
#  theme_bw()+
#  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
#        strip.background = element_rect(colour=NA,fill="grey"),
#        strip.text.x = element_text(face="bold"))          
#
#prare2 <- ggrarecurve(obj=rareres,
#                      factorNames="Status",
#                      shadow=FALSE,
#                     indexNames=c("Observe", "Chao1", "ACE")
#) +
#  scale_color_manual(values = c("Healthy" = "green", "Recovered" = "blue",
#                                "Diseased" = "pink"))+
#  theme_bw()+
#  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
#        strip.background = element_rect(colour=NA,fill="grey"),
#        strip.text.x = element_text(face="bold"))

#p_rare / prare1 / prare2

MPSE <- as.MPSE(physeq_merged_Species)

cols = c("Control" = "green",  "Infection" = "pink")


MPSE %<>%
  mp_rrarefy(.abundance=Abundance) %>%
  mp_cal_rarecurve(.abundance=RareAbundance, chunks=500)
p_rare <- MPSE %>%
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve,
    .alpha = c(Observe, Chao1, ACE),
  ) +
  theme(
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.spacing.y = unit(0.01,"cm"),
    legend.text = element_text(size=4)
  )
prare1 <- MPSE %>%
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve,
    .alpha = c(Observe, Chao1, ACE),
    .group = Treatment
  ) +
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  theme_bw()+
  theme(
    axis.text=element_text(size=8), panel.grid=element_blank(),
    strip.background = element_rect(colour=NA,fill="grey"),
    strip.text.x = element_text(face="bold")
  )
prare2 <- MPSE %>%
  mp_plot_rarecurve(
    .rare = RareAbundanceRarecurve,
    .alpha = c(Observe, Chao1, ACE),
    .group = Treatment,
    plot.group = TRUE
  ) +
  scale_color_manual(values = cols)+
  scale_fill_manual(values = cols) +
  theme_bw()+
  theme(
    axis.text=element_text(size=8), panel.grid=element_blank(),
    strip.background = element_rect(colour=NA,fill="grey"),
    strip.text.x = element_text(face="bold")
  )
(p_rare / prare1 / prare2) + patchwork::plot_annotation(tag_levels="A")


ggsave("rarefaction.png", width = 12, height = 14, dpi = 600)



save.image("drosophila.Rdata")

load("drosophila.Rdata")
