# install R to get R to newest version
install.packages("installr")
library(installr)
# new menu item appears, where you can update R to the newest version



# in new R version
# install phyloseq package
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

# Check version
packageVersion('phyloseq')

install.packages("vegan")
install.packages("tidyverse")

# set working directory
getwd()
setwd("D:/")
getwd()

# read otu and sample data into R
otus <- read.csv("janine.csv", sep = ";", row.names = 1)
metafile <- read.csv("janine_meta.csv", sep = ";", row.names = 1)

# remove column named "Summe"
otus <- subset(otus, select = -Summe)

# remove sample Temora 4
names <- "Tem4"
# names <- c("Tem3", "Tem5")
otus <- subset(otus, !(rownames(otus) %in% names))
# get data into phyloseq
require(phyloseq)
metafile2 <- sample_data(metafile)
otu2 <- otu_table(otus, taxa_are_rows = FALSE)

								
# merge meta data and OTU table as phyloseq object
janine_ps <- merge_phyloseq(otu2, metafile2)


# this is the function we call to split our data into different subsets
get_sample_subsets <- function(ps, Habitat, threshold){
	sample_subset <- sample_data(ps)[ which(sample_data(ps)$Habitat == Habitat),]
	phy_subset <- merge_phyloseq(otu_table(ps),
								 sample_subset)
	phy_subset2 <- filter_taxa(phy_subset, function (x) {sum(x > 0) > threshold}, prune = TRUE)
	return(phy_subset2)
}

# these are the parameters passed to function
ps <- janine_ps
habitats <- c("Acartia", "Temora", "Water", "Particle")
threshold <- 1

# this is the for loop which calls the subsetting function 
# for each Habitat and OTU threshold
div_subset_list <- list() 

if(length(div_subset_list) == 0) {

			for (habitat in habitats) {
				print(paste0("habitat is ", habitat))
				tmp <-	get_sample_subsets(ps = ps, 
									   Habitat = habitat, 
			otus						   threshold = threshold)
				div_subset_list[[paste(habitat,  
										  "min_reads_per_OTU", 
										  threshold, 
										  sep = "_")]] <- tmp
			
		
	}
print(div_subset_list)
} else {
	print("list is not empty, abort to prevend appending...")
}

# how to get data from list?

# adress by index
div_subset_list[1]  # as list element
div_subset_list[[1]]  # as data.frame (or whatever the original class was)

# adress by name
names(div_subset_list)
div_subset_list[["Acartia_min_reads_per_OTU_1"]] 

# for example
otu_table(div_subset_list[["Acartia_min_reads_per_OTU_1"]])

# estimate diversity indeces for all subsamples
richness_list <-
lapply(div_subset_list, estimate_richness, measures = c("Observed", 
													 "Chao1", 
													 "ACE", 
												     "Shannon", 
												     "Simpson", 
												     "InvSimpson", 
												     "Fisher"))
													 
# append to one dataframe and save for Excel
richness_all_habitats <- do.call(rbind, richness_list)
write.csv(richness_all_habitats, file = "richness_janine.csv")

# plot richness
plot_richness(div_subset_list[[1]])
plot_richness(div_subset_list[[2]])
plot_richness(div_subset_list[[3]])
plot_richness(div_subset_list[[4]])

# save diversity data of a specific habitat
richness_acartia <- 
estimate_richness(div_subset_list[[1]], measures = c("Observed", 
													 "Chao1", 
													 "ACE", 
												     "Shannon", 
												     "Simpson", 
												     "InvSimpson", 
												     "Fisher"))
richness_acartia