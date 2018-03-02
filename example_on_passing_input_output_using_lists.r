#see also http://faculty.nps.edu/sebuttre/home/R/lists.html for background on lists
#see also http://faculty.nps.edu/sebuttre/home/R/apply.html for functions which can be applied to lists

############load example data
# example data is a silvaNGS count table for OTUs (operational taxonomic unit) with separated taxonomic levels
otu_test<-read.csv(file.choose(),sep=";") #if from excel
head(otu_test)
#create a subset of interest e.g. genera that start with "P"; we only get the names, not the data because "$genus"
selection_of_interest<-subset(otu_test, grepl("P", genus))$genus
selection_of_interest
class(selection_of_interest)
#some names are duplicates, also the level information from the original data is still included
selection_of_interest<-sort(selection_of_interest)
selection_of_interest<-sub(' .*$','', selection_of_interest)
#we use "sub" to adjust the names , but it also changes the class from factor to character (what in this case is good)!
class(selection_of_interest)
#remove duplicate entries
selection_of_interest<-unique(selection_of_interest)
selection_of_interest
#store this information in a list ("container") and giving it a name, so it is accessible via reference
container<-list("OTUs_starting_with_P"=selection_of_interest)
container
str(container)
names(container)

############now use this data as input for the actual subsetting process
#define a function (not mandatory, but fancy)
OTU.subset<-function(dataset,column,otu){
	output=subset(dataset,grepl(otu,column))
	return(output)
	}
	
#test the function
OTU.subset(otu_test,otu_test$genus,"Paenisporosarcina")

#use a for loop (or something like *apply) to process all items in the "selection of interest"
#check loop first
for (OTU in selection_of_interest){
	print(OTU)
}

#now store the selected OTU subsets in the list called "container" by creating a new "place" in there, called "P_subsets"
for (OTU in selection_of_interest){
container$P_subsets[[OTU]] = OTU.subset(otu_test,otu_test$genus,OTU)
}
#how to access data in a list
names(container)
container
container$P_subsets
container$P_subsets$Paludibacter
container$P_subsets[["Paludibacter"]]

#now we use the output from the subsetting as input for the next function, 
#in this case we just check how many rows each subset has
for (OTU in selection_of_interest){
print(nrow(container$P_subsets[[OTU]]))
}
