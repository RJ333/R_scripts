# set working directory
setwd("D:/1_Fachliches/gene_richness_selected_genes")

# you can make a function out of this
import.multiple.table.files <- function(mypath, mypattern, sep, header_read, ...)
{
tmp.list.1 <- list.files(mypath, pattern = mypattern, recursive = TRUE)
tmp.list.2 <- list(length = length(tmp.list.1))
  for (i in 1:length(tmp.list.1)){tmp.list.2[[i]]<-read.table(tmp.list.1[i], sep = sep, header = header_read, ...)}
names(tmp.list.2)<-basename(list.files(mypath, pattern = mypattern, recursive = TRUE))
tmp.list.2
}

# use it like this
list.gene_richness <- import.multiple.table.files("D:/1_Fachliches/gene_richness_selected_genes", "counts.cov$", sep = " ", header_read = FALSE )
list.gene_length <- import.multiple.table.files("D:/1_Fachliches/gene_richness_selected_genes", "trimmed.cov$", sep = " ", header_read = FALSE )


str(list.gene_richness)
str(list.gene_length)

head(list.gene_length)

# awesome head function for e.g. data frames inside a list
# thanks to https://gist.github.com/pimentel/256fc8c9b5191da63819

head.list <- function(obj, n = 6L, ...)
{
    stopifnot(length(n) == 1L)
    origN <- n
    n <- if (n < 0L)
        max(length(obj) + n, 0L)
    else min(n, length(obj))
    lapply(obj[seq_len(n)], function(x)
           {
               tryCatch({
                   head(x, origN, ...)
               }, error = function(e) {
                   x
               })
           })
}
environment(head.list) <- asNamespace('utils')

head.list(list.gene_richness)
head.list(list.gene_length)

# rename column headers
colnames_richness <- c("contig_id", "base_pos", "coverage", "sample", "gene")
list.gene_richness <- lapply(list.gene_richness, setNames, colnames_richness)

colnames_length <- c("contig_id", "gene_length", "gene", "contig_length")
list.gene_length <- lapply(list.gene_length, setNames, colnames_length)

head.list(list.gene_richness)
head.list(list.gene_length)

# remove everything before "="
# this approach only works if you modify existing columns
head.list(list.gene_length)

list.gene_length <- lapply(rapply(list.gene_length, function(x) 
  sub('.*\\=', '', x), how = "list"), 
  as.data.frame)
  
head.list(list.gene_length)
# add new column with contig_ids without appended fraction number
# this approach is suited to ADD a new column
head.list(list.gene_richness)

list.gene_richness <- lapply(list.gene_richness, function(x) 
{x$full_contig_id <- sub('\\..*', '', x$contig_id);
x})

head.list(list.gene_richness)

# merge coverage and length information
# without SIMPLIFY it creates a mess of data inside the new object
list.gene_merged <- list()
list.gene_merged <- mapply(function(x,y) {
  as.data.frame(merge(x, y, by.x = "full_contig_id", by.y = "contig_id"))
  }, x = list.gene_richness, y = list.gene_length, SIMPLIFY = FALSE)

# check the outcome, [10] is phnM here, don't forget the head.list function  
  
head.list(list.gene_merged)
lapply(list.gene_richness,nrow) 
lapply(list.gene_richness[10],nrow)  # 218253
lapply(list.gene_length[10],nrow)  # 43
sapply(list.gene_merged[10],nrow)  # 251580
sapply(list.gene_merged,str)

# apply (type coercion?) 
# lapply = list apply 
# sapply = simple apply returns vector
# rapply = recursive apply 
# mapply = multiple apply
# tapply = apply for ragged arrays  (different length for vectors in array)
