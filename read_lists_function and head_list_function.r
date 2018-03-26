# a function to load multiple files of a type from a folder into a list
import.multiple.table.files <- function(mypath, mypattern, sep, header_read, ...)
{
tmp.list.1 <- list.files(mypath, pattern = mypattern, recursive = TRUE)
tmp.list.2 <- list(length = length(tmp.list.1))
  for (i in 1:length(tmp.list.1)){tmp.list.2[[i]]<-read.table(tmp.list.1[i], sep = sep, header = header_read, ...)}
names(tmp.list.2)<-basename(list.files(mypath, pattern = mypattern, recursive = TRUE))
tmp.list.2
}

# use it like this
my_list_of_count_data.frames <- import.multiple.table.files("D:/1_Fachliches/gene_richness_selected_genes", "counts.cov$", sep = " ", header_read = FALSE )


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

head.list(mylist)

# example of aggregate and apply
StandSpieltag <- data.frame(apply(punkte,2, FUN=function(x) {aggregate(x,by=list(punkte[,1]),sum)[,2]}))

# apply
# lapply = list apply 
# sapply = simple apply
# rapply = recursive apply 
# mapply = multiple apply
# tapply = apply for ragged arrays
