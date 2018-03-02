######functions? first try:
##read required files
kallisto_prokka_concoct2_metaxa_checkm<-read.csv(file.choose(),row.names=1,sep=";")
meta_omics<-read.csv(file.choose(),row.names=1,sep=";")




###subset data
gene.subset<-function(dataset,column,gene){
	genestart<-paste0(gene,"_[0-999]")
	subset(dataset,grepl(genestart,column))
	}
	
dataset<-kallisto_prokka_concoct2_metaxa_checkm
column<-kallisto_prokka_concoct2_metaxa_checkm$genes
gene<-"gyrA"

gene_sub<-gene.subset(kallisto_prokka_concoct2_metaxa_checkm,kallisto_prokka_concoct2_metaxa_checkm$genes,gene)

###check subset
head(gene_sub)
names(gene_sub)
nrow(gene_sub)

###select data for further investigation
gene.select<-function(subset_table,column_list){
	subset_table[,columnlist]}

columnlist<-c(1,5,6,7,9,11,13,15,17,19,21,23,25,26,27,40,41)
#testing
gene_sel<-gene.select(gene_sub,columnlist)


###melt data of gene subset
melt.selected.gene<-function(gene_sel,columns_to_keep){
	require(reshape)
	melt(gene_sel,id=columns_to_keep)}

columns_to_keep<-c("contig_id","ec_number","genes","bin_gt1000","bin_nocutup","marker_lineage","completeness")
#testing
molten_gene<-melt.selected.gene(gene_sel,columns_to_keep)

###check molten subset
head(molten_gene)
nrow(molten_gene)

###merge with meta data
merge.molten.gene<-function(molten_gene,metatable,...){
	merge(molten_gene,metatable,...)}
#testing	
merged_gene<-merge.molten.gene(molten_gene,meta_omics,by.x="variable",by.y="row.names",all.x=T)

###use only non-zero entries
drop.value<-function(table,column_condition){
	droplevels(subset(table,column_condition))
}
#testing
dropped_gene<-drop.value(merged_gene,merged_gene$value!=0.0000000)

###check output
head(dropped_gene)
dropped_gene$value


###get genes and contig_ids from dropped genes per sample

###create object to collect data
gene_contigs<-vector() 
gene_contigs<-table((subset(dropped_gene,time=="6"&treatment=="glyph"))$contig_id)
gene_contigs<-rbind(gene_contigs,table((subset(dropped_gene,time=="7"&treatment=="glyph"))$contig_id))
gene_contigs<-rbind(gene_contigs,table((subset(dropped_gene,time=="8"&treatment=="glyph"))$contig_id))
gene_contigs<-rbind(gene_contigs,table((subset(dropped_gene,time=="10"&treatment=="glyph"))$contig_id))
gene_contigs<-rbind(gene_contigs,table((subset(dropped_gene,time=="12"&treatment=="glyph"))$contig_id))
gene_contigs<-rbind(gene_contigs,table((subset(dropped_gene,time=="15"&treatment=="glyph"))$contig_id))
gene_contigs<-rbind(gene_contigs,table((subset(dropped_gene,time=="19"&treatment=="glyph"))$contig_id))
gene_contigs<-rbind(gene_contigs,table((subset(dropped_gene,time=="6"&treatment=="control"))$contig_id))
gene_contigs<-rbind(gene_contigs,table((subset(dropped_gene,time=="12"&treatment=="control"))$contig_id))
gene_contigs<-rbind(gene_contigs,table((subset(dropped_gene,time=="19"&treatment=="control"))$contig_id))
###add sample names to the rows
row.names(gene_contigs)<-samples_omic
###run loop to count unique contigs per sample
unique_contigs<-vector()
for (i in 1:10){
unique_contigs<-rbind(unique_contigs,count.unique.contigs(i))
}
###bind the unique contigs to the subset table
shape.table<-function(table_gene_contigs,table_unique_contigs,...){
	table_gene_contigs<-cbind(table_gene_contigs,table_unique_contigs)
	table_gene_contigs<-as.data.frame(table_gene_contigs)
	names(table_gene_contigs)[ncol(table_gene_contigs)] <- "unique_contigs"
	table_gene_contigs<-as.data.frame(rbind(table_gene_contigs,colSums(table_gene_contigs)))
	row.names(table_gene_contigs)[nrow(table_gene_contigs)]<-"sum_per_contig"
	return(table_gene_contigs)
}
#testing
gene_contigs<-shape.table(gene_contigs,unique_contigs)

###prepare contig data for plotting
prepare.contig.table<-function(shaped_table,meta_table,...){
	shaped_table<-merge(shaped_table,meta_table,...)
	shaped_table_plot<-subset(shaped_table,treatment != "NA")
	row.names(shaped_table_plot)<-samples_omic
	shaped_table_plot<-shaped_table_plot[,-1]
	return(shaped_table_plot)
}

#testing
gene_contig_plot<-prepare.contig.table(gene_contigs,meta_omics_twodigit, by="row.names",all.x=T)


##############################################do the same stuff with genes instead of contigs

#testing
dropped_gene<-drop.value(merged_gene,merged_gene$value!=0.0000000)

###empty object before running
gene_genes<-vector() 
gene_genes<-table((subset(dropped_gene,time=="6"&treatment=="glyph"))$genes)
gene_genes<-rbind(gene_genes,table((subset(dropped_gene,time=="7"&treatment=="glyph"))$genes))
gene_genes<-rbind(gene_genes,table((subset(dropped_gene,time=="8"&treatment=="glyph"))$genes))
gene_genes<-rbind(gene_genes,table((subset(dropped_gene,time=="10"&treatment=="glyph"))$genes))
gene_genes<-rbind(gene_genes,table((subset(dropped_gene,time=="12"&treatment=="glyph"))$genes))
gene_genes<-rbind(gene_genes,table((subset(dropped_gene,time=="15"&treatment=="glyph"))$genes))
gene_genes<-rbind(gene_genes,table((subset(dropped_gene,time=="19"&treatment=="glyph"))$genes))
gene_genes<-rbind(gene_genes,table((subset(dropped_gene,time=="6"&treatment=="control"))$genes))
gene_genes<-rbind(gene_genes,table((subset(dropped_gene,time=="12"&treatment=="control"))$genes))
gene_genes<-rbind(gene_genes,table((subset(dropped_gene,time=="19"&treatment=="control"))$genes))
###add sample names to the rows
row.names(gene_genes)<-samples_omic
###run loop to count unique isoforms per sample
###set function
count.gene.isoforms <- function(row,gene){ 
	sum(gene_genes[row,] > 0)
}
###run function with loop
unique_gene_isoforms<-vector()
for (i in 1:10){
unique_gene_isoforms<-rbind(unique_gene_isoforms,count.gene.isoforms(i))
}
###bind the unique genes to the subset table
gene_genes<-shape.table(gene_genes,unique_gene_isoforms)

###prepare gene data for plotting
gene_genes_plot<-prepare.contig.table(gene_genes,meta_omics_twodigit, by="row.names",all.x=T)

#merge with meta data
columns_to_keep<-c("time","day","treatment","cells_ml","glyph_mg_L")
#testing
molten_gene_gene<-melt.selected.gene(gene_genes_plot,columns_to_keep)

###generate factor for cell number ratio to other values
find.cellfactor<-function(column){
	max(column)
}

cellfactor<-find.cellfactor(gene_contig_plot$unique_contigs)


###adjust functions
count.unique.contigs <- function(row,gene){ sum(gene_contigs[row,] > 0)}
