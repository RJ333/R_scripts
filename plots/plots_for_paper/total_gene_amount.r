gene_amount<-read.csv(file.choose(),sep=";") #contains multiple entries for same name because of nmt_1, nmt_2 and so on
gene_amount2<-subset(gene_amount,!duplicated(contig_ec_gene)) #or unique(gene_amount)

nrow(gene_amount)
#[1] 247797
nrow(gene_amount2)
#[1] 247684

gene_amount_molten<-melt(gene_amount2)
head(gene_amount_molten)

##getting rid of zero values
gene_amount_molten_nonzero<-subset(gene_amount_molten,value!=0.0000000)
nrow(gene_amount_molten)
#2476840
nrow(gene_amount_molten_nonzero)
#1691258

write.csv(gene_amount_molten_nonzero,file="gene_amount_molten_nonzero.csv")

#count rows per variable
nrow(subset(gene_amount_molten_nonzero,variable=="A1_tpm"))  #slow
table(gene_amount_molten_nonzero$variable) #quick :-)

# A1_tpm  A2_tpm  A3_tpm  A4_tpm  A5_tpm  A6_tpm  A7_tpm  B8_tpm  B9_tpm B10_tpm 
# 153005  144253  161262  194177  170604  174176  182224  174750  158323  178484

str(all_unique2)
all_unique3<-all_unique2
all_unique3$total_gene_richness<-as.factor(all_unique3$total_gene_richness)
all_unique3$total_contig_richness<-as.factor(all_unique3$total_contig_richness)
all_unique3$gene_amount<-as.factor(all_unique3$gene_amount)
all_unique3$unique_contigs<-as.factor(all_unique3$unique_contigs)
all_unique3$unique_gene<-as.factor(all_unique3$unique_gene)
str(all_unique3)