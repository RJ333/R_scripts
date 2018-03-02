#####dca vegan mit omics

#counttable: otus in columns, samples in row.names -> e.g.
tomics_for_vegan2_nonzero_tpm
#metadata: metadata as dataframe in columns, samples in row.names, same order 
#angabe der zu berücksichtigenden paramater mit ~...+...+...
#informationen in summary: x und y werte für samples und species, welche achse erklärt welchen faktor, wo liegen die schwerpunkte

head(tomics_for_vegan2_nonzero_tpm)
names(tomics_for_vegan2_nonzero_tpm)

#metafile erstellen
meta_omics_tpm<-read.csv(file.choose(),sep=";",row.names=1)
head(meta_omics_tpm)
str(meta_omics_tpm)
#meta_omics_tpm$day<-as.factor(meta_omics_tpm$day)


#dca
dca_omics_nonzero_tpm<-decorana(tomics_for_vegan2_nonzero_tpm)
plot(dca_omics_nonzero_tpm,type="t",display="sites",main="dca_omics_nonzero_tpm",xlim=c(-4,4))
plot(dca_omics_nonzero_tpm,type="p",display="species",main="dca_omics_nonzero_tpm",xlim=c(-4,4))
text(dca_omics_nonzero_tpm,labels=row.names(tomics_for_vegan2_nonzero_tpm),cex=0.8)

summary(dca_omics_nonzero_tpm)