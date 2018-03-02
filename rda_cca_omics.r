####cca + rda omics

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


#cca
cca_omics_nonzero_tpm<-cca(tomics_for_vegan2_nonzero_tpm ~ duration,data=meta_omics_tpm)
plot(cca_omics_nonzero_tpm,type="t",display="sites",main="cca_omics_nonzero_tpm",xlim=c(-4,4))
plot(cca_omics_nonzero_tpm,type="p",display="species",main="cca_omics_nonzero_tpm",xlim=c(-4,4))
text(cca_omics_nonzero_tpm,labels=row.names(tomics_for_vegan2_nonzero_tpm),cex=0.8)

summary(cca_omics_nonzero_tpm)

#rda
rda_omics_nonzero_tpm<-rda(tomics_for_vegan2_nonzero_tpm ~ duration,data=meta_omics_tpm)
plot(rda_omics_nonzero_tpm,type="t",display="sites",main="rda_omics_nonzero_tpm",xlim=c(-4,4))
plot(rda_omics_nonzero_tpm,type="p",display="species",main="rda_omics_nonzero_tpm",xlim=c(-4,4))
text(rda_omics_nonzero_tpm,labels=row.names(tomics_for_vegan2_nonzero_tpm),cex=0.8)

summary(rda_omics_nonzero_tpm)