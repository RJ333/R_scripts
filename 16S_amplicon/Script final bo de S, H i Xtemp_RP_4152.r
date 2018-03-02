Y<-read.csv(file.choose(),sep=',')
summary(Y)

#remove lables
dataY<-Y[c(2:55)]     # here I have OTUs in rows and samples in colums
summary(dataY)

#transpose
tY<-t(dataY)    # I will have here OTUs in colums and samples in row. This is how you should have your file to start with, but check how it goes with the labels (make sure they are not consider a #variable)

P<-min(tY)
P<-6500 # a


T<-100 #number of iterations
Sp<-length(tY[1,])
Site<-length(tY[,1])
Tot<-rowSums(tY)
DataR<-matrix(0, nrow=Site, ncol=Sp)
XtempTot<-matrix(0, nrow=Site, ncol=Sp)  # this one is the one that you shuould keep if you only want a subsampled OTU table
H_sim<-matrix(0, nrow=T, ncol=Site)      # the next 6 lines are here to calculate the means for richness, chao, diversity (shannon) out of the 100 permutations if you don't need this values you can skip #these lines and continue at SR_sim...... The problem is that the subsampled table has decimal numbers and then you will have to round it to calculate richness. I personally think that if you are #going to check on Richness and diversity you should keep them. This steps, of course make the whole thing much longer.
S.obs_sim<-matrix(0, nrow=T, ncol=Site)
S.chao1_sim<-matrix(0, nrow=T, ncol=Site)
se.chao1_sim<-matrix(0, nrow=T, ncol=Site)
S.ACE_sim<-matrix(0, nrow=T, ncol=Site)
se.ACE_sim<-matrix(0, nrow=T, ncol=Site)
SR_sim<-matrix(0, nrow=Site, ncol=T)
for (s in 1:Site) {                    # here is where the actual randomitzation starts, it actually says that for every site it should calculate everything that comes afterwards and do that 100 times
    DataR[s,]<-tY[s,]/Tot[s]
	}
DataN<-tY
DataCum<-matrix(0, nrow=Site, ncol=Sp)
DataCum[,1]<-DataR[,1]
for (s in 1:Site) {
    for (s1 in 2:Sp) {
        DataCum[s,s1]<-DataR[s,s1]+DataCum[s,s1-1]
	}} 
tt<-0
for (st in 1:T){
        st
        tt<-tt+1
        Xtemp<-matrix(0, nrow=Site, ncol=Sp)
        for (s in 1:Site) {
        Random<-runif(P)
	#Plocka ut P individer slumpmässig art relativt förekomst
		for (n in 1:P) { 
		    In<-findInterval(Random[n],DataCum[s,])
		    In<-In+1
		    Xtemp[s,In]<-Xtemp[s,In]+1
                rm(In)
		} }
	estR<-estimateR(Xtemp)   # once again, if you are not interested in the richness or diversity you can skip the next 8 lines. Keep the XtempTot (this is related to the subsampled OTU table)
	S.obs_sim[tt,]<-estR[1,]
	S.chao1_sim[tt,]<-estR[2,]
	se.chao1_sim[tt,]<-estR[3,]
	S.ACE_sim[tt,]<-estR[4,]
	se.ACE_sim[tt,]<-estR[5,]
	Div<-diversity(Xtemp)
	H_sim[tt,]<-Div
	XtempTot<-XtempTot+Xtemp
	rm(Xtemp)
	}
XtempAve<-XtempTot/100
Smean<-colMeans(S.obs_sim)  # you can skip all the rest if you are not interested in Richness and diversity
Chao1mean<-colMeans(S.chao1_sim)
seChao1mean<-colMeans(se.chao1_sim)
ACEmean<-colMeans(S.ACE_sim)
seAcemean<-colMeans(se.ACE_sim)
Hmean<-colMeans(H_sim)

outData2<-rbind(Hmean)

outData<-rbind(Smean, Chao1mean, seChao1mean, ACEmean, seAcemean)