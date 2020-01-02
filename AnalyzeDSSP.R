#Analyze DSSP
library(data.table)
library(ggplot2)
library(reshape2)
load("~/GitHub/DENV/DENV/fitnesstable_All.Rdata")
Wtable[M>=8 | L>=7,HC:="HC"]
table(Wtable$HC)
Wtable[!(M>=8 | L>=7),HC:="LC"]
Wt<-Wtable[HC=="HC"&muttype=="NonSyn"&passage=="MF.1",]
Wt[,aggBin:=cut(scale(Aggregation),breaks = c(-100,-.5,.5,100))]

#DEFINE Max accessibility per residue
maxAcc<-read.delim("~/Google_Drive/DengueDataAndAnalysis/MaxAccessibility.txt",header = T)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

INDIR<-"~/Research/Structures/"
files<-list.files(INDIR,pattern = "dssp.csv",full.names = F)
files
dsspList<-lapply(files,FUN = function(f){
  input<-data.table::fread(file = paste(INDIR,f,sep = ""))
  print(input)
  input[,label:=f]
  return(input)
})
allDSSP<-rbindlist(dsspList)

sd<-dcast.data.table(allDSSP,formula = pos+aa~.,value.var="acc",fun.aggregate = sd)
acc<-dcast.data.table(allDSSP,formula = pos+aa~.,value.var="acc",fun.aggregate = min)
phi<-dcast.data.table(allDSSP,formula = pos+aa~.,value.var="phi",fun.aggregate = median)
psi<-dcast.data.table(allDSSP,formula = pos+aa~.,value.var="psi",fun.aggregate = median)
NH<-dcast.data.table(allDSSP,formula = pos+aa~.,value.var="acc",fun.aggregate = length)

allDSSP[,h_ohn2_E:=as.numeric(limma::strsplit2(h_ohn2,split = ",")[,2])]
allDSSP[,h_nho1_E:=as.numeric(limma::strsplit2(h_nho1,split = ",")[,2])]
allDSSP[,h_ohn1_E:=as.numeric(limma::strsplit2(h_ohn1,split = ",")[,2])]
allDSSP[,h_nho2_E:=as.numeric(limma::strsplit2(h_nho2,split = ",")[,2])]
allDSSP$index=1:nrow(allDSSP)
allDSSP[,h_E:=sum(h_ohn1_E,h_ohn2_E,h_nho1_E,h_nho2_E),by='index']

E<-dcast.data.table(allDSSP,formula = pos+aa~.,value.var="h_E",fun.aggregate = mean)

bound<-cbind(acc,N=NH$.,SD=sd$.,phi=phi$.,psi=psi$.,E=E$.)

merged<-merge(Wt,bound,by = 'pos',allow.cartesian = T)
filteredM<-merged[aa==wtRes,]

mergedM<-merge(filteredM,maxAcc,by.x = "wtRes",by.y="AA")

mergedM[,rSASA:=./Emp]
mergedM[,zSASA:=scale(rSASA)]

ggplot(mergedM)+geom_boxplot(aes(reorder(wtRes,rSASA),rSASA))

mergedM[,phiBin:=cut(phi,breaks = seq(-180,180,10))]
mergedM[,psiBin:=cut(psi,breaks = seq(-180,180,10))]
mergedM[,accBin:=cut(scale(rSASA),breaks = c(seq(-3.5,0.5,1),10))]
mergedM[,accBin:=cut(scale(SASA),breaks = c(seq(-3.5,0.5,1),10))]

ggplot(mergedM[reg=="NS3"])+
  geom_boxplot(aes(accBin,wrel),varwidth = T)+ylim(0,3.5)

ggplot(mergedM[reg=="NS5"])+
  stat_summary(aes(status,rSASA),fun.data = function(X){
    y=median(X)
    ymin=y-(mad(X))
    ymax=y+(mad(X))
    return(data.frame(y,ymin,ymax))
  })

ggplot(mergedM)+facet_grid(~reg,margins = T)+
  stat_summary(aes(reorder(status,wrel,median),rSASA),fun.data = mean_cl_boot
  )

ggplot(mergedM)+facet_grid(~reg,margins = T)+
  geom_boxplot(varwidth = T,aes(reorder(status,wrel,median),zSASA))

ggplot(mergedM)+facet_grid(~reg,margins = T)+
  geom_boxplot(aes(rSASA>0.5,wrel))+coord_cartesian(ylim = c(0,4))

#geom_smooth(aes(E,wrel),method = "gam")



(mergedM$wrel~mergedM$zSASA>0)

ggplot(mergedM)+geom_point(aes())


binnedPhiPsi<-dcast(mergedM[],phiBin+psiBin+accBin~.,value.var = "wrel",fun.aggregate = median)
countedPhiPsi<-dcast(mergedM[],phiBin+psiBin+accBin~.,value.var = "wrel",fun.aggregate = length)
binnedPhiPsi$N<-countedPhiPsi$.
ggplot(binnedPhiPsi[binnedPhiPsi$N>10,])+geom_point(aes(phiBin,psiBin,size=N,col=.))+scale_color_viridis_c()+facet_wrap(~accBin)#Dedicated to Mauricio.
ggplot(mergedM)+geom_boxplot(varwidth=T,aes(as.factor(accBin),wrel))+ylim(0,3)+facet_wrap(~hydrophobic)
ggplot(mergedM)+geom_boxplot(varwidth=T,aes(as.factor(accBin),wrel))+ylim(0,3)+facet_wrap(~charged)
ggplot(mergedM)+geom_boxplot(varwidth=T,aes(as.factor(accBin),wrel))+ylim(0,3)+facet_wrap(~acidic)
ggplot(mergedM)+geom_boxplot(varwidth=T,aes(as.factor(accBin),wrel))+ylim(0,3)+facet_wrap(~basic)
ggplot(mergedM[reg%in%c("E","NS3","NS5","NS1")])+geom_boxplot(varwidth = T,aes(as.factor(accBin),wrel))+ylim(0,3)+facet_grid(~reg)

ggplot(mergedM)+geom_boxplot(varwidth = T,aes(as.factor(TMstat),wrel))+ylim(0,3)+facet_wrap(~accBin)

ggplot(Wt)+geom_boxplot(varwidth=T,aes(as.factor(TMstat),wrel))+ylim(0,3)
ggplot(Wt)+geom_boxplot(varwidth=T,aes(as.factor(feature),wrel))+ylim(0,3)
ggplot(Wt)+geom_boxplot(varwidth=T,aes(as.factor(feature),wrel))+ylim(0,3)
ggplot(Wt)+geom_bin2d(varwidth=T,aes(as.factor(Aggregation>0),status))

binnedAgg<-dcast(mergedM[reg=="NS5"],Aggregation~.,value.var = "wrel",fun.aggregate = median)

hist(scale(Wt$Aggregation))
table(dcast(mergedM[reg=="NS5"],pos+Class~.)$Class)

acc<-by(data = Wt$acc,INDICES = Wt$wtRes,FUN = function(X){median(X)})

plot(acc)

Wt[,alladjAcc:=acc/max(acc)]

ggplot(Wt)+geom_histogram(binwidth = 0.1,aes(wrel,fill=(alladjAcc>1)))+coord_cartesian(xlim = c(0,5))
ggplot(Wt)+geom_violin(draw_quantiles = c(0.5),aes(wtRes,wrel,fill=(acc==0)))+coord_cartesian(ylim = c(0,5))

N<-acast(Wt[wtRes!="X"&mutRes!="X",],wtRes~mutRes,value.var = "wrel",fun.aggregate = median)

N=-(N-1)
N>1=1
N<0=0
for(i in 1:20){
  N[i,i]<-0
}
image(N)
imputeN<-ape::ultrametric(N)
heatmap(imputeN,labCol= rownames(N),labRow = rownames(N),scale="none",col=viridis::viridis(100))
scaleAA<-cmdscale(as.dist(imputeN,upper = T))
pcAA<-princomp(imputeN)
plot(pcAA$loadings)
text(pcAA$loadings,labels=rownames(N))
merges<-data.frame(scaleAA,aaMap,N)
ggplot(merges)+
  geom_point(aes(x = X1,X2),size=5)+
  geom_text(aes(x = X1,X2,label=let.1),color="white")
text(scaleAA,labels=rownames(N))

ggplot(mergedM)+ ggpubr::theme_pubr()+
  geom_linerange(aes(index,ymin=0.5,ymax=rSASA,col=SS),stat="identity")+
  geom_text(show.legend = F,aes(index,-0.1,label=AA),cex=1.9)+
  xlab("Position")+ylab("Accessibility")+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
  scale_color_brewer(palette = "Paired")

ggplot(mergedM)+ ggpubr::theme_pubr()+
  rollup(mergedM,by="pos",j = zSASA)
  #geom_text(show.legend = F,aes(index,-0.1,label=AA),cex=1.9)+
  xlab("Position")+ylab("Accessibility")+scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1.0))+
  scale_color_brewer(palette = "Paired")

ggplot(data.table(zoo::rollmean(fill=NA,mergedM[order(num),mean(zSASA),by=num],k=21)))+geom_linerange(aes(num,ymax=V1,ymin=0))

ggplot(m[,])+geom_point(aes(pos,y=rSASA,col=binnedW.x))+scale_color_viridis_c()
