library(Rpdb)
library(bio3d)
library(seqinr)
library(data.table)

load("/Users/ptdolan/Google Drive/DenguePanelsAndOutput4-18/fitnesstable_All.Rdata")
correction=0

maxReturn<-function(X,Y){return(X[which(Y==max(Y))])}
meanReturn<-function(X,Y){return(mean(Y))}

WtableMap<-Wtable[muttype=="NonSyn",]
WtableMap[,Wrep:=meanReturn(wrel,wrel.ciLower),by=c("regPos","host","set")]
WtableMap[,Wrep:=maxReturn(wrel,wrel.ciLower),by=c("regPos","host","set")]

WtableMap<- dcast(Wtable,value.var = "wrel.ciLower",formula = regPos+host+set+reg+wtRes~.,fun.aggregate = max,name='Wrep')
WtableMap[,aa3:=aa123(wtRes)]

PDBmap<-function(file,region,correction=0){  
  PDB=bio3d::read.pdb(file)
  cPDB<-clean.pdb(PDB)
  atomTable<-data.table(cPDB$atom)
  for(S in unique(WtableMap$set)){
    for(H in unique(WtableMap$host)){
      Wt<-WtableMap[(WtableMap$set==S)&(WtableMap$host==H),]
      atomTable[,b:=1]
      #HERE
      atomTable[,b:=Wt$Wrep[((Wt$regPos+correction)==resno)&(resid==Wt$aa3)&(Wt$reg==region)], by=eleno]
      atomTable[,matchR:=Wt$wtRes[((Wt$regPos+correction)==resno)&resid==(Wt$aa3)&(Wt$reg==region)],by=eleno]
      print(atomTable)
      cPDB$atom<-atomTable
      cPDB$R<-matchR
      bio3d::write.pdb(pdb = cPDB,file = paste("~/Research/Structures/Fitness_",H,"_",S,"_",file,".pdb",sep=""))
    }}}

PDBmap("1ok8","E",0)
PDBmap("4O6B","NS1",0)
PDBmap("3EVG","NS5",0)



