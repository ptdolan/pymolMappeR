library(Rpdb)
library(bio3d)
library(seqinr)
library(data.table)

load("/Users/ptdolan/Google Drive/DenguePanelsAndOutput4-18/fitnesstable_All.Rdata")
correction=0

maxReturn<-function(X,Y){return(X[which(Y==max(Y))])}

WtableMap<-Wtable[muttype=="NonSyn",]
WtableMap[,Wrep:=maxReturn(wrel,wrel.ciLower),by=c("regPos","host","set")]
WtableMap[,aa3:=aa123(wtRes)]

PDBmap<-function(file,region,correction=0){  
  PDB=bio3d::read.pdb(file)
  cPDB<-clean.pdb(PDB)
  atomTable<-data.table(cPDB$atom)
  for(S in unique(WtableMap$set)){
    for(H in unique(WtableMap$host)){
      atomTable[,b:=1]
      
      atomTable[,b:=     WtableMap$Wrep [((WtableMap$regPos+correction)==resno)&aa321(resid)==(WtableMap$wtRes)&(WtableMap$set==S)&(WtableMap$host==H)&(WtableMap$reg==region)],by=resno]
      atomTable[,matchR:=WtableMap$wtRes[((WtableMap$regPos+correction)==resno)&resid==(WtableMap$aa3)&(WtableMap$set==S)&(WtableMap$host==H)&(WtableMap$reg==region)],by=resno]
      print(atomTable)
      cPDB$atom<-atomTable
      bio3d::write.pdb(pdb = cPDB,file = paste("~/Research/Structures/Fitness_",H,"_",S,"_",file,".pdb",sep=""))
    }}}


PDBmap("1ok8","E",0)
PDBmap("4O6B","NS1",0)
PDBmap("3EVG","NS5",0)



