library(Rpdb)
library(data.table)
library(muscle)
library(Biostrings)
library(seqinr)
library(lettercase)


#Full Length sequence (experimental data sequence)
templateSeq<-AAString("MNNQRKKAKNTPFNMLKRERNRVSTVQQLTKRFSLGMLQGRGPLKLFMALVAFLRFLTIPPTAGILKRWGTIKKSKAINVLRGFRKEIGRMLNILNRRRRSAGMIIMLIPTVMAFHLTTRNGEPHMIVSRQEKGKSLLFKTEDGVNMCTLMAMDLGELCEDTITYKCPLLRQNEPEDIDCWCNSTSTWVTYGTCTTMGEHRREKRSVALVPHVGMGLETRTETWMSSEGAWKHVQRIETWILRHPGFTMMAAILAYTIGTTHFQRALIFILLTAVTPSMTMRCIGMSNRDFVEGVSGGSWVDIVLEHGSCVTTMAKNKPTLDFELIKTEAKQPATLRKYCIEAKLTNTTTESRCPTQGEPSLNEEQDKRFVCKHSMVDRGWGNGCGLFGKGGIVTCAMFRCKKNMEGKVVQPENLEYTIVITPHSGEEHAVGNDTGKHGKEIKITPQSSITEAELTGYGTVTMECSPRTGLDFNEMVLLQMENKAWLVHRQWFLDLPLPWLPGADTQGSNWIQKETLVTFKNPHAKKQDVVVLGSQEGAMHTALTGATEIQMSSGNLLFTGHLKCRLRMDKLQLKGMSYSMCTGKFKVVKEIAETQHGTIVIRVQYEGDGSPCKIPFEIMDLEKRHVLGRLITVNPIVTEKDSPVNIEAEPPFGDSYIIIGVEPGQLKLNWFKKGSSIGQMFETTMRGAKRMAILGDTAWDFGSLGGVFTSIGKALHQVFGAIYGAAFSGVSWTMKILIGVIITWIGMNSRSTSLSVTLVLVGIVTLYLGVMVQADSGCVVSWKNKELKCGSGIFITDNVHTWTEQYKFQPESPSKLASAIQKAHEEGICGIRSVTRLENLMWKQITPELNHILSENEVKLTIMTGDIKGIMQAGKRSLRPQPTELKYSWKTWGKAKMLSTESHNQTFLIDGPETAECPNTNRAWNSLEVEDYGFGVFTTNIWLKLKEKQDVFCDSKLMSAAIKDNRAVHADMGYWIESALNDTWKIEKASFIEVKNCHWPKSHTLWSNGVLESEMIIPKNLAGPVSQHNYRPGYHTQITGPWHLGKLEMDFDFCDGTTVVVTEDCGNRGPSLRTTTASGKLITEWCCRSCTLPPLRYRGEDGCWYGMEIRPLKEKEENLVNSLVTAGHGQVDNFSLGVLGMALFLEEMLRTRVGTKHAILLVAVSFVTLITGNMSFRDLGRVMVMVGATMTDDIGMGVTYLALLAAFKVRPTFAAGLLLRKLTSKELMMTTIGIVLLSQSTIPETILELTDALALGMMVLKMVRNMEKYQLAVTIMAILCVPNAVILQNAWKVSCTILAVVSVSPLLLTSSQQKTDWIPLALTIKGLNPTAIFLTTLSRTSKKRSWPLNEAIMAVGMVSILASSLLKNDIPMTGPLVAGGLLTVCYVLTGRSADLELERAADVKWEDQAEISGSSPILSITISEDGSMSIKNEEEEQTLTILIRTGLLVISGLFPVSIPITAAAWYLWEVKKQRAGVLWDVPSPPPMGKAELEDGAYRIKQKGILGYSQIGAGVYKEGTFHTMWHVTRGAVLMHKGKRIEPSWADVKKDLISYGGGWKLEGEWKEGEEVQVLALEPGKNPRAVQTKPGLFKTNAGTIGAVSLDFSPGTSGSPIIDKKGKVVGLYGNGVVTRSGAYVSAIAQTEKSIEDNPEIEDDIFRKRRLTIMDLHPGAGKTKRYLPAIVREAIKRGLRTLILAPTRVVAAEMEEALRGLPIRYQTPAIRAEHTGREIVDLMCHATFTMRLLSPVRVPNYNLIIMDEAHFTDPASIAARGYISTRVEMGEAAGIFMTATPPGSRDPFPQSNAPIIDEEREIPERSWNSGHEWVTDFKGKTVWFVPSIKAGNDIAACLRKNGKKVIQLSRKTFDSEYVKTRTNDWDFVVTTDISEMGANFKAERVIDPRRCMKPVILTDGEERVILAGPMPVTHSSAAQRRGRIGRNPKNENDQYIYMGEPLENDEDCAHWKEAKMLLDNINTPEGIIPSMFEPEREKVDAIDGEYRLRGEARKTFVDLMRRGDLPVWLAYRVAAEGINYADRRWCFDGVKNNQILEENVEVEIWTKEGERKKLKPRWLDARIYSDPLALKEFKEFAAGRKSLTLNLITEMGRLPTFMTQKARDALDNLAVLHTAEAGGRAYNHALSELPETLETLLLLTLLATVTGGIFLFLMSGRGIGKMTLGMCCIITASILLWYAQIQPHWIAASIILEFFLIVLLIPEPEKQRTPQDNQLTYVVIAILTVVAATMANEMGFLEKTKKDLGLGSIATQQPESNILDIDLRPASAWTLYAVATTFVTPMLRHSIENSSVNVSLTAIANQATVLMGLGKGWPLSKMDIGVPLLAIGCYSQVNPITLTAALFLLVAHYAIIGPGLQAKATREAQKRAAAGIMKNPTVDGITVIDLDPIPYDPKFEKQLGQVMLLVLCVTQVLMMRTTWALCEALTLATGPISTLWEGNPGRFWNTTIAVSMANIFRGSYLAGAGLLFSIMKNTTNTRRGTGNIGETLGEKWKSRLNALGKSEFQIYKKSGIQEVDRTLAKEGIKRGETDHHAVSRGSAKLRWFVERNMVTPEGKVVDLGCGRGGWSYYCGGLKNVREVKGLTKGGPGHEEPIPMSTYGWNLVRLQSGVDVFFIPPEKCDTLLCDIGESSPNPTVEAGRTLRVLNLVENWLNNNTQFCIKVLNPYMPSVIEKMEALQRKYGGALVRNPLSRNSTHEMYWVSNASGNIVSSVNMISRMLINRFTMRYKKATYEPDVDLGSGTRNIGIESEIPNLDIIGKRIEKIKQEHETSWHYDQDHPYKTWAYHGSYETKQTGSASSMVNGVVRLLTKPWDVVPMVTQMAMTDTTPFGQQRVFKEKVDTRTQEPKEGTKKLMKITAEWLWKELGKKKTPRMCTREEFTRKVRSNAALGAIFTDENKWKSAREAVEDSRFWELVDKERNLHLEGKCETCVYNMMGKREKKLGEFGKAKGSRAIWYMWLGARFLEFEALGFLNEDHWFSRENSLSGVEGEGLHKLGYILRDVSKKEGGAMYADDTAGWDTRITLEDLKNEEMVTNHMEGEHKKLAEAIFKLTYQNKVVRVQRPTPRGTVMDIISRRDQRGSGQVGTYGLNTFTNMEAQLIRQMEGEGVFKSIQHLTITEEIAVQNWLARVGRERLSRMAISGDDCVVKPLDDRFASALTALNDMGKIRKDIQQWEPSRGWNDWTQVPFCSHHFHELIMKDGRVLVVPCRNQDELIGRARISQGAGWSLRETACLGKSYAQMWSLMYFHRRDLRLAANAICSAVPSHWVPTSRTTWSIHAKHEWMTTEDMLTVWNRVWIQENPWMEDKTPVESWEEIPYLGKREDQWCGSLIGLTSRATWAKNIQAAINQVRSLIGNEEYTDYMPSMKRFRREEEEAGVLW")

#load("/Users/ptdolan/Google Drive/DenguePanelsAndOutput4-18/fitnesstable_All.Rdata")

WtableMap<-Wtable[muttype=="NonSyn",]

dispFuncFitness<-function(X){
  if(length(X)>0){
    return(max(X))
  }
  else{return(-1)}
}

#Control function to edit for specific analysis
pdbMapFitness<-function(pdbF,Wt=WtableMap,aligner="water"){
  for (S in c("A","B")){
    for (H in c("Mosquito","Human")){
      inputWs<-Wt[(Wt$set==S)&(Wt$host==H),]
      Aligned<-pdbAligner(pdbF,templateSeq,aligner = aligner)
      pdbO<-Aligned[[1]]
      mapped<-Aligned[[2]]
      atom<-data.table(pdbO$atoms)
      atom[,temp:=-1]
      atom[recname=="ATOM"&mapped=="mapped", temp:=dispFunc(inputWs[(inputWs$pos==resid),]), by=c("chainid","resid")]
      #pdbO$atoms$resid<-atom$resid
      pdbO$atoms$temp<- (atom$temp+1.05)
      Rpdb::write.pdb(x = pdbO, file = paste("~/Research/Structures/Fitness_",H,"_",S,"_",pdbF,".pdb",sep=""))
    }
  }
}


Wtable[,FE:=(1-wrel)*freq]

pdbMapFitness<-function(pdbF,Wt=WtableMap,aligner="water"){
  for (S in c("A","B")){
    for (H in c("Mosquito","Human")){
      inputWs<-Wt[(Wt$set==S)&(Wt$host==H),]
      Aligned<-pdbAligner(pdbF,templateSeq,aligner = aligner)
      pdbO<-Aligned[[1]]
      mapped<-Aligned[[2]]
      atom<-data.table(pdbO$atoms)
      atom[,temp:=-1]
      atom[recname=="ATOM"&mapped=="mapped",temp:=dispFunc(inputWs$wrel.ciLower[(inputWs$pos==(resid))]),by=c("chainid","resid")]
      #pdbO$atoms$resid<-atom$resid
      pdbO$atoms$temp<- (atom$temp+1.05)
      Rpdb::write.pdb(x = pdbO, file = paste("~/Research/Structures/Fitness_",H,"_",S,"_",pdbF,".pdb",sep=""))
    }
  }
}


pdbAligner<-function(pdbF, templateSeq, aligner="water"){
  pdbInput<-curl::curl_download(paste("https://files.rcsb.org/view/",pdbF,".pdb",sep=""),"curl.pdb")
  pdbO<- Rpdb::read.pdb("curl.pdb")
  atomTable<-data.table(pdbO$atoms)
  orig<-pdbO$atoms$resid
  print(pdbO$title)
  #generate muscle alignment of chain
  for(chain in unique(atomTable$chainid)){
    print(paste("Aligning chain ",chain,"...",sep = ""))
    seqtable<-dcast.data.table(atomTable[chainid==chain],chainid+resid~.,value.var = "resname",fun.aggregate = function(X){X[1]})
    idtable<-dcast.data.table(atomTable[chainid==chain],chainid+resid~.,value.var = "resid",fun.aggregate = function(X){X[1]})
    structureSeqString<-seqtable$. %>% str_lowercase() %>% str_title_case() %>% a() %>% na.omit() %>% paste(collapse = "") %>% Biostrings::AAString()
    #print(structureSeqString)

    if (length(structureSeqString)<3){
      print("...skipping short chain.")
    }
    else{
      alnDF<-data.frame()
      if (aligner=="water"){
          print("Aligning with Water...")
          water<-pairwiseAlignment(type = 'overlap',pattern = structureSeqString,subject=templateSeq) #End-free alignment "overlap" (no end gap penalties)
          
          boundaries <- water %>% Views() %>% ranges() %>% as.data.frame()
          indels <- water %>% subject() %>% indel() %>% as.data.frame()
          
          alignment <- AAStringSet(list(AAString(toString(pattern(water))),AAString(toString(subject(water)))))
          print(paste(pdbF,chain,".phy",sep=""))
          write.phylip(alignment,filepath = paste(pdbF,chain,".phy",sep="")) # Output alignment for QC
          alnMat <- alignment %>% Biostrings::as.matrix() %>% t()
          print(boundaries)
          rangeSub<-boundaries$start:boundaries$end
          if (boundaries$width>1){
            if(length(rangeSub)!=nrow(alnMat)){
              print("Removing gaps from reference...")
              for (row in 1:nrow(indels)){
                for(gap in 1:indels$width[row]){
                  rangeSub<-append(rangeSub,values = indels$start[row],after=indels$start[row])
                }
              }
            }
          alnDF<-data.frame(alnMat,pos=rangeSub)
          }
        }
      else{
        print("Aligning with Muscle...")
        alignment <- Biostrings::AAStringSet(list(structureSeqString,templateSeq)) %>% muscle(gapextend=-4.0,gapopen=-12)
        print(paste(pdbF,chain,".phy",sep=""))
        write.phylip(alignment,filepath = paste(pdbF,chain,".phy",sep="")) # Output alignment for QC
        alnMat <- alignment %>% Biostrings::as.matrix() %>% t()
        alnDF<-data.frame(alnMat,pos=1:nrow(alnMat))
        }

      #clean up gaps and assign residue positions based on full-length protein
      if (nrow(alnDF)>0){
        gapCleanedDF<-alnDF[alnDF[,1]!="-",] #remove gaps in query.
        filteredSeq<-seqtable[!is.na(seqtable$. %>% str_lowercase() %>% str_title_case() %>% a())]
        jointDF<-data.frame(filteredSeq, gapCleanedDF)
        atomTable[chainid==chain&recname=="ATOM",resid:=jointDF[(jointDF$chainid==chain)&(jointDF$resid==resid),"pos"],by=eleid]
        atomTable[chainid==chain&recname=="ATOM",mapped:="mapped",by=eleid]
        }
    }
  } 
  pdbO$atoms<-with(atomTable,atoms(recname, eleid, elename, alt, resname, chainid, resid,
                                   insert, x1, x2, x3, occ, temp, segid, basis='xyz'))
  return(list(pdbO,atomTable$mapped))
  
}

pdbMap("4M9F")
pdbMap("2JSF")
pdbMap("1OAM")
pdbMap("1OK8")
pdbMap("4UTC")
pdbMap("4O6B")
pdbMap("3EVG")
pdbMap("2R29")
pdbMap("2R69")
pdbMap("3C5X")
pdbMap("3C6E")

input<-read.csv("~/Research/ProteinBioinformatics/inputPDBs.txt", sep = "\t", header = F)

for (PDB in input$V2){ pdbMap(PDB) }
