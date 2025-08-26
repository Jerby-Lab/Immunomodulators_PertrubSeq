# GitHub Code, Akana, Yoe et al., 2025; Module 2
# Figure 5 related analyses: In situ Perturb-seq application to A375 xenografts in NSG mice.

insituPerturbSeqA375_main<-function(){
  files1<-list.files(get.file(""),pattern = "InsituPerturbSeq_A375_",full.names = T)
  files1<-files1[grepl("_mal",files1)]
  rA<-lapply(files1, readRDS)
  names(rA)<-paste0("tumor",1:4)
  overwrite<<-T
  R<-insituPerturbSeqA375_DEGs(rA)
  Rcv<-insituPerturbSeqA375_LOOCV(R$DEGs_perTumor,rA,min.n.genes = 5,q1 = 1000,sig ="BH.sig",pCut = )
  return()
}

insituPerturbSeqA375_DEGs<-function(rA){
  file0<-get.file("InsituPerturbSeq_A375_DEGs_perTumor.rds")
  file1<-get.file("InsituPerturbSeq_A375_DEGs_HLM.rds")
  file2<-get.file("InsituPerturbSeq_A375_DEGs_MAST.rds")
  
  if(!overwrite&file.exists(file1)){
    print(paste("Loading",file1))
    return(list(DEGs_perTumor = readRDS(file0),HLM = readRDS(file1),MAST = readRDS(file2)))
  }
 
  R<-lapply(rA, function(r) insituPerturbSeqA375_DEGs1(r))
  names(R)<-names(rA)
  Rh<-lapply(R, function(rslts) return(rslts$HLM))
  Rm<-lapply(R, function(rslts) return(rslts$MAST))
  rsltsH<-combine.GAs(Rh,q1 = 250,pCut = 1e-4)
  rsltsM<-combine.GAs(Rm,q1 = 250,pCut = 1e-4)
  saveRDS(R,file = file0)
  saveRDS(rsltsH,file = file1)
  saveRDS(rsltsM,file = file2)
  return(list(DEGs_perTumor = R,HLM = rsltsH,MAST = rsltsH))
}

insituPerturbSeqA375_DEGs1<-function(r){
  rslts2<-PerturbSeq_DEGs_MAST(r)
  rslts1<-insituPerturb_DEG_HLM(r = r)
  return(list(name = r$name,HLM = rslts1,MAST = rslts2))
}

insituPerturbSeqA375_LOOCV<-function(R,rA,min.n.genes = 5,q1 = 1000,sig ="BH.sig",pCut = 0.05){
  savefile<-get.file("InsituPerturbSeq_A375_LOOCV.rds")
  if(file.exists(savefile)){return(readRDS(savefile))}
  R2<-list()
  Rh<-lapply(R, function(rslts) return(rslts$HLM))
  Rm<-lapply(R, function(rslts) return(rslts$MAST))
  rsltsH<-combine.GAs(Rh,q1 = q1)
  rsltsM<-combine.GAs(Rm,q1 = q1)
  sigH<-rsltsH[[sig]]
  sigM<-rsltsM[[sig]]
  sigM<-setdiff.lists(sigM,rsltsH$sig)
  
  targets<-names(R[[1]]$HLM$full)
  for(x in names(R)){
    print(x)
    r<-rA[[x]]
    if(length(R)==2){
      rsltsH<-R[[setdiff(names(R),x)]]$HLM
      rsltsM<-R[[setdiff(names(R),x)]]$MAST
    }else{
      Rh<-lapply(R[setdiff(names(R),x)], function(rslts) return(rslts$HLM))
      Rm<-lapply(R[setdiff(names(R),x)], function(rslts) return(rslts$MAST))
      rsltsH<-combine.GAs(Rh,q1 = q1,pCut = pCut)
      rsltsM<-combine.GAs(Rm,q1 = q1,pCut = pCut)
    }
    
    S1<-cbind.data.frame(OE.fix.scores(get.OE(r,rsltsH[[sig]],min.n.genes = min.n.genes)),add = NA)
    S2<-cbind.data.frame(OE.fix.scores(get.OE(r,setdiff.lists(rsltsM[[sig]],rsltsH$sig),min.n.genes = min.n.genes)),add = NA)
    S3<-OE.fix.scores(get.OE(r,sigH,min.n.genes = min.n.genes))
    S4<-OE.fix.scores(get.OE(r,sigM,min.n.genes = min.n.genes))
    idx1<-match(targets,colnames(S1));idx1[is.na(idx1)]<-ncol(S1);S1<-S1[,idx1];colnames(S1)<-targets
    idx2<-match(targets,colnames(S2));idx2[is.na(idx2)]<-ncol(S2);S2<-S2[,idx2];colnames(S2)<-targets
    S3n<-center.matrix(S3,dim = 2,sd.flag = T)
    S4n<-center.matrix(S4,dim = 2,sd.flag = T)
    if(x==names(R)[1]){
      R2$S1<-S1;R2$S2<-S2
      R2$S3<-S3;R2$S4<-S4
      R2$S3n<-S3n;R2$S4n<-S4n
      R2$targets<-r$targets
    }else{
      R2$S1<-rbind(R2$S1,S1);R2$S2<-rbind(R2$S2,S2)
      R2$S3<-rbind(R2$S3,S3);R2$S4<-rbind(R2$S4,S4)
      R2$S3n<-rbind(R2$S3n,S3n);R2$S4n<-rbind(R2$S4n,S4n)
      R2$targets<-c(R2$targets,r$targets)
    }
    boxplot(OE.pred(r,CV = F,scores = S1)[,1],
            OE.pred(r,CV = F,scores = S2)[,1],
            OE.pred(r,CV = F,scores = S3)[,1],
            OE.pred(r,CV = F,scores = S4)[,1])
    print(OE.pred(r,CV = F,scores = S1))
  }
  
  R2$AUROCh<-OE.pred(R2,CV = F,scores = R2$S1)
  R2$AUROCm<-OE.pred(R2,CV = F,scores = R2$S2)
  R2$AUROCh_all<-OE.pred(R2,CV = F,scores = R2$S3)
  R2$AUROCm_all<-OE.pred(R2,CV = F,scores = R2$S4)
  R2$AUROCh_allN<-OE.pred(R2,CV = F,scores = R2$S3n)
  R2$AUROCm_allN<-OE.pred(R2,CV = F,scores = R2$S4n)
  X<-union.multiple.mats(R2[8:13])
  X<-X[,!grepl("ttest",colnames(X))]
  R2$sumAUC<-X
  saveRDS(R2,file = savefile)
  return(R2)
}

insituPerturbSeq_LOOCV_ROCs<-function(r1,S1,S2,x,names12 = c("Train","Test")){
  b1<-is.element(r1$targets,c(x,"NTC","Empty"))&!is.na(S1[,x])&!is.na(S2[,x])
  y<-is.element(r1$targets[b1],x)
  if(sum(y)<10){return()}
  P<-list(DEGs1 = S1[b1,x],DEGs2 = S2[b1,x])
  Y<-list(DEGs1 = y,DEGs2 = y)
  names(P)<-names12
  names(Y)<-names12
  plot.3aucs.mat(P = P,Y = Y,main = paste0(x," DEGs"))
  return()
}

combine.GAs<-function(R1,q1 = 250,pCut = 0.05){
  targets<-names(R1[[1]]$full)
  Z<-union.multiple.mats(lapply(R1, function(rslts) return(rslts$zscores)))
  
  P1<-get.pval.from.zscores(Z)
  idx<-get.strsplit(colnames(P1),".",2)
  P1g<-t(laply(targets,function(x) fisher.combine(P1[,idx==x])))
  colnames(P1g)<-paste0(targets,".up")
  
  P2<-get.pval.from.zscores(-Z)
  idx<-get.strsplit(colnames(P2),".",2)
  P2g<-t(laply(targets,function(x) fisher.combine(P2[,idx==x])))
  colnames(P2g)<-paste0(targets,".down")
  
  sig<-get.top.elements(cbind.data.frame(P1g,P2g),q = q1,min.ci = pCut)
  rslts<-list(P.up = P1g, P.down = P2g,sig = sig)
  rslts$BH<-p.adjust.mat(cbind.data.frame(P1g,P2g))
  rslts$BH.sig<-get.top.elements(rslts$BH,q = 5000,min.ci = 0.05)
  return(rslts)
}


