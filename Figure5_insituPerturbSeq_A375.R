# GitHub Code, Akana, Yoe et al., 2025
# Figure 5: In situ Perturb-seq application to A375 xenografts in NSG mice.

#************** Figure 5: In situ Perturb-seq applied to A375 xenografts **************#

Fig5_insituPerturbSeqA375<-function(){
  Rcv<-readRDS(get.file("InsituPerturbSeq_A375_LOOCV.rds"))
  rslts<-readRDS(get.file("InsituPerturbSeq_A375_DEGs_HLM.rds"))
  rF<-readRDS(get.file("InsituPerturbSeq_A375_Slide2.rds"))
  r<-readRDS(get.file("InsituPerturbSeq_A375_Slide2_tumor3_mal.rds"))
  r<-prep4OE(r,n.cat = 20)
  r$scores<-OE.fix.scores(get.OE(r,rslts$sig,min.n.genes = 5))
  Fig5bc_insituPerturbSeqA375_spatialMaps(r,rF)
  Fig5d_insituPerturbSeqA375_LOOCV(Rcv)
  Fig5e_insituPerturbSeqA375_GAsig(rslts$sig)
}

Fig5bc_insituPerturbSeqA375_spatialMaps<-function(r,rF){
  rF<-set.list(rF,rF$samples=="tumor3",name = "tumor3")
  b1<-rF$coor$x_global_px<=max(r$coor$x_global_px)
  b2<-rF$coor$y_global_px<=max(r$coor$y_global_px)
  rF<-set.list(rF,b1&b2)
  
  pdf(get.file(paste0("Fig5b_insituPerturbSeq_A375.pdf")),width = 15,height = 10)
  par(mfrow=c(2,3),oma = c(0, 0, 0, 0),xpd = T)
  insituPlot(rF,r,target = "MYC",plotIdx = c(T,F,T,F))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  
  pdf(get.file(paste0("Fig5c_insituPerturbSeq_A375.pdf")),width = 15,height = 10)
  par(mfrow=c(2,3),oma = c(0, 0, 0, 0),xpd = T)
  insituPlot(rF,r,target = "CD44",plotIdx = c(T,F,F,T))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

Fig5d_insituPerturbSeqA375_LOOCV<-function(R1){
  idx<-c("CD44","VAV1","MYC")
  pdf(get.file("Fig5d_insitu_A375_GA_CV.pdf"),width = 3.5,height = 10)
  par(mfrow=c(3,1),oma = c(0, 0, 0, 0),xpd = T)
  lapply(idx, function(x) insituPerturbSeq_LOOCV_ROCs(R1,R1$S1,R1$S2,x,names12 = c("Spatial GA (test)","Non-spatial GA (test)")))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
}

Fig5e_insituPerturbSeqA375_GAsig<-function(sig){
  pdf(get.file("Fig5e_insituPerturbSeqA375_GAsig.pdf"),width = 6,height = 6)
  summarize.program(sig,plot.flag = T,horiz = F)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
}

