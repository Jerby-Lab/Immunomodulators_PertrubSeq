# GitHub Code, Akana, Yoe et al., 2025
# Figure 6 related analyses: In situ Perturb-seq application to B16 tumors in C57BL/6 mice.

#************** B16 in situ Perturb-Seq data analyses **************#

insituPSeq_B16_main<-function(){
  overwrite<<-T
  rF<-readRDS(get.file("InsituPerturbSeq_B16.rds"))
  r<-set.list(rF,grepl("pool",rF$samples),name = "InsituPSeq_B16_pools_FINAL")
  
  # Map perturbation impact on the perturbed cancer cells
  r_mal<-readRDS(get.file("InsituPerturbSeq_B16_mal.rds"))
  rslts<-insituPSeq_B16_malignant_DEG(r_mal)
  
  # Map perturbation impact on nearby immune and stromal cells
  rslts<-insituPSeq_B16_pMCPs(r)
}

insituPSeq_B16_malignant_DEG<-function(r){
  file1<-"InsituPerturbSeq_B16_DEGs.rds"
  if(!overwrite&file.exists((file1))){
    print(paste("Loading",file1))
    return(readRDS(file1))
  }
  rslts<-insituPerturb_DEG_HLM(r,outputfile = file1)
  rslts$colocalize<-insituPS_colocalizeTest(r)
  saveRDS(rslts,file = file1)
  return(rslts)
}

insituPSeq_B16_pMCPs<-function(r){
  file1<-get.file("InsituPerturbSeq_B16_pMCPs")
  if(!overwrite&file.exists((file1))){
    print(paste("Loading",file1))
    return(readRDS(file1))
  }
  rslts<-list()
  rslts$Macrophages<-insituPS_prtOtherCells(r1 = r,cell.type = "Macrophage",full = T,partial.cor.flag = partial.cor.flag)
  rslts$CAFs<-insituPS_prtOtherCells(r1 = r,cell.type = "CAF",full = T,partial.cor.flag = partial.cor.flag)
  rslts$Endothelial<-insituPS_prtOtherCells(r1 = r,cell.type = "Endothelial",full = T,partial.cor.flag = partial.cor.flag)
  rslts$Tcell<-insituPS_prtOtherCells(r1 = r,cell.type = "T.cell",full = T,partial.cor.flag = partial.cor.flag)
  sig<-lapply(rslts, function(x) return(x$sig))
  rslts$sig<-unlist(sig,recursive = F)
  rslts$cell2env_Lifr<-insituPSeq_B16_Lifr(r)
  saveRDS(rslts,file = file1)
  return(rslts)
}

insituPSeq_B16_Lifr<-function(r){
  r$fovs<-r$frames
  r1<-get.frames(r,1000)
  X1<-average.mat.rows(t(r1$tpm),r1$frames)
  r1<-set.list(r1,!is.element(r1$targets,c("UD","None"))&is.element(r1$cell.types,"Malignant"))
  table(r1$cell.types,r1$targets)
  r1$Clcf1<-X1[r1$frames,"Clcf1"]>median(X1[r1$frames,"Clcf1"])
  r1$Osm<-X1[r1$frames,"Osm"]>median(X1[r1$frames,"Osm"])
  r1$Lif<-X1[r1$frames,"Lif"]>median(X1[r1$frames,"Lif"])
  r1$Clcf1_lifr<-r1$Clcf1&(r1$targets=="Lifr")
  r1$Osm_lifr<-r1$Osm&(r1$targets=="Lifr")
  r1$Lif_lifr<-r1$Lif&(r1$targets=="Lifr")
  
  r2<-set.list(r1,is.element(r1$targets,c("Lifr","NTC","Empty")));print(table(r2$targets))
  r2$comp<-r2$comp/max(r2$comp)
  rslts<-list(name = "Lifr_GAsigs")
  rslts$HLM<-apply.formula.all.HLM(r = r2,X = r2$targets=="Lifr",Y = r2$tpm,
                                   MARGIN = 1,formula = "y ~ (1 | fovs) + comp + x + Clcf1 + Osm + Lif + Clcf1_lifr + Osm_lifr + Lif_lifr")$HLM
  colnames(rslts$HLM)<-gsub("xTRUE","Z",colnames(rslts$HLM))
  rslts$sig<-get.top.cor(rslts$HLM1[,c("Z","Clcf1_lifrTRUE","Osm_lifrTRUE","Lif_lifrTRUE")],q = 200,min.ci = 2)
  return(rslts)
}
