# GitHub Code, Akana, Yoe et al., 2025
# Figure 6: In situ Perturb-seq application to B16 tumors in C57BL/6 mice.

#************** Figure 6: In situ Perturb-Seq applied to B16 tumors **************#

Fig6_insituPerturbSeqB16<-function(){
  r<-readRDS(get.file("InsituPerturbSeq_B16_referenceMap.rds"))
  r1<-readRDS(get.file("InsituPerturbSeq_B16_slide2_FOVs2plot.rds"))
  rslts1<-readRDS(get.file("InsituPerturbSeq_B16_DEGs.rds"))
  rslts2<-readRDS(get.file("InsituPerturbSeq_B16_pMCPs"))
  
  Fig6a_insituPerturbSeqB16_spatialMaps(r1)
  Fig6b_insituPerturbSeqB16_UMAPs(r)
  Fig6e_insituPerturbSeqB16_GA.pMCPs_barplot(rslts1,rslts2)
  Fig6f_insituPertrubSeqB16_pMCPs_Tcells(rslts2)
  return()
}

Fig6a_insituPerturbSeqB16_spatialMaps<-function(r1){
  cell2rgb <- list(CAF = col2rgb('#DE8C00'),Endothelial = col2rgb('#B79F00'),
                   Malignant = col2rgb("#00B8E7"),Macrophage = col2rgb('#F8766D'),
                   T.cell = col2rgb('#F564E3'),UD = col2rgb("grey90"))
  idx<-c("365","363","371","358")
  seg_path<-get.file("InsituPerturbSeq_B16_slide2_FOV")
  lapply(idx, function(x) Fig6a_insituPerturbSeqB16_spatialMaps1(seg_path,cell2rgb,r1 = r1,x = x))
  return()
}

Fig6a_insituPerturbSeqB16_spatialMaps1<-function(seg_path,cell2rgb,r1,x){
  segfile<-paste0(seg_path,x,".csv")
  b<-r1$frames==paste0("F00",x)
  celltypes<-r1$cell.types[b];
  names(celltypes)<-get.strsplit(r1$cells[b],"_",2)
  outfile<-paste0("Fig6a_",gsub(".csv",".jpg",basename(segfile)))
  samplename<-gsub(".csv","",basename(segfile))
  v<-inSituPlotSegmentation1(seg_path = segfile,celltypes = celltypes,
                             trans.flag = T,cell2rgb = cell2rgb,samplename = samplename,
                             outpath = dirname(get.file(outfile)),outfile = outfile)
  return()
}

Fig6b_insituPerturbSeqB16_UMAPs<-function(r){
  D1<-r$seurat
  p1 <- DimPlot(D1, reduction = "umap", group.by = "celltype", label = TRUE,label.size = 6,
                repel = TRUE) + xlab("UMAP1") + ylab("UMAP2") +
    NoLegend() + ggtitle("Reference Map: B16 in situ Perturb-Seq")
  p2 <- DimPlot(D1, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
                label.size = 6, repel = TRUE) + xlab("UMAP1") + ylab("UMAP2") +
    NoLegend() + ggtitle("Reference Map: B16 in situ Perturb-Seq")
  colnames(r$scores)<-paste(colnames(r$scores),"signature expression")
  l<-c(list(p1,p2),umap.ggplot(r$umap,labels = r$scores,size = 0.001))
  
  pdf(get.file("Fig6b_insituPerturbSeqB16_UMAPs.pdf"),width = 22,height = 16)
  v<-call_multiplot(l,nplots = 6,cols = 3)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

Fig6e_insituPerturbSeqB16_GA.pMCPs_barplot<-function(rslts1,rslts2){
  
  rslts2<-rslts2[c("Macrophages","CAFs","Endothelial","Tcell")]
  cell.types<-c("Macrophages","Fibroblasts","Endothelial cells","T cells")
  names(rslts2)<-cell.types
  sigA<-lapply(rslts2,function(x) return(x$sig))
  
  file1<-get.file("Fig6e_insituPertrubSeqB16_GA_pMCP_barplots.pdf")
  
  pdf(file1,width = 8,height = 8)
  par(mfrow=c(2,2),oma = c(0, 0, 0, 0),xpd = T)
  X1<-summarize.program(rslts1$sig,plot.flag = T,main = "GA signatures, B16 in situ Pertrub-Seq")
  X<-GA.MCP_barplot(sigA,main = "GA MCPs, B16 in situ Pertrub-Seq")
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

Fig6f_insituPertrubSeqB16_pMCPs_Tcells<-function(rslts){
  file1<-get.file("Fig6f_insituPertrubSeqB16_pMCPs_Tcells.pdf")
  rslts<-rslts[c("Macrophages","CAFs","Endothelial","Tcell")]
  cell.types<-c("Macrophages","Fibroblasts","Endothelial cells","T cells")
  names(rslts)<-cell.types
  pdf(file1,width = 8,height = 4)
  v<-call_multiplot(GA.MCP_boxplot(rslts$`T cells`,cell.type = "T cells",x = c("Cd274","Ccl22","Wnt3"),q1 = 0.1),nplots = 3,cols = 3)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
}


