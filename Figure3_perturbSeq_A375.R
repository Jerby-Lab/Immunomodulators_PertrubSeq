# GitHub Code, Akana, Yoe et al., 2025; Module 1
# Figure 3: Perturb-seq application to A375 melanoma cells in monoculture and coculture with NY-ESO-1 TCR T cells.

#************** Figure 3: Perturb-seq applied to A375 cells in monoculture and coculture with T cells **************#

Fig3_PerturbSeq<-function(){
  r<-readRDS(get.file("PerturbA375_processed_wSeurat.rds"))
  rNTC<-readRDS(get.file("PerturbA375_NTCs.rds"))
  r1<-readRDS(get.file("PerturbA375_monoculture.rds"))
  X<-readRDS(get.file("PerturbA375_GAsig_vs_hits.rds"))
  sig<-readRDS(get.file("PerturbA375_GA_sig.rds"))
  sigCo<-readRDS(get.file("PerturbA375_covsmono_sig.rds"))
  sigHits<-readRDS(get.file("ImmuneHits.rds"))
  sgRNAQC<-readRDS(get.file("PerturbA375_sgRNA_zscoreBasedSimilarities.rds"))
  
  Fig3b_PerturbSeq_sgRNA.zscores(sgRNAQC)
  Fig3c_PerturbSeq_umaps(r)
  Fig3d_PerturbSeq_cocultureSigHeatmap(rNTC,sigCo = sigCo,sigHits = sigHits)
  Fig3e_PerturbSeq_GA.vs.hits(X)
  Fig3f_PerturbSeq_coRegHeatmap(sig,sigHits)
  return()
}

Fig3b_PerturbSeq_sgRNA.zscores<-function(sgRNAQC){
  pdf(get.file("Fig3b_sgRNA_zscorePearson.pdf"),width = 10,height = 4)
  call_multiplot(list(sgRNAQC$monoculture$plot,sgRNAQC$coculture$plot),nplots = 2,cols = 2);
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
}

Fig3c_PerturbSeq_umaps<-function(r,recompute.embedding = F){
  if(recompute.embedding){
    r<-seuratW_get.embedding(r,no.genes = 5000,n.pcs = 20,cd.flag = T,umap.flag = T,
                             norm.flag = T,full.flag = F,return.model = F)
  }
  X<-cbind.data.frame(Conditions = paste(r$labels,ifelse(r$ccB[,"Melanoma_cell_cycle_Tirosh"]=="cycling","cycling","non-cycling")),
                      MYC.CRISPRa = ifelse(r$targets=="MYC","MYC","Other"),
                      Cell.cycle = r$cc.scores[,1],
                      BIRC3 = r$tpm["BIRC3",])
  colnames(X)<-gsub("."," ",colnames(X),fixed = T)
  pdf(get.file("Fig3c_umaps.pdf"),width = 58,height = 14)
  call_multiplot(umap.ggplot(r$umap,labels = X,size = 1),nplots = 4,cols = 4)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

Fig3d_PerturbSeq_cocultureSigHeatmap<-function(r,sigCo,sigHits){
  hits<-c("BIRC3",unlist(sigHits[c("Res","Sen")]))
  g<-intersect(sigCo$edgeR500.up,hits)
  D1<-r$seurat
  r$cc<-ifelse(r$ccB[,"Melanoma_cell_cycle_Tirosh"]=="cycling","cycling","non-cycling")
  idx<-c('Monoculture, cycling','Monoculture, non-cycling','Coculture, non-cycling','Coculture, cycling')
  Idents(D1)<-factor(paste(r$labels,r$cc,sep = ", "),levels = idx)
  p<-DoHeatmap(D1,features = g,group.colors = c("#00BFC4","#C77CFF","#7CAE00","#F8766D"))
  pdf(get.file("Fig3d_heatmapCoVsMono.pdf"),width = 7,height = 10)
  print(p)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
}

Fig3e_PerturbSeq_GA.vs.hits<-function(Xf,main = "GA signatures vs. CRISPRa hits",MYC.trim = T,cex.names = 1){
  Xp<-Xf[order(Xf$N.hits,decreasing = F),-1]
  if(MYC.trim){
    Xp<-Xp[!is.element(rownames(Xp),c("genesA","MYC")),]
    rownames(Xp)<-gsub("MYCtop","MYC (top 500)",rownames(Xp))
  }else{
    Xp<-Xp[!is.element(rownames(Xp),c("genesA","MYCtop")),]
  }
  idx<-c("up.N.res","up.N.sen","down.N.res","down.N.sen","N","p.res","p.sen")
  idx<-intersect(idx,colnames(Xp))
  Xp<-Xp[Xp$N.hits>0,idx]
  v<-Xp[,"N"];names(v)<-rownames(Xp)
  
  pdf(get.file("Fig3e.pdf"),width = 14)
  par(mfrow=c(1,4),oma = c(0, 0, 0, 0),xpd = T)
  call_barplot(Xp[,1:4],beside = F,horiz = T,cex.names = cex.names,ylab = "GA signature",
               xlab = "No. of genes",legend.loc = "bottomright",las = 1,
               main = main,legend.names = c("Up, resistance","Up, sensitizing","Down, resistance","Down, sensitizing"))
  barplot(v,beside = F,horiz = T,cex.names = cex.names,ylab = "GA signature",
          las = 1,main = main,xlab = "No. of genes")
  if(!is.element("p.res",colnames(Xp))){return()}
  v1<-(-log10(Xp[,"p.res"]));names(v1)<-rownames(Xp)
  v2<-(-log10(Xp[,"p.sen"]));names(v2)<-rownames(Xp)
  barplot(v1,beside = F,horiz = T,cex.names = cex.names,ylab = "GA signature",xlab = "-log10(p-value)",
          las = 1,main = paste(main,"resistance hits enrichment",sep = "\n"),col = ifelse(v1>2,'#F8766D',"grey"))
  barplot(v2,beside = F,horiz = T,cex.names = cex.names,las = 1,ylab = "GA signature",xlab = "-log10(p-value)",
          main = paste(main,"sensitizing hits enrichment",sep = "\n"),col = ifelse(v2>2,'#F8766D',"grey"))
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

Fig3f_PerturbSeq_coRegHeatmap<-function(sig,sigHits,cytoscape.output = F){
  sig$MYC.up<-sig$MYC.top.up
  sig$MYC.down<-sig$MYC.top.down
  sig$MYC.top.up<-NULL;sig$MYC.top.down<-NULL
  sigHits<-setdiff.list1(sigHits,intersect(sigHits$Res,sigHits$Sen))
  X<-melt(sig)
  X<-cbind.data.frame(X,list.2.ids(X$value,sigHits[c("Res","Sen")]))
  X<-X[X$Anno!="",]
  sigCoReg<-split(X$L1,X$value)
  P<-GO.enrichment.lapply.v1(sigCoReg,unique(unlist(sigCoReg)),sigCoReg,valuType = "ob")
  X1<-list.2.ids(rownames(P),sigHits[c("Res","Sen")])
  pdf(get.file("Fig3f_coRegHeatmap.pdf"))
  call_heatmap(P,cluster.flag = "both",row.labels = X1,col.labels = X1,scale = "row",
               m.value = "Co-regulation score",k = 5,legend.flag = F)
  dev.off();par(font.axis = 2);par(font.lab = 2);par(font = 2)
  return()
}

ED_Figs_4_5_onTargetGeneActivaiton<-function(r1){
  D1<-r1$seurat
  targets<-sort(unique(r1$targets))
  prts<-sort(unique(r1$prt))
  g<-intersect(r1$genes,unique(r1$targets))
  Idents(D1)<-factor(r1$targets,levels = targets)
  p1<-call.dotPlotSeurat(D1,sort(g,decreasing = T),idents = targets,main = "On target activation per target",flip.flag = T)
  Idents(D1)<-factor(r1$prt,levels = prts)
  p2<-call.dotPlotSeurat(D1,sort(g,decreasing = T),idents = prts,main = "On target activation per sgRNA",flip.flag = F)
  pdf(get.file("ED_Fig4a.pdf"),width = 20,height = 20);print(p1);dev.off();par(font.axis = 2)
  pdf(get.file("ED_Fig4b_5.pdf"),width = 20,height = 60);print(p2);dev.off();par(font.axis = 2)
  return()
}

call.dotPlotSeurat<-function(D1,features,idents,main,flip.flag = T){
  p<-DotPlot(D1,features = features,idents = idents,cols="RdBu",col.max = 2)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  p<-p+labs(x = "Genes",y = "Perturbation")
  p<-p+theme(legend.text = element_text(size=10),
             legend.title = element_text(size=12),
             axis.text=element_text(size=10),
             axis.title=element_text(size=12,face="bold"))+ggtitle(main)
  if(flip.flag){p<-p+coord_flip()}
  return(p)
}


