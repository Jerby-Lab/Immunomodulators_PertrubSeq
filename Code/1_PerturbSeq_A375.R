# GitHub Code, Akana, Yoe et al., 2025
# Figure 3 related analyses: Perturb-seq application to A375 melanoma cells in monoculture and coculture with NY-ESO-1 TCR T cells.

#************** Analysis of A375 Perturb-seq data **************#

PerturbA375_main<-function(){
  # Download the Perturb-seq data
  r<-readRDS(get.file("PerturbA375_processed.rds"))
  r1<-readRDS(get.file("PerturbA375_monoculture_wSeurat.rds"))
  r2<-readRDS(get.file("PerturbA375_coculture.rds"))
  rNTC<-readRDS(get.file("PerturbA375_NTCs.rds"))
  sig<-readRDS(get.file("PerturbA375_GA_sig.rds"))
  overwrite<<-F
  
  # Visualize on-target gene activation per gene and per sgRNA
  ED_Figs_4_5_onTargetGeneActivaiton(r1)
  
  # Test for consistency across different sgRNAs targeting the same gene
  rsltsZscores<-PerturbA375_sgRNA_zscoreBasedSimilarities_compute(r1,r2,sig)
  
  # Identify differential expressed genes (DEGs) in coculture compared to monoculture
  sigCo<-PerturbA375_DEGs_co.vs.mono(rNTC)
  
  # Identify GA signatures with and without filtering of lowly expressed genes
  PerturbA375_DEGs_edgeR(r = r,r1 = r1,r2 = r2,filter.flag = T)
  PerturbA375_DEGs_edgeR(r = r,r1 = r1,r2 = r2,filter.flag = F)
  R<-PerturbA375_DEGs_edgeR_load(filter.flag = T)
  R1<-PerturbA375_DEGs_edgeR_load(filter.flag = F)
  sig<-PerturbA375_saveFinalGAsignatures(R = R,R1 = R1)
  
  # Construct a Perturb-seq-based regulatory signaling network
  lig2rec<-readRDS(get.file("lig2rec.rds"))
  sigHits <- readRDS(get.file("CRISPRaA375_hits.rds"))
  PerturbSeqNet<-PerturbA375_convergingTargets(sig = sig,sigHits = sigHits,lig2rec = lig2rec,MYC.trim = T)
  
  return()
}

# **************** Differential Gene Expression Analyses, A375 **************** #

PerturbA375_DEGs_co.vs.mono<-function(r){
  file1<-get.file("PerturbA375_covsmono_sig.rds")
  file2<-get.file("PerturbA375_covsmono.rds")
  if(file.exists(file1)&&!overwrite){return(readRDS(file1))}
  r$samples<-gsub("_","",paste0(get.strsplit(r$cells,"_",1),r$prt))
  r$samples<-gsub("Monoculture","",r$samples)
  r$samples<-gsub("Coculture","",r$samples)
  rslts<-scRNA_DEGs_edgeR(r1 = r,batch = r$samples,labels1 = r$labels,ref1 = "Monoculture")
  saveRDS(rslts$sig,file = file1)
  saveRDS(rslts,file = file2)
  return(rslts$sig)
}

PerturbA375_DEGs_edgeR<-function(r,r1,r2,filter.flag = T){
  # All
  rsltsAt<-PerturbSeq_DEGs_edgeR(r,per.target = T,filter.flag = filter.flag)
  # Monoculture
  rslts1t<-PerturbSeq_DEGs_edgeR(r1,per.target = T,filter.flag = filter.flag)
  rslts1cv<-PerturbSeq_DEGs_edgeR_sgRNA_CV(r1,filter.flag = filter.flag)
  # Coculture
  rslts2t<-PerturbSeq_DEGs_edgeR(r2,per.target = T,filter.flag = filter.flag)
  rslts2cv<-PerturbSeq_DEGs_edgeR_sgRNA_CV(r2,filter.flag = filter.flag)
  return()
}

PerturbA375_DEGs_edgeR_load<-function(filter.flag){
  dir1<-get.file("PerturbA375_DEGs_")
  str1<-ifelse(filter.flag,"_filter.rds",".rds")
  R<-list(name = ifelse(filter.flag,"","_withFilter"))
  R$target1<-readRDS(paste0(dir1,"monoculture_Target",str1))
  R$target2<-readRDS(paste0(dir1,"coculture_Target",str1))
  R$targetA<-readRDS(paste0(dir1,"_Target",str1))
  R$cv1<-readRDS(paste0(dir1,"monoculture_sgRNA_CV",str1))
  R$cv2<-readRDS(paste0(dir1,"coculture_sgRNA_CV",str1))
  return(R)
}

PerturbA375_saveFinalGAsignatures<-function(R,R1){
  g<-names(R$targetA$full)
  sigOntarget<-union.lists(union.lists(R1$target1$sig,R1$target2$sig),R1$targetA$sig)
  sigOntarget<-intersect.list1(sigOntarget,g)
  X<-melt(sigOntarget);X$prt<-get.strsplit(X$L1,".",1);X<-X[grepl("up",X$L1),]
  X$same<-X$value==X$prt
  
  sig1<-R$target1$sig
  sig2<-R$target2$sig
  sigA<-R$targetA$sig
  sigA<-union.lists(sigA,sigOntarget)
  sigF<-union.lists(union.lists(sig1,sig2),sigA)
  sigF$genesA<-sort(unique(c(rownames(R$target1$zscores),
                             rownames(R$target2$zscores),
                             rownames(R$targetA$zscores))))
  sigA$genesA<-rownames(R$targetA$zscores)
  sigMYC<-get.top.cor(R$targetA$zscores,q = 250,idx = "MYC")
  sigA$MYC.top.up<-sigMYC$MYC.up
  sigA$MYC.top.down<-sigMYC$MYC.down
  summary(sigF);summary(sigA)
  barplot_sigs(sigF[1:122],sigA[1:122])
  saveRDS(sigF,file = get.file("PerturbA375_GA_sigFull.rds"))
  saveRDS(sigA,file = get.file("PerturbA375_GA_sig.rds"))
  return(sigA)
}

# **************** Differential Gene Expression Analyses - Generic Code **************** #

PerturbSeq_DEGs_edgeR<-function(r1,per.target = T,n1 = 10,filter.flag = F){
  if(per.target){
    targets<-setdiff(get.abundant(r1$targets,20,decreasing = T),"NTC")
    name<-paste0(r1$name,"_Target")
  }else{
    targets<-setdiff(get.abundant(r1$prt,20,decreasing = T),r1$prt[r1$targets=="NTC"])
    name<-paste0(r1$name,"_sgRNA")
  }
  if(filter.flag){name<-paste0(name,"_filter")}else{n1 = 5;print("No filter")}
  savefile<-get.file(paste0("PerturbA375_DEGs_",name,".rds"))
  if(!overwrite&file.exists(savefile)){return(readRDS(savefile))}
  
  merged <- SingleCellExperiment(assays = list(counts = r1$cd),mainExpName = "counts")
  r1$samples<-gsub("_",".",r1$samples,fixed = T)
  ids<-paste(r1$samples,r1$prt,sep = "_")
  summed <- aggregateAcrossCells(merged, id = ids)
  X<-lapply(targets, function(x){
    X1<-PerturbSeq_DEGs_edgeR1(summed, target = x,n1 = n1,filter.flag = filter.flag)
    return(X1)})
  names(X)<-targets
  X<-X[laply(X,is.list)]
  names(X)[names(X)=="KIAA0368"]<-"ECPAS"
  X0<-union.multiple.mats(lapply(X,function(X1){
    d1<-ifelse(X1$comparison=="factor(target)non",-1,1)
    X1$table$logFC<-d1*X1$table$logFC
    X1<-X1$table
    X1$zscoresPval<-get.cor.zscores(c = X1$logFC,p = X1$PValue)
    X1$zscoresFDR<-get.cor.zscores(c = X1$logFC,p = X1$FDR)
    return(X1)
  }))
  X1<-X0[,grepl("zscoresPval",colnames(X0))]
  X2<-X0[,grepl("zscoresFDR",colnames(X0))]
  colnames(X1)<-gsub(".zscoresPval","",colnames(X1))
  colnames(X2)<-gsub(".zscoresFDR","",colnames(X2))
  idx<-match(colnames(X1),rownames(X1))
  if(all(is.na(idx))){
    idx<-match(get.strsplit(colnames(X1),"_",1),rownames(X1))
  }
  X3<-cbind.data.frame(target = colnames(X1),diag(as.matrix(X1)[idx,]))
  rslts<-list(name = name,full = X,
              zscoresPval = X1,zscores = X2,onTarget = X3)
  
  rslts$sig<-get.top.cor(rslts$zscores,q = 10000,min.ci = 2)
  if(!per.target){
    rslts$sum<-PerturbSeq_sgRNA_DEG.consistency(sig = rslts$sig,jac.flag = 2,
                                                g = rownames(rslts$zscores),
                                                name = rslts$name)
  }
  rslts$cells<-r1$cells
  saveRDS(rslts,file = savefile)
  return(rslts)
}

PerturbSeq_DEGs_edgeR1<-function(summed,target = "MYC",n1 = 10,filter.flag = F){
  b2<-grepl(paste0("_",target,"_"),summed$ids)
  current <- summed[,grepl("non_targeting",summed$ids)|b2]
  y <- DGEList(counts(current), samples=colData(current))
  discarded <- current$ncells < n1
  y <- y[,!discarded]
  summary(discarded)
  Xdesign<-cbind.data.frame(batch = get.strsplit(y$samples$ids,"_",1),
                            target = get.strsplit(y$samples$ids,"_",2))
  if(length(unique(Xdesign$target))<2){print("Only 1 group.");return(NA)}
  if(filter.flag){
    keep <- filterByExpr(y, group=Xdesign$target)
    y <- y[keep,]
    summary(keep)
  }
  
  y <- calcNormFactors(y)
  design<-model.matrix(~factor(batch) + factor(target),Xdesign)
  
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust=TRUE)
  res <- glmQLFTest(fit, coef=ncol(design))
  X<-topTags(res,n = nrow(res$coefficients))
  return(X)
}

PerturbSeq_DEGs_edgeR_sgRNA_CV<-function(r,rslts,overwrite = T,filter.flag){
  # Identify GA signature for each target while holding out one of the 3 sgRNAs.
  # Test the predictive performances of the GA signature on the cells with the held out sgRNA.
  # Predictive performances are determined based on OE of the signatures and AUROC when considering all the cells
  # or only the control and target cells.
  library(scater)
  library(edgeR)
  prts<-sort(get.abundant(r$prt[r$targets!="NTC"],abn.c = 10))
  if(filter.flag){r$name<-paste0(r$name,"_filter")}
  savefile<-paste0("/Volumes/ljerby/PerturbSeq/Results/DEGs/PerturbA375/DEGs_edgeR_",r$name,"_sgRNA_CV.rds")
  if(!overwrite&file.exists(savefile)){return(readRDS(savefile))}
  
  if((!missing(rslts)|file.exists(savefile))&!overwrite){
    if(missing(rslts)){rslts<-readRDS(savefile)}
  }else{
    merged <- SingleCellExperiment(assays = list(counts = r$cd),mainExpName = "counts")
    summed <- aggregateAcrossCells(merged, id = paste(r$samples,r$prt,sep = "_"))
    X<-lapply(prts, function(x) PerturbSeq_DEGs_edgeR_sgRNA_CV1(summed,prt = x,filter.flag = filter.flag))
    names(X)<-prts
    X<-X[laply(X,is.list)]
    names(X)[names(X)=="KIAA0368"]<-"ECPAS"
    X0<-union.multiple.mats(lapply(X,function(X1){
      d1<-ifelse(X1$comparison=="factor(target)non",-1,1)
      X1$table$logFC<-d1*X1$table$logFC
      X1<-X1$table
      X1$zscoresPval<-get.cor.zscores(c = X1$logFC,p = X1$PValue)
      X1$zscoresFDR<-get.cor.zscores(c = X1$logFC,p = X1$FDR)
      return(X1)
    }))
    X1<-X0[,grepl("zscoresPval",colnames(X0))]
    X2<-X0[,grepl("zscoresFDR",colnames(X0))]
    colnames(X1)<-gsub(".zscoresPval","",colnames(X1))
    colnames(X2)<-gsub(".zscoresFDR","",colnames(X2))
    idx<-match(colnames(X1),rownames(X1))
    if(all(is.na(idx))){
      idx<-match(get.strsplit(colnames(X1),"_",1),rownames(X1))
    }
    X3<-cbind.data.frame(target = colnames(X1),diag(as.matrix(X1)[idx,]))
    rslts<-list(name = r$name,full = X,
                zscoresPval = X1,zscoresFDR = X2,onTarget = X3)
    rslts$sig<-get.top.cor(rslts$zscoresFDR,q = 1000,min.ci = 2)
    saveRDS(rslts,file = savefile)
  }
  
  sig<-rslts$sig[laply(rslts$sig,length)>1]
  r<-prep4OE(r[setdiff(names(r),"seurat")],n.cat = 100)
  r$scores<-OE.fix.scores(get.OE(r,sig))
  rslts$scores<-cbind.data.frame(cells = r$cells,targets = r$targets,
                                 prts = r$prt,scores = r$scores)
  prts<-intersect(prts,colnames(r$scores))
  rslts$sgRNA.pred<-get.mat(prts,c("auc.vs.NTC","auc.vs.all","ttest.vs.NTC","ttest.vs.all"))
  for(x1 in prts){
    x<-get.strsplit(x1,"_",1)
    b1<-is.element(r$prt,x1)|r$b.ctrl
    b2<-is.element(r$prt,x1)|r$targets!=x
    rslts$sgRNA.pred[x1,]<-c(get.auc(r$scores[b1,x1],r$targets[b1]==x),
                             get.auc(r$scores[b2,x1],r$targets[b2]==x),
                             t.test.labels(r$scores[b1,x1],r$targets[b1]==x,alternative = "greater"),
                             t.test.labels(r$scores[b2,x1],r$targets[b2]==x,alternative = "greater"))
  }
  rslts$sgRNA.pred<-cbind.data.frame(N = table(r$prt)[prts],rslts$sgRNA.pred)
  colnames(rslts$sgRNA.pred)[1:2]<-c("sgRNA","N.cells")
  rslts$sgRNA.pred$targets<-get.strsplit(rslts$sgRNA.pred$sgRNA,"_",1)
  print(paste0(round(100*mean(rslts$sgRNA.pred$auc.vs.NTC>0.5)),"% of sgRNAs show AUC > 0.5."))
  print(paste0(round(100*mean(rslts$sgRNA.pred$auc.vs.NTC>0.6)),"% of sgRNAs show AUC > 0.6."))
  print(paste0(round(100*mean(rslts$sgRNA.pred$auc.vs.NTC>0.7)),"% of sgRNAs show AUC > 0.7."))
  print(paste0(round(100*mean(rslts$sgRNA.pred$ttest.vs.NTC<0.05)),"% of sgRNAs show consistent effects."))
  print(paste0(round(100*mean(rslts$sgRNA.pred$ttest.vs.all<0.05)),"% of sgRNAs show consistent effects."))
  print("Saving results...")
  saveRDS(rslts,file = savefile)
  return(rslts)
}

PerturbSeq_DEGs_edgeR_sgRNA_CV1<-function(summed,prt,filter.flag){
  target<-get.strsplit(prt,"_",1)
  b1<-grepl(paste0("_",target),summed$ids)
  b2<-!grepl(paste0("_",prt),summed$ids)
  b2<-b1&b2
  
  current <- summed[,grepl("non_targeting",summed$ids)|b2]
  y <- DGEList(counts(current), samples=colData(current))
  discarded <- current$ncells < 10
  y <- y[,!discarded]
  
  Xdesign<-cbind.data.frame(batch = get.strsplit(y$samples$ids,"_",1),
                            target = get.strsplit(y$samples$ids,"_",2))
  if(length(unique(Xdesign$target))<2){print("Only 1 group.");return(NA)}
  design<-model.matrix(~factor(batch) + factor(target),Xdesign)
  if(filter.flag){
    keep <- filterByExpr(y, group=Xdesign$target)
    y <- y[keep,]
    summary(keep)
  }
  
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust=TRUE)
  res <- glmQLFTest(fit, coef=ncol(design))
  X<-topTags(res,n = nrow(res$coefficients))
  return(X)
}

PerturbSeq_DEGs_edgeR_per_sgRNA<-function(r,rslts,overwrite = T){
  # Identify GA signature for each target while holding out one of the 3 sgRNAs.
  # Test the predictive performances of the GA signature on the cells with the held out sgRNA.
  # Predictive performances are determined based on OE of the signatures and AUROC when considering all the cells
  # or only the control and target cells.
  library(scater)
  library(edgeR)
  prts<-sort(get.abundant(r$prt[r$targets!="NTC"],abn.c = 10))
  savefile<-paste0("/Volumes/ljerby/PerturbSeq/Results/DEGs/PerturbA375/DEGs_edgeR_",r$name,"_per_sgRNA.rds")
  if(!overwrite&file.exists(savefile)){return(readRDS(savefile))}
  
  
  merged <- SingleCellExperiment(assays = list(counts = r$cd),mainExpName = "counts")
  summed <- aggregateAcrossCells(merged, id = paste(r$samples,r$prt,sep = "_"))
  X<-lapply(prts, function(x) PerturbSeq_DEGs_edgeR_per_sgRNA1(summed,prt = x))
  names(X)<-prts
  X<-X[laply(X,is.list)]
  names(X)[names(X)=="KIAA0368"]<-"ECPAS"
  X0<-union.multiple.mats(lapply(X,function(X1){
    d1<-ifelse(X1$comparison=="factor(target)non",-1,1)
    X1$table$logFC<-d1*X1$table$logFC
    X1<-X1$table
    X1$zscoresPval<-get.cor.zscores(c = X1$logFC,p = X1$PValue)
    X1$zscoresFDR<-get.cor.zscores(c = X1$logFC,p = X1$FDR)
    return(X1)
  }))
  X1<-X0[,grepl("zscoresPval",colnames(X0))]
  X2<-X0[,grepl("zscoresFDR",colnames(X0))]
  colnames(X1)<-gsub(".zscoresPval","",colnames(X1))
  colnames(X2)<-gsub(".zscoresFDR","",colnames(X2))
  idx<-match(colnames(X1),rownames(X1))
  if(all(is.na(idx))){
    idx<-match(get.strsplit(colnames(X1),"_",1),rownames(X1))
  }
  X3<-cbind.data.frame(target = colnames(X1),diag(as.matrix(X1)[idx,]))
  rslts<-list(name = r$name,full = X,
              zscoresPval = X1,zscoresFDR = X2,onTarget = X3)
  rslts$sig<-get.top.cor(rslts$zscoresFDR,q = 1000,min.ci = 2)
  saveRDS(rslts,file = savefile)
  
  sig<-rslts$sig[laply(rslts$sig,length)>1]
  r<-prep4OE(r[setdiff(names(r),"seurat")],n.cat = 100)
  r$scores<-OE.fix.scores(get.OE(r,sig))
  rslts$scores<-cbind.data.frame(cells = r$cells,targets = r$targets,
                                 prts = r$prt,scores = r$scores)
  prts<-intersect(prts,colnames(r$scores))
  rslts$sgRNA.pred<-get.mat(prts,c("auc.vs.NTC","auc.vs.all","ttest.vs.NTC","ttest.vs.all"))
  for(x1 in prts){
    x<-get.strsplit(x1,"_",1)
    b1<-is.element(r$prt,x1)|r$b.ctrl
    b2<-is.element(r$prt,x1)|r$targets!=x
    rslts$sgRNA.pred[x1,]<-c(get.auc(r$scores[b1,x1],r$targets[b1]==x),
                             get.auc(r$scores[b2,x1],r$targets[b2]==x),
                             t.test.labels(r$scores[b1,x1],r$targets[b1]==x,alternative = "greater"),
                             t.test.labels(r$scores[b2,x1],r$targets[b2]==x,alternative = "greater"))
  }
  rslts$sgRNA.pred<-cbind.data.frame(N = table(r$prt)[prts],rslts$sgRNA.pred)
  rslts$sgRNA.pred$targets<-get.strsplit(rslts$sgRNA.pred$N.Var1,"_",1)
  print(paste0(round(100*mean(rslts$sgRNA.pred$auc.vs.NTC>0.5)),"% of sgRNAs show AUC > 0.5."))
  print(paste0(round(100*mean(rslts$sgRNA.pred$auc.vs.NTC>0.6)),"% of sgRNAs show AUC > 0.6."))
  print(paste0(round(100*mean(rslts$sgRNA.pred$auc.vs.NTC>0.7)),"% of sgRNAs show AUC > 0.7."))
  print(paste0(round(100*mean(rslts$sgRNA.pred$ttest.vs.NTC<0.05)),"% of sgRNAs show consistent effects."))
  print(paste0(round(100*mean(rslts$sgRNA.pred$ttest.vs.all<0.05)),"% of sgRNAs show consistent effects."))
  print("Saving results...")
  saveRDS(rslts,file = savefile)
  return(rslts)
}

scRNA_DEGs_edgeR<-function(r1,batch,labels1,ref1){
  library(scater);library(edgeR)
  ids<-paste(batch,labels1,sep = "_")
  merged <- SingleCellExperiment(assays = list(counts = r1$cd),mainExpName = "counts")
  current <- aggregateAcrossCells(merged, id = ids)
  y <- DGEList(counts(current), samples=colData(current))
  discarded <- current$ncells < 10
  y <- y[,!discarded]
  summary(discarded)
  Xdesign<-cbind.data.frame(batch = get.strsplit(y$samples$ids,"_",1),
                            labels = get.strsplit(y$samples$ids,"_",2))
  # keep <- filterByExpr(y, group=Xdesign$target)
  # y <- y[keep,]
  y <- calcNormFactors(y)
  # design<-model.matrix(~factor(batch) + factor(labels),Xdesign)
  design<-model.matrix(~factor(labels),Xdesign)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust=TRUE)
  res <- glmQLFTest(fit, coef=ncol(design))
  
  rslts<-topTags(res,n = nrow(res$coefficients))
  View(rslts$table[order(rslts$table$PValue,decreasing = F),])
  d1<-ifelse(rslts$comparison==paste0("factor(labels)",ref1),-1,1)
  
  X1<-rslts$table
  X1$logFC<-d1*X1$logFC
  X1$zscoresPval<-get.cor.zscores(c = X1$logFC,p = X1$PValue)
  X1$zscores<-get.cor.zscores(c = X1$logFC,p = X1$FDR)
  rslts$table<-X1
  
  sig1<-get.top.cor(rslts$table[abs(rslts$table$logFC)>0.25,],q = Inf,min.ci = 2,idx = "zscores")
  sig1t<-get.top.cor(rslts$table,q = 500,min.ci = 2,idx = "zscores")
  names(sig1)<-c("edgeR.up","edgeR.down")
  names(sig1t)<-c("edgeR500.up","edgeR500.down")
  rslts$sig<-c(sig1,sig1t)
  return(rslts)
}

# **************** QC **************** #

PerturbA375_sgRNA_zscoreBasedSimilarities_compute<-function(r1,r2,sig){
  savefile<-get.file("PerturbA375_sgRNA_zscoreBasedSimilarities.rds")
  if(!overwrite&file.exists(savefile)){return(readRDS(savefile))}
  r1<-PerturbSeq_getZscores(r1)
  r2<-PerturbSeq_getZscores(r2)
  rslts1<-PerturbSeq_sgRNA_zscoreBasedSimilarities(r1,sig,gene.selection = "de.level",deg.level = 20,main = "Monooculture",n1.genes = 1000, n2.genes = 5)
  rslts2<-PerturbSeq_sgRNA_zscoreBasedSimilarities(r2,sig,gene.selection = "de.level",deg.level = 20,main = "Coculture",n1.genes = 1000, n2.genes = 5)
  rslts<-list(monoculture = rslts1,coculture = rslts2,main = "sgRNA_zscoreBasedSimilarities")
  saveRDS(rslts,file = savefile)
  return(rslts)
}

PerturbSeq_getZscores<-function(r1){
  batches<-unique(r1$samples)
  r1$NTC_avg<-t(average.mat.rows(t(r1$tpm[,r1$b.ctrl]),r1$samples[r1$b.ctrl],f = colMeans))
  r1$NTC_sd<-t(average.mat.rows(t(r1$tpm[,r1$b.ctrl]),r1$samples[r1$b.ctrl],f = colSD))
  Z<-r1$tpm;Z[]<-NA
  r1$zscores<-Z
  r1$zscoresNorm<-Z
  for(x in batches){
    b<-is.element(r1$samples,x)
    r1$zscores[,b]<-sweep(r1$tpm[,b],MARGIN = 1,STATS = r1$NTC_avg[,x],FUN = "-")
  }
  r1$ids<-paste(r1$samples,r1$prt,sep = "_")
  r1$zscoresAv<-t(average.mat.rows(t(r1$zscores),r1$ids,f = colMeans))
  ids<-colnames(r1$zscoresAv)
  for(x in batches){
    ids<-gsub(paste0(x,"_"),"",ids)
  }
  f0<-function(X){
    Xa<-t(average.mat.rows(t(X),ids = ids,f = colMeans))
    Xa[Xa>10&!is.na(Xa)]<-10
    Xa[Xa<(-10)&!is.na(Xa)]<-(-10)
    return(Xa)
  }
  r1$zscoresAvF<-f0(r1$zscoresAv)
  return(r1)
}

PerturbSeq_sgRNA_zscoreBasedSimilarities<-function(r1,sig,main = "A375, monoculture",plot.flag = F,gene.selection = "de.level",
                                                   deg.level = 20,n1.genes = 1000, n2.genes = 4){
  X_GA.size<-PerturbSeq_get.GA.size(sig)
  batches<-unique(r1$samples)
  Z<-r1$zscoresAvF
  sig<-sig[laply(sig,length)>deg.level]
  v<-rowMeans(r1$NTC_avg)
  b<-rowSums(r1$NTC_sd>0)==length(batches)
  b1<-b&v>quantile(v[b],1-(n1.genes/sum(b)))
  b2<-b&is.element(rownames(Z),get.abundant(unlist(sig),abn.c = n2.genes))
  if(gene.selection == "both"){b<-b1|b2}
  if(gene.selection == "exp.level"){b<-b1}
  if(gene.selection == "de.level"){b<-b2}
  genes<-rownames(Z)[b]
  X<-cor(Z[b,])
  targets<-unique(get.strsplit(names(sig),sep = ".",idx = 1))
  f1<-function(x){
    b<-is.element(get.strsplit(colnames(X),"_",1),x)
    X1<-X[b,b]
    return(X1[upper.tri(X1)])
  }
  f2<-function(x){
    b<-is.element(get.strsplit(colnames(X),"_",1),x)
    X1<-X[b,!b]
    return(X1)
  }
  f.sum.per.sgRNA<-function(x1){
    b<-is.element(get.strsplit(colnames(X),"_",1),get.strsplit(x1,"_",1))
    b<-b&!is.element(colnames(X),x1)
    v<-sort(X[x1,b],decreasing = T)
    if(length(v)<3){v<-c(v,rep(NA,3-length(v)))}
    if(length(v)>10){v<-v[1:3]}
    return(v)
  }
  
  cor.same<-lapply(targets, f1)
  cor.different<-lapply(targets, f2)
  cor.NTCs<-f1("non")
  print(paste("No. of genes used =",length(genes)))
  rslts<-list(gene.selection = gene.selection,genes = genes,
              deg.level = deg.level,n1.genes = n1.genes, n2.genes = n2.genes,
              cor = X,cor.same = cor.same,cor.different = cor.different,cor.NTCs = cor.NTCs)
  rslts$sum<-cbind.data.frame(cor = c(unlist(cor.same),unlist(cor.different),cor.NTCs),
                              labels = c(ifelse(unlist(cor.same)<2,"Target-pairs",""),
                                         ifelse(unlist(cor.different)<2,"Different-target pairs"),
                                         ifelse(cor.NTCs<2,"NTC-pairs","")))
  b1<-rslts$sum$labels=="NTC-pairs"
  b2<-rslts$sum$labels=="Different-target pairs"
  b3<-rslts$sum$labels=="Target-pairs"
  rslts$pvals<-c(same.vs.other.targets = ranksum.test(rslts$sum$cor[b1],rslts$sum$cor[b3]),
                 same.vs.NTC = ranksum.test(rslts$sum$cor[b2],rslts$sum$cor[b3]))
  rslts$medians<-c(NTCs = median(rslts$sum$cor[b1]),
                   different.targets = median(rslts$sum$cor[b2]),
                   same.target = median(rslts$sum$cor[b3]))
  X1<-laply(rownames(X), f.sum.per.sgRNA)
  rownames(X1)<-rownames(X)
  colnames(X1)<-c("sgRNA1","sgRNA2","sgRNA3")
  X1<-as.data.frame(X1)
  X1$Target<-get.strsplit(rownames(X1),"_",1)
  X1$GA.size<-X_GA.size[match(X1$Target,names(X_GA.size))]
  X1$GA.size[is.na(X1$GA.size)]<-0
  X1$N.cells<-table(r1$prt)[rownames(X1)]
  rslts$per.sgRNA_sum<-X1
  rslts$plot<-gg.densityplot(rslts$sum$cor,rslts$sum$labels,xlab = "Pearson correlation coefficient (r)",
                             title = paste0(main,"\n",call_format.pval(rslts$pvals)))
  return(rslts)
}

PerturbSeq_get.GA.size<-function(sig){
  X1<-cbind.data.frame(sig = names(sig),n1 = laply(sig,length))
  X1$target<-get.strsplit(X1$sig,".",1)
  X2<-X1[,c("n1","n1")]
  X2<-as.matrix(X2)
  rownames(X2)<-NULL
  X3<-average.mat.rows(X2,ids = X1$target,f = colSums)[,1]
  return(X3)
}

# **************** Perturb-seq network analyses **************** # 

PerturbA375_convergingTargets<-function(sig,sigHits,go.env,lig2rec,MYC.trim = F){
  genes<-sig$genesA
  savefile<-get.file("PerturbA375_Network.rds")
  if(MYC.trim){
    sig$MYC.up<-sig$MYC.top.up
    sig$MYC.down<-sig$MYC.top.down
    savefile<-get.file("PerturbA375_Network_MYCtrim.rds")
  }
  sig<-sig[setdiff(names(sig),c("genesA","MYC.top.up","MYC.top.down"))]
  sig_noMYC<-sig[!grepl("MYC",names(sig))]
  
  rslts1<-PerturbSeq_convergingTargets(sig = sig,hitRes = sigHits$Res,hitSen = sigHits$Sen,
                                       genes = genes,overlap.removal = "strict",
                                       n1 = 3,n2 = 2,sig2 = sig_noMYC,name = "final")
  
  rslts1$lig2rec_perturb<-sig.lig2rec(lig2rec,sig = rslts1$sigT[c("Sen","Res")])
  rslts1$lig2rec_hits<-sig.lig2rec(lig2rec,sig = sigHits[c("Sen","Res")])
  genes_ligRec<-unique(intersect(c(lig2rec$Ligand,lig2rec$Receptor),genes))
  rslts1$ligRec_enrich<-GO.enrichment.lapply(rslts1$lig2rec_hits$sig,genes_ligRec,rslts1$lig2rec_perturb$sig)
  rslts1$ligRec_both<-intersect.lists(rslts1$lig2rec_hits$sig,rslts1$lig2rec_perturb$sig)
  rslts1$ligRec_both<-lapply.unique.sort(rslts1$ligRec_both)
  
  # Genes that are both resistance hits and interact with resistance (CRISPRa and Perturb-Seq) hits
  sigRes<-list(CRISPRa = sigHits$Res,PerturbSeq = rslts1$sigT$Res)
  g<-intersect(rslts1$ligRec_both$res,sigHits$Res)
  X1<-rslts1$lig2rec_hits$Res.net[is.element(rslts1$lig2rec_hits$Res.net$`Ligand/receptor`,g),]
  X2<-rslts1$lig2rec_perturb$Res.net[is.element(rslts1$lig2rec_perturb$Res.net$`Ligand/receptor`,g),]
  X<-rbind(X1,X2);X<-unique(X)
  X<-cbind.data.frame(X,g1 = list.2.ids(X$Res.gene,sigRes),g2 = list.2.ids(X$`Ligand/receptor`,sigRes))
  colnames(X)<-c("G1","G2","G1.type","G2.type")
  X$G1.type[grepl("&",X$G1.type)]<-"Both"
  X$G2.type[grepl("&",X$G2.type)]<-"Both"
  X$G1.type<-paste0("Res.",X$G1.type)
  X$G2.type<-paste0("Res.",X$G2.type)
  X$edge<-"LigandReceptor"
  idx1<-c(paste0(sigHits$Res,".up"),paste0(sigHits$Sen,".down")) # Resistance
  X2<-melt(intersect.list1(sig[idx1],unique(X$G1)))
  colnames(X2)<-c("G1","G2")
  X2$edge<-paste0("Reg.",get.strsplit(X2$G2,".",2))
  X2$G2<-get.strsplit(X2$G2,".",1)
  X2$G1.type<-X$G1.type[match(X2$G1,X$G1)]
  X2$G2.type<-ifelse(is.element(X2$G2,sigHits$Res),"Res","Other")
  X2$G2.type[is.element(X2$G2,sigHits$Sen)]<-"Sen"
  X2<-X2[X2$G1!=X2$G2,]
  X<-rbind(X,X2[,colnames(X)])
  rslts1$lig2rec_subNetworkRes<-X
  saveRDS(rslts1,file = savefile)
  return(rslts1)
}

PerturbSeq_convergingTargets<-function(sig,n1,hitRes,hitSen,genes,sig2,n2,overlap.removal="strict",score.flag = F,name = ""){
  idx1<-c(paste0(hitRes,".up"),paste0(hitSen,".down")) # Resistance
  idx2<-c(paste0(hitSen,".up"),paste0(hitRes,".down")) # Sensitizing
  t1<-get.abundant(unlist(sig[idx1]),abn.c = n1)
  t2<-get.abundant(unlist(sig[idx2]),abn.c = n1)
  g<-unique(sort(unlist(sig)))
  X<-cbind.data.frame(v1 = table(unlist(sig[idx1]))[g],v2 = table(unlist(sig[idx2]))[g])
  X$v1.Freq[is.na(X$v1.Freq)]<-0
  X$v2.Freq[is.na(X$v2.Freq)]<-0
  X$scores<-X$v1.Freq-X$v2.Freq
  rownames(X)<-g
  X.scores<-cbind.data.frame(X,list.2.ids(rownames(X),list(Res = hitRes,Sen = hitSen)))
  if(score.flag){
    t1<-g[X$scores<(-n1)]
    t2<-g[X$scores>n1]
  }
  if(!missing(sig2)){
    t1<-unique(c(t1,get.abundant(unlist(sig2[idx1]),abn.c = n2)))
    t2<-unique(c(t2,get.abundant(unlist(sig2[idx2]),abn.c = n2)))
  }
  t12<-intersect(t1,t2)
  
  if(overlap.removal=="permissive"){
    t1<-setdiff(t1,t12)
    t2<-setdiff(t2,t12)
  }
  if(overlap.removal=="strict"){
    t1<-setdiff(t1,unlist(sig[idx2]))
    t2<-setdiff(t2,unlist(sig[idx1]))
  }
  
  sigT<-list(ResA = t1,SenA = t2, Both = t12)
  sigT$Res<-setdiff(sigT$ResA,hitSen)
  sigT$Sen<-setdiff(sigT$SenA,hitRes)
  
  # Supplementary Tables and input for Cytoscape
  X<-rbind(cbind.data.frame(melt(intersect.list1(sig,intersect(sigT$Sen,hitSen))),target = "Sen"),
           cbind.data.frame(melt(intersect.list1(sig,intersect(sigT$Res,hitRes))),target = "Res"))
  X$prt<-get.strsplit(X$L1,".",1)
  X$d<-get.strsplit(X$L1,".",2)
  X$prt.type<-ifelse(is.element(X$prt,hitRes),"Res","Sen")
  table(X$prt.type,is.element(X$prt,hitSen))
  X.node.anno<-melt(list(Sen = hitSen, Res = hitRes))
  
  rslts<-list(sigT = sigT,
              network = X,
              network_anno = X.node.anno,
              overlap.removal = overlap.removal,score.flag = score.flag, n1 = n1,
              hitRes = hitRes, hitSen = hitSen, sig = sig,
              sum = X.scores)
  if(!missing(sig2)){rslts$sig2<-sig2;rslts$n2<-n2}
  return(rslts)
}

sig.lig2rec<-function(lig2rec,sig,sigIDs,view.flag = F){
  X1a<-lig2rec[is.element(lig2rec[,2],sig$Sen),c(2,1,3)]
  X1b<-lig2rec[is.element(lig2rec[,1],sig$Sen),]
  colnames(X1a)<-c("Sen.gene","Ligand/receptor","Pair")
  colnames(X1b)<-c("Sen.gene","Ligand/receptor","Pair")
  X1<-unique(rbind.data.frame(X1a,X1b)[,1:2])
  
  X2a<-lig2rec[is.element(lig2rec[,2],sig$Res),c(2,1,3)]
  X2b<-lig2rec[is.element(lig2rec[,1],sig$Res),]
  colnames(X2a)<-c("Res.gene","Ligand/receptor","Pair")
  colnames(X2b)<-c("Res.gene","Ligand/receptor","Pair")
  X2<-unique(rbind.data.frame(X2a,X2b)[,1:2])
  
  if(!missing(sigIDs)){
    X1<-cbind.data.frame(X1,list.2.ids(X1$`Ligand/receptor`,sigIDs))
    X2<-cbind.data.frame(X2,list.2.ids(X2$`Ligand/receptor`,sigIDs))
  }
  if(view.flag){
    View(X1,title = "Sensitizing ligand-receptor")
    View(X2,title = "Resistance ligand-receptor")
  }
  sig<-list(sen = unique(X1$`Ligand/receptor`),
            res = unique(X2$`Ligand/receptor`))
  return(list(Sen.net = X1,Res.net = X2,sig = sig))
}
