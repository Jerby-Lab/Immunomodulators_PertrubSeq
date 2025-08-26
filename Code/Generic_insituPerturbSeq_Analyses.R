# GitHub Code, Akana, Yoe et al., 2025
# Generic code for Perturb-seq and in situ Perturb-seq data analyses.

GA.combine_HLM<-function(rslts1,rslts2,x){
  if(missing(x)){
    rslts<-lapply(names(rslts1$full), function(x) GA.combine_HLM(rslts1 = rslts1,rslts2 = rslts2,x = x))
    names(rslts)<-names(rslts1$full)
    sig<-unlist(lapply(rslts, function(x) x$sig),recursive = F)
    summary(sig)
    rslts$sig<-sig
    return(rslts)
  }
  
  X<-union.multiple.mats(list(X1 = rslts1$full[[x]][,c("up","down")],
                              X2 = rslts2$full[[x]][,c("up","down")]))
  X<-cbind.data.frame(up = p.adjust(fisher.combine(X[,c("X1.up","X2.up")]),method = "BH"),
                      down = p.adjust(fisher.combine(X[,c("X1.down","X2.down")]),method = "BH"),
                      X)
  sig<-get.top.elements(X[,c("up","down")],q = 250,min.ci = 0.05,sort.flag = F)
  return(list(X = X,sig = sig))
}

GA.combine_MAST<-function(rslts1,rslts2,x){
  if(missing(x)){
    f<-function(X) return(cbind.data.frame(up = get.pval.from.zscores(X$zscores),down = get.pval.from.zscores(-X$zscores),X))
    rslts1$full<-lapply(rslts1$full,f)
    rslts2$full<-lapply(rslts2$full,f)
    rslts<-lapply(names(rslts1$full), function(x) GA.combine_MAST(rslts1 = rslts1,rslts2 = rslts2,x = x))
    names(rslts)<-names(rslts1$full)
    sig<-unlist(lapply(rslts, function(x) x$sig),recursive = F)
    summary(sig)
    rslts$sig<-sig
    return(rslts)
  }
  
  X<-union.multiple.mats(list(X1 = rslts1$full[[x]][,c("up","down")],
                              X2 = rslts2$full[[x]][,c("up","down")]))
  X<-cbind.data.frame(up = p.adjust(fisher.combine(X[,c("X1.up","X2.up")]),method = "BH"),
                      down = p.adjust(fisher.combine(X[,c("X1.down","X2.down")]),method = "BH"),X)
  sig<-get.top.elements(X[,c("up","down")],q = 1000,min.ci = 0.1,sort.flag = F)
  return(list(X = X,sig = sig))
}

PerturbSeq_DEGs_MAST<-function(r = NULL,abn.c = 20,latent.vars,outputfile,
                                  per.guide = F,rslts,return.flag = T){
  name1<-r$name
  if(!missing(latent.vars)){name1<-paste0(r$name,"_wBatchControl")}
  if(per.guide){name1<-paste0(name1,"_perGuide")}else{name1<-paste0(name1,"_perTarget")}
  if(!missing(outputfile)&&!overwrite){
    if(file.exists(outputfile)){
      print(paste("Loading",outputfile))
      rslts<-readRDS(outputfile)
      if(return.flag){return(rslts)}
    }
  }
  
  if(is.null(r$seurat)){
    r<-seuratW_get.embedding(r,no.genes = 5000,n.pcs = 10,cd.flag = T,umap.flag = T,
                             norm.flag = T,full.flag = T)
  }
  
  if(per.guide){
    print("Running per target.")
    v<-r$guides
    v[r$b.ctrl]<-"NTC"
    r$targets<-v
  }
  D1<-r$seurat
  Idents(D1) <- r$targets
  print(table(Idents(D1)))
  
  r$genes.dr<-t(average.mat.rows(t(r$tpm),ids = r$targets))
  targets<-setdiff(get.abundant(r$targets,abn.c,decreasing = F),c("NTC"))
  X.zscores<-get.mat(r$genes,targets)
  X.logFC<-get.mat(r$genes,targets)
  X.delta<-get.mat(r$genes,targets)
  X<-list()
  if(!missing(rslts)){
    X.zscores<-rslts$zscores
    X.logFC<-rslts$logFC
    X<-rslts$full
  }
  n1<-1
  for(x in setdiff(targets,names(X))){
    print(x)
    if(sum(!is.na(X.zscores[,x]))>0){next()}
    if(!missing(latent.vars)){
      print(paste("Using",latent.vars,"as confounders."))
      X1<-FindMarkers(D1,ident.1 = x,ident.2 = "NTC",test.use = "MAST",logfc.threshold = 0,latent.vars = latent.vars)
    }else{
      X1<-FindMarkers(D1,ident.1 = x,ident.2 = "NTC",test.use = "MAST",logfc.threshold = 0)
    }
    
    X1$fdr<-p.adjust(X1$p_val,method = "BH")
    X1$fdr[X1$fdr==0]<-1e-100
    X1$zscores<-get.cor.zscores(X1$avg_log2FC,X1$fdr)
    X1<-cbind.data.frame(X1,DR = r$genes.dr[rownames(X1),c("NTC",x)])
    X1<-cbind.data.frame(X1,delta = X1$DR.NTC-X1[,paste0("DR.",x)])
    X[[x]]<-X1
    
    idx<-match(r$genes,rownames(X1))
    X1<-X1[idx,]
    X.zscores[,x]<-X1[,"zscores"]
    X.logFC[,x]<-X1[,"avg_log2FC"]
    X.delta[,x]<-r$genes.dr[,x]-r$genes.dr[,"NTC"]
    n1<-n1+1
    if((n1%%10)==0){
      rslts<-list(name = name1,full = X, zscores = X.zscores,
                  logFC = X.logFC,delta = X.delta)
    }
  }
  
  sig<-get.top.cor(X.zscores,q = 1000,min.ci = 1)
  rslts<-list(name = name1,full = X, zscores = X.zscores,
              logFC = X.logFC,delta = X.delta,sig = sig)
  r$seurat<-NULL
  r$targetsFC<-scFoldChange(r = r,targets = r$targets)
  r<-prep4OE(r,n.cat = 100)
  sig<-sig[laply(sig,length)>10]
  r$scores<-OE.fix.scores(get.OE(r,sig))
  if(!is.null(r$guides)){
    rslts$scores<-cbind.data.frame(name = name1,cells = r$cells,targets = r$targets,guides = r$guides,scores = r$scores)
  }else{
    rslts$scores<-cbind.data.frame(name = name1,cells = r$cells,targets = r$targets,scores = r$scores)
  }
  
  rslts<-scFoldChange_fix(logFC_new = r$targetsFC$logFC,rslts = rslts)
  if(per.guide){
    rslts$sum<-get_sgRNA_GA_consistency(sig = rslts$sig,jac.flag = 2,
                                        g = rownames(rslts$zscores),name = rslts$name)
  }
  if(missing(outputfile)){return(rslts)}
  print("Saving results...")
  saveRDS(rslts,file = outputfile)
  return(rslts)
}

insituPerturb_DEG_HLM<-function(r,outputfile,formula = "y ~ (1 | fovs) + x + comp"){
  if(!missing(outputfile)&&file.exists(outputfile)&&!overwrite){
    print(paste("Loading",outputfile))
    return(readRDS(outputfile))}
  targets<-setdiff(sort(unique(r$targets)),c("NTC","Empty"))
  if(sum(is.element(c("NTC","Empty"),r$targets))==0){return("Error - Empty ORF not found.")}
  R<-list()
  zscoresF<-get.mat(m.rows = r$genes,m.cols = targets)
  zscores<-get.mat(m.rows = r$genes,m.cols = targets)
  for(x in targets){
    r1<-set.list(r,is.element(r$targets,c(x,"NTC","Empty")))
    r1$comp<-r1$comp/max(r1$comp)
    X<-apply.formula.HLM(r = r1,X = r1$targets==x,Y = r1$tpm,MARGIN = 1,formula = formula)
    X<-cbind.data.frame(up = get.pval.from.zscores(X[,"Z"]),
                        down = get.pval.from.zscores(-X[,"Z"]),X)
    zscoresF[,x]<-X[r$genes,"Z"]
    zscores[,x]<-X[r$genes,"Z"]
    zscores[X[r$genes,"singular"]==1,x]<-NA
    R[[x]]<-X
  }
  
  rslts<-list(name = r$name,full = R,
              zscores = zscores,
              zscoresF = zscoresF,
              formula = "y ~ (1 | fovs) + x + comp",
              sig = get.top.cor(zscores,q = 500,min.ci = 2))
  
  rslts$BH.up<-p.adjust.mat(get.pval.from.zscores(rslts$zscores))
  rslts$BH.down<-p.adjust.mat(get.pval.from.zscores(-rslts$zscores))
  sig1<-get.top.elements(rslts$BH.up,min.ci = 0.1,q = 500);names(sig1)<-paste0(names(sig1),".up")
  sig2<-get.top.elements(rslts$BH.down,min.ci = 0.1,q = 500);names(sig2)<-paste0(names(sig2),".down")
  rslts$BH.sig<-c(sig1,sig2)
  if(missing(outputfile)){return(rslts)}
  saveRDS(rslts,file = outputfile)
  return(rslts)
}

OE.pred<-function(r,CV = F,scores = r$scores){
  if(CV){
    X1<-OE.pred(r,CV = F,scores = r$scores1)
    X2<-OE.pred(r,CV = F,scores = r$scores2)
    X<-union.multiple.mats(list(train = X1,test = X2))
    return(X)
  }
  
  targets<-colnames(scores)
  X<-laply(targets,function(x) OE.pred1(r,x,scores = scores))
  rownames(X)<-targets
  colnames(X)<-c("AUROC","ttest")
  return(X)
}

OE.pred1<-function(r,x,scores = r$scores){
  b<-is.element(r$targets,c(x,"NTC","Empty"))&!is.na(scores[,x])
  # print(table(is.na(scores[b,x])))
  # print(table(is.na(b)))
  if(sum(b)<10){return(c(NA,NA))}
  p<-c(get.auc(scores[b,x],r$targets[b]==x),t.test.labels(scores[b,x],r$targets[b]==x))
  return(p)
}

scFoldChange<-function(r,targets = r$targets,full.flag = F){
  r$tpmAv1<-t(average.mat.rows(t(r$tpm),ids = targets))
  r$cdAv<-t(average.mat.rows(t(r$cd),ids = targets))
  r$comp.readsAv<-colSums(r$cdAv)
  r$cdAvNorm<-(1000000*(sweep(r$cdAv,2,r$comp.readsAv,FUN = '/')))+1
  r$tpmAv2<-log2(r$cdAvNorm)
  if(!is.element("NTC",colnames(r$cdAvNorm))){
    if(is.null(r$targetsFC$cdAvNorm)){print("Need target level information.");return(r)}
    r$logFC<-log2(sweep(r$cdAvNorm,1,r$targetsFC$cdAvNorm[,"NTC"],FUN = '/'))
  }else{
    r$logFC<-log2(sweep(r$cdAvNorm,1,r$cdAvNorm[,"NTC"],FUN = '/'))
  }
  idx<-c('tpmAv1','cdAv','comp.readsAv','cdAvNorm','tpmAv2','logFC')
  if(full.flag){return(r)}
  return(r[idx])
}

scFoldChange_fix<-function(logFC_new,rslts,r){
  rslts$old_zscores<-rslts$zscores
  rslts$old_logFC<-rslts$logFC
  rslts$old_sig<-get.top.cor(rslts$zscores,q = 10000,min.ci = 2)
  Z <- rslts$zscores
  if(missing(logFC_new)){
    r<-scFoldChange(r)
    logFC_new<-r$logFC
  }
  rslts$logFC<-logFC_new[rownames(Z),colnames(Z)]
  print(table(rslts$logFC<0,Z<0,abs(Z)>3))
  Z<-abs(Z)*sign(rslts$logFC)
  print(table(rslts$logFC<0,Z<0,abs(Z)>3))
  rslts$zscores<-Z
  rslts$sig<-get.top.cor(rslts$zscores,q = 10000,min.ci = 2)
  return(rslts)
}





