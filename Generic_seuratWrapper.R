# GitHub Code, Akana, Yoe et al., 2025
# Seurat wrapper for clustering, embedding, and analyzing single-cell data.

#************** Seurat Wrapper **************#

seuratW_get.embedding<-function(r,no.genes = 2000,n.pcs = 10,cd.flag = T,PCs = NULL,find.var = T,
                               resolution = 0.4,norm.flag = T,tsne.flag = F,cluster.flag = T,
                               return.model = F,umap.flag = T,full.flag = F,do.scale = T,do.center = T,seed = 1234){
  
  D1<-seuratW_process(r = r,
                     no.genes = no.genes,
                     scores = NULL,
                     plot.flag = T,
                     n.pcs = n.pcs,
                     find.var = find.var,
                     cd.flag = cd.flag,
                     PCA.approx = T,
                     resolution = resolution,
                     cluster.flag = cluster.flag,
                     norm.flag = norm.flag,
                     tsne.flag = tsne.flag,
                     umap.flag = umap.flag,
                     PCs = PCs,
                     do.scale = do.scale,
                     do.center = do.center,seed = seed,
                     return.model = return.model)
  
  idx<-call_match(r$cells,rownames(D1@meta.data));table(is.na(idx))
  r$pca<-Embeddings(object = D1, reduction = "pca")[idx,]
  r$pca.load<-Loadings(object = D1,reduction = "pca")
  
  if(cluster.flag){r$clusters<-paste0("C",FetchData(D1,vars = "seurat_clusters")[idx,])}
  if(tsne.flag){r$tsne<-Embeddings(object = D1, reduction = "tsne")[idx,]}
  if(umap.flag){r$umap<-Embeddings(object = D1, reduction = "umap")[idx,]}
  if(full.flag){r$seurat<-D1}
  return(r)
}

seuratW_process<-function(r,no.genes = 1000, n.pcs = 10,tsne.flag = T,umap.flag = T,find.var = T,do.scale = T,do.center = T,
                         cluster.iter = F,cd.flag = F,scores = NULL,plot.flag = F,fileName = NULL,cluster.flag = T,
                         D1,resolution = 0.4,norm.flag = T,tsne.method = "Rtsne",PCA.approx = T,PCs,seed = 1234,
                         return.model = F){
  set.seed(seed)
  print(paste(n.pcs,"PCs used."))
  if(length(n.pcs)==1){n.pcs<-1:n.pcs}
  if(missing(D1)){
    if(cd.flag){D1.data<-round(r$cd,2)}else{D1.data<-round(r$tpm,2)}
    D1<-seuratW_createObject(r$genes,D1.data,r$sampleName,norm.flag = norm.flag,
                             find.var = find.var,n.genes = no.genes,
                             do.scale = do.scale,do.center = do.center)
  }
  
  if(find.var){
    if(is.element("meta.data",slotNames(D1@assays$RNA))){
      X.var<-D1@assays$RNA@meta.data
      hvg.D1 <-unique(X.var$var.features[!is.na(X.var$var.features)])
    }else{
      hvg.D1<-D1@assays$RNA@var.features # changed on 03.31.2024 
    }
  }else{
    hvg.D1<-rownames(D1.data)
  }
  
  D1 <- RunPCA(object = D1,features = hvg.D1,npcs = max(n.pcs),approx = PCA.approx)
  X<-Embeddings(object = D1, reduction = "pca")[,n.pcs]
  # v<-Loadings(object = D1, reduction = "pca")
  
  if(!missing(PCs)&!is.null(PCs)){
    print("Using previous PCs!")
    idx.pca<-match(rownames(D1@reductions$pca@cell.embeddings),
                   rownames(PCs$cell.embeddings))
    D1@reductions$pca@feature.loadings<-PCs$feature.loadings
    D1@reductions$pca@cell.embeddings<-PCs$cell.embeddings[idx.pca,]
  }
  
  b<-duplicated(X)
  if(any(b)){
    print(paste("Removing",sum(b),"cells."))
    D1<-subset(D1,cells = rownames(X)[!b])
    # D1<-subset(D1,cells = rownames(D1@meta.data)[!b])
  }
  
  if(tsne.flag){D1 <- RunTSNE(object = D1, reduction.use = "pca",dims = n.pcs, do.fast = TRUE)}
  if(umap.flag){D1<-RunUMAP(D1,dims = n.pcs,return.model = return.model)}
  
  if(cluster.flag){
    D1 <- FindNeighbors(object = D1,reduction = "pca",dims = n.pcs)
    D1 <- FindClusters(D1,resolution = resolution,return.model = return.model)
    # reduction.type = "pca",
    # dims.use = 1:n.pcs, save.SNN = T, resolution = resolution)
  }
  
  if(cluster.iter){
    while(length(unique(D1@ident))==1){
      resolution<-resolution+0.05
      print(paste("Clustering with resolution",resolution))
      D1 <- FindClusters(D1,resolution = resolution)
    }
    D1@resolution<-resolution
  }
  
  idx<-call_match(rownames(D1@meta.data),colnames(D1.data))
  if(!is.null(r$samples)){D1@meta.data<-cbind.data.frame(D1@meta.data,samples = r$samples[idx])}
  if(!is.null(scores)){D1@meta.data<-cbind.data.frame(D1@meta.data,scores[idx,])}
  return(D1)
}

seuratW_createObject<-function(genes,D1.data,D1.name,find.var = T,norm.flag = T,n.genes,do.scale = T,do.center = T){
  print(paste("Processing",D1.name,"..."))
  # D1.data<-round(D1.data,2)
  b<-!duplicated(colnames(D1.data))
  if(any(!b)){
    D1.data<-D1.data[genes,!duplicated(colnames(D1.data))];dim(D1.data)
  }else{
    D1.data<-D1.data[genes,];dim(D1.data)
  }
  
  D1 <- CreateSeuratObject(counts = D1.data)
  if(norm.flag){
    D1 <- NormalizeData(object = D1)
    D1 <- ScaleData(object = D1,do.scale = do.scale,do.center = do.center)
  }
  
  if(find.var){D1 <- FindVariableFeatures(object = D1, do.plot = FALSE,nfeatures = n.genes)}
  # D1@meta.data$source <- D1.name
  return(D1)
}