# GitHub Code, Akana, Yoe et al., 2025
# Generic code for plotting in situ Perturb-Seq data.

#************** General plotting of in situ Perturb-Seq data **************#

insitu_plotLayers<-function(X,labels,b1,b2,cex = 0.4, main = "",col1,cex.axis = 1,cex.main = 3,
                            ylab = "Spatial Coordinate 1",xlab = "Spatial Coordinate 2",set.flag = F){
  if(missing(col1)){
    col1<-labels.2.colors(labels,color.spec = "rgb")
  }
  call_plot(X[,1],X[,2],cex = cex,icol = col1,pch = 16,labels = labels,
            main = main,cex.axis = cex.axis,cex.main = cex.main,
            ylab = ylab,xlab = xlab,set.flag = set.flag)
  points(X[b1,1],X[b1,2],cex = cex,col = col1[b1],pch = 16)
  points(X[b2,1],X[b2,2],cex = cex,col = col1[b2],pch = 16)
  return()
}

insituPlot<-function(r,r1,target = "CD44",nplots = 2,width = 15,q1 = 0.9,q2 = 0.9,cex = 0.3,plotIdx){
  b<-r$targets==target
  col1<-ifelse(b,"red",ifelse(r$targets=="None","deepskyblue","black"))
  if(plotIdx[1]){
    insitu_plotLayers(r$coor,labels = paste(ifelse(b,"","No"),target,"ORF"),cex = cex,
                      main = paste(target,"ORF detection"),
                      b1 = r$targets!="None",b2 = b,col1 = col1,
                      cex.axis = 1,cex.main = 2,
                      ylab = "Spatial Coordinate 1",xlab = "Spatial Coordinate 2")
  }
  if(plotIdx[2]){
    plot.layers(r1$coor,labels = r1$scores[,target],cex.main = 2,
                b.top = b,cex = cex,legend.flag = F,main = paste(target,"GA"),cex.axis = 1e-10)
  }
  if(plotIdx[3]){
    plot.layers(r1$coor,labels = r1$scores[,target]>quantile(r1$scores[,target],q1),cex.main = 2,
                b.top = b,cex = cex,legend.flag = F,cex.axis = 1e-10,main = paste(target,"GA"),
                ylab = "Spatial Coordinate 1",xlab = "Spatial Coordinate 2")
  }
  if(plotIdx[4]){
    plot.layers(r1$coor,labels = r1$tpm[target,]>quantile(r1$tpm[target,],q2),cex.main = 2,
                b.top = b,cex = cex,legend.flag = F,main = target,cex.axis = 1e-10,
                ylab = "Spatial Coordinate 1",xlab = "Spatial Coordinate 2")
  }
  return()
}


#************** MCP plots **************#

GA.MCP_barplot<-function(sigA,main = ""){
  idx<-sort(unique(unlist(lapply(sigA,function(x) return(names(x))))))
  idx<-unique(get.strsplit(idx,".",1))
  X<-union.multiple.mats(lapply(sigA, function(x) summarize.program(x,idx)))
  X<-X[order(rowSums(X),decreasing = T),]
  colnames(X)<-paste0(gsub("."," (",colnames(X),fixed = T),")")
  call_barplot(X,beside = F,xlab = "Perturbation",
               ylab = "No. DEGs",main = main)
  return(X)
}

GA.MCP_boxplot<-function(rslts,cell.type,x,q1 = 0.1){
  if(length(x)>1){
    p<-lapply(x,function(x1) GA.MCP_boxplot(rslts,cell.type = cell.type,x = x1,q1 = q1))
    return(p)
  }
  y1<-discretize.3.labels(rslts$frames.mal2[,x],q = q1)
  p<-call_boxplot(rslts$scores2[,x],y1,labels = y1,title.size = 10,
                  lab.size = 10,ylab = paste("GA expression,",cell.type),
                  main = x,unique.x = c("Low","Moderate","High"),
                  xlab = paste0(x,"\nabundance"),remove.legend = T)
  p<-p+stat_compare_means(comparisons = list(c("Low","Moderate"),c("Moderate","High"),c("Low","High")),
                          label = "p.signif",method = "t.test")+theme(legend.position="none")
  # p<-p+theme(aspect.ratio = 0.4)
  return(p)
}

summarize.program<-function(sig,idx,plot.flag = F,main = "",horiz = F){
  if(missing(idx)){idx<-sort(unique(get.strsplit(names(sig),".",1)))}
  X<-cbind.data.frame(up = laply(sig[paste0(idx,".up")],length),
                      down = laply(sig[paste0(idx,".down")],length))
  rownames(X)<-idx
  if(!plot.flag){return(X)}
  X<-X[order(rowSums(X),decreasing = !horiz),]
  colnames(X)<-c("Up regulated","Down regulated")
  call_barplot(X,beside = F,xlab = "Perturbation",ylab = "No. DEGs",
               icol = c("gray6","cadetblue"),main = main,horiz = horiz)
  return(X)
}

#************** Segmentation in situ plots **************#

inSituPlotSegmentation1 <- function (seg_path,celltypes,cell2rgb,samplename,outpath = "~/",outfile = "out.jpg",
                                     background = "black",cont_field = "",low_qc_color = 0,contvals = NULL,QC.flag = F,trans.flag = T){
  cellseg = read.csv(seg_path)
  
  if(QC.flag){
    v<-unique(unlist(cellseg))
    print(paste(length(setdiff(paste0("c",v),names(celltypes))),"cells not in the cell type vector"))
  }
  colnames(cellseg) <- unlist(lapply(colnames(cellseg), function(x) {
    strsplit(x, split = "X0.")[[1]][2]
  }))
  colnames(cellseg)[1] <- "0"
  if(trans.flag){
    cellseg<-cellseg[seq(nrow(cellseg),1,-1),seq(ncol(cellseg),1,-1)]
    cellseg<-cellseg[seq(nrow(cellseg),1,-1),seq(ncol(cellseg),1,-1)]
    print(dim(cellseg))
    cellseg<-t(cellseg)
  }
  
  cellmask <- cellseg
  cellmask[cellmask != 0] <- low_qc_color
  if (background == "white") {
    cellmask[cellmask == 0] <- 1
  }
  else {
    cellmask[cellmask == 0] <- 0
  }
  cellmask <- EBImage::Image(as.matrix(cellmask))
  cellmask <- EBImage::channel(cellmask, "rgb")
  levels <- setdiff(unique(celltypes), cont_field)
  for (type in levels) {
    cellids = names(celltypes[celltypes == type])
    cellidx = unlist(lapply(cellids, function(x) as.integer(strsplit(x,
                                                                     split = "c")[[1]][2])))
    celltype_mask <- t(apply(cellseg, 1, function(x) {
      x %in% cellidx
    }))
    celltype_rgb_1 <- cellmask@.Data[, , 1]
    celltype_rgb_1[celltype_mask] <- as.numeric(cell2rgb[[type]][1])/255
    cellmask[, , 1] <- celltype_rgb_1
    celltype_rgb_2 <- cellmask@.Data[, , 2]
    celltype_rgb_2[celltype_mask] <- as.numeric(cell2rgb[[type]][2])/255
    cellmask[, , 2] <- celltype_rgb_2
    celltype_rgb_3 <- cellmask@.Data[, , 3]
    celltype_rgb_3[celltype_mask] <- as.numeric(cell2rgb[[type]][3])/255
    cellmask[, , 3] <- celltype_rgb_3
  }
  if (cont_field != "") {
    cellids = names(celltypes[celltypes == cont_field])
    cellidx = unlist(lapply(cellids, function(x) as.integer(strsplit(x,
                                                                     split = "c")[[1]][2])))
    celltype_mask <- t(apply(cellseg, 1, function(x) {
      x %in% cellidx
    }))
    cellids <- paste0(samplename, "_c", cellseg[celltype_mask])
    colormap <- t(sapply(unique(contvals), simplify = T,
                         grDevices::col2rgb))/255
    celltype_rgb_1 <- cellmask@.Data[, , 1]
    celltype_rgb_1[celltype_mask] <- colormap[contvals[cellids],
                                              1]
    cellmask[, , 1] <- celltype_rgb_1
    celltype_rgb_2 <- cellmask@.Data[, , 2]
    celltype_rgb_2[celltype_mask] <- colormap[contvals[cellids],
                                              2]
    cellmask[, , 2] <- celltype_rgb_2
    celltype_rgb_3 <- cellmask@.Data[, , 3]
    celltype_rgb_3[celltype_mask] <- colormap[contvals[cellids],
                                              3]
    cellmask[, , 3] <- celltype_rgb_3
  }
  
  EBImage::writeImage(cellmask, files = paste0(outpath,"/",outfile))
}


