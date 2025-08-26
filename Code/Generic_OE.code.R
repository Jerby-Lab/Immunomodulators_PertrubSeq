get.OE<-function(r,sig,method2= F,min.n.genes){
  if(!missing(min.n.genes)){sig<-sig[laply(sig,length)>min.n.genes]}
  if(method2){
    scores<-get.OE1.method2(r,sig)
  }else{
    scores<-get.OE1(r,sig)
  }
  names(sig)<-gsub(" ",".",names(sig))
  two.sided<-unique(gsub(".up","",gsub(".down","",names(sig))))
  b<-is.element(paste0(two.sided,".up"),names(sig))&
    is.element(paste0(two.sided,".down"),names(sig))
  if(any(b)){
    two.sided<-two.sided[b]
    scores2<-as.matrix(scores[,paste0(two.sided,".up")]-scores[,paste0(two.sided,".down")])
    colnames(scores2)<-two.sided
    scores<-cbind(scores2,scores)
  }
  
  if(!is.null(r$cells)){
    rownames(scores)<-r$cells
  }else{
    if(!is.null(r$samples)){
      rownames(scores)<-r$samples 
    }
  }
  return(scores)
}

prep4OE<-function(r,n.cat = 50){
  r$zscores<-center.matrix(r$tpm,dim = 1,sd.flag = T)
  X<-10*((2^r$tpm)-1)
  r$genes.dist<-log2(rowMeans(X,na.rm = T)+1)
  r$genes.dist.q<-discretize.prvt(r$genes.dist,n.cat = n.cat)
  b<-rowSums(is.na(r$zscores))==0
  if(any(!b)){r<-set.list(r,b)}
  r$binZ<-average.mat.rows(r$zscores,r$genes.dist.q,f = colMeans)
  return(r)
}

get.OE1 <- function(r,sig){
  if(is.list(sig)){
    scores<-t(laply(sig, function(g) get.OE1(r,g)))
    rownames(scores)<-r$cells
    colnames(scores)<-names(sig)
    return(scores)
  }
  g<-sig
  b<-is.element(r$genes,g)
  n1<-average.mat.rows(1*cbind(b,b),paste0("Bin",r$genes.dist.q),f = colSums)
  rand.scores<-t(r$binZ)%*%n1[,1]
  if(sum(b)==1){
    raw.scores<-r$zscores[b,]
  }else{
    raw.scores<-colSums(r$zscores[b,])
  }
  scores<-(raw.scores-rand.scores)/sum(b)
  return(scores)
}

OE.fix.scores<-function(X0){
  b2<-!grepl(".",colnames(X0),fixed = T)
  X<-subset(x = as.matrix(X0),select = which(b2))
  b1<-!is.element(get.strsplit(colnames(X0),".",1),colnames(X))
  X<-cbind.data.frame(X0[,b2],X0[,b1])
  colnames(X)<-c(colnames(X0)[b2],colnames(X0)[b1])
  X[,grepl("down",colnames(X))]<-(-X[,grepl("down",colnames(X))])
  colnames(X)<-get.strsplit(colnames(X),".",1)
  return(X)
}

sig.remove.target.genes<-function(sig,sep = "_",idx = 1){
  targets<-get.strsplit(names(sig),sep = sep,idx = idx)
  sig1<-lapply(1:length(sig), function(x){return(setdiff(sig[[x]],targets[x]))})
  names(sig1)<-names(sig)
  print(setdiff.lists.by.idx(sig,sig1))
  return(sig1)
}
