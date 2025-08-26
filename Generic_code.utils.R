# GitHub Code, Akana, Yoe et al., 2025
# Generic code needed to regenerate the figures and re-run the analyses.

#************** Generic code **************#

umap.ggplot<-function(umapX,labels,labels.name = "",main = "",size = 0.2,xlim1, ylim1,reorder.flag = F){
  if(is.matrix(labels)|is.data.frame(labels)){
    p<-lapply(colnames(labels),function(x){
      print(x)
      p1<-umap.ggplot(umapX,labels = labels[,x],labels.name = labels.name,main = x)
      return(p1)
    })
    names(p)<-colnames(labels)
    return(p)
  }
  if(reorder.flag){
    print("reordering")
    idx<-order(labels)
    labels<-labels[idx]
    umapX<-umapX[idx,]
  }
  
  xylabs <- colnames(umapX)
  colnames(umapX)<-c("UMAP1","UMAP2")
  X <- cbind.data.frame(umapX,col = labels)
  p <- ggplot(X, aes(x = UMAP1, y = UMAP2, color = col)) + geom_point(size = size) + labs(color = labels.name)
  p <- p + ggtitle(main)
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    xlab(xylabs[1]) + ylab(xylabs[2])
  if(!missing(xlim1)){
    xlim(xlim1)
    ylim(ylim1)
  }
  
  if(!is.numeric(labels)){
    return(p)
  }
  # mid<-median(labels)
  # p<-p+scale_color_gradient2(midpoint=0, low="blue", mid="gray",high="red", space ="Lab")
  p<-p+scale_color_gradient2(midpoint=mean(labels),
                             low="blue", mid="gray",high="red", space ="Lab")
  
  return(p)
  
  # p+geom_point(colour="black",pch=21, size=0.3)
  # p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #           panel.background = element_rect(fill = "darkgray"),
  #           axis.line = element_line(colour = "black"))
  
  
}

center.matrix<-function(m,dim = 1,sd.flag = F){
  if(dim == 1){
    zscores<-sweep(m,1,rowMeans(m,na.rm = T),FUN = '-')
  }else{
    zscores<-sweep(m,2,colMeans(m,na.rm = T),FUN = '-')
  }
  if(sd.flag){
    zscores<-sweep(zscores,dim,apply(m,dim,function(x) (sd(x,na.rm = T))),FUN = '/')
  }
  return(zscores)
}

discretize.prvt<-function(v,n.cat,q1){
  q1<-quantile(v,seq(from = (1/n.cat),to = 1,by = (1/n.cat)),na.rm = T)
  u<-matrix(data = 1,nrow = length(v))
  for(i in 2:n.cat){
    u[(v>=q1[i-1])&(v<=q1[i])]<-i
  }
  u<-paste0("Q",u)
  return(u)
}

set.list<-function (r,b,sampleName){
  set.field<-function (v,b){
    d <- dim(v)
    d.b<-length(b)
    if(!is.null(d)){
      if(d[1]==d.b){v <- subset(v,subset = b)}
      if(d[2]==d.b){v <- v[,b]}
    }else{if(length(v)==d.b){v <- v[b]}}
    return(v)
  }
  rn<-lapply(r, set.field, b = b)
  if(!missing(sampleName)){rn$sampleName<-sampleName}
  return(rn)
}

average.mat.rows<-function(m,ids,f = colMeans){
  ids.u<-sort(unique(ids))
  m1<-get.mat(ids.u,colnames(m))
  for(x in ids.u){
    b<-is.element(ids,x)
    if(sum(b)==1){
      m1[x,]<-m[b,]
    }else{
      m1[x,]<-f(m[b,])
    }
  }
  return(m1)
}

get.abundant<-function(v,abn.c = 2,boolean.flag = F,top,decreasing = T){
  m<-as.matrix(table(v))
  m<-as.matrix(m[order(m,decreasing = decreasing),])
  if(!missing(top)){
    abn.c<-m[top]
  }
  m<-m[m>=abn.c,]
  abn.names<-names(m)
  if(boolean.flag){
    b<-is.element(v,abn.names)
    return(b)
  }
  return(abn.names)
}

labels.mat.2.logical.mat<-function(M,filter.flag = F){
  M<-as.matrix(M)
  B<-NULL
  for(i in 1:ncol(M)){
    Bi<-labels.2.mat(M[,i])
    li<-colnames(Bi)
    if(filter.flag){
      Bi<-as.matrix(Bi[,-1])
      colnames(Bi)<-li[-1]  
    }
    B<-cbind(B,Bi)
  }
  return(B)
}

labels.2.mat<-function(x){
  B<-apply(as.matrix(unique(x)),1,function(xi) is.element(x,xi))  
  colnames(B)<-unique(x)
  return(B)
}

get.hyper.p.value<-function(b1,b2,full.flag = T){
  p1<-NA;p2<-NA;e<-0;
  if(any(b1)&&any(b2)){
    p1<-max(1-phyper(sum(b1&b2)-1, sum(b1), sum(!b1), sum(b2)),1e-17)
    e<-sum(b2)*(sum(b1)/length(b2))
    p2<-sum(b1&b2)/e
  }
  if (full.flag){
    p<-c(p1,p2,sum(b1&b2),e)
    names(p)<-c('hyper.p.value','ob.vs.exp','ob','exp')
    return(p)  
  }
  return(p1)
}

get.hyper.p.value.mat<-function(m1,m2,vcut = 1e-3,full.flag = FALSE){
  if(is.null(colnames(m1))){colnames(m1)<-1:ncol(m1)}
  if(is.null(colnames(m2))){colnames(m2)<-1:ncol(m2)}
  P<-get.mat(colnames(m1),colnames(m2))
  ob.vs.ex<-P;ob<-P;ex<-P;frc<-P
  for(i in 1:ncol(m1)){
    for(j in 1:ncol(m2)){
      p <- get.hyper.p.value(m1[,i],m2[,j],vcut)
      P[i,j] <- p[1]
      frc[i,j] <- sum(m1[,i]&m2[,j])/sum(m1[,i])
      ob.vs.ex[i,j] <- p[2]
      ob[i,j] <- p[3]
      ex[i,j] <- p[4]
    }
  }
  if(!full.flag){return(P)}
  results<-list(p = P, frc = frc,ob.vs.ex = ob.vs.ex,ob=ob, ex=ex)
  results$full<-cbind(matrix(P),matrix(ob.vs.ex),
                      matrix(ob),matrix(ex))
  return(results)  
  
}

get.mat<-function(m.rows,m.cols,data = NA){
  m<-matrix(data = data, nrow = length(m.rows),ncol = length(m.cols),
            dimnames = list(m.rows,m.cols))
  return(m)
}

scRNA_get.tpm<-function(r){
  r$cd<-as.matrix(r$cd)
  r$tpm<-1000000*sweep(r$cd,2,r$comp.reads,FUN = '/')
  r$tpm<-log2((r$tpm/10)+1)
  colnames(r$tpm)<-r$cells
  return(r)
}

intersect.list1<-function(l,g,n1=0){
  l1<-lapply(l, function(x) intersect(x,g))
  l1<-l1[laply(l1,length)>n1]
  return(l1)
}

ranksum.test.mat<-function(m,b,zscores.flag = T,two.sided=F){
  if(two.sided){
    p<-as.matrix(apply(m,1,function(x) ranksum.test(x[b],x[!b])))
  }else{
    p<-t(apply(m,1,function(x) c(ranksum.test(x[b],x[!b],alternative = 'greater'),
                                 ranksum.test(x[b],x[!b],alternative = 'less'))))
    colnames(p)<-c('more','less')
    if(zscores.flag){
      p<-cbind(p,get.p.zscores(p))
      colnames(p)[3]<-"zscores"
    }
  }
  return(p)
}

ranksum.test<-function(v1,v2,alternative="two.sided"){
  p <- NA
  if (sum(!is.na(v1))==0|sum(!is.na(v2))==0){return(p)}else{
    p<-wilcox.test(v1,v2,alternative = alternative)$p.value  
  }
  return(p)
}

get.p.zscores<-function(p){
  b<-p[,1]>0.5
  b[is.na(b)]<-F
  zscores<-(-log10(p[,1]))
  zscores[b]<-log10(p[b,2])
  # signficiant in p[,1] will be positive
  # signficiant in p[,2] will be negative
  return(zscores)
}

select.new.cells<-function(cells,samples,unused.cells,n2){
  b.use<-is.element(cells,unused.cells)
  u.samples<-unique(samples)
  b.final<-rep(F,length(cells))
  for(x in u.samples){
    b1<-is.element(samples,x)
    b.use1<-b1&b.use
    idx1<-which(b.use1)
    if(sum(b.use1)>n2){
      idx<-sample.wrapper(which(b.use1),n2)
    }else{
      idx<-unique(c(which(b.use1),sample.wrapper(which(b1&!b.use),n2-sum(b.use1))))
    }
    b.final[idx]<-T
  }
  return(b.final)
}

sample.wrapper<-function(x,n1){
  if(length(x)<n1){return(x)}
  return(sample(x = x,size = n1))
}

sample.per.label<-function(labels,size,boolean.flag = T,v,remove.flag = F,lb = 0,seed = 1234){
  if(missing(v)){v<-1:length(labels)}
  set.seed(seed)
  ul<-unique(labels)
  vr<-lapply(ul, function(x){
    b<-is.element(labels,x)
    if(size<1){
      vr<-sample(v[b],size = min(max(sum(b)*size,lb),sum(b)))
    }else{
      vr<-sample(v[b],size = min(size,sum(b)))
    }
    return(vr)
  })
  vr<-sort(unlist(vr))
  b<-is.element(v,vr)
  if(remove.flag){
    b.small<-!is.element(labels,get.abundant(labels[b],size-1))
    b[b.small]<-F
    vr<-v[b]
  }
  if(boolean.flag){return(b)}
  return(vr)
}

t.test.mat<-function(m,b,two.sided=F,rankf = F,fold.changeF = F){
  if(length(b)!=ncol(m)){
    print("Error. Inconsistent no. of samples.")
    return()
  }
  if(sum(b)<2||sum(!b)<2){
    return(get.mat(rownames(m),c('more','less',"zscores")))
  }
  if(two.sided){
    p<-as.matrix(apply(m,1,function(x) t.test(x[b],x[!b])$p.value))
  }else{
    p<-t(apply(m,1,function(x) c(t.test(x[b],x[!b],alternative = 'greater')$p.value,
                                 t.test(x[b],x[!b],alternative = 'less')$p.value)))
    colnames(p)<-c('more','less')
    p<-cbind(p,get.p.zscores(p))
    colnames(p)[3]<-"zscores"
  }
  if(rankf){
    p<-cbind(p,rank(p[,1]),rank(p[,2]))
    colnames(p)[4:5]<-c("rank.more","rank.less")
  }
  if(fold.changeF){
    p<-cbind.data.frame(p,pos.mean = rowMeans(m[,b]),neg.mean = rowMeans(m[,!b]))
    p$logFC<-log2(p$pos.mean/p$neg.mean)
  }
  
  return(p)
}

get.top.elements<-function (m,q = 100,min.ci = NULL,main = "",sort.flag = T){
  top.l<-list()
  v<-rownames(m)
  for (i in 1:ncol(m)){
    mi<-m[,i];mi<-mi[!is.na(mi)]
    idx<-order(mi,decreasing = F)
    ci <- mi[idx[min(q,length(mi))]]
    ci <- min(ci,min.ci)
    b <- m[,i]<=ci
    b[is.na(m[,i])]<-F
    if(sort.flag){
      top.l[[i]]<-sort(v[b])
    }else{
      top.l[[i]]<-v[b][order(m[b,i])]
    }
    
  }
  if(main!=""){main<-paste0(main,".")}
  names(top.l)<-paste0(main,colnames(m))
  return(top.l)
}

call_match<-function(v1,v2){
  rmv<-c('_',"-",'.'," ",":","+")
  v1<-casefold(call_gsub(pattern = rmv,replacement = '_',x = v1))
  v2<-casefold(call_gsub(pattern = rmv,replacement = '_',x = v2))
  idx<-match(v1,v2)
  return(idx)
}

call_gsub<-function(pattern,replacement = '',x){
  for(i in 1:length(pattern)){
    x<-gsub(pattern = pattern[i],replacement = replacement ,x = x,fixed = T)
  }
  return(x)
}

call_plot<-function(x, y = NULL,labels,regression.flag = F,icol = NULL,set.flag = F,cor.flag = F,
                    pch=16,cex=0.5,main="",ylab = "UMAP2",xlab = "UMAP1", cex.axis = 0.6,add.N = F,cex.main = 1){
  main<-capitalize(main)
  if(add.N&length(unique(labels))<30){
    labels<-add.n.of.samples(labels)
  }
  if(set.flag){
    par(mar=c(8, 7, 4.1, 12.1), xpd=TRUE)
  }
  if(is.null(icol)){
    icol<-labels.2.colors(labels)
  }
  if(is.null(y)){
    if(missing(xlab)){xlab<-colnames(x)[1]}
    if(missing(ylab)){ylab<-colnames(x)[2]}
    y<-x[,2];x<-x[,1]
  }
  
  if(cor.flag){
    xy.cor<-spearman.cor(y,x)
    main <- paste(main, "\nR =",format(xy.cor[1],digits = 2),"P =",format(xy.cor[2],scientific = T,digits = 2))
  }
  plot(x,y,col=icol,pch=pch,cex=cex,main=main,ylab=ylab,xlab = xlab,cex.axis = cex.axis,cex.main = cex.main)  
  
  labels<-gsub(" ","_",labels)
  l<-(max(x,na.rm = T)-min(x,na.rm = T))/20
  if(length(unique(labels))<30){
    if(length(pch)==length(labels)){
      map<-unique(paste(labels,icol,pch))
      labels.n<-as.matrix(table(labels))
      idx<-match(get.strsplit(map,' ',1),names(labels.n))
      map[,1]<-paste0(map[,1]," (N = ",m[idx],")")
      legend(x = max(x,na.rm = T)+l,
             y = max(y,na.rm = T),
             legend = get.strsplit(map,' ',1), 
             col = get.strsplit(map,' ',2),
             inset=c(-0.5,0),
             bty = "n",lty= NA, lwd = 0,cex = 0.7,pch = pch)
    }else{
      map<-unique(paste(labels,icol,pch))
      legend(x = max(x,na.rm = T)+l,
             y = max(y,na.rm = T),inset = c(-0.5,0),
             legend = gsub("_"," ",get.strsplit(map,' ',1)), 
             col = get.strsplit(map,' ',2),xpd = T,
             bty = "n",lty= NA, lwd = 0,cex = 0.7,pch = pch)
    }
    
  }
  if(regression.flag ==1){
    b<-!is.na(x)&!is.na(y)
    v<-lowess(x[b],y[b])
    lines(v)
    return(v)
  }
  if(regression.flag ==2){
    b<-!is.na(x)&!is.na(y)
    ulabels<-unique(labels)
    for(i in ulabels){
      bi<-b&labels==i
      v<-lowess(x[bi],y[bi])
      lines(v)
    }
    
  }
  
}

labels.2.colors<-function(x.class,x = NULL,number.flag = F,color.spec = "hsv"){
  palette("default")
  call_col<-c("black","red",setdiff(c(palette(),"cadetblue","gray","darkgreen","darkorange","darkviolet","gold3",
                                      "lightpink","deeppink2","deepskyblue",rainbow(20)),c("black","red")))
  call_col<-c("black","red","cadetblue","gray","darkgreen","darkorange","darkviolet","gold3",
              "lightpink","deeppink2","deepskyblue",palette(),rainbow(20))
  
  no.classes<-length(unique(x.class))
  # no.classes>length(call_col)
  if(number.flag){
    call_col<-match(x.class,sort(unique(x.class)))
  }else{
    if(is.numeric(x.class[1])){
      call_col<-plotrix::color.scale(x.class,c(0,10),0.8,0.8,color.spec = color.spec,
                                     xrange = c(min(x.class,na.rm = T)-1,max(x.class,na.rm = T)+1))
      if(color.spec=="hsv"){
        call_col<-plotrix::color.scale(x.class,c(0,10),0.8,0.8,color.spec = "hsv")
      }else{
        call_col<-plotrix::color.scale(x.class,c(0,10),extremes = c("green","red"),color.spec = "rgb")
      }
      # xrange = c(min(x.class,na.rm = T)-1,max(x.class,na.rm = T)+1))
      # icol<-plotrix::color.scale(x.class,c(0,10),0.8,0.8,extremes = c("green","red"),
      #                              color.spec = "hsv")
    }else{
      call_col<-call_col[match(x.class,sort(unique(x.class)))]
    }
  }
  if(!is.null(x)){
    names(call_col)<-x  
  }
  return(call_col)
}

get.strsplit<-function(v,sep,idx){
  v<-as.character(v)
  vi<-laply(strsplit(v,split = sep,fixed = T),function(x) x[idx])
  return(vi)
}

merge.datasets<-function(r1,r2,add.sample.id = F, features.list = NULL, sampleName = paste0(r1$sampleName,'+',r2$sampleName)){
  cd.flag<-!is.null(r1$cd)&!is.null(r2$cd)
  r<-list()
  genes<-intersect(rownames(r1$tpm),rownames(r2$tpm))  
  
  idx1<-match(genes,r1$genes);idx2<-match(genes,r2$genes)
  if(!is.null(r1$tpm)){
    r$tpm<-cbind(r1$tpm[idx1,],r2$tpm[idx2,])
  }
  if(cd.flag){
    r$cd<-cbind(r1$cd[idx1,],r2$cd[idx2,])
  }
  if(!is.null(r1$samples)){
    r$samples<-c(r1$samples,r2$samples)
  }
  
  r$genes<-genes
  if(add.sample.id){
    r$cells<-c(paste(r1$sampleName,r1$cells,sep = '_'),
               paste(r2$sampleName,r2$cells,sep = '_'))  
    colnames(r$cd)<-r$cells
    colnames(r$tpm)<-r$cells
  }else{
    if(!is.null(dim(r1$cells))){
      r$cells<-rbind(r1$cells,r2$cells)
    }else{
      r$cells<-c(r1$cells,r2$cells)
    }
  }
  
  r$comp<-c(r1$comp,r2$comp)
  r$comp.reads<-c(r1$comp.reads,r2$comp.reads)
  r$sampleName<-sampleName
  if(is.null(r1$labels)){r1$labels<-rep(r1$sampleName,length(r1$cells))}
  if(is.null(r2$labels)){r2$labels<-rep(r2$sampleName,length(r2$cells))}
  if(!is.null(features.list)){
    extra.fields<-lapply(features.list, function(x) rbind(as.matrix(r1[[x]]),as.matrix(r2[[x]])))
    names(extra.fields)<-features.list
    r<-c(r,extra.fields)
  }
  r$labels<-c(r1$labels,r2$labels)
  r$work.dir<-r1$work.dir
  return(r)
}

cap.mat<-function(M,cap = 0.01,MARGIN = 1){
  Z<-apply(M,MARGIN = MARGIN,function(x){
    q9<-quantile(x,1-cap)
    q1<-quantile(x,cap)
    x[x>q9]<-q9;x[x<q1]<-q1
    return(x)
  })
  if(MARGIN==1){Z<-t(Z)}
  return(Z)
}

call_plot.multilabels<-function(X,labels,main = NULL,xlab = "UMAP1",ylab="UMAP2",add.N = F,
                                pch = 16, cex = 0.3, cex.axis = 0.6,set.flag = F){
  laply(1:ncol(labels),function(i){
    call_plot(X,labels = labels[,i],
              main = ifelse(is.null(main),colnames(labels)[i],
                            paste(main,colnames(labels)[i],sep = ":")),
              xlab = xlab,ylab = ylab,pch = pch,cex.axis = cex.axis,cex = cex,
              set.flag = set.flag,add.N = add.N)
    return(i)})
}

call_view<-function(m,v,return.flag = F,order.flag = T,alternative = "two.sided",main = "m"){
  if(length(v)==1){
    b<-grepl(v,rownames(m),fixed = T)
  }else{
    b<-is.element(rownames(m),v)
  }
  if(is.numeric(m[1,1])){
    print(ranksum.test.labels(m[,1],b,alternative = alternative))  
    m<-m[b,]
    m<-m[order(m[,1],decreasing = alternative=="greater"),]
  }else{
    m<-m[b,]
  }
  View(m,main)
  if(return.flag){
    return(m)  
  }
}

ranksum.test.labels<-function(v,b,alternative="two.sided"){
  p<-NA
  v<-v[!is.na(b)]
  b<-b[!is.na(b)]
  v1<-v[b];v2<-v[!b]
  if (sum(!is.na(v1))==0|sum(!is.na(v2))==0){
    return(p)
  }else{
    p<-wilcox.test(v1,v2,alternative = alternative)$p.value  
  }
  return(p)
}

call_boxplot<-function (y,x,unique.x, f = median,
                        ylab = '',xlab = '',main = '',labels='black',
                        legend.name = "",add.anova = F,blank.flag = T,
                        cex = 0.7,order.flag = T,b.ref = NULL,p.val.show = is.null(b.ref)){
  b<-is.infinite(y)|is.na(y)
  if(length(labels)==length(x)){
    labels<-labels[!b] 
  }
  y<-y[!b];x<-x[!b]
  if(p.val.show){
    a<-aov(y ~ as.factor(x))
    kt<-kruskal.test(y~as.factor(x))
    p.kt = format(kt$p.value,digits=2)
    p.anova=format(unlist(summary(a))['Pr(>F)1'],digits=2)
    if(add.anova){
      main = paste(main,'\n(ANOVA p-value = ',p.anova,'\nKruskal p-value =',p.kt,')',sep = '')
    }
  }
  if(!is.null(b.ref)){
    main <- paste0(main,"\n(",call_format.pval(t.test.labels(y,b.ref,alternative = "greater")),", ",
                   "AUC = ",round(get.auc(y,b.ref),2),")")
  }
  
  if(missing(unique.x)){
    unique.x<-unique(x)
    x.med<-apply(as.matrix(unique.x),1,function(xi) f(y[is.element(x,xi)]))
    names(x.med)<-unique.x
    if(order.flag){
      unique.x<-unique.x[order(x.med)]  
    }
    #labels<-labels[unique.x]
  }
  
  #idx<-order(y.med[y])
  x <- data.frame(name = x, val = y)
  x$name <- factor(x$name, levels = unique.x)
  p <- ggplot(x, aes(x = name, y = val,color = labels)) + #geom_point(stat="identity") +
    labs(y = ylab, title = main, x = xlab) + 
    geom_boxplot() + theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size = rel(cex)),
          axis.text.y = element_text(size = rel(cex)),
          axis.title.y = element_text(size = rel(cex)),
          plot.title = element_text(size = rel(1)))
  p <- p + labs(color = legend.name)
  if(blank.flag){
    p <- p+ theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())
  }
  return(p)
  #labs(title = paste(screen.name,'ANOVA p-value = ',p.anova)) + 
}

call.dotPlot<-function(X,cex = 12,midpoint = 0,color.name = 'Effect size'){
  ggplot(data = X,aes(x=cell.type, y = Gene, color = Estimate, size = Z)) + 
    geom_point() + 
    cowplot::theme_cowplot() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = cex),
          axis.text.y = element_text(size=cex)) +
    scale_colour_gradient2(low = "blue",high = "red",midpoint = midpoint,
                           oob = scales::squish, name = color.name)+
    geom_point(shape = 1,colour = "black")
}

convert.to.vector<-function(str){
  return(paste0("'",gsub(", ","','",str),"'"))
}

plot.layers<-function(x, y = NULL,labels,b.top,red.top = F,regression.flag = F,icol = NULL,set.flag = F,cor.flag = F,
                      pch=16,cex=0.3,main="",ylab = "tSNE2",xlab = "tSNE1", cex.axis = 0.6,cex.main = 1,
                      add.N = F,grey.zeros = F,legend.flag = T,color.spec = "rgb"){
  
  regl<-call_plot(x = x,y = y,labels,regression.flag,icol = icol,
                  set.flag = set.flag,cor.flag = cor.flag,
                  pch = pch,cex = cex,main = main,ylab = ylab,xlab = xlab,
                  cex.axis = cex.axis,cex.main = cex.main,
                  add.N = add.N,legend.flag = legend.flag,color.spec = color.spec)
  if(is.null(y)){
    v<-colnames(x)
    if(xlab==""){xlab<-v[1]}
    if(ylab==""){ylab<-v[2]}
    y<-x[,2];x<-x[,1]
  }
  if(red.top){
    points(x[b.top],y[b.top],cex = cex,col = "red",pch = 1)
  }else{
    if(is.null(icol)){
      if(grey.zeros){
        icol<-rep("grey",length(labels))
        icol[labels>0]<-labels.2.colors(labels[labels>0])
      }else{
        icol <- labels.2.colors(labels)
      }
    }
    points(x[b.top],y[b.top],cex = cex,col = icol[b.top],pch = 16)
  }
  return(regl)
  
}

grepl.list<-function(sig,x){
  sig<-lapply(sig, function(v) v[grepl(x,v)])
  sig<-sig[laply(sig,length)>0]
  return(sig)
}

lapply.unique.sort<-function(l){
  l<-lapply(l, function(x) unique(sort(x)))
  return(l)
}

discretize.3.labels<-function(X,q = 0.1,verbose = F){
  if(q>0.5){q<-(1-q)}
  if(verbose){
    print(paste("Low quantile<=",q))
    print(paste("High quantile>=",1-q))
  }
  
  f<-function(v){
    b.low<-v<=quantile(v,q,na.rm = T)
    b.high<-v>=quantile(v,1-q,na.rm = T)
    labels<-ifelse(b.high,"High",ifelse(b.low,"Low","Moderate"))
    if(!any(is.na(labels))){
      labels<-factor(labels,levels = c("High","Moderate","Low"))
    }
    return(labels)
  }
  if(!is.matrix(X)){return(f(X))}
  B<-apply(X,2,f)
  return(B)
}

intersect.list1<-function(l,g,n1=0,HG.universe = NULL,prf = ""){
  l1<-lapply(l, function(x) intersect(x,g))
  l1<-l1[laply(l1,length)>n1]
  if(prf!=""){
    names(l1)<-paste(prf,names(l1),sep = ".")
  }
  if(!is.null(HG.universe)){
    p<-GO.enrichment.lapply(l[names(l1)],genes = HG.universe,list(g))
    names(l1)<-paste(names(l1),format(p,scientific = T,digits= 3),sep = "P = ")
  }
  return(l1)
}

setdiff.list1<-function(l,g,n1=0,prf = ""){
  l1<-lapply(l, function(x) setdiff(x,g))
  l1<-l1[laply(l1,length)>n1]
  if(prf!=""){
    names(l1)<-paste(prf,names(l1),sep = ".")
  }
  return(l1)
}

add.n.of.samples<-function(l,n.flag = T,sep = " "){
  num.samples<-table(l)
  # print(num.samples)
  idx<-match(l,names(num.samples))
  if(n.flag){
    l<-paste0(l,sep,"(n = ",num.samples[idx],")")
  }else{
    l<-paste0(l,sep,"(",num.samples[idx],")")
  }
  
  # print(table(l))
  return(l)
}

discretize.mat.q<-function(X,q1 = 0.9,strict = F){
  qv<-t(apply(X,2,function(x) quantile(x,q1,na.rm = T)))
  if(strict){
    B<-sweep(X,2,qv,FUN = "-")>0
  }else{
    B<-sweep(X,2,qv,FUN = "-")>=0
  }
  rownames(B)<-rownames(X)
  return(B)
}

gg.densityplot<-function(x,labels,title="",subtitle="",xlab = "Scores",legend.name="",caption="",alpha = 0.4){
  theme_set(theme_classic())
  p<-t.test.labels(x,labels==labels[1])
  subtitle<-paste0(subtitle," (",call_format.pval(p),")")
  mpg <- cbind.data.frame(cty = x,cyl = labels)
  # Plot
  g <- ggplot(mpg, aes(cty))
  # g <- ggplot(mpg, aes(cty,color = cyl))
  
  g <- g + geom_density(aes(fill=factor(cyl)), alpha = alpha) + 
    labs(title=title, 
         x=xlab,
         subtitle=subtitle,
         fill=legend.name)
  # g<-g+ xlim(min(x),max(x))
  
  g<-g+ xlim(min(x)-(sd(x)*3),max(x)+(sd(x)*3))
  
  # g <- g + theme_classic()
  # 
  # caption=caption,
  # labs(title="Density plot", 
  #      subtitle="City Mileage Grouped by Number of cylinders",
  #      caption="Source: mpg",
  #      x="City Mileage",
  #      fill="# Cylinders")
  return(g)
}

call_multiplot<-function(plotlist,nplots = 4,cols = 2){
  flag<-F
  while(!is.null(plotlist)&!flag){
    idx<-1:min(nplots,length(plotlist))
    print(multiplot(plotlist = plotlist[idx],cols = cols))
    flag<-(min(nplots,length(plotlist))+1)>length(plotlist)
    if(max(idx)>=length(plotlist)){return(plotlist)}
    plotlist<-plotlist[(max(idx)+1):length(plotlist)]
  }
  return(plotlist)
}

intersect.lists<-function(l1,l2){
  idx<-call_match(names(l1),names(l2))
  names(l1)<-names(l2)[idx]
  L<-lapply(names(l1), function(x) intersect(l1[[x]],l2[[x]]))
  names(L)<-names(l1)
  return(L)
}

intersect.lists.by.idx<-function(l1,l2,remove.empty = F){
  L<-lapply(1:length(l1), function(x) intersect(l1[[x]],l2[[x]]))
  names(L)<-names(l1)
  if(!remove.empty){return(L)}
  L<-L[laply(L,length)>0]
  return(L)
}

union.lists<-function(l1,l2,unique.flag = T,disregard.names = F){
  if(disregard.names){
    names(l2)<-names(l1)
  }else{
    idx<-call_match(names(l1),names(l2))
  }
  # names(l1)<-names(l2)[idx]
  if(unique.flag){
    L<-lapply(names(l1), function(x) unique(sort(c(l1[[x]],l2[[x]]))))
  }else{
    L<-lapply(names(l1), function(x) c(l1[[x]],l2[[x]]))
  }
  
  names(L)<-names(l1)
  return(L)
}

setdiff.lists<-function(l1,l2){
  idx<-call_match(names(l1),names(l2))
  names(l1)<-names(l2)[idx]
  L<-lapply(names(l1), function(x) setdiff(l1[[x]],l2[[x]]))
  names(L)<-names(l1)
  # print(summary(L))
  return(L)
}

setdiff.lists.by.idx<-function(l1,l2){
  L<-lapply(1:length(l1), function(x) setdiff(l1[[x]],l2[[x]]))
  names(L)<-names(l1)
  print(summary(L))
  return(L)
}

p.adjust.mat<-function(m,method = "BH"){
  P<-apply(m,2,function(x) p.adjust(x,method = method))
  return(P)
}

call_format.pval<-function(p,prnt.flag = F,d = "="){
  if(length(p)>1){
    P<-laply(p,call_format.pval)
    P<-gsub("p = ","",P)
    # P<-paste0("P",1:length(P)," = ",P)
    P<-paste("p =",paste(P,collapse = ", "))
    return(P)
  }
  if(abs(p)>1){
    p<-10^(-abs(p))
  }
  if(p>0.05){p<-paste("p",d,round(p,3));return(p)}
  if(p==0){p<-"p < 1 x 10-20";return(p)}
  p<-gsub("e"," x 10",paste("p",d,format(p,scientific = T,digits= 3)))
  p<-gsub("-0","-",p)
  
  if(prnt.flag){
    p<-paste0("(",p,")")
  }
  return(p)
}

convert.to.vector.of.names<-function(v,R.flag = T){
  if(R.flag){
    return(paste0("c('",paste(v,collapse = "','"),"')"))
  }else{
    return(paste(v,collapse = ", "))
  }
  
}

cap.mat<-function(M,cap = 0.01,MARGIN = 1){
  Z<-apply(M,MARGIN = MARGIN,function(x){
    q9<-quantile(x,1-cap)
    q1<-quantile(x,cap)
    x[x>q9]<-q9;x[x<q1]<-q1
    return(x)
  })
  if(MARGIN==1){Z<-t(Z)}
  return(Z)
}

get.auc<-function(p1,y1){
  pr <- prediction(p1, y1)
  auc <- performance(pr, measure = "auc")
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  # prf <- performance(pr, measure = "prec", x.measure = "rec")
  # plot(prf,ylim = c(0,1))
  # abline(h = mean(y1))
  return(auc)
}

comparative.pval.analysis<-function(p,p.cut,genes,go.env,rank.cut = 2){
  v<-t(apply(p,1,rank))
  v[p>p.cut]<-100
  rs<-list()
  rs$m<-v
  rs$sig<-get.top.elements(v,100,rank.cut)
  rs$go<-GO.enrichment.lapply(go.env,genes,rs$sig)
  return(rs)
}

get.top.elements.comb<-function(m,q = 100,min.ci = NULL){
  b<-rowSums(is.na(m))==0
  m<-m[b,]
  mc<-cbind(-rowMin(m),rowMax(m))
  colnames(mc)<-c("more","less")
  return(get.top.elements(mc,q = q,min.ci = min.ci))
}

list.2.boolean.mat<-function(l,ids = NULL){
  if(is.null(ids)){
    ids<-sort(unique(unlist(l)))
  }
  B<-t(laply(l,function(x) is.element(ids,x)))
  colnames(B)<-names(l)
  rownames(B)<-ids
  return(B)
}

call_plot.multilabels<-function(X,labels,main = NULL,xlab = "UMAP1",ylab="UMAP2",add.N = F,
                                pch = 16, cex = 0.3, cex.axis = 0.6,set.flag = F){
  laply(1:ncol(labels),function(i){
    call_plot(X,labels = labels[,i],
              main = ifelse(is.null(main),colnames(labels)[i],
                            paste(main,colnames(labels)[i],sep = ":")),
              xlab = xlab,ylab = ylab,pch = pch,cex.axis = cex.axis,cex = cex,
              set.flag = set.flag,add.N = add.N)
    return(i)})
}

get.top.cor<-function(m,q = 100,min.ci = 0,idx = NULL, add.prefix ="",sort.flag = T){
  m<-as.matrix(m)
  if(is.null(colnames(m))){colnames(m)<-1:ncol(m)}
  m.pos<-(-m);m.neg<-m
  
  colnames(m.pos)<-paste0(colnames(m.pos),".up")
  colnames(m.neg)<-paste0(colnames(m.neg),".down")
  v<-get.top.elements(cbind(m.pos,m.neg),q,min.ci = (-abs(min.ci)),sort.flag = sort.flag)
  names(v)<-c(colnames(m.pos),colnames(m.neg))
  if(!is.null(idx)){
    v<-v[paste(idx,c("up","down"),sep = ".")]
  }
  names(v)<-paste0(add.prefix,names(v))
  # v<-v[order(names(v),decreasing = T)]
  return(v)
}

list.2.mat<-function(l){
  if(length(l)<2){
    m<-as.matrix(l[[1]])
    colnames(m)<-names(l)
    return(m)
  }
  n1<-max(laply(l,length))
  m<-t(laply(l,function(x) c(x,matrix(data = "",nrow = n1-length(x)+1))))
  m<-m[1:n1,]
  if(n1==1){
    names(m)<-names(l)
  }else{
    colnames(m)<-names(l) 
  }
  return(m)
}

get.mat<-function(m.rows,m.cols,data = NA){
  m<-matrix(data = data, nrow = length(m.rows),ncol = length(m.cols),
            dimnames = list(m.rows,m.cols))
  return(m)
}

get.residuals<-function(X,g,MARGIN = 1){
  if(MARGIN == 2){return(t(get.residuals(t(X),g,MARGIN = 1)))}
  g<-as.matrix(g)
  b<-rowSums(is.na(g))==0
  g<-as.matrix(g[b,])
  residuals<-matrix(nrow = nrow(X),ncol = ncol(X))
  rownames(residuals)<-rownames(X)
  colnames(residuals)<-colnames(X)
  f<-function(y){
    b1<-!is.na(y)
    y<-y[b1];g<-g[b1,]
    v<-matrix(nrow = length(b1))
    v[b1]<-lm(y~.,data = as.data.frame(g))$residuals
    return(v)
  }
  residuals[,b]<-t(apply(X[,b],1,f))
  return(residuals)
}

logic.mat.2.labels<-function(B){
  v<-colnames(B)
  labels<-apply(B,1,function(x){
    if(sum(x,na.rm = T)>1){return("Umbi.")}
    if(sum(sum(x,na.rm = T)==1)){return(v[x])}
    return("?")})
  return(labels)
}

union.multiple.mats<-function(mat.list){
  g<-sort(as.character(unique(unlist(lapply(mat.list, rownames)))))
  for(i in 1:length(mat.list)){
    m1<-mat.list[[i]]
    colnames(m1)<-paste(names(mat.list)[i],colnames(m1),sep = ".")
    idx<-match(g,rownames(m1))
    if(i==1){m<-m1[idx,]}else{m<-cbind(m,m1[idx,])}
  }
  rownames(m)<-g
  return(m)
}

ranksum.test.binary.mat<-function(v,B){
  r<-t(apply(B,1,function(b) {
    c(ranksum.test.labels(v,b,alternative = "greater"),
      ranksum.test.labels(v,b,alternative = "less"))}))
  return(r)
}

get.pval.from.zscores<-function(z){
  p<-10^(-abs(z))
  b<-z<0&!is.na(z)
  p[b]<-1-p[b]
  return(p)
}

sample.per.label <-function(labels,size,boolean.flag = T,v,remove.flag = F,lb = 0,seed = 1234){
  set.seed(seed)
  if(missing(v)){
    v<-1:length(labels)
  }
  ul<-unique(labels)
  vr<-NULL
  for(i in 1:length(ul)){
    b<-is.element(labels,ul[i])
    if(size<1){
      vr<-c(vr,sample(v[b],size = min(max(sum(b)*size,lb),sum(b))))
    }else{
      vr<-c(vr,sample(v[b],size = min(size,sum(b))))
    }
  }
  b<-is.element(v,vr)
  if(remove.flag){
    b.small<-!is.element(labels,get.abundant(labels[b],size-1))
    b[b.small]<-F
    vr<-v[b]
  }
  if(boolean.flag){return(b)}
  return(vr)
}

boxplot.test<-function(vals,labels,las = 1,main = "",cex = 1.1,ylab = NULL,alternative = "two.sided",boxwex = 0.5,
                       t.test.flag = T,dots.flag = F,legend.flag = T,ref.label = labels[1],xlab = ""){
  no.labels<-length(unique(labels))
  if(no.labels==2){col = c("cadetblue","gray")}
  if(no.labels==3){col = c("cadetblue","gray70","brown4")}
  if(no.labels>3){col = c("cadetblue","lightblue","gray70","brown4","orange")}
  if(t.test.flag==T){
    p<-t.test.labels(vals,is.element(labels,ref.label),alternative = alternative)
    if(alternative=="greater"){
      auc <- get.auc(vals,is.element(labels,ref.label))
    }else{
      auc <- get.auc(vals,!is.element(labels,ref.label))
    }
    auc<-round(auc,digits = 2)
  }else{
    p<-ranksum.test.labels(vals,is.element(labels,ref.label),alternative = alternative)
  }
  main<-set.titles(gsub("_"," ",main))
  if(is.null(ylab)){ylab <- main}
  if(t.test.flag!="none"){
    main <- gsub("_"," ",paste0(main,"\nt-test ",call_format.pval(p),#format(p,scientific = T,digits = 2),
                                "\nAUROC = ",auc))
  }
  # if(violin.flag){
  #   violin_plot_labels(v = vals,G = labels,main = main,col = col,ylab = ylab)
  # }else{
  boxplot(vals~labels,las=las,cex.axis = 1,col = col,ylab = ylab,main = main,cex.main = cex,
          cex.lab = cex,xlab = xlab,boxwex = boxwex)
  l<-as.matrix(table(labels))
  l<-paste0(rownames(l)," (n = ",l,")")
  if(legend.flag){legend("topleft",legend = l,fill = col,cex = 0.8)}
  if(dots.flag){
    stripchart(vals~labels, vertical = TRUE,method = "jitter", 
               add = TRUE, pch = 20, col = 'gray30')
  }
  return(p)
}

center.matrix<-function(m,dim = 1,sd.flag = F){
  if(dim == 1){
    zscores<-sweep(m,1,rowMeans(m,na.rm = T),FUN = '-')
  }else{
    zscores<-sweep(m,2,colMeans(m,na.rm = T),FUN = '-')
  }
  if(sd.flag){
    zscores<-sweep(zscores,dim,apply(m,dim,function(x) (sd(x,na.rm = T))),FUN = '/')
  }
  return(zscores)
}

plot.bimodal.distribution<-function(y,density.flag = F,xlab = "",main = "",
                                    pos.label = "pos",neg.label = "neg",seed = 1234){
  set.seed(seed)
  mixmdl = normalmixEM(y)
  main = paste(main,paste("loglikelihood =",round(mixmdl$loglik)),sep = "\n")
  plot(mixmdl,which=2,xlab2=xlab,main2 = main)
  lines(density(y), lty=2, lwd=2)
  mixmdl$rank<-rank(y)
  idx1<-order(mixmdl$mu)[1]
  mixmdl$labels<-mixmdl$rank>round(mixmdl$lambda[idx1]*length(y))
  if(density.flag){
    call_densityplot(y,mixmdl$labels)
  }
  if(length(mixmdl$all.loglik)==1001){
    mixmdl$labels[]<-NA  
  }
  mixmdl$labelsF<-ifelse(mixmdl$labels,pos.label,neg.label)
  y.pos<-y[mixmdl$labels]
  y.neg<-y[!mixmdl$labels]
  mixmdl$labelsF[(y<(mean(y.pos)-sd(y.pos)))&mixmdl$labels]<-paste0(pos.label,"?")
  mixmdl$labelsF[(y>(mean(y.neg)+sd(y.neg)))&!mixmdl$labels]<-paste0(neg.label,"?")
  return(mixmdl)
}

inter.vs.intra.corr<-function(m,b,B = NULL,main = '',lab1 = 'Within',lab2 = 'Between'){
  diag(m)<-NA
  m<-m[b,];m1<-m[,b];m2<-m[,!b]
  v<-c(matrix(m1),matrix(m2))
  l<-c(ifelse(matrix(m1)>(-100),1,0),ifelse(matrix(m2)>(-100),0,1))
  b<-!is.na(v)&!is.na(l)
  #multiple.hist(v[b],l[b])
  p<-ranksum.test.labels(v,l==1,alternative = 'greater')
  print(p)
  print(paste('Within:',round(mean(v[l==1],na.rm = T),2),round(median(v[l==1],na.rm = T),2)))
  print(paste('Between:',round(mean(v[l==0],na.rm = T),2),round(median(v[l==0],na.rm = T),2)))
  l<-c(ifelse(matrix(m1)>(-100),lab1,0),ifelse(matrix(m2)>(-100),lab2,1))
  boxplot(v~l,main = paste(main,'\n P =',format(p,scientific=T,digits = 2)),ylab = 'Spearman correlation')
  if(!is.null(B)){
    L<-list(within<-matrix(m1),between<-matrix(m2))
    for(i in 1:ncol(B)){
      L[[i+2]]<-m[,B[,i]]  
    }
    names(L)<-c('Within','Between',colnames(B))
    boxplot(L,las=2)
    
  }
  
}

divid.for.one.cv<-function(features,train.f = 0.7){
  features<-apply(features,1,function(x) paste(x,collapse = '_'))
  unique.features<-unique(features)
  b.train<-matrix(data = F,nrow = length(features))
  rownames(b.train)<-rownames(features)
  for(i in 1:length(unique.features)){
    idx<-which(is.element(features,unique.features[i]))
    n.train<-round(length(idx)*train.f)
    idx.train<-sample(idx,size = n.train)
    b.train[idx.train]<-T
  }
  return(b.train)
}

GO.enrichment<-function(go.env,genes,selected.genes,add.genes = F,sort.flag = F){
  b<-is.element(genes,selected.genes)
  p<-laply(go.env,function(x) get.hyper.p.value(is.element(genes,x),b))
  rownames(p)<-names(go.env)
  p<-as.data.frame(p)
  if(add.genes){
    p<-cbind.data.frame(p,
                        Genes = laply(go.env,function(x) paste(intersect(selected.genes,x),collapse = ", ")))
  }
  if(sort.flag){
    p<-p[order(p[,1]),]
  }
  p$p.adj<-p.adjust(p$hyper.p.value,method = "BH")
  return(p)
}

GO.enrichment.lapply<-function(go.env,genes,sig){
  f<-function(x1){
    b1<-is.element(genes,x1)
    s.b1<-sum(b1)
    s.n.b1<-sum(!b1)
    p<-laply(go.env,function(x2){
      b2<-is.element(genes,x2)
      s.b2<-sum(b2)
      if(s.b2==0){return(0)}
      return(phyper(sum(b1&b2)-1, s.b1, s.n.b1, s.b2))
    })
  }
  m<-1-t(laply(sig,f))
  m[m==0]<-1e-17
  rownames(m)<-names(go.env)
  colnames(m)<-names(sig)
  return(m)
}

GO.enrichment.lapply.v1<-function(go.env,genes,sig,valuType = 1){
  m<-t(laply(sig,function(x) GO.enrichment(go.env,genes,x)[,valuType]))
  colnames(m)<-names(sig)
  rownames(m)<-names(go.env)
  return(m)
}

GO.enrichment.lapplyF<-function(go.env,genes,sig){
  sig.go<-lapply(sig,function(x) GO.enrichment(go.env = go.env,genes = genes,selected.genes = x,add.genes = T))
  names(sig.go)
  sig.go$sum<-t(laply(sig.go,function(X) return(X[,1])))
  colnames(sig.go$sum)<-names(sig)
  rownames(sig.go$sum)<-rownames(sig.go[[1]])
  return(sig.go)
}

call_plot<-function(x, y = NULL,labels,regression.flag = F,icol = NULL,set.flag = F,cor.flag = F,legend.flag = T,
                    pch=16,cex=0.5,main="",ylab = "UMAP2",xlab = "UMAP1", cex.axis = 0.6,add.N = F,cex.main = 1,
                    color.spec = "rgb"){
  main<-capitalize(main)
  if(add.N&length(unique(labels))<30){
    labels<-add.n.of.samples(labels)
  }
  if(set.flag){
    par(mar=c(8, 7, 4.1, 12.1), xpd=TRUE)
  }
  if(is.null(icol)){
    icol<-labels.2.colors(labels,color.spec = color.spec)
  }
  if(is.null(y)){
    if(missing(xlab)){xlab<-colnames(x)[1]}
    if(missing(ylab)){ylab<-colnames(x)[2]}
    y<-x[,2];x<-x[,1]
  }
  
  if(cor.flag){
    xy.cor<-spearman.cor(y,x)
    main <- paste(main, "\nR =",format(xy.cor[1],digits = 2),"P =",format(xy.cor[2],scientific = T,digits = 2))
  }
  plot(x,y,col=icol,pch=pch,cex=cex,main=main,ylab=ylab,xlab = xlab,cex.axis = cex.axis,cex.main = cex.main)  
  
  labels<-gsub(" ","_",labels)
  l<-(max(x,na.rm = T)-min(x,na.rm = T))/20
  if(length(unique(labels))<30&legend.flag){
    if(length(pch)==length(labels)){
      map<-unique(paste(labels,icol,pch))
      labels.n<-as.matrix(table(labels))
      idx<-match(get.strsplit(map,' ',1),names(labels.n))
      map[,1]<-paste0(map[,1]," (N = ",m[idx],")")
      print(as.integer(get.strsplit(map,' ',3)))
      legend(x = max(x,na.rm = T)+l,
             y = max(y,na.rm = T),
             legend = get.strsplit(map,' ',1), 
             col = get.strsplit(map,' ',2),
             inset=c(-0.5,0),
             bty = "n",lty= NA, lwd = 0,cex = 0.7,pch = pch)
    }else{
      map<-unique(paste(labels,icol,pch))
      legend(x = max(x,na.rm = T)+l,
             y = max(y,na.rm = T),inset = c(-0.5,0),
             legend = gsub("_"," ",get.strsplit(map,' ',1)), 
             col = get.strsplit(map,' ',2),xpd = T,
             bty = "n",lty= NA, lwd = 0,cex = 0.7,pch = pch)
    }
    
  }
  if(regression.flag ==1){
    b<-!is.na(x)&!is.na(y)
    v<-lowess(x[b],y[b])
    lines(v)
    return(v)
  }
  if(regression.flag ==2){
    b<-!is.na(x)&!is.na(y)
    ulabels<-unique(labels)
    for(i in ulabels){
      bi<-b&labels==i
      v<-lowess(x[bi],y[bi])
      lines(v)
    }
    
  }
  
  
}

get.strsplit<-function(v,sep,idx){
  v<-as.character(v)
  vi<-laply(strsplit(v,split = sep,fixed = T),function(x) x[idx])
  return(vi)
}

get.abundant<-function(v,abn.c = 2,boolean.flag = F,top,decreasing = T){
  m<-as.matrix(table(v))
  m<-as.matrix(m[order(m,decreasing = decreasing),])
  if(!missing(top)){
    abn.c<-m[top]
  }
  m<-m[m>=abn.c,]
  abn.names<-names(m)
  if(boolean.flag){
    b<-is.element(v,abn.names)
    return(b)
  }
  return(abn.names)
}

colSD <- function(m){
  m<-apply(m,2,function(x) sd(x,na.rm = T))
  return(m)
}

rowMax<-function(X){
  y<-apply(X,1,function(x) max(x,na.rm = T))
  return(y)
}

colMax<- function(X){
  return(t(rowMax(t(X))))
}

labels.2.mat<-function(x){
  B<-apply(as.matrix(unique(x)),1,function(xi) is.element(x,xi))  
  colnames(B)<-unique(x)
  return(B)
}

call_barplot<-function(X,names.arg = rownames(X),legend.names=colnames(X),sort.idx = 0,
                       main ='',xlab='',ylab='',icol=NULL,beside=T,ymax = max(X),legend.loc = "topright",
                       cex.names = 0.8,las = 2,horiz = F,ylim = NULL){
  # par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  if(is.null(icol)){
    icol<-labels.2.colors(colnames(X))
  }
  if(is.null(names.arg)){
    names.arg<-1:nrow(X)
  }
  if(sort.idx>0){
    X<-X[order(X[,sort.idx]),]
    names.arg<-rownames(X)
  }
  v<-barplot(t(X),width
             =0.1, 
             main=main,horiz = horiz,
             names.arg=names.arg, ylab=ylab, 
             xlab=xlab, col=icol,ylim = ylim,
             cex.names=cex.names,beside = beside,las=las,border = 0.01)
  
  if(!is.null(legend.names)){
    legend(legend.loc,horiz = F,
           legend = legend.names,#inset=c(-0.8,0),
           col = icol,bty='n',lty= 0, # line style
           lwd = 5,pch = 15,xpd = T)  
  }
  
}

get.hyper.p.value<-function(b1,b2,full.flag = T){
  p1<-NA;p2<-NA;e<-0;
  if(any(b1)&&any(b2)){
    p1<-max(1-phyper(sum(b1&b2)-1, sum(b1), sum(!b1), sum(b2)),1e-17)
    e<-sum(b2)*(sum(b1)/length(b2))
    p2<-sum(b1&b2)/e
  }
  if (full.flag){
    p<-c(p1,p2,sum(b1&b2),e)
    names(p)<-c('hyper.p.value','ob.vs.exp','ob','exp')
    return(p)  
  }
  return(p1)
}

get.hyper.p.value.mat<-function(m1,m2,vcut = 1e-3,diag.flag = FALSE,full.flag = FALSE){
  if (diag.flag){
    P=matrix(nrow = ncol(m1),ncol = 1)
    ob.vs.ex=P
    rownames(P)<-colnames(m1)
    for(i in 1:ncol(m1)){
      p=get.hyper.p.value(m1[,i],m2[,i],vcut)
      P[i]=p[1]
      ob.vs.ex[i]=p[2]
    }
  }else{
    if(is.null(colnames(m1))){colnames(m1)<-1:ncol(m1)}
    if(is.null(colnames(m2))){colnames(m2)<-1:ncol(m2)}
    P<-get.mat(colnames(m1),colnames(m2))
    ob.vs.ex<-P;ob<-P;ex<-P
    for(i in 1:ncol(m1)){
      for(j in 1:ncol(m2)){
        p=get.hyper.p.value(m1[,i],m2[,j],vcut)
        P[i,j]=p[1]
        ob.vs.ex[i,j]=p[2]
        ob[i,j]=p[3]
        ex[i,j]=p[4]
      }
    }
  }
  if (full.flag){
    results<-list()
    results$p<-P
    results$ob.vs.ex=ob.vs.ex
    results$ob=ob
    results$ex=ex
    results$full<-cbind(matrix(P),matrix(ob.vs.ex),
                        matrix(ob),matrix(ex))
    return(results)  
  }else{
    return(P)
  }
  
}

semi.random <-function(v,labels){
  ul<-unique(labels)
  vr<-v;vr[]<-NA
  for(i in 1:length(ul)){
    b<-is.element(labels,ul[i])
    vr[b]<-sample(v[b])
  }
  return(vr)
}

pcor.mat<-function(v1,v2 = NULL,v3, method = 'spearman',
                   use = "pairwise.complete.obs",match.flag = F,
                   alternative = "two.sided",upper.tri.flag = F){
  if(is.null(v2)){
    v2<-v1
  }
  if(!is.matrix(v1)){v1<-as.matrix(v1)}
  if(!is.matrix(v2)){v2<-as.matrix(v2)}
  if(!is.matrix(v3)){v3<-as.matrix(v3)}
  if(match.flag){
    n=ncol(v1)
    results<-matrix(data = NA,nrow = n,ncol = 2)
    rownames(results)<-colnames(v1)
    for(i in 1:ncol(v1)){
      b<-!is.na(v1[,i])&!is.na(v2[,i])
      v1i<-v1[b,i];v2i<-v2[b,i]
      c.i <- tryCatch(pcor.test(v1i,v2i,v3[b,],method = method),
                      error = function(err){return(NA)})
      if(is.list(c.i)){
        results[i,1] <- c.i$estimate
        results[i,2] <- c.i$p.value  
      }
    }
    colnames(results)<-c("R","P")
  }else{
    n1=ncol(v1)
    m<-matrix(nrow = n1,ncol = ncol(v2))
    rownames(m)<-colnames(v1)
    colnames(m)<-colnames(v2)
    results<-list(cor = m, p = m)
    for(i in 1:n1){
      f<-function(x){
        b<-!is.na(v1[,i])&!is.na(x)
        if(ncol(v3)>1){
          c.i<-tryCatch(pcor.test(v1[b,i],x[b],v3[b,],method = method),
                        error = function(err){return(NA)})}
        else{
          c.i<-tryCatch(pcor.test(v1[b,i],x[b],v3[b],method = method),
                        error = function(err){return(NA)})
        }
        if(is.list(c.i)){
          return(c(c.i$estimate,c.i$p.value))
        }else{
          return(c(NA,NA))
        }
      }
      c.i <- apply(v2,2,f)
      results$cor[i,] <- c.i[1,]
      results$p[i,] <- c.i[2,]
    }
    if(ncol(v2)==1){
      results<-cbind(results$cor,results$p)
      colnames(results)<-c('R','P')
    }
  }
  if(upper.tri.flag){
    results$up <- cbind(results$cor[upper.tri(results$cor)],
                        results$p[upper.tri(results$p)])
  }
  return(results)
}

spearman.cor<-function(v1,v2 = NULL,method = 'spearman',
                       use = "pairwise.complete.obs",match.flag = F,
                       alternative = "two.sided",upper.tri.flag = F){
  if(is.null(v2)){
    v2<-v1
  }
  if(!is.matrix(v1)){v1<-as.matrix(v1)}
  if(!is.matrix(v2)){v2<-as.matrix(v2)}
  if(match.flag){
    n=ncol(v1)
    if(is.null(colnames(v1))){colnames(v1)<-1:ncol(v1)}
    results<-get.mat(m.cols = c("R","P"),m.rows = colnames(v1))
    for(i in 1:ncol(v1)){
      c.i <- cor.test(v1[,i],v2[,i],method = method,use = use, alternative = alternative)
      results[i,1] <- c.i$estimate
      results[i,2] <- c.i$p.value
    }
  }else{
    n1=ncol(v1)
    m<-matrix(nrow = n1,ncol = ncol(v2))
    rownames(m)<-colnames(v1)
    colnames(m)<-colnames(v2)
    results<-list(cor = m, p = m)
    for(i in 1:n1){
      f<-function(x){
        c.i<-cor.test(v1[,i],x,method = method,use = use, alternative = alternative);
        c(c.i$estimate,c.i$p.value)}
      c.i <- apply(v2,2,f)
      results$cor[i,] <- c.i[1,]
      results$p[i,] <- c.i[2,]
    }
    if(ncol(v2)==1){
      results<-cbind(results$cor,results$p)
      colnames(results)<-c('R','P')
    }
  }
  if(upper.tri.flag){
    results$up <- cbind(results$cor[upper.tri(results$cor)],
                        results$p[upper.tri(results$p)])
  }
  return(results)
}

fisher.combine <- function(p){
  p.f<-apply(p,1,get.fisher.p.value)
  return(p.f)
}

get.fisher.p.value<-function(p){
  p<-p[!is.na(p)]
  if(length(p)==1){
    p.fisher=p
  }else{
    p.fisher<- 1 - pchisq(-2*sum(log(p),na.rm = T), 2*sum(!is.na(p)))    
  }
  return(p.fisher)
}

get.onesided.p.value <- function(c,p){
  p[p==0] <-1e-17
  p.one.side <- p
  p.one.side[] <- NA
  b<-c>0&!is.na(c)
  p.one.side[b]=p[b]/2
  b<-c<=0&!is.na(c)
  p.one.side[b]=1-(p[b]/2)
  return(p.one.side)
}

get.cor.zscores<-function(c,p){
  v<-cbind(get.onesided.p.value(c,p),get.onesided.p.value(-c,p))
  z<-get.p.zscores(v)
  return(z)
}

t.test.mat<-function(m,b,two.sided=F,rankf = F,fold.changeF = F){
  if(length(b)!=ncol(m)){
    print("Error. Inconsistent no. of samples.")
    return()
  }
  if(sum(b)<2||sum(!b)<2){
    return(get.mat(rownames(m),c('more','less',"zscores")))
  }
  if(two.sided){
    p<-as.matrix(apply(m,1,function(x) t.test(x[b],x[!b])$p.value))
  }else{
    p<-t(apply(m,1,function(x) c(t.test(x[b],x[!b],alternative = 'greater')$p.value,
                                 t.test(x[b],x[!b],alternative = 'less')$p.value)))
    colnames(p)<-c('more','less')
    p<-cbind(p,get.p.zscores(p))
    colnames(p)[3]<-"zscores"
  }
  if(rankf){
    p<-cbind(p,rank(p[,1]),rank(p[,2]))
    colnames(p)[4:5]<-c("rank.more","rank.less")
  }
  if(fold.changeF){
    p<-cbind.data.frame(p,pos.mean = rowMeans(m[,b]),neg.mean = rowMeans(m[,!b]))
    p$logFC<-log2(p$pos.mean/p$neg.mean)
    p$more.padj<-p.adjust(p$more)
    p$less.padj<-p.adjust(p$less)
  }
  rownames(p)<-rownames(m)
  return(p)
}

ranksum.test.mat<-function(m,b,zscores.flag = T,two.sided=F){
  if(two.sided){
    p<-as.matrix(apply(m,1,function(x) ranksum.test(x[b],x[!b])))
  }else{
    p<-t(apply(m,1,function(x) c(ranksum.test(x[b],x[!b],alternative = 'greater'),
                                 ranksum.test(x[b],x[!b],alternative = 'less'))))
    colnames(p)<-c('more','less')
    if(zscores.flag){
      p<-cbind(p,get.p.zscores(p))
      colnames(p)[3]<-"zscores"
    }
  }
  rownames(p)<-rownames(m)
  return(p)
}

ranksum.test<-function(v1,v2,alternative="two.sided"){
  p=NA
  if (sum(!is.na(v1))==0|sum(!is.na(v2))==0){
    return(p)
  }else{
    p<-wilcox.test(v1,v2,alternative = alternative)$p.value  
  }
  return(p)
}

ranksum.test.labels<-function(v,b,alternative="two.sided"){
  p=NA
  v<-v[!is.na(b)]
  b<-b[!is.na(b)]
  v1<-v[b];v2<-v[!b]
  if (sum(!is.na(v1))==0|sum(!is.na(v2))==0){
    return(p)
  }else{
    p<-wilcox.test(v1,v2,alternative = alternative)$p.value  
  }
  return(p)
}

t.test.labels<-function(v,b,alternative="two.sided"){
  p<-t.test(v[b],v[!b],alternative = alternative)$p.value
  return(p)
}

labels.2.colors<-function(x.class,x = NULL,number.flag = F,color.spec = "hsv"){
  palette("default")
  icol<-c("black","red","cadetblue","gray","darkgreen","darkorange","darkviolet","gold3",
          "lightpink","deeppink2","deepskyblue",palette(),rainbow(50))
  no.classes<-length(unique(x.class))
  if(number.flag){
    icol<-match(x.class,sort(unique(x.class)))
  }else{
    if(is.numeric(x.class[1])){
      icol<-plotrix::color.scale(x.class,c(0,10),0.8,0.8,color.spec = color.spec)
    }else{
      icol<-icol[match(x.class,sort(unique(x.class)))]
    }
  }
  if(!is.null(x)){
    names(icol)<-x  
  }
  return(icol)
}

call_boxplot <- function (y,x,unique.x, f = median,lab.size = 16, title.size = 20,remove.legend = F,
                          ylab = '',xlab = '',main = '',labels=NULL,
                          legend.name = "",add.anova = F,blank.flag = T,
                          cex = 0.7,order.flag = T,b.ref = NULL,p.val.show = is.null(b.ref)){
  b<-is.infinite(y)|is.na(y)
  if(length(labels)==length(x)){
    labels<-labels[!b] 
  }
  y<-y[!b];x<-x[!b]
  if(p.val.show){
    a<-aov(y ~ as.factor(x))
    kt<-kruskal.test(y~as.factor(x))
    p.kt = format(kt$p.value,digits=2)
    p.anova=format(unlist(summary(a))['Pr(>F)1'],digits=2)
    if(add.anova){
      main = paste(main,'\n(ANOVA p-value = ',p.anova,'\nKruskal p-value =',p.kt,')',sep = '')
    }
  }
  if(!is.null(b.ref)){
    main <- paste0(main,"\n(",call_format.pval(t.test.labels(y,b.ref,alternative = "greater")),", ",
                   "AUC = ",round(get.auc(y,b.ref),2),")")
  }
  
  if(missing(unique.x)){
    unique.x<-unique(x)
    x.med<-apply(as.matrix(unique.x),1,function(xi) f(y[is.element(x,xi)]))
    names(x.med)<-unique.x
    if(order.flag){
      unique.x<-unique.x[order(x.med)]  
    }
    #labels<-labels[unique.x]
  }
  
  #idx<-order(y.med[y])
  x <- data.frame(name = x, val = y)
  x$name <- factor(x$name, levels = unique.x)
  p <- ggplot(x, aes(x = name, y = val,fill = labels)) + #geom_point(stat="identity") +
    labs(y = ylab, title = main, x = xlab) + 
    geom_boxplot() + theme_bw()
  p <- p + labs(color = legend.name) + scale_fill_discrete(name = legend.name)
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1,size = lab.size),
                 axis.title=element_text(size = lab.size),
                 title = element_text(size = title.size),
                 legend.title=element_text(size = lab.size), 
                 legend.text=element_text(size = lab.size))
  if(blank.flag){
    p <- p+ theme(panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank())
  }
  if(remove.legend){
    p<-p+theme(legend.position = "none")
  }
  return(p)
  #labs(title = paste(screen.name,'ANOVA p-value = ',p.anova)) + 
}

call_heatmap <- function(m,main = '',col.labels = NULL,row.labels = NULL,k = 3,filter.na = T,
                         cexCol = ifelse(ncol(m)>70,0.00001,1),
                         cexRow = ifelse(nrow(m)>70,0.00001,1),
                         m.value = '',scale = "none",
                         cluster.flag = "none",sym.flag = F,
                         xlab = "",ylab = "",row.col = NULL,col.col = NULL,legend.flag = T,
                         method = 'euclidean',symm = F, palette = "redblue"){
  if(filter.na){
    b.row<-!is.na(rowSums(m))
    m<-m[b.row,]
    if(!is.null(row.labels)){row.labels<-subset(row.labels,b.row)}
    b.col<-!is.na(colSums(m))
    m<-m[,b.col]
    if(!is.null(col.labels)){col.labels<-subset(col.labels,b.col)}
  }
  
  if(cluster.flag!="none"){
    if(method =="cor"){
      hc <- hclust(as.dist(2-cor(m)), method="complete");
      hr <- hclust(as.dist(2-cor(t(m))), method="complete");
    }else{
      hc <- hclust(dist(t(m),method = method), method="complete");
      hr <- hclust(dist(m,method =method), method="complete");
    }
    if(sym.flag){
      hc<-hr
    }
    Rowv <- as.dendrogram(hr)
    Colv <- as.dendrogram(hc)
    if(!is.null(col.labels)){
      col.labels<-cbind.data.frame(col.labels,clusters = paste0("C",cutree(hc, k = k)))
    }else{
      col.labels<-cbind.data.frame(clusters = paste0("C",cutree(hc, k = k)))
    }
    
  }else{
    Rowv <- NA;Colv <- NA
    hc<-T;hr<-T
  }
  
  col.col <- NULL
  if(!is.null(col.labels)&is.null(col.col)){
    col.labels<-as.data.frame(col.labels)
    col.col<-(apply(col.labels,2,labels.2.colors))
    colnames(col.col)<-colnames(col.labels)
  }
  if(!is.null(row.labels)&is.null(row.col)){
    row.labels<-as.data.frame(row.labels)
    row.col<-t(as.matrix(t(laply(1:ncol(row.labels),function(i) return(labels.2.colors(row.labels[,i]))))))
    if(ncol(row.labels)==1){row.col<-t(row.col)}
    row.labels<-t(as.matrix(row.labels))
    rownames(row.col)<-rownames(row.labels)
  }
  # myheatcol <- hcl.colors(100, palette = palette, alpha = NULL, rev = T, fixup = T)
  
  myheatcol <- redblue(50)
  # myheatcol <- colorpanel(50, "green", "white", "pink")
  
  myheatcol <- c(rep(myheatcol[1],10),myheatcol,rep(myheatcol[length(myheatcol)],10))
  myheatcol<-myheatcol[seq(length(myheatcol),1,-1)]
  call_heatmap1(m,main,Rowv=Rowv, Colv=Colv,
                m.value,cexRow,cexCol,myheatcol,scale,col.col,row.col,
                cluster.flag = cluster.flag,
                xlab = xlab,ylab = ylab,symm = symm)
  if(!legend.flag||(is.null(col.labels)&&is.null(row.labels))){return(col.labels)}
  legend.place<-rep(c("topright","right","bottomright","left",
                      "bottomleft","topleft"),5)
  legend.place<-rep(c("topright","right","bottomright","left",
                      "bottomleft","topleft"),5)
  if(!is.null(row.labels)&&identical(col.labels,t(row.labels))){
    row.labels<-NULL
  }
  
  if(!is.null(col.labels)){
    coltitles<-colnames(col.labels)
    col.col<-as.matrix(col.col[,!duplicated(t(col.labels))])
    col.labels<-as.matrix(col.labels[,!duplicated(t(col.labels))])
    n1<-ifelse(is.matrix(col.labels),ncol(col.labels),1)
    
    for(i in 1:n1){
      v<-unique(paste(col.col[,i],col.labels[,i],sep = "?"))
      # legend(legend.place[i],legend = get.strsplit(v,"?",2),col=get.strsplit(v,"?",1),
      #        pch = 15,cex = 0.4,title = coltitles[i],title.adj = 0.1)
      legend(x = (1+i)*0.1,y = 0.99,legend = get.strsplit(v,"?",2),col=get.strsplit(v,"?",1),
             pch = 15,cex = 0.4,title = coltitles[i],title.adj = 0.1)
    }
    legend.place<-legend.place[(i+1):length(legend.place)]
  }else{n1<-0}
  if(!is.null(row.labels)){
    rowtitles<-rownames(row.labels)
    n2<-ifelse(is.matrix(row.labels),nrow(row.labels),1)
    for(i in 1:n2){
      v<-unique(paste(row.col[i,],row.labels[i,],sep = "?"))
      # legend(legend.place[i],legend = get.strsplit(v,"?",2),
      #        col=get.strsplit(v,"?",1),pch = 15,cex = 0.4,
      #        title = rowtitles[i],title.adj = 0.1)
      
      legend(x = (1+n1+i)*0.1,y = 0.99,legend = get.strsplit(v,"?",2),
             col=get.strsplit(v,"?",1),pch = 15,cex = 0.4,
             title = rowtitles[i],title.adj = 0.1)
    }
  }
  return(list(hr = hr,col.labels = col.labels))
}

call_heatmap1<-function(m,main,Rowv,Colv,m.value,cexRow,cexCol,
                        myheatcol,scale,col.col = NULL,row.col = NULL,
                        cluster.flag = "none",
                        xlab = "",ylab = "",symm = F){
  ColSideColorsSize = 2
  if(cluster.flag=="none"){Rowv = NA;Colv = NA}
  if(cluster.flag=="row"){Colv = NA}
  if(cluster.flag=="col"){Rowv = NA}
  if(is.null(col.col)&is.null(row.col)){
    heatmap.3(m, main = main,  Rowv = Rowv,Colv = Colv,
              col = myheatcol, density.info="none",margins=c(10,10),
              ColSideColorsSize = ColSideColorsSize,
              RowSideColorsSize = 0.5,scale = scale,#RowAxisColors = 1,
              key=TRUE,
              KeyValueName = m.value,keysize = 1,
              symm = symm,cexRow=cexRow,cexCol=cexCol,
              dendrogram = cluster.flag,
              xlab = xlab,ylab = ylab,symbreaks = T)
    return()
  }
  if(is.null(col.col)){col.col<-as.matrix(rep("black",ncol(m)));ColSideColorsSize = 0.5}
  if(is.null(row.col)){row.col<-t(as.matrix(rep("black",nrow(m))));RowSideColorsSize = 0.5}
  heatmap.3(m, main = main,  Rowv = Rowv,Colv = Colv,
            col=myheatcol, density.info="none",margins=c(10,10),
            ColSideColorsSize = ColSideColorsSize,
            RowSideColorsSize = ColSideColorsSize,scale = scale,#RowAxisColors = 1,
            key=TRUE , KeyValueName = m.value,symm = symm,cexRow=cexRow,cexCol=cexCol,
            ColSideColors=col.col,RowSideColors = row.col,dendrogram = cluster.flag,
            xlab = xlab,ylab = ylab)
}

get.top.elements <- function (m,q = 100,min.ci = NULL,main = "",sort.flag = T){
  top.l<-list()
  v<-rownames(m)
  for (i in 1:ncol(m)){
    mi<-m[,i];mi<-mi[!is.na(mi)]
    idx<-order(mi,decreasing = F)
    ci <- mi[idx[min(q,length(mi))]]
    ci <- min(ci,min.ci)
    b <- m[,i]<=ci
    b[is.na(m[,i])]<-F
    if(sort.flag){
      top.l[[i]]<-sort(v[b])
    }else{
      top.l[[i]]<-v[b][order(m[b,i])]
    }
    
  }
  if(main!=""){main<-paste0(main,".")}
  names(top.l)<-paste0(main,colnames(m))
  return(top.l)
}

set.list<-function (r,b,name){
  set.field<-function (v,b){
    d <- dim(v)
    d.b<-length(b)
    if(!is.null(d)){
      if(d[1]==d.b){v <- v[b,]} #subset(v,subset = b)}
      if(d[2]==d.b){v <- v[,b]}
    }else{if(length(v)==d.b){v <- v[b]}}
    return(v)
  }
  rn<-lapply(r, set.field, b = b)
  if(!missing(name)){rn$name<-name}
  return(rn)
}

set.field <- function (v,b){
  d <- dim(v)
  d.b<-length(b)
  if(!is.null(d)){
    if(d[1]==d.b){v <- v[b,]}
    if(d[2]==d.b){v <- v[,b]}
  }else{if(length(v)==d.b){v <- v[b]}}
  return(v)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

plot.3aucs<-function(p,y,b,names,main){
  p1<-plot.auc(p,y,plot.flag = F,subplotF = F)
  p2<-plot.auc(p[b],y[b],plot.flag = F,subplotF = F)
  p3<-plot.auc(p[!b],y[!b],plot.flag = F,subplotF = F)
  a<-c(get.auc(p,y),
       get.auc(p[b],y[b]),
       get.auc(p[!b],y[!b]))
  a<-round(a,2)
  plot(p1,ylim = c(0,1),main = main)
  plot(p2,ylim = c(0,1),col = "red",add = T)
  plot(p3,ylim = c(0,1),col = "blue",add = T)
  abline(a = 0,b = 1,col = "gray30")
  legend(x = 0.4,y = 0.3,col = c("black","red","blue"),
         legend = paste0(names," (AUC = ",a,")"),lty = 1,lwd = 3,
         cex = 0.8)
}

plot.3aucs.mat<-function(P,Y,names1 = names(Y),main = ""){
  p1<-plot.auc(P[[1]],Y[[1]],plot.flag = F,subplotF = F)
  p2<-plot.auc(P[[2]],Y[[2]],plot.flag = F,subplotF = F)
  a<-c(get.auc(P[[1]],Y[[1]]),
       get.auc(P[[2]],Y[[2]]))
  a<-round(a,2)
  plot(p1,ylim = c(0,1),main = main)
  plot(p2,ylim = c(0,1),col = "red",add = T)
  col<-c("black","red")
  if(length(P)>2){
    p3<-plot.auc(P[[3]],Y[[3]],plot.flag = F,subplotF = F)
    a<-c(a,get.auc(P[[3]],Y[[3]]))
    plot(p3,ylim = c(0,1),col = "blue",add = T)
    col<-c("black","red","blue")
  }
  
  abline(a = 0,b = 1,col = "gray30")
  legend(x = 0.4,y = 0.3,col = col,
         legend = paste0(names1," (AUC = ",a,")"),lty = 1,lwd = 3,
         cex = 0.8)
}

plot.auc<-function(p1,y1,main = "",subplotF = T,precF = F,add = F,col = "black",plot.flag = T,cex = 20){
  if(subplotF){
    par(mfrow=c(1,2),oma = c(0, 0, 3, 0))
  }
  pr <- ROCR::prediction(p1, y1)
  auc <- performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  if(!plot.flag){return(prf)}
  plot(prf,main = paste0(main,'\nAUC =',round(auc,digits = 2)),
       add = add,col = col)
  abline(a = 0,b = 1,col = "gray30")
  if(precF){
    prf <- performance(pr, measure = "prec", x.measure = "rec")
    plot(prf,ylim = c(0,1))
    abline(h = mean(y1))
  }
  
  return(auc)
}



