apply.generic.mixed.logistic<-function(r,x,Y = r$tpm,MARGIN = 1){
  m<-apply(Y,MARGIN = MARGIN,function(y){generic.mixed.logistic(y>median(y),x,r)})
  if(MARGIN==1){
    m<-t(m)
  }
  colnames(m)<-c("Estimate","P")
  m<-cbind.data.frame(Z = get.cor.zscores(m[,"Estimate"],m[,"P"]),m)
  return(m)
}

apply.generic.HLM<-function(r,x,Y = r$tpm,MARGIN = 1,remove.dropouts = F,bulk.flag = F){
  if(remove.dropouts){
    m<-t(apply(Y,MARGIN = MARGIN,function(y){generic.HLM(y,x,r,b = y>0)}))
  }else{
    if(bulk.flag){
      m<-t(apply(Y,MARGIN = MARGIN,function(y){generic.HLM.bulk(y,x,r)}))
    }else{
      m<-t(apply(Y,MARGIN = MARGIN,function(y){generic.HLM(y,x,r)}))
    }
    
  }
  colnames(m)<-c("Estimate","P")
  m<-cbind.data.frame(Z = get.cor.zscores(m[,"Estimate"],m[,"P"]),m)
  return(m)
}

generic.HLM<-function(y,x,r0,b = NULL){
  r0$x<-x;r0$y<-y
  if(!is.null(b)){
    r0<-list(comp.reads = r$comp.reads[b], samples = r$samples[b],
             x = x[b],y = y[b])
  }
  f<-function(r0){
    M1 <- with(r0, lmer (y ~ comp.reads + (1 + comp.reads| samples) + x, mtcars))
    if(is.numeric(x)){
      c1<-summary(M1)$coef["x",]
    }else{
      c1<-summary(M1)$coef["xTRUE",]
    }
    idx<-match(c("Estimate","Pr(>|t|)"),names(c1))
    c1<-c1[idx]
    return(c1)
  }
  c1<-tryCatch({f(r0)},
               error = function(err){return(c(NA,NA))})
  return(c1)
}

generic.HLM.bulk<-function(y,x,r0,b = NULL, xc = ""){
  r0$x<-x;r0$y<-y
  if(!is.null(b)){
    r0<-list(samples = r$samples[b],
             x = x[b],y = y[b],xc = r[[xc]][b])
  }
  f<-function(r0){
    if(!is.null(r[[xc]])){
      M1 <- with(r0, lmer (y ~ (1 | samples) + x + xc, mtcars)) 
    }else{
      M1 <- with(r0, lmer (y ~ (1 | samples) + x, mtcars))
    }
    if(is.numeric(x)){
      c1<-summary(M1)$coef["x",]
    }else{
      c1<-summary(M1)$coef["xTRUE",]
    }
    idx<-match(c("Estimate","Pr(>|t|)"),names(c1))
    c1<-c1[idx]
    return(c1)
  }
  c1<-tryCatch({f(r0)},
               error = function(err){return(c(NA,NA))})
  return(c1)
}

generic.mixed.logistic<-function(y,x,r0,b = NULL){
  r0$x<-x;r0$y<-y
  if(!is.null(b)){
    r0<-list(comp.reads = r$comp.reads[b], samples = r$samples[b],
             x = x[b],y = y[b])
  }
  f<-function(r0){
    M1 <- with(r0, glmer (y ~ comp.reads + (1 | samples) + x, family = binomial(link="logit"),
                          control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=20000))))
    if(is.numeric(x)){
      c1<-summary(M1)$coef["x",]
    }else{
      c1<-summary(M1)$coef["xTRUE",]
    }
    idx<-match(c("Estimate","Pr(>|z|)"),names(c1))
    c1<-c1[idx]
    error.code<-summary(M1)$optinfo$conv$lme4$code
    print(error.code)
    if(!is.null(error.code)){
      c1<-c(NA,NA)
    }
    return(c1)
  }
  c1<-tryCatch({f(r0)},
               error = function(err){return(c(NA,NA))})
  return(c1)
}

apply.formula.mixed.logistic<-function(r,X,Y,MARGIN = 2,formula = "y ~ comp.reads + (1 | samples) + x",all.flag = F){
  if(is.matrix(Y)){
    m<-t(apply(Y,MARGIN = MARGIN,FUN = function(y){formula.mixed.logistic(r,y,X,formula = formula,all.flag = all.flag)}))
  }else{
    m<-t(apply(X,MARGIN = MARGIN,FUN = function(x){formula.mixed.logistic(r,Y,x,formula = formula,all.flag = all.flag)}))
  }
  if(!all.flag){
    colnames(m)<-c("Estimate","P")
    m<-cbind.data.frame(Z = get.cor.zscores(m[,"Estimate"],m[,"P"]),m)
  }
  return(m)
}

formula.mixed.logistic<-function(r,y,x,formula = "y ~ comp.reads + (1 | samples) + x",all.flag = F){
  r$x<-x;r$y<-y
  f<-function(r){
    M0 <- with(r, glmer (formula = formula, family = binomial(link="logit"),
                         control=glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=20000))))
    M1<-summary(M0)$coef[,c("Estimate","Pr(>|z|)")]
    if(!all.flag){
      if(is.numeric(x)){
        c1<-M1["x",]
      }else{
        c1<-M1["xTRUE",]
      }
    }else{
      c1<-get.cor.zscores(M1[,"Estimate"],M1[,"Pr(>|z|)"])
    }
    error.code<-summary(M0)$optinfo$conv$lme4$code
    print(error.code)
    if(!is.null(error.code)){
      # c1<-c(NA,NA)
    }
    print(c1)
    return(c1)
  }
  c1<-tryCatch({f(r)},
               error = function(err){return(c(NA,NA))})
  return(c1)
}

log.comp<-function(r){
  if(max(r$comp.reads)>2){
    print("Transforming the complexity values")
    r$max.comp<-c(max(r$comp.reads),max(r$comp))
    r$comp.reads<-log(r$comp.reads)
    r$comp.reads<-r$comp.reads/max(r$comp.reads)
    r$comp<-r$comp/max(r$comp)
  }
  return(r)
}

apply.formula.HLM<-function(r,X,Y,MARGIN = 1,formula = "y ~ (1 | samples) + x",ttest.flag = F){
  if(is.matrix(Y)){
    if(ttest.flag){
      m1<-t.test.mat(Y,X)
      b<-rowSums(p.adjust.mat(m1[,1:2])<0.1,na.rm = T)>0
      m1<-m1[b,];Y<-Y[b,]
      print(paste("Testing",sum(b),"genes that show a signal."))
    }
    m<-t(apply(Y,MARGIN = MARGIN,function(y){formula.HLM(y,X,r,formula = formula)}))
  }else{
    m<-t(apply(X,MARGIN = MARGIN,function(x){formula.HLM(Y,x,r,formula = formula)}))
  }
  colnames(m)<-c("Estimate","P","singular")
  m<-cbind.data.frame(Z = get.cor.zscores(m[,"Estimate"],m[,"P"]),m)
  if(ttest.flag){
    m<-cbind.data.frame(m,ttest = m1)
  }
  m<-cbind.data.frame(m,padj = p.adjust(m[,"P"],method = "BH"))
  return(m)
}

apply.formula.all.HLM<-function(r,X,Y,MARGIN = 1,formula = "y ~ (1 | samples) + x",ttest.flag = F){
  if(is.matrix(Y)){
    m<-t(apply(Y,MARGIN = MARGIN,function(y){
      P<-formula.HLM(y,X,r,formula = formula,return.all = T)
      Z<-c(get.cor.zscores(P[,"Estimate"],P[,"Pr(>|t|)"]),P[1,"singular"])
      names(Z)<-c(rownames(P),"singular")
      return(Z)
    }))
  }else{
    m<-t(apply(X,MARGIN = MARGIN,function(x){
      P<-formula.HLM(Y,x,r,formula = formula,return.all = T)
      Z<-get.cor.zscores(P[,"Estimate"],P[,"Pr(>|t|)"])
      return(Z)
      return(c(Z,Singular = P[1,"singular"]))
    }))
  }
  return(list(HLM = m,formula = formula))
}

formula.HLM<-function(y,x,r0, formula = "y ~ (1 | samples) + x",
                      val = ifelse(is.numeric(x),"","TRUE"),return.all = F){
  r0$x<-x;r0$y<-y
  f<-function(r0){
    M1 <- with(r0, lmer (formula = formula))
    if(return.all){
      c1<-summary(M1)$coef[,c("Estimate","Pr(>|t|)")]
      c1<-cbind.data.frame(c1,singular = isSingular(M1))
    }else{
      c1<-summary(M1)$coef[paste0("x",val),]
      idx<-match(c("Estimate","Pr(>|t|)"),names(c1))
      c1<-c1[idx]
      return(c(c1,singular = isSingular(M1)))
    }
    return(c1)
  }
  c1<-tryCatch({f(r0)},
               error = function(err){return(c(NA,NA,NA))})
  # print(c1)
  return(c1)
}

formula.HLM.multiple.genes<-function(y,prd.genes,r, formula = "y ~ (1 | samples) + comp.reads + "){
  X<-lapply(prd.genes, function(x) r$tpm[x,])
  names(X)<-prd.genes
  r0<-c(r,X)
  r0$y<-y
  formula<- paste0(formula,paste(prd.genes,collapse = " + "))
  M1 <- with(r0, lmer (formula = formula))
  # print(paste("Convergence:",M1@optinfo$conv$lme4$code))
  # print(paste("Convergence:",M1@optinfo$conv$lme4$messages))
  v<-summary(M1)$coef
  v<-cbind.data.frame(Z = get.cor.zscores(v[,"Estimate"],v[,"Pr(>|t|)"]),
                      v[,c("Estimate","Pr(>|t|)")])
  if(!is.null(M1@optinfo$conv$lme4$code)){
    print("Did not converge!")
    v$Z<-NA 
  }
  return(v)
}

formula.hurdle<-function(){
  
  library(pscl)
  mod.hurdle <- hurdle(visits ~ ., data = nmes)
  mod.hurdle <- hurdle(visits ~ ., data = nmes)
  
  
  mod.hurdle <- hurdle(y ~ x + r$comp.reads, random = )
  
  with(r0,hurdle(formula = "y ~ samples + x + comp.reads"))
  hurdle(formula = "y ~ (1 | samples) + x + comp.reads",data = r0)
  de3<-apply.formula.HLM(r,Y = r$tpm[rownames(de1),],X = r$treated,
                         MARGIN = 1,formula = "y ~ (1 | samples) + x + comp.reads")
  
  library(GLMMadaptive)
  r0<-set.list(r,sample.per.label(r0$cells,r0$samples,size = 10,boolean.flag = T))
  r0$y<-as.matrix(r0$cd["MT2A",])
  r0$x<-as.matrix(r0$treated)
  r0$samples<-as.matrix(r0$samples)
  r0$comp.reads<-as.matrix(r0$comp.reads)
  rownames(r0$x)<-rownames(r0$y)
  rownames(r0$samples)<-rownames(r0$y)
  rownames(r0$comp.reads)<-rownames(r0$y)
  x<-r0$x
  y<-r0$y
  km1 <- mixed_model(fixed = y ~ x + comp.reads, random = ~ 1 | samples,
                     zi_fixed = y ~ x + comp.reads,zi_random = ~ 1 | samples,
                     data = r0,
                     family = hurdle.lognormal(), n_phis = 1)
  
  r$x<-r$tpm["XIST",]
  r<-log.comp(r)
  
  mixed_model(y ~ x + comp.reads +gender, random = ~ 1 | sample, data = as.data.frame(r),
              family = zi.negative.binomial(),
              zi_fixed = ~ sex)
  g<-sig$Z.up[1]
  r0<-cbind.data.frame(y = r$tpm[g,],x = r$clusters=="C0",
                       comp.reads = r$comp.reads,samples = r$samples)
  v<-mixed_model(y ~ x + comp.reads, random = ~ 1 | samples, data = r0,
                 family = zi.poisson(),
                 zi_fixed = ~ x + comp.reads,
                 zi_random = ~ 1 | samples)
  v<-summary(v)
  v$coef_table
  v$coef_table_zi
  
  
  with(r0,mixed_model(fixed = y ~ x + comp.reads, random = ~ 1 | samples,
                      zi_fixed = y ~ x + comp.reads,zi_random = ~ 1 | samples,
                      data = r0,
                      family = hurdle.lognormal(), n_phis = 1))
  
  mixed_model(fixed = )
}



