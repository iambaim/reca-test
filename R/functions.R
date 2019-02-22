#################################################################            
#' @export   
compress.mat<-function(caa.mat,weight)
{
  char.mat<-NULL
  for(i in 1:ncol(caa.mat))
    char.mat<-paste(char.mat,caa.mat[,i])
  sum.wt<-tapply(weight,char.mat,sum)
  row.no<-match(names(sum.wt),char.mat)
  if(length(row.no)>1)
    new.mat<-cbind(caa.mat[row.no,],sum.wt)
  else
    new.mat<-matrix(c(caa.mat[row.no,],sum.wt),nrow=1)
  new.mat<-new.mat[sum.wt>0,]
  new.mat
}
####################################################################
#' @export   
make.constr.fit<-function(nage,cov.lev,i.obs)
{
  A<-make.constr.matrix(c(nage,cov.lev))
  i.obs2<-i.obs
  if(nage>1)i.obs2<-rep(i.obs,times=nage)+
    prod(cov.lev)*(rep(1:nage,each=length(i.obs))-1)
  n.o2 =length(i.obs2)
  A.o = A[,i.obs2]
  V.o = diag(1,n.o2)-t(A.o)%*%solve(A%*%t(A),A.o)
  e = eigen(V.o,symmetric=TRUE)
  r = sum(abs(e$val)>0.0000001)
  n.constr = n.o2-r
  if(n.constr>0)
  {
    e$val[(r+1):n.o2] = e$val[r]
    CC = t(e$vec[,(r+1):n.o2])
  }
  else
  {
    CC = matrix(1)
  }
  x<-e$val
  Winv = e$vec%*%diag(1/x,nrow=length(x))%*%t(e$vec)
  
  return(list(Sigma=Winv,constr=CC,n.constr.cell=as.integer(n.constr)))
}
####################################################################
#' @export   
make.constr.predict<-function(nage,cov.lev,i.obs,i.unobs){
  if(length(i.unobs)==0)
    return(list(E.u=matrix(0),C.u=matrix(0)))
  #Make constraint matrix
  A <- make.constr.matrix(c(nage,cov.lev))
  ncell<-prod(cov.lev)
  i.obs2<-i.obs
  i.unobs2<-i.unobs         #Unobserved cells of interest
  i.Unobs2<-c(1:ncell)[-i.obs]  #All unobserved cells
  #Assumes observed and unobserved indices are repeated for each age group
  if(nage>1)
  {
    i.obs2<-rep(i.obs2,times=nage)+
      prod(cov.lev)*(rep(1:nage,each=length(i.obs2))-1)
    i.unobs2<-rep(i.unobs2,times=nage)+
      ncell*(rep(1:nage,each=length(i.unobs2))-1)
    i.Unobs2<-rep(i.Unobs2,times=nage)+
      ncell*(rep(1:nage,each=length(i.Unobs2))-1)
  }
  A.o <- A[,i.obs2]
  #First consider all unobserved cells
  A.u <- A[,-i.obs2,drop=F]
  B <- solve(A.u%*%t(A.u),A.u)
  E.u <- -t(B)%*%A.o
  V.u <- diag(1,ncol(A.u))-t(A.u)%*%B
  #Pick out relevant ones
  ind <- match(i.unobs2,i.Unobs2)
  E.u <- E.u[ind,]
  V.u <- V.u[ind,ind]
  #Check effective rank of V.u
  V.u.eig <- eigen(V.u,symmetric=TRUE)
  r <- sum(V.u.eig$val> 1e-12)
  if(r==0)
  {#Singular V.u, cells deterministically given
    C.u <- matrix(0)
  }
  else
  {
    x<-sqrt(V.u.eig$val[1:r])
    C.u <- V.u.eig$vec[,1:r,drop=TRUE]%*%diag(x,nrow=length(x))
  }
  return(list(E.u=E.u,C.u=C.u))
}

####################################################################
#' @export   
make.constr.predict.old<-function(nage,cov.lev,i.obs,i.unobs)
{
  A<-make.constr.matrix(c(nage,cov.lev))
  i.obs2<-i.obs
  if(nage>1)
    i.obs2<-rep(i.obs,times=nage)+
    prod(cov.lev)*(rep(1:nage,each=length(i.obs))-1)
  n.o2 =length(i.obs2)
  A.o = A[,i.obs2]
  ncell<-prod(cov.lev)
  n.u<-length(i.unobs)
  i.unobs2<-i.unobs
  if(nage>1)i.unobs2<-rep(i.unobs,times=nage)+
    ncell*(rep(1:nage,each=length(i.unobs))-1)
  n.u2 = length(i.unobs2)
  if(n.u2>0)
  {
    A.u = A[,i.unobs2]
    B = t(A.u)%*%ginv(A.u%*%t(A.u))
    E.u = -B%*%A.o
    V.u = diag(1,n.u2)-B%*%A.u
    V.u.eig = eigen(V.u,symmetric=TRUE)
    r = sum(V.u.eig$val> (dim(V.u)[1]*V.u.eig$val[1]*1e-10))
    if(r>0){
      x<-sqrt(V.u.eig$val[1:r])
      C.u = V.u.eig$vec[,1:r]%*%diag(x,nrow=length(x))}
    else {C.u<-matrix(0,nrow=length(i.unobs2),ncol=1)
    V.u<-matrix(0,nrow=length(i.unobs2),ncol=length(i.unobs2))}
  }
  else
  {
    E.u = matrix(1)
    C.u = matrix(1)
  }
  return(list(E.u=E.u,C.u=C.u))
}
####################################################################
#' @export   
make.constr.matrix<-function(n.lev){
  d = n.lev[n.lev>1]
  ld = length(d)
  if(ld==2)
    A = make.constr.matrix2(d)
  else if(ld==3)
    A = make.constr.matrix3(d)
  else if(ld==4)
    A = make.constr.matrix4(d)
  else if(ld==5)
    A = make.constr.matrix5(d)
  else
    cat("Not valid number of factors for use of cell effects")
  
  A
}
####################################################################
#' @export   
enumerate.cell = function(cov,cov.lev)
{
  if(is.null(dim(cov)))cov<-matrix(cov,nrow=1)
  pr = 1
  cell = 0
  for(i in 1:length(cov.lev)){
    if(cov.lev[i]>1)
    {
      cell = cell + pr * (cov[,i]-1)
      pr = pr*cov.lev[i]
    }
  }
  cell = cell+1
  cell
}
####################################################################
#' @export   
make.predict.cell = function(ncat,biotic,landings)
{
  in.landings<-biotic$info[,"in.landings"]==1
  interaction<-biotic$info[,"interaction"]==1&in.landings
  if(sum(interaction)>1){
    landings.names<-rownames(biotic$info)[in.landings]
    interaction.names<-rownames(biotic$info)[interaction]
    predict.int.cov<-landings[,interaction.names]
    biotic.int.cov<-biotic$CovariateMatrix[,interaction.names]
    nlev.int<-biotic$info[,"nlev"][interaction]
    biotic.cell<-unique(enumerate.cell(biotic.int.cov,nlev.int))
    predict.cell<-enumerate.cell(predict.int.cov,nlev.int)
    ## predict.cell used directly - corresponds to fitted cells
    ## res and n.u, n.o not used - but will be removed later
    
    res<-list(E.u=matrix(1),V.u=matrix(1),C.u=matrix(1))
    n.u<-n.o<-0
    if(0){ 
      ### CHANGED 12/1-18
      ### Use all combination of cells - from enumerate.cell - same cell number as in fitting
      
    unobs.cell<-setdiff(predict.cell,biotic.cell)
    n.o<-length(biotic.cell)
    n.u<-length(unobs.cell)
    obs.cell<-is.element(predict.cell,biotic.cell)
    if(n.u>0)
    {
      res<-make.constr.predict(ncat,nlev.int,biotic.cell,unobs.cell)
    }
    else
    {
      ind.u = 0
      res<-list(E.u=matrix(1),V.u=matrix(1),C.u=matrix(1))
    }
    predict.cell[obs.cell]<-pmatch(predict.cell[obs.cell],sort(biotic.cell),duplicates.ok=T)
    predict.cell[!obs.cell]<-length(unique(biotic.cell))+
      as.integer(as.factor(predict.cell[!obs.cell]))
    }
  }
  else {predict.cell<-NULL
  res<-list(E.u=matrix(1),V.u=matrix(1),C.u=matrix(1))
  n.u<-n.o<-0}
  list(predict.cell=predict.cell,
       cell.u.dist=list(n.u=n.u*ncat,n.o=n.o*ncat,
                        E=res$E.u,C=res$C.u))
}
################################################################            
#' @export   
make.mod.ind<-function(model){
  ind<-c(model$year,2*model$seas,3*model$gear,4*model$area)
  if(sum(ind)>0) ind<-ind[ind>0]
  else ind<-NULL
  ind
}
#################################################################            
#' @export   
make.tot.factors <- function(nAges,cov.lev,real.cell,predict.cov,model.list,inc.haulsize){
  null.dist<-list(n.u=0,n.o=0,ind.u=0,E=matrix(1),C=matrix(1))
  cell.u.dist =list(age=list(Int=null.dist),
                    lga=list(Int=null.dist,Slp=null.dist),
                    wgl=list(Int=null.dist,Slp=null.dist),
                    hsz=list(Int=null.dist))
  fac<-list(age=list(Int=-1),
            lga=list(Int=-1,Slp=-1),
            wgl=list(Int=-1,Slp=-1),
            hsz=list(Int=-1))
  
  cell.mat<-NULL
  cell.ind<-ncol(predict.cov)+1
  
  ind1<-c(1,2,2,3,3,4)#age,lga,wgl,hsz
  ind2<-c(1,1,2,1,2,1)#int,slp
  model.list[[4]]<-model.list[[2]]
  cov.lev[[4]]<-cov.lev[[2]]
  
  for(i in 1:6){
    model1<-model.list[[ind1[i]]]
    n.model<-ind2[i]
    if(length(model1)>=n.model){
      use.model<-model1[[n.model]]
      if(!is.null(use.model)){
        ind<-make.mod.ind(use.model)
        if(use.model$cell&!is.null(ind)){
          real.cell[[4]]<-real.cell[[1]]  
          ind<-ind+4*(ind1[i]==3)
          if(ind1[i]==1)ncat<-nAges
          else ncat<-1
          foo<-make.predict.cell(ncat,cov.lev[[ind1[i]]][[ind2[i]]][2:(length(ind)+1)],predict.cov[,ind])
          cell.mat<-cbind(cell.mat,foo$predict.cell)
          ind<-c(ind,cell.ind)
          cell.ind<-cell.ind+1
          cell.u.dist[[ind1[i]]][[ind2[i]]]<-foo$cell.u.dist
        }}
      fac[[ind1[i]]][[ind2[i]]]<-c(1,ind+1)
    }
    else fac[[ind1[i]]][[ind2[i]]]<--1
  }
  
  
  tot.factors<-cbind(1,predict.cov,cell.mat,-1)
  
  ind.haul = dim(tot.factors)[2]
  for(i in 1:6){
    if(i<6)model1<-model.list[[ind1[i]]]
    else model1<-model.list[[ind1[2]]]
    n.model<-ind2[i]
    if(length(model1)>=n.model){
      use.model<-model1[[n.model]]
      if(!is.null(use.model)){
        if(i==1&use.model$haul&inc.haulsize)
          fac[[ind1[i]]][[ind2[i]]]<-c(fac[[ind1[i]]][[ind2[i]]],1)
        if(use.model$haul){
          fac[[ind1[i]]][[ind2[i]]]<-c(fac[[ind1[i]]][[ind2[i]]],ind.haul)
        }
        if(use.model$boat&i!=6){
          fac[[ind1[i]]][[ind2[i]]]<-c(fac[[ind1[i]]][[ind2[i]]],ind.haul)
        }
        
      }}}
  
  list(factors=tot.factors,fac=fac,cell.u.dist=cell.u.dist)
}
#################################################################            
#' @export   
make.predict.cov<-function(caa.vars,covtab,weight,model){
  year.ind<-match(caa.vars$year,changeyears(covtab$Year$Realname))
  seas.ind<-match(caa.vars$seas,covtab$Season$Realname)
  gear.ind<-match(caa.vars$gear,covtab$Gear$Realname)
  area.ind<-match(caa.vars$area,covtab$Area$Realname)
  cov<-NULL
  if(model$age$Int$year)cov<-cbind(cov,covtab$Year$index[year.ind])
  else cov<-cbind(cov,rep(1,length(caa.vars$year)))
  if(model$age$Int$seas)cov<-cbind(cov,covtab$Season$index[seas.ind])
  else cov<-cbind(cov,rep(1,length(caa.vars$year)))
  if(model$age$Int$gear)cov<-cbind(cov,covtab$Gear$index[gear.ind])
  else cov<-cbind(cov,rep(1,length(caa.vars$year)))
  if(model$age$Int$area)cov<-cbind(cov,covtab$Area$index[area.ind])
  else cov<-cbind(cov,rep(1,length(caa.vars$year)))
  
  if(model$wgl$Int$year)cov<-cbind(cov,covtab$Year$index[year.ind])
  else cov<-cbind(cov,rep(1,length(caa.vars$year)))
  if(model$wgl$Int$seas)cov<-cbind(cov,covtab$Season$wglindex[seas.ind])
  else cov<-cbind(cov,rep(1,length(caa.vars$year)))
  if(model$wgl$Int$gear)cov<-cbind(cov,covtab$Gear$wglindex[gear.ind])
  else cov<-cbind(cov,rep(1,length(caa.vars$year)))
  if(model$wgl$Int$area)cov<-cbind(cov,covtab$Area$wglindex[area.ind])
  else cov<-cbind(cov,rep(1,length(caa.vars$year)))
  cov<-cbind(cov,caa.vars$seas)
  cov<-compress.mat(cov,weight)
  if(is.null(nrow(cov)))cov<-matrix(cov,nrow=1)
  
  
  list(predict.cov=matrix(cov[,1:8],ncol=8),season=cov[,9],catch=cov[,10])
}
#################################################################            
#' @export   
make.caa.params<-function(stoxdata,common)
{
  length.int<-common$lengthresCM
  maxlength<-common$maxlength
  catch<-stoxdata$Landings$LiveWeightKG
  qq<-(1:floor(maxlength/length.int))*length.int
  l.int<-log(qq)
  N.lint = length(l.int)
  nMCvec<-ceiling(common$npredictMC*catch/sum(catch))
  inc.cts.var<-sum(stoxdata$AgeLength$info[,"continuous"]*(1-stoxdata$AgeLength$info[,"in.landings"]))>0
  
  list(int=as.integer(length.int),
       l.int=as.double(l.int),
       N.lint=as.integer(N.lint),
       burnin=as.integer(common$caa.burnin),
       nMCvec=as.integer(nMCvec),
       inputfolder=common$inputfolder,
       filename.mcmc1=common$filename.mcmc1,
       filename.mcmc2=common$filename.mcmc2,
       filename.hsz.mcmc2=common$filename.hsz.mcmc2,
       filename.hsz.it=common$filename.hsz.it,
       filename.predict=common$filename.predict,
       inc.haulsize=as.integer(inc.cts.var)
  )
}

####################################################################
#' @export   
new.make.cov <- function(model,factors,ncat,rand.seas,rand.gear,Covtabs,inc.haulsize=F,
                         old.cell=NULL,wgl=F)
{
  cov<-cbind(factors$year,factors$seas,factors$gear,factors$area)
  nlev<-NULL
  nlev<-c(nlev,length(unique(factors$year)))
  if(rand.seas){
    if(!wgl)nlev<-c(nlev,max(Covtabs$Season$index))
    if(wgl)nlev<-c(nlev,max(Covtabs$Season$wglindex))}
  else nlev<-c(nlev,length(unique(factors$seas)))
  if(rand.gear){
    if(!wgl)nlev<-c(nlev,max(Covtabs$Gear$index))
    if(wgl)nlev<-c(nlev,max(Covtabs$Gear$wglindex))}
  else nlev<-c(nlev,length(unique(factors$gear)))
  if(!wgl)nlev<-c(nlev,max(Covtabs$Area$index))
  if(wgl)nlev<-c(nlev,max(Covtabs$Area$wglindex))
  ind<-make.mod.ind(model)
  cov<-cbind(1,cov[,ind])
  fix<-1
  if(model$year)fix<-c(fix,1)
  if(model$seas&!rand.seas)fix<-c(fix,1)
  if(model$seas&rand.seas)fix<-c(fix,0)
  if(model$gear&!rand.gear)fix<-c(fix,1)
  if(model$gear&rand.gear)fix<-c(fix,0)
  fix<-c(fix,rep(0,sum(c(model$area,model$boat,model$cell))))
  if(model$area)
    ispat<-1+sum(c(model$year,model$seas,model$gear,model$area))
  else ispat<--1
  if(!is.null(ind))nFac<-c(1,nlev[ind])
  else nFac<-1
  #  for(i in 1:ncol(cov))
  #    nFac<-c(nFac,length(unique(cov[,i])))
  #  if(model$area)
  #    nFac[ispat]<-nArea
  if(model$cell)
  {real.cell<-enumerate.cell(cov[,-1],nFac[-1])
  if(is.null(old.cell))cell<-as.integer(as.factor(real.cell))
  else cell<-make.new.cell(real.cell,old.cell)
  res = make.constr.fit(ncat,nFac,sort(unique(cell)))
  cov = cbind(cov,cell)
  icell<-ncol(cov)
  nFac<-c(nFac,max(cell))
  }
  else
  {
    res = list(cell=0,Sigma=matrix(1),constr=matrix(1),n.constr.cell=0)
    icell = -1
    real.cell<-NULL
  }
  if(model$boat){cov<-cbind(cov,factors$boat)
  iboat<-ncol(cov)
  nFac<-c(nFac,length(unique(factors$boat)))}
  else iboat<--1
  if(inc.haulsize){
    if(model$haul){cov<-cbind(cov,-9)
    ihaulsize<-ncol(cov)
    nFac<-c(nFac,1)
    fix<-c(fix,1)}
    else ihaulsize<--1
  }
  else ihaulsize<--1
  
  if(model$haul)
  {
    nHaul<-nrow(cov)
    cov = cbind(cov,1:nHaul)
    ihaul = ncol(cov)
    nFac<-c(nFac,nHaul)
    fix<-c(fix,0)
  }
  else
    ihaul = -1
  list(cov=cov,ncov=ncol(cov),ispat=ispat,nFac=nFac,fix=fix,icell=icell,ihaul=ihaul,
       iboat=iboat,ihaulsize=ihaulsize,
       Sigma.cell=res$Sigma,constr.cell=res$constr,
       n.constr.cell=res$n.constr.cell,real.cell=real.cell)
}
####################################################################
#' @export   
make.new.cell<-function(new.cell,old.cell){
  in.old<-is.element(new.cell,old.cell)
  new.cell[in.old]<-as.integer(as.factor(new.cell[in.old]))
  n.old<-max(old.cell)
  new.cell[!in.old]<-as.integer(as.factor(new.cell[!in.old]))+n.old
  new.cell
}
#######################################################################
#' @export   
new.make.cov.info<-function(model,factors,nAges,rand.seas,
                            rand.gear,Covtabs,inc.haulsize=F,old.cell=NULL,wgl=F)
{
  real.cell<-NULL
  if(!is.null(model))
  {
    mcov = new.make.cov(model,factors,nAges,rand.seas,rand.gear,Covtabs,inc.haulsize,old.cell,wgl)
    mcov$cov = as.matrix(mcov$cov)
    if(!is.null(mcov$real.cell))
      real.cell<-sort(unique(mcov$real.cell))
    cov.lev<-mcov$nFac
  }
  else
  {
    mcov = list(ncov=0,ispat=-1,icell=-1,ihaul=-1,iboat=-1,nFac=0,fix=0,cov=0,
                Sigma.cell=matrix(1),constr.cell=matrix(1),n.constr.cell=0)
    real.cell<-cov.lev<-NULL
  }
  list(mcov=mcov,real.cell=real.cell,cov.lev=cov.lev)
}
#######################################################################
#' @export   
make.constr.matrix2 = function(d)
{
  ncell = prod(d)
  df = prod(d-1)
  nconstr = ncell - df
  A = matrix(0,nrow=nconstr,ncol=ncell)
  cc = 0
  ind = 1:d[2]
  for(j in 1:d[1])
  {
    cc = cc+1
    A[cc,(j-1)*d[2]+ind] = 1
  }
  ind = (1:d[1])*d[2]-d[2]+1
  for(j in 1:(d[2]-1))
  {
    cc = cc+1
    A[cc,j-1+ind] = 1
  }
  A
}
####################################################################
#' @export   
make.constr.matrix3 = function(d)
{
  ncell = prod(d)
  df = prod(d-1)
  nconstr = ncell - df
  A = matrix(0,nrow=nconstr,ncol=ncell)
  cc = 0
  ind = 1:d[3]
  for(j in 1:d[1])
    for(k in 1:d[2])
    {
      cc = cc+1
      A[cc,(j-1)*d[2]*d[3]+(k-1)*d[3]+ind] = 1
    }
  ind = (1:d[2])*d[3]-d[3]+1
  for(j in 1:d[1])
    for(k in 1:(d[3]-1))
    {
      cc = cc+1
      A[cc,(j-1)*d[2]*d[3]+k-1+ind] = 1
    }
  ind = (1:d[1])*d[2]*d[3]-d[2]*d[3]+1
  for(j in 1:(d[2]-1))
    for(k in 1:(d[3]-1))
    {
      cc = cc+1
      A[cc,(j-1)*d[3]+k-1+ind] = 1
    }
  A
}
####################################################################
#' @export   
make.constr.matrix4 = function(d)
{
  ncell = prod(d)
  df = prod(d-1)
  nconstr = ncell - df
  A = matrix(0,nrow=nconstr,ncol=ncell)
  cc = 0
  ind = 1:d[4]
  for(j in 1:d[1])
    for(k in 1:d[2])
      for(l in 1:d[3])
      {
        cc = cc+1
        A[cc,(j-1)*d[2]*d[3]*d[4]+(k-1)*d[3]*d[4]+(l-1)*d[4]+ind] = 1
      }
  ind = (1:d[3])*d[4]-d[4]+1
  for(j in 1:d[1])
    for(k in 1:d[2])
      for(l in 1:(d[4]-1))
      {
        cc = cc+1
        A[cc,(j-1)*d[2]*d[3]*d[4]+(k-1)*d[3]*d[4]+l-1+ind] = 1
      }
  ind = (1:d[2])*d[3]*d[4]-d[3]*d[4]+1
  for(j in 1:d[1])
    for(k in 1:(d[3]-1))
      for(l in 1:(d[4]-1))
      {
        cc = cc+1
        A[cc,(j-1)*d[2]*d[3]*d[4]+(k-1)*d[4]+l-1+ind] = 1
      }
  ind = (1:d[1])*d[2]*d[3]*d[4]-d[2]*d[3]*d[4]+1
  for(j in 1:(d[2]-1))
    for(k in 1:(d[3]-1))
      for(l in 1:(d[4]-1))
      {
        cc = cc+1
        A[cc,(j-1)*d[3]*d[4]+(k-1)*d[4]+l-1+ind] = 1
      }
  A
}
####################################################################
#' @export   
make.constr.matrix5 = function(d)
{
  ncell = prod(d)
  df = prod(d-1)
  nconstr = ncell - df
  A = matrix(0,nrow=nconstr,ncol=ncell)
  cc = 0
  ind = 1:d[5]
  for(j in 1:d[1])
    for(k in 1:d[2])
      for(l in 1:d[3])
        for(m in 1:d[4])
        {
          cc = cc+1
          A[cc,(j-1)*d[2]*d[3]*d[4]*d[5]+(k-1)*d[3]*d[4]*d[5]+(l-1)*d[4]*d[5]+(m-1)*d[5]+ind] = 1
        }
  ind = (1:d[4])*d[5]-d[5]+1
  for(j in 1:d[1])
    for(k in 1:d[2])
      for(l in 1:d[3])
        for(m in 1:(d[5]-1))
        {
          cc = cc+1
          A[cc,(j-1)*d[2]*d[3]*d[4]*d[5]+(k-1)*d[3]*d[4]*d[5]+(l-1)*d[4]*d[5]+m-1+ind] = 1
        }
  ind = (1:d[3])*d[4]*d[5]-d[4]*d[5]+1
  for(j in 1:d[1])
    for(k in 1:d[2])
      for(l in 1:(d[4]-1))
        for(m in 1:(d[5]-1))
        {
          cc = cc+1
          A[cc,(j-1)*d[2]*d[3]*d[4]*d[5]+(k-1)*d[3]*d[4]*d[5]+(l-1)*d[5]+m-1+ind] = 1
        }
  ind = (1:d[2])*d[3]*d[4]*d[5]-d[3]*d[4]*d[5]+1
  for(j in 1:d[1])
    for(k in 1:(d[3]-1))
      for(l in 1:(d[4]-1))
        for(m in 1:(d[5]-1))
        {
          cc = cc+1
          A[cc,(j-1)*d[2]*d[3]*d[4]*d[5]+(k-1)*d[4]*d[5]+(l-1)*d[5]+m-1+ind] = 1
        }
  ind = (1:d[1])*d[2]*d[3]*d[4]*d[5]-d[2]*d[3]*d[4]*d[5]+1
  for(j in 1:(d[2]-1))
    for(k in 1:(d[3]-1))
      for(l in 1:(d[4]-1))
        for(m in 1:(d[5]-1))
        {
          cc = cc+1
          A[cc,(j-1)*d[3]*d[4]*d[5]+(k-1)*d[4]*d[5]+(l-1)*d[5]+m-1+ind] = 1
        }
  A
}
####################################################################
#' @export   
make.kappa<-function(caadata){
  use<-!is.na(caadata$data$ALDER)
  haulstart<-match(unique(caadata$data$SERIENR[use]),caadata$data$SERIENR)
  x<-caadata$data$sgopt[haulstart]
  y<-caadata$data$ggopt[haulstart]
  z<-caadata$data$agopt[haulstart]
  x1<-length(unique(x))>1
  y1<-length(unique(y))>1
  z1<-length(unique(z))>1
  if(x1&y1&z1)m<-model.matrix(~x+y+z-1)
  if(x1&y1&!z1)m<-model.matrix(~x+y-1)
  if(x1&!y1&z1)m<-model.matrix(~x+z-1)
  if(!x1&y1&z1)m<-model.matrix(~y+z-1)
  tab<-NULL
  k<-0
  if(x1+y1+z1>1){
    tab<-table(x,y,z)
    k<-kappa(m)
  }
  list(k=k,tab=tab)
}

