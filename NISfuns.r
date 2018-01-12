"gen.data"<-function(funind,n=NULL,p=NULL){
if(is.null(n)) n=400
if(is.null(p)) p=1000
y0=NULL
if(funind<=4){
s=3*2^(funind-1)
u=matrix(rnorm(n*(p-50)),n,p-50)
v=matrix(rnorm(n*50),n,50)
w=u[,1:s]
wcoef=(-1)^(1:s+1)/5
x0=w%*%wcoef
x2=as.vector(x0)+sqrt(25-s)/5*v
x=cbind(u,x2)
betatrue=(-1)^(1:s+1)
y0=w%*%betatrue
y=y0+rnorm(n)*sqrt(3)
}  else if(funind==5){
s=3
x=matrix(rnorm(n*p),n,p)
x[,2]=-1/3*x[,1]^3+rnorm(n)
y0= x[,1]+x[,2]+x[,3]
y=y0+rnorm(n)*sqrt(3)
}
#if(funind<=4){
#s=3*2^(funind-1)
#u=matrix(rnorm(n*(p-50)),n,p-50)
#v=matrix(rnorm(n*50),n,50)
#w=u[,1:s]
#wcoef=(-1)^(1:s+1)/5
#x0=w%*%wcoef
#x2=as.vector(x0)+sqrt(25-s)/5*v
#x=cbind(u,x2)
#betatrue=(-1)^(1:s+1)
#y=w%*%betatrue+rnorm(n)
#}  else if(funind==5){
#s=3
#x=matrix(rnorm(n*p),n,p)
#x[,2]=-1/3*x[,1]^3+rnorm(n)
#y=x[,1]+x[,2]+x[,3]+rnorm(n)
#}
#
else if(funind==6){
s=4
x=matrix(runif(n*p),n,p)
y0=5*f1(x[,1])+3*f2(x[,2])+4*f3(x[,3])+6*f4(x[,4])
y=y0+rnorm(n)*sqrt(1.74)
} else if(funind==7){
s=4
w=matrix(runif(n*p),n,p)
u=runif(n)
x=(w+u)/2
y0=5*f1(x[,1])+3*f2(x[,2])+4*f3(x[,3])+6*f4(x[,4])
y=y0+rnorm(n)*sqrt(1.74)
} else if(funind==8){
s=12
x=matrix(runif(n*p),n,p)
y0=(f1(x[,1])+f2(x[,2])+f3(x[,3])+f4(x[,4])+
1.5*f1(x[,5])+1.5*f2(x[,6])+1.5*f3(x[,7])+1.5*f4(x[,8])+
2*f1(x[,9])+2*f2(x[,10])+2*f3(x[,11])+2*f4(x[,12]))
y=y0+rnorm(n)*sqrt(0.5184)
} else if(funind==9){
s=12
w=matrix(runif(n*p),n,p)
u=runif(n)
x=(w+u)/2
y0=(f1(x[,1])+f2(x[,2])+f3(x[,3])+f4(x[,4])+
1.5*f1(x[,5])+1.5*f2(x[,6])+1.5*f3(x[,7])+1.5*f4(x[,8])+
2*f1(x[,9])+2*f2(x[,10])+2*f3(x[,11])+2*f4(x[,12]))
y=y0+rnorm(n)*sqrt(0.5184)
} else if(funind==10){
s=4
x=matrix(rnorm(p*n, mean=0, sd=1), n, p)
betatrue <- c(2,2,2,-3*sqrt(2))
truerho=0.5
corrmat=diag(rep(1-truerho, p))+matrix(truerho, p, p)
corrmat[,4]=sqrt(truerho)
corrmat[4, ]=sqrt(truerho)
corrmat[4,4]=1
cholmat=chol(corrmat)
x=x%*%cholmat
y0 <- x[,1:s]%*%betatrue
y=y0+rnorm(n)
} else if(funind==11){
s=12
x=matrix(runif(n*p),n,p)
y=(f1(x[,1])+f2(x[,2])+f3(x[,3])+f4(x[,4])+
1.5*f1(x[,5])+1.5*f2(x[,6])+1.5*f3(x[,7])+1.5*f4(x[,8])+
2*f1(x[,9])+2*f2(x[,10])+2*f3(x[,11])+2*f4(x[,12]))+0.5*rnorm(n)*sqrt(0.5184)
} else if(funind==12){
s=12
w=matrix(runif(n*p),n,p)
u=runif(n)
x=(w+u)/2
y=(f1(x[,1])+f2(x[,2])+f3(x[,3])+f4(x[,4])+
1.5*f1(x[,5])+1.5*f2(x[,6])+1.5*f3(x[,7])+1.5*f4(x[,8])+
2*f1(x[,9])+2*f2(x[,10])+2*f3(x[,11])+2*f4(x[,12]))+0.5*rnorm(n)*sqrt(0.5184)
}  else if(funind>=18 & funind<=22)
{
s=4
x=matrix(runif(n*p),n,p)
y0=3*f1(x[,1])+3*f2(x[,2])+2*f3(x[,3])+2*f4(x[,4])
y=y0+rnorm(n)*sqrt(3.3843)*sqrt(2^(20-funind))
}   else if(funind>=28 & funind<=32)
{
s=4
w=matrix(runif(n*p),n,p)
u=runif(n)
x=(w+u)/2
y0=3*f1(x[,1])+3*f2(x[,2])+2*f3(x[,3])+2*f4(x[,4])
y=y0+rnorm(n)*sqrt(3.3843)*sqrt(2^(30-funind))
}
return(data=list(x=x,y=y,n=n,p=p,s=s,y0=y0))

}

f1<-function(x){x}
f2<-function(x){(2*x-1)^2}
f3<-function(x){sin(2*pi*x)/(2-sin(2*pi*x))}
f4<-function(x){0.1*sin(2*pi*x)+0.2*cos(2*pi*x)+0.3*sin(2*pi*x)^2+0.4*cos(2*pi*x)^3+0.5*sin(2*pi*x)^3}


simuNIS<-function(funind=NULL, repgroup=NULL, n=NULL, p=NULL, trace.it=NULL, knots=NULL, maxloop=NULL, eps0=NULL, randSeed = 0)
{
  if(is.null(funind)) stop('No given data generating scheme.')
  if(is.null(trace.it)) trace.it=FALSE
  if(is.null(maxloop)) maxloop=10
  if(is.null(eps0)) eps0=1e-6
  if(is.null(n)) n=400
  if(is.null(p)) p=1000
  if(is.null(knots)) knots=ceiling(n^0.2)

  trace.it=as.logical(trace.it)
  set.seed(randSeed)
    cat('Random seed=', randSeed, '...\n')
    data=gen.data(funind,n,p)
    testdata=gen.data(funind,n/2,p)
  
    greedINIS.fit = greedINIS(data, testdata,  knots=knots, eps0=eps0, DOISIS=TRUE,trace=trace.it)
    
    aINIS.fit = adaptINIS(data, testdata,  knots=knots, eps0=eps0, DOISIS=TRUE, maxloop=maxloop,trace=trace.it)
    

  return(list(greedyINIS = greedINIS.fit, INIS = aINIS.fit))
  #filename=paste(newpath,'/','funind=',funind,'_knots=',knots,'_r0=',r0,'_r1=',r1,'.RData',sep='')
  #save(aINIS.fit, greedINIS.fit, penGAM.fit,lmsis.fit,results,file=filename)
}


###########################################################################
#greedINIS: g-INIS algorithm, using knots=n^0.2
###########################################################################
greedINIS <- function(data, testdata=NULL, lambda.pen.list=NULL, folds=NULL, quant=NULL, gnum=1, kfold=NULL, knots=NULL,  eps0=1e-6, DOISIS=TRUE, maxloop=20, trace=FALSE, detailed=FALSE){
  t0=proc.time()[1]
  cat('starting greedINIS, g-INIS-penGAM algorithm, adatively choose number of variables\n')
  x=data$x
  y=data$y
  n <- nrow(x)
  p <- ncol(x)


  #if(is.null(nsis)) nsis=min(floor(n/log(n)),p-1)
  if(is.null(knots)) knots=ceiling(n^0.2)
  if(is.null(folds)) {
    temp= sample(1:n, n, replace = FALSE)
    if(is.null(kfold)) kfold=5
    for(i in 1:kfold){
      folds[[i]]=setdiff(1:n, temp[seq(i, n, kfold)])
    }
  }
  if(is.null(quant)) quant=1
  df0 <- knots+1

  xbs=matrix(0,n,df0*p)


  for(i in 1:p)
    xbs[,(i-1)*(df0)+(1:df0)]=ns(x[,i],df=df0)

  tempresi <- rep(0,p)

  curloop=1
  for(i in 1:p)
  {

    tempfit <-lm.fit (x=cbind(1,xbs[,(i-1)*df0+1:df0]),y=y)
    tempresi[i] <- sum(tempfit$residuals^2)
  }


  used.list<- tempresi

  used.sort <- sort(used.list, method= "sh", index=TRUE, decreasing=FALSE)
  initRANKorder <- used.sort$ix


  mindex <- sample(1:n)
  mresi=NULL
  for(i in 1:p){
    tempfit <-lm.fit (x=cbind(1,xbs[,(i-1)*df0+1:df0]),y=y[mindex])
    mresi[i] <- sum(tempfit$residuals^2)
  }
  resi.thres = quantile(mresi,1-quant)
  nsis <- min(sum(used.list<resi.thres), floor(n/df0/2),gnum)

  SISind <- sort(initRANKorder[1:nsis])
  if(!DOISIS) return (list(initRANKorder=initRANKorder, SISind=SISind, nsis=nsis))

  cat('loop ', curloop, '...SISind ', SISind,'\n')
  pick.ind=initRANKorder[1:nsis]


  xnew=x[,pick.ind]

  fit.tmp <- cv.penGAM(data=list(x=xnew,y=y), folds=folds, lambda.pen.list=NULL, knots=knots, eps0=eps0,cv.trace=trace)
  ISISind <- sort(pick.ind[fit.tmp$penGAMind])
  cat('loop ', curloop, '...ISISind ', ISISind,'\n')
  test= 1
  normal.exit = 1

  detail.pickind=NULL
  detail.ISISind=NULL
  detail.pickind=as.list(detail.pickind)
  detail.ISISind=as.list(detail.ISISind)

  detail.pickind[[curloop]]=pick.ind
  detail.ISISind[[curloop]]=ISISind

  if(length(ISISind)==0) {
    warnings('No variable was selected by this variant of ISIS.')
    normal.exit = 0
    if(detailed){
      return(list(fit=fit.tmp,initRANKorder=initRANKorder, detail.pickind=detail.pickind, detail.ISISind=detail.ISISind,
                  SISind=SISind, ISISind=ISISind,  nsis=nsis, normal.exit=normal.exit))
    } else
    {
      return(list(fit=fit.tmp,SISind=SISind, ISISind=ISISind))
    }
  }


  while(test){
    oldISISind=ISISind
    curloop=curloop+1
    remind <- setdiff(1:p, ISISind)

    tempresi=rep(0,length(remind))


    oldInd = rep(1:df0,length(ISISind))+df0*rep(ISISind-1,each=df0)

    for(i in 1:length(remind))
    {
      tempfit <-lm.fit (x=cbind(1,xbs[,c(oldInd,(remind[i]-1)*df0+1:df0)]),y=y)
      tempresi[i]<-sum(tempfit$residuals^2)
    }

    used.list<- tempresi

    used.sort <- sort(used.list, method= "sh", index=TRUE, decreasing=FALSE)


    mindex <- sample(1:n)
    mresi=NULL
    for(i in 1:length(remind))
    {
      tempfit <-lm.fit (x=cbind(1,xbs[,oldInd],xbs[mindex,(remind[i]-1)*df0+1:df0]),y=y)
      mresi[i] <- sum(tempfit$residuals^2)
    }
    resi.thres = quantile(mresi,1-quant)
    newind <- min(sum(used.list<resi.thres), floor(n/df0/2)-length(ISISind),gnum)

    if(newind==0)
    {
      test=0
      break
    }
    new.pickind <- sort(remind[used.sort$ix[1:newind]])

    pick.ind=c(ISISind, new.pickind)

    cat('loop ', curloop, '...SISind ', sort(pick.ind),'\n')
    xnew=x[,pick.ind]

    fit.tmp <- cv.penGAM(data=list(x=xnew,y=y), folds=folds, lambda.pen.list=NULL, knots=knots, eps0=eps0 ,cv.trace=trace)
    ISISind <- sort(pick.ind[fit.tmp$penGAMind])
    cat('loop ', curloop, '... ISISind ', ISISind,'\n')

    detail.pickind[[curloop]]=pick.ind
    detail.ISISind[[curloop]]=ISISind

    if(setequal(oldISISind, ISISind)) test=0
    if(curloop>=maxloop) {
      test=0
      normal.exit=0
    }

  }
  final.fit<-fit.tmp$fit
  cv.error<-fit.tmp$cv.error
  if(!is.null(testdata)){
    testx=testdata$x
    testy=testdata$y
    pred.error<-mean((predict(final.fit,as.matrix(testx[,pick.ind]))-testy)^2)
  } else pred.error=NULL


  ISISind=sort(ISISind)
  SISind=sort(SISind)
  ptime = proc.time()[1]-t0
  cat('finishing adaptINIS...\n')
  if(detailed){
    return(list(fit=fit.tmp, initRANKorder=initRANKorder, detail.pickind=detail.pickind, detail.ISISind=detail.ISISind,
                SISind=SISind, ISISind=ISISind,  nsis=nsis, cv.error=cv.error, pred.error=pred.error, ptime = ptime, normal.exit=normal.exit))
  } else
  {
    return(list(fit=fit.tmp, ISISind=ISISind, cv.error=cv.error, pred.error=pred.error, ptime = ptime, normal.exit=normal.exit))
  }
}

###########################################################################
#adaptINIS: INIS-penGAM algorithm, using knots=n^0.2
###########################################################################
adaptINIS <- function(data, testdata=NULL, lambda.pen.list=NULL, folds=NULL, quant=NULL, kfold=NULL, knots=NULL,  eps0=1e-6, DOISIS=TRUE, maxloop=10, trace=FALSE, detailed=FALSE){
  t0=proc.time()[1]
  cat('starting adaptINIS, INIS-penGAM algorithm, adatively choose number of variables\n')
  x=data$x
  y=data$y
  n <- nrow(x)
  p <- ncol(x)
  
  
  #if(is.null(nsis)) nsis=min(floor(n/log(n)),p-1)
  if(is.null(knots)) knots=ceiling(n^0.2)
  if(is.null(folds)) {
    temp= sample(1:n, n, replace = FALSE)
    if(is.null(kfold)) kfold=5
    for(i in 1:kfold){
      folds[[i]]=setdiff(1:n, temp[seq(i, n, kfold)])
    }
  }
  if(is.null(quant)) quant=1
  df0 <- knots+1
  
  xbs=matrix(0,n,df0*p)
  
  
  for(i in 1:p)
    xbs[,(i-1)*(df0)+(1:df0)]=ns(x[,i],df=df0)
  
  tempresi <- rep(0,p)
  
  curloop=1
  for(i in 1:p)
  {
    
    tempfit <-lm.fit (x=cbind(1,xbs[,(i-1)*df0+1:df0]),y=y)
    tempresi[i] <- sum(tempfit$residuals^2)
  }
  
  
  used.list<- tempresi
  
  used.sort <- sort(used.list, method= "sh", index=TRUE, decreasing=FALSE)
  initRANKorder <- used.sort$ix
  
  
  mindex <- sample(1:n)
  mresi=NULL
  for(i in 1:p){
    tempfit <-lm.fit (x=cbind(1,xbs[,(i-1)*df0+1:df0]),y=y[mindex])
    mresi[i] <- sum(tempfit$residuals^2)
  }
  resi.thres = quantile(mresi,1-quant)
  nsis <- max(min(sum(used.list<resi.thres), floor(n/df0/3)),2)
  
  
  SISind <- sort(initRANKorder[1:nsis])
  if(!DOISIS) return (list(initRANKorder=initRANKorder, SISind=SISind, nsis=nsis))
  
  cat('loop ', curloop, '...SISind ', SISind,'\n')
  pick.ind=initRANKorder[1:nsis]
  
  
  xnew=x[,pick.ind]
  
  fit.tmp <- cv.penGAM(data=list(x=xnew,y=y), folds=folds, lambda.pen.list=NULL, knots=knots, eps0=eps0,cv.trace=trace)
  ISISind <- sort(pick.ind[fit.tmp$penGAMind])
  cat('loop ', curloop, '...ISISind ', ISISind,'\n')
  test= 1
  normal.exit = 1
  
  detail.pickind=NULL
  detail.ISISind=NULL
  detail.pickind=as.list(detail.pickind)
  detail.ISISind=as.list(detail.ISISind)
  
  detail.pickind[[curloop]]=pick.ind
  detail.ISISind[[curloop]]=ISISind
  
  if(length(ISISind)==0) {
    warnings('No variable was selected by this variant of ISIS.')
    normal.exit = 0
    if(detailed){
      return(list(fit=fit.tmp,initRANKorder=initRANKorder, detail.pickind=detail.pickind, detail.ISISind=detail.ISISind,
                  SISind=SISind, ISISind=ISISind,  nsis=nsis, normal.exit=normal.exit))
    } else
    {
      return(list(fit=fit.tmp,SISind=SISind, ISISind=ISISind))
    }
  }
  
  
  while(test){
    oldISISind=ISISind
    curloop=curloop+1
    remind <- setdiff(1:p, ISISind)
    
    tempresi=rep(0,length(remind))
    
    
    oldInd = rep(1:df0,length(ISISind))+df0*rep(ISISind-1,each=df0)
    
    for(i in 1:length(remind))
    {
      tempfit <-lm.fit (x=cbind(1,xbs[,c(oldInd,(remind[i]-1)*df0+1:df0)]),y=y)
      tempresi[i]<-sum(tempfit$residuals^2)
    }
    
    used.list<- tempresi
    
    used.sort <- sort(used.list, method= "sh", index=TRUE, decreasing=FALSE)
    
    
    mindex <- sample(1:n)
    mresi=NULL
    for(i in 1:length(remind))
    {
      tempfit <-lm.fit (x=cbind(1,xbs[,oldInd],xbs[mindex,(remind[i]-1)*df0+1:df0]),y=y)
      mresi[i] <- sum(tempfit$residuals^2)
    }
    resi.thres = quantile(mresi,1-quant)
    newind <- min(sum(used.list<resi.thres), floor(n/df0/3)-length(ISISind))
    
    if(newind==0)
    {
      test=0
      break
    }
    new.pickind <- sort(remind[used.sort$ix[1:newind]])
    
    pick.ind=c(ISISind, new.pickind)
    
    cat('loop ', curloop, '...SISind ', sort(pick.ind),'\n')
    xnew=x[,pick.ind]
    
    fit.tmp <- cv.penGAM(data=list(x=xnew,y=y), folds=folds, lambda.pen.list=NULL, knots=knots, eps0=eps0 ,cv.trace=trace)
    ISISind <- sort(pick.ind[fit.tmp$penGAMind])
    cat('loop ', curloop, '... ISISind ', ISISind,'\n')
    
    detail.pickind[[curloop]]=pick.ind
    detail.ISISind[[curloop]]=ISISind
    for(testind in (curloop-1):1){
      if(setequal(detail.ISISind[[testind]], ISISind)){
        test=0
        break
      }
    }
    #if(setequal(oldISISind, ISISind)) test=0
    if(curloop>=maxloop) {
      test=0
      normal.exit=0
    }
    
  }
  final.fit<-fit.tmp$fit
  cv.error<-fit.tmp$cv.error
  if(!is.null(testdata)){
    testx=testdata$x
    testy=testdata$y
    pred.error<-mean((predict(final.fit,testx[,pick.ind])-testy)^2)
  } else pred.error=NULL
  
  
  ISISind=sort(ISISind)
  SISind=sort(SISind)
  ptime = proc.time()[1]-t0
  cat('finishing adaptINIS...\n')
  if(detailed){
    return(list(fit=fit.tmp, initRANKorder=initRANKorder, detail.pickind=detail.pickind, detail.ISISind=detail.ISISind,
                SISind=SISind, ISISind=ISISind,  nsis=nsis, cv.error=cv.error, pred.error=pred.error, ptime = ptime, normal.exit=normal.exit))
  } else
  {
    return(list(fit=fit.tmp, ISISind=ISISind, cv.error=cv.error, pred.error=pred.error, ptime = ptime, normal.exit=normal.exit))
  }
}


#########################################################################################################
compModelSize<-function(funind=NULL, n=NULL, p=NULL, repgroup=NULL,  r0=NULL, r1=NULL,  knots=NULL, maxloop=NULL, eps0=NULL, path=NULL)
{
  if(is.null(funind)) stop('No given data generating scheme.')

  if(is.null(eps0)) eps0=1e-6

  if(is.null(n)) n=400
  if(is.null(p)) p=1000
  if(is.null(knots)) knots=ceiling(n^0.2)
  if(is.null(r0)&is.null(r1)){
    if(is.null(repgroup))
    {
      r0=r1=1
    }
    else{
      r0=(repgroup-1)*10+1
      r1=(repgroup)*10
    }
  }

  if(is.null(path)) path='../'

  ########################################
  ####loading packages
  ########################################
  lib.loc=paste(path, 'packages', sep='')
  library(zoo, lib.loc=lib.loc)
  library(fda, lib.loc=lib.loc)
  library(grplasso, lib.loc=lib.loc)
  library(penGAM, lib.loc=lib.loc)
  library(SIS, lib.loc=lib.loc)

  library(mgcv)

  source('gen.datanew.r')

  newpath=paste('result',funind,sep='')
  if(!file.exists(newpath))
    dir.create(newpath)

  s0=gen.data(funind,n,p)$s
  truemodel<-1:s0

  sizelist=matrix(0,r1,3)
  outputfilename=paste('results/result', funind,'/',funind,'modelsize.txt',sep='')
  for(repind in r0:r1){
    set.seed(repind)
    cat('Random seed=', repind, '...\n')
    data=gen.data(funind,n,p)
    sizelist[repind,]=c(getSISSize(data,truemodel=truemodel),getPenGamSize(data,truemodel=truemodel),getLmSISSize(data,truemodel=truemodel))
    cat(repind,sizelist[repind,],'\n',sep=' ',file=outputfilename, append=T)
  }

  filename=paste('result', funind,'/',funind,'_r0=',r0,'_r1=',r1,'modelsize.RData',sep='')

  save(truemodel,sizelist,file=filename)
}
#######################################################################################################
#######################################################################################################
cv.penGAM<- function(data, testdata=NULL, knots=NULL, folds=NULL, kfold=NULL, lambda.pen.list=NULL, model=LinReg(), control = grpl.control(trace=0),
                     eps0=1e-6 ,cv.trace=FALSE, detailed=FALSE){
  t0=proc.time()[1]
  x=data$x
  y=data$y
  x=as.matrix(x)
  n=nrow(x)

  if(is.null(folds)) {
    temp= sample(1:n, n, replace = FALSE)
    if(is.null(kfold)) kfold=5
    for(i in 1:kfold){
      folds[[i]]=setdiff(1:n, temp[seq(i, n, kfold)])
    }
  }
  kfold=length(folds)
  if(is.null(knots)) knots=ceiling(n^0.2)

  if(is.null(lambda.pen.list)){
    lambda.pen.list=seq(1,0.02,-0.02)
  }
  error <- rep(0, length(lambda.pen.list))
  if(cv.trace) cat('Beginning Cross Validation.....\n')

  for(i in 1:kfold){
    if(cv.trace) cat('Cross Validation...', 'fold ', i,'...\n')
    train <- folds[[i]]
    test <- setdiff(1:n, folds[[i]])
    fit <- penGAM(as.matrix(x[train,]), y[train], lambda.pen=lambda.pen.list, lambda.curv=0,
                  knots=knots, model=LinReg(), control = grpl.control(trace=0))
    fit.pred <- predict(fit, as.matrix(x[test,]))
    fit.mat<-t(fit.pred[,1,])
    error <- error+apply((fit.mat-y[test])^2,2,mean)
  }

  if(cv.trace) cat('Ending Cross Validation.....\n')
  best.lambda.pen.ind <- (1:length(lambda.pen.list))[which.min(error)]


  fit <- penGAM(x, y, lambda.pen=lambda.pen.list[best.lambda.pen.ind], lambda.curv=0,
                knots=knots, model=LinReg(), control = grpl.control(trace=0))

  cv.error <- mean((predict(fit,x)-y)^2)

  coef.mat <- matrix(fit$coef[1,1,][-1],nrow=knots+2,ncol=ncol(x))

  coef.norm <- apply(coef.mat^2, 2, sum)

  non.zero.ind <- which(coef.norm>eps0)

  if(!is.null(testdata)){
    testx=as.matrix(testdata$x)
    testy=testdata$y
    pred.error<-mean((predict(fit,testx)-testy)^2)
  } else pred.error=NULL

  ptime = proc.time()[1]-t0
  if(detailed){
    return (list=list(fit=fit, best.lambda.pen.ind=best.lambda.pen.ind, best.lambda=lambda.pen.list[best.lambda.pen.ind],
                      penGAMind=non.zero.ind, cv.error=cv.error, pred.error=pred.error, penGAMfit=fit, ptime = ptime))
  }
  else
  {
    return (list=list(fit=fit, penGAMind=non.zero.ind, cv.error=cv.error, pred.error=pred.error, ptime = ptime))
  }

}
########################################################################################
getSISSize <- function(data, truemodel,  knots=NULL){
  x = data$x
  y = data$y
  n = nrow(x)
  p = ncol(x)

  if(is.null(knots)) knots=ceiling(n^0.2)
  df0 <- knots+1

  xbs=matrix(0,n,df0*p)


  for(i in 1:p)
    xbs[,(i-1)*(df0)+(1:df0)]=ns(x[,i],df=df0)

  tempresi <- rep(0,p)

  curloop=1
  for(i in 1:p)
  {

    tempfit <-lm.fit (x=cbind(1,xbs[,(i-1)*df0+1:df0]),y=y)
    tempresi[i] <- sum(tempfit$residuals^2)
  }


  used.list<- tempresi

  used.sort <- sort(used.list, method= "sh", index=TRUE, decreasing=FALSE)
  amsis.ind <- used.sort$ix


  for(i in 1:p)
  {
    if(length(setdiff(truemodel,amsis.ind[1:i]))==0) {
      amsis.size=i
      break
    }
  }

  return(amsis.size)

}
#######################################################################################
########################################################################################
getPenGamSize <- function(data,  truemodel, lambda.pen.list=NULL, knots=NULL, eps0=1e-6){
  x = data$x
  y = data$y
  n = nrow(x)
  p = ncol(x)
  if(is.null(knots)) knots=ceiling(n^0.2)
  if(is.null(lambda.pen.list)){
    lambda.max=0.99
    lambda.pen.list=lambda.max*(2^seq(0,-10,-0.1))
  }
  lam.len=length(lambda.pen.list)

  fit <- penGAM(x, y, lambda.pen=lambda.pen.list, lambda.curv=0,
                knots=knots, model=LinReg(), control = grpl.control(trace=0))

  coef.mat <- array(fit$coef[,1,-1],c(lam.len,knots+2,p))

  coef.norm <- apply(coef.mat^2, c(1,3), sum)


  nonzeros = coef.norm[,truemodel]>eps0;
  lambdas = apply(nonzeros, 1, all);
  ##str(fit.a);

  pickone = max(which(lambdas==F))+1;

  if (pickone==lam.len+1){
    mlasso <- p;
  }
  if (pickone<lam.len+1){
    mlasso <- sum(coef.norm[pickone,]>eps0);
  }
  return(mlasso)

}

########################################################################################
getLmSISSize <- function(data, truemodel){
  x = data$x
  y = data$y
  n = nrow(x)
  p = ncol(x)


  lmsis.corr <- NULL
  for(i in 1:p)
  {
    lmsis.corr[i] <- -abs(cor(x[,i],y))
  }
  lmsis.sort <- sort(lmsis.corr, method= "sh", index=TRUE, decreasing=FALSE)
  lmsis.ind <- lmsis.sort$ix

  for(i in 1:p)
  {
    if(length(setdiff(truemodel,lmsis.ind[1:i]))==0) {
      lmsis.size=i
      break
    }
  }
  return(lmsis.size)

}
########################################################################################
getAmSISSize <- function(data, truemodel, rank.method='corr'){
  x = data$x
  y = data$y
  n = nrow(x)
  p = ncol(x)
  q <- ceiling(n^0.2)
  amsis.corr <- NULL
  amsis.resi <- NULL
  amsis.norm <- NULL
  for(i in 1:p)
  {
    #if(i%%(round(p/10))==0) cat('First Step:', round(i*100/p), 'percent complete\n')
    amsis.fit <-gam (y~s(x[,i],k=q,fx=TRUE, bs="cr"))
    amsis.corr[i] <- mean(fitted(amsis.fit)*y)
    amsis.resi[i] <- amsis.fit$deviance
    amsis.norm[i] <- sum(fitted(amsis.fit)^2)
  }
  amsis.list<-switch(rank.method, corr=-amsis.corr, resi=amsis.resi, norm=-amsis.norm)
  amsis.sort <- sort(amsis.list, method= "sh", index=TRUE, decreasing=FALSE)
  amsis.ind <- amsis.sort$ix
  lmsis.corr <- NULL
  for(i in 1:p)
  {
    lmsis.corr[i] <- -abs(cor(x[,i],y))
  }
  lmsis.sort <- sort(lmsis.corr, method= "sh", index=TRUE, decreasing=FALSE)
  lmsis.ind <- lmsis.sort$ix

  for(i in 1:p)
  {
    if(length(setdiff(truemodel,amsis.ind[1:i]))==0) {
      amsis.size=i
      break
    }
  }
  for(i in 1:p)
  {
    if(length(setdiff(truemodel,lmsis.ind[1:i]))==0) {
      lmsis.size=i
      break
    }
  }
  return(size=list(amsis.size=amsis.size,lmsis.size=lmsis.size))

}
####################

####################################################
#######################################################################################
##################################################
"gen"<-function(n,p,funind,t=0){
  if(funind==4){
    w=matrix(runif(n*p),n,p)
    u=runif(n)
    v=runif(n)
  } else {
    w=truncnorm.gen(n,p)
    u=as.vector(truncnorm.gen(n,1))
    v=as.vector(truncnorm.gen(n,1))
  }
  x=matrix(0,n,p)
  x[,1:4]=(w[,1:4]+t*u)/(1+t)
  x[,5:p]=(w[,5:p]+t*v)/(1+t)
  epsi=rnorm(n)
  s=3
  y=switch(funind, -3+3*x[,1]+4*x[,2]-8*x[,3]+5*x[,4]+epsi, -7+8*x[,1]-3*x[,2]+10*x[,3]^3-6*x[,4]*(x[,4]-1)+epsi , -5+8*x[,1]^3+10*x[,2]*(1-x[,2])-10*x[,3]^5-8*x[,4]^2+epsi, -4+4*x[,1]+cos(2*pi*x[,2])-8*x[,3]^3+sqrt(x[,4]*(1-x[,4]))*sin(2*pi*(1+2^((9-4*s)/5))/(x[,4]+2^((9-4*s)/5)))+epsi)
  data=list(x=x,y=y)
}

############################################################################################################

"truncnorm.gen"<-function(n,p){
  norm.trunc<-rnorm(n*p)
  ind= abs(norm.trunc)>2
  test=sum(ind)
  count=0
  while(test>0){
    tmp=rnorm(test)
    test=sum(abs(tmp)>2)
    norm.trunc[ind]= tmp
    ind=abs(norm.trunc)>2
  }
  abs(matrix(norm.trunc,n,p))/2
}

