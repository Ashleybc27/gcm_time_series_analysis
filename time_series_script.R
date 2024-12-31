### final time series analysis script 

#libraries 
library(tidyverse)
library(WaveletComp)
library(SPEI)

#processed file location
setwd("~/Documents/projects/plumas.gcms/")
annualfiles <- list.files(pattern="_annual_pr.csv|_annual_tmax.csv")
#dayfiles <- list.files(pattern="dailypr|_mean_dailytasmax")
monthfiles <- list.files(pattern="_mean_monthlytasmax|_mean_monthlypr")
#monthfiles <- append(monthfiles,list.files(path="monthly/",pattern="monthlypr"))
allcharacteristics <- c() ##adding results to this 

###trend and detrend 
trendanalysis <- function(files){ ##load the prdtrend and tmpdrend functions
  prtrend <- c()
  tmptrend <-c()
  for (file in files){ ##apply lm and output residuals for analysis
    name <- file
    splits <- strsplit(name,c('_'))
    newname <- splits[[1]][1]
    print(newname)
    
    if (splits[[1]][3] == "monthlypr.csv"){
    mod <- prdtrend(name)
    residal <- as.data.frame(mod$residuals)
    residal$month <- rep(seq(1:12),times=94)
    residal$year <- rep(2006:2099, each=12)
    #residal$date <- seq(as.Date("2006/01/01"), length.out=nrow(residal), by = "month")
    #residal$year <- (seq(from = 2006,length.out=nrow(residal)))
    write.csv(residal, file=paste0("detrended/",newname,"_detrendmonthlypr.csv"))
    toadd <- cbind(mod$coefficients[2],newname) ##this returns the trend
    prtrend <- rbind(prtrend,toadd)}
    
    else if (splits[[1]][3] == "monthlytasmax.csv"){
    mod <- tmpdtrend(name)
    residal <- as.data.frame(mod$residuals)
    residal$month <- rep(seq(1:12),times=94)
    residal$year <- rep(2006:2099, each=12)
    #residal$date <- seq(as.Date("2006/01/01"), length.out=nrow(residal), by = "month")
    #residal$year <- (seq(from = 2006,length.out=nrow(residal)))
    write.csv(residal, file=paste0("detrended/",newname,"_detrendmonthlymaxtmp.csv"))
    toadd <- cbind(mod$coefficients[2],newname) ##this returns the trend
    tmptrend <- rbind(tmptrend,toadd)}
  }
  colnames(prtrend) <- c("prtrend","model")
  colnames(tmptrend) <- c("tmptrend","model")
  alltrend <- merge(prtrend,tmptrend, by="model")
  return(alltrend)
}

test <- c()
test <- trendanalysis(monthfiles)

allcharacteristics <- test


###variance 
setwd("~/Documents/projects/plumas.gcms/detrended/")
annualde <- list.files(pattern="detrendannualypr|detrendannualmaxtmp")
##interannual 
findintervar <- function(x){
  intervariancepr <- c()
  intervariancetmp <- c()
  for (j in x){
    dat <- read.csv(j,header=T)
    print(head(dat))
    name <- strsplit(j,c('_'))
    if (name[[1]][2] == "detrendannualypr.csv"){
    vr <- mad(dat$mod.residuals,center=mean(dat$mod.residuals))
    datadd <- cbind(vr,name[[1]][1])
    intervariancepr <- rbind(intervariancepr,datadd)}
    
    else if (name[[1]][2] == "detrendannualmaxtmp.csv"){
    vr <- mad(dat$mod.residuals,center=mean(dat$mod.residuals))
    datadd <- cbind(vr,name[[1]][1])
    intervariancetmp <- rbind(intervariancetmp,datadd)}
  }
  colnames(intervariancepr) <- c('interannualpr','model')
  colnames(intervariancetmp) <- c('interannualmaxtmp','model')
  interann <- merge(intervariancepr,intervariancetmp,by="model")
  return(interann)
}
intervarout <-c()
intervarout <- findintervar(annualde)

allcharacteristicsbak <- allcharacteristics
allcharacteristics <- merge(allcharacteristics, intervarout, by="model")

###intra-annual variance 
monthlyde <- list.files(pattern="detrendmonthlypr|detrendmonthlymaxtmp")

findintravar <- function(x){
  intravariancepr <- c()
  intravariancetmp <- c()
  for (k in x){
    splits <- strsplit(k,c('_')) #get namebase
    newname <- splits[[1]][1]
    dat <-as.data.frame(read.csv(k,header=T))
    print(k)
    print(head(dat))
    print(nrow(dat))
    #dat$year <- year(dat$date)
    vr <- dat %>% group_by(year) %>% summarise(anvar = mad(mod.residuals,center=mean(mod.residuals),na.rm=T))
    vr$model <- paste0(newname)
    mn <- mean(vr$anvar, na.rm=T)
    mod1 <- lm(vr$anvar~vr$year)
    slp <- mod1$coefficients[2]
    toadd <- cbind(newname,mn,slp)
    if (splits[[1]][2] == "detrendmonthlypr.csv"){
      intravariancepr <- rbind(intravariancepr,toadd)
    }
    else if (splits[[1]][2] == "detrendmonthlymaxtmp.csv"){
      intravariancetmp <- rbind(intravariancetmp,toadd)
    }
  }
  colnames(intravariancepr) <- c('model','intraannualpr','intraannualprtrend')
  colnames(intravariancetmp) <- c('model','intraannualmaxtmp','intraannualmaxtmptrend')
  intraann <- merge(intravariancepr,intravariancetmp,by="model")
  return(intraann)
}
intraannout <- c()
intraannout <- findintravar(monthlyde)
write.csv(intraannout,"../intrafix.csv")
allcharacteristicsbak <- allcharacteristics
allcharacteristics <- merge(allcharacteristics, intraannout, by="model")

##Periodicity
setwd("/Volumes/HananLab/Ash_rawdata/regionalMACA/")

periodicityfunc <- function(annualfiles){
  tempper <- c()
  prper <- c()
  for (l in annualfiles){
    splits <- strsplit(l,c('_')) #get namebase
    newname <- splits[[1]][1]
    #print(newname)
    dat <- read.csv(l, header=T)
    print(head(dat))
    if (splits[[1]][3] == "pr.csv"){
      result <- analyze.wavelet(dat,my.series = 'pr')
      out <- as.data.frame(cbind(result$Period,result$Power.avg.pval))
      print(max(out$V2))
      out$tc <- trunc(out$V1)
      ag <- out %>% group_by(tc) %>% summarise(pval=max(V2))
      toadd <- cbind(min(ag$tc[ag$pval>=.95]),max(ag$tc[ag$pval>=.95]),newname)
      prper <- rbind(prper,toadd)
    }
    else if (splits[[1]][3] == "tmax.csv"){
      result <- analyze.wavelet(dat,my.series = 'tmax')
      out <- as.data.frame(cbind(result$Period,result$Power.avg.pval))
      print(max(out$V2))
      out$tc <- trunc(out$V1)
      ag <- out %>% group_by(tc) %>% summarise(pval=max(V2))
      toadd <- cbind(min(ag$tc[ag$pval>=.95]),max(ag$tc[ag$pval>=.95]),newname)
      tempper <- rbind(tempper,toadd)
    }
  }
  colnames(prper) <- c('minprper','maxprper','model')
  colnames(tempper) <- c('mintempper','maxtempper','model')
  periodicity<- merge(prper,tempper,by="model")
  return(periodicity)
}
periodanaly <- c()
periodanaly <- periodicityfunc(annualfiles)
periodanaly <- as.data.frame(periodanaly)

allcharacteristicsbak <- allcharacteristics
allcharacteristics <- merge(allcharacteristics, periodanaly, by="model")

###coherency 
prfiles <- list.files(pattern = "_annual_pr.csv")

coherfunc <- function(prfiles){
  coher <- c()
  for (g in prfiles){
    splits <- strsplit(g,c('_')) #get namebase
    newname <- splits[[1]][1]
    tmpname <- paste0(newname,"_annual_tmax.csv")
    print(newname)
    pr <- read.csv(g, header=T)
    tmp <- read.csv(tmpname,header=T)
    coherdf <- as.data.frame(cbind(pr$year,pr$pr,tmp$tmax))#coherdf <- cbind(prdf, tasmindf$mod.residuals)
    colnames(coherdf) <- c('year','pr','tasmax')
    print(head(coherdf))
    result <- analyze.coherency(coherdf,my.pair=c('pr','tasmax'))
    out <- as.data.frame(cbind(result$Period,result$Coherence.avg.pval))
    print(head(out))
    out$tc <- trunc(out$V1)
    ag <- out %>% group_by(tc) %>% summarise(pval=max(V2))
    toadd <- cbind(min(ag$tc[ag$pval>=.95]),max(ag$tc[ag$pval>=.95]),newname)
    coher <- rbind(coher,toadd)
  }
  colnames(coher) <- c("mincoher","maxcoher","model")
  return(coher)
}
co <- c()
co <- coherfunc(prfiles)
co <- as.data.frame(co)
allcharacteristicsbak <- allcharacteristics
allcharacteristics <- merge(allcharacteristics, co, by="model")

##SPI 
monthpr <- list.files(pattern="_monthly_pr.csv")
historic <- "~/Documents/projects/all_gridmet/sn_hist_monthpr.csv"
###make sure to load the peak and len functions below 

spifunc <- function(historic,monthpr){
  hist <- read.csv(historic,header = T)
  histts <- ts(hist, start=c(1979,1),frequency=12)
  droughtout <- c()
  for (h in monthpr){
    splits <- strsplit(h,c('_'))
    futdat <- read.csv(h,header=T)
    if (nrow(futdat) == 1128){
      futclip <- subset(futdat, futdat[,2] >=2017)
      futclip2 <- subset(futdat,futdat[,2] == 2016 & futdat[,3] >=10)
      futneed <- rbind(futclip2, futclip)
      futal <- rbind(histts,futneed)
      futallts <- ts(futal, start=c(1979,10), frequency=12)
      futspi <- spi(ts(futallts[,4], start=c(1979,1), frequency=12),12, ref.start=c(1979,10),ref.end=c(2021,12),na.rm=T)
      peaks <- as.data.frame(peak1(futspi,splits[[1]][1]))
      print(head(peaks))
      droughts <- as.data.frame(droughtlen(peaks))
      toadd <- cbind(mean(droughts$lnth),max(droughts$lnth),splits[[1]][1])
      droughtout <- rbind(droughtout,toadd)
    }
    else{
      print(paste0('bad model is',splits))
    }
  }
  colnames(droughtout) <- c("averagedrought","maxdrought","model")
  return(droughtout)
}
drought <- c()
droughts <- spifunc(historic,monthpr)
droughts <- as.data.frame(droughts)
allcharacteristicsbak <- allcharacteristics
allcharacteristics <- merge(allcharacteristics, droughts, by="model")

write.csv(allcharacteristics,"finaltimeseries_sn.csv")









