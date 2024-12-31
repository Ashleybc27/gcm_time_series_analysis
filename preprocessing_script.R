##preprocessing

##take raw netcdfs and aggregate to a single time series 

##libraries 
library(exactextractr)
library(sf)
library(ncdf4)
library(raster)
library(tidyverse)
##location of rawfiles 
setwd("/Users/ashley/Documents/projects/sagehen_rh/rawgcms/")

##area to clip 
basin <- st_read("~/Documents/projects/sagehen_rh/sage_shp/basin_area.shp")

##function to turn each day of the netcdf array 
#into a raster and then extract area weighted daily means

####NOTE#### 
##this function flips data based on MACA setup 
##gridMET needs to be flipped differently so always check orientation of data 
dayfunc<-function(d){
  rv<-raster(apply(d,1,rev),crs=crs(meanprRas),
             xmn=corners[1], xmx=corners[2], ymn=corners[3], ymx=corners[4])
  #rv2 <- flip(rv,direction='y')
  zs<-exact_extract(rv,basin,'mean') ##"fast" zonal statistics
  return(zs) 
}

##file names 
files <- list.files(pattern = "agg_macav2metdata_pr_*")

##processing function
##NOTE this is very slow 
aggdatpr <- function(files){
  dates <- seq(as.Date("2006/1/1"), as.Date("2099/12/31"), by = "day")
  nc.data.one <- nc_open(files[1])
  pr_array <- ncvar_get(nc.data.one, "precipitation")
  fillvalue <- ncatt_get(nc.data.one,"precipitation", "_FillValue")
  pr_array[pr_array==fillvalue$value] <- NA
  lon <- ncvar_get(nc.data.one,'lon')
  lat <- ncvar_get(nc.data.one,'lat')
  meanpr_array<-apply(pr_array,c(1,2),mean)
  meanprRas<-raster(meanpr_array)
  #plot(meanprRas)
  #str(meanprRas)
  meanprRas<- raster(apply(meanpr_array,1,rev))
  crs(meanprRas) = "+proj=longlat +datum=WGS84 +no_defs"
  corners<-c(min(lon)-360-0.0417/2, max(lon)-360+0.0417/2, min(lat)-0.0417/2, max(lat)+0.0417/2)
  extent(meanprRas) = extent(corners)
  nc_close(nc.data)
  for (k in files){
    name <- strsplit(k,"_")
    nc.data <- nc_open(k)
    pr_array <- ncvar_get(nc.data, "precipitation")
    fillvalue <- ncatt_get(nc.data,"precipitation", "_FillValue")
    pr_array[pr_array==fillvalue$value] <- NA
    nc_close(nc.data)
    pr<-apply(pr_array,3,dayfunc)
    zsdf<-as.data.frame(pr)
    zsdf$date<-dates
    zsdf$year<-year(zsdf$date)
    zsdf$month <- month(zsdf$date)
    annualpr<-zsdf%>% group_by(year) %>% summarize(pr=sum(pr))
    monthpr<-zsdf%>% group_by(year,month) %>% summarize(pr=sum(pr))
    head(zsdf)
    head(annualpr)
    write.csv(zsdf, paste0("../newprocessed/",name[[1]][4],"_dailypr.csv"))
    write.csv(annualpr, paste0("../newprocessed/",name[[1]][4],"_annualpr.csv"))
    write.csv(monthpr, paste0("../newprocessed/",name[[1]][4],"_monthlypr.csv"))
  }
}

test <- aggdatpr(files)



####extra 

nc.data <- nc_open("agg_macav2metdata_pr_BNU-ESM_r1i1p1_rcp85_2006_2099_CONUS_daily.nc")
strsplit("agg_macav2metdata_pr_BNU-ESM_r1i1p1_rcp85_2006_2099_CONUS_daily.nc", "_")
print(nc.data)

dates <- seq(as.Date("2006/1/1"), as.Date("2099/12/31"), by = "day")
pr_array <- ncvar_get(nc.data, "precipitation")
fillvalue <- ncatt_get(nc.data,"precipitation", "_FillValue")
pr_array[pr_array==fillvalue$value] <- NA
lon <- ncvar_get(nc.data,'lon')
lat <- ncvar_get(nc.data,'lat')
nc_close(nc.data)

meanpr_array<-apply(pr_array,c(1,2),mean)
meanprRas<-raster(meanpr_array)
plot(meanprRas)
str(meanprRas)
meanprRas<- raster(apply(meanpr_array,1,rev))
crs(meanprRas) = "+proj=longlat +datum=WGS84 +no_defs"
corners<-c(min(lon)-360-0.0417/2, max(lon)-360+0.0417/2, min(lat)-0.0417/2, max(lat)+0.0417/2)
extent(meanprRas) = extent(corners)

plot(meanprRas)
#writeRaster(meanTmaxRas,"AveTmax_2018_2046_HADGEMES_RCP85.tif", overwrite=TRUE)
pr<-apply(pr_array,3,dayfunc)
zsdf<-as.data.frame(pr)
zsdf$date<-dates
zsdf$year<-year(zsdf$date)
zsdf$month <- month(zsdf$date)
annualpr<-zsdf%>% group_by(year) %>% summarize(pr=sum(pr))
monthpr<-zsdf%>% group_by(year,month) %>% summarize(pr=sum(pr))
head(zsdf)
head(annualpr)
write.csv(zsdf, "BNU-ESM_dailypr.csv")
write.csv(annualpr, "BNU-ESM_mean_annualpr.csv")
write.csv(monthpr, "BNU-ESM_mean_monthlypr.csv")

##function to put netcdf into the right orientation and average cells to single value 
dayfunc<-function(d){
  rv<-raster(apply(d,1,rev),crs=crs(meantasmaxRas),
             xmn=corners[1], xmx=corners[2], ymn=corners[3], ymx=corners[4])
  #rv2 <- flip(rv,direction='y')
  zs<-exact_extract(rv,basin,'mean') ##"fast" zonal statistics
  return(zs) 
}

