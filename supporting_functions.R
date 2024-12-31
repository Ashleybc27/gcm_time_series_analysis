###supporting functions

prdtrend <- function(file){ ##function to get the linear model 
  dt <- read.csv(file,header=T)
  print(tail(dt))
  dt$idx <- (1:nrow(dt))
  mod <- lm(dt$pr~dt$idx)
  return(mod)
}

tmpdtrend <- function(file){ ##function to get the linear model 
  dt <- read.csv(file,header=T)
  print(tail(dt))
  dt$idx <- (1:nrow(dt))
  mod <- lm(dt$tasmax~dt$idx)
  return(mod)
}

peak1 <- function(futspi,newname){ ##list of spi outputs
  boolout <- c()
  for (r in futspi$fitted){
    if (is.na(r)){
      print(paste0("is bad",r))
    }
    else{
      if (r <= -1){
        nd <- cbind(r,newname,1)
        boolout <- rbind(boolout,nd)
      }
      if (r >= -1){
        nt <- cbind(r,newname,0)
        boolout <- rbind(boolout,nt)
      }
    }
  }
  return(boolout)
}

droughtlen <- function(x){ ###pass it the db 
  x$r <- as.numeric(x$r)
  masteridx <- 0
  lnthdf <- c()
  tot <- (nrow(x))
  for(i in 1:nrow(x)){
    if (i <= masteridx){
      print('no')
    }
    else{
      if (x[i,3] == 1){
        print('yes')
        n <- i
        lnth <- 1
        repeat{
          n <- n+1
          lnth <- lnth+1
          print(lnth)
          #print(tstdat[n,3])
          if (x[n,3] == 0 | is.na(x[n,3])){
            print('breaking')
            stidx <- n - (lnth-1)
            stidx1 <- stidx -1
            stdidx2 <- stidx-2
            if (x[stidx1,1]-x[stidx,1] >= 0 & x[stdidx2,1]-x[stidx1,1] >=0){##both positive (downward curve)          
              lnth <- lnth + 2
            }
            if (x[stidx1,1]-x[stidx,1] >= 0 & x[stdidx2,1]-x[stidx1,1] <=0){## positive then negative (valley)
              lnth <- lnth + 1
            } 
            else{}
            ndidx <- n 
            ndisx1 <- n + 1
            ndidx2 <- n + 2
            if (ndidx2 > tot){
              print('end')
            }
            else{
              if (x[ndisx1,1]-x[ndidx,1] >= 0 & x[ndidx2,1]-x[ndisx1,1] >=0){##both positive (downward curve)          
                lnth <- lnth + 2
                print('add 2')
              }
              if (x[ndisx1,1]-x[ndidx,1] >= 0 & x[ndidx2,1]-x[ndisx1,1] <=0){## positive then negative (valley)
                lnth <- lnth + 1
                print('add 1')
              } 
              else{
                print('add 0')
              }
            }
            fnlnth <- cbind(stidx, lnth)
            print(fnlnth)
            masteridx <- ndidx
            print(masteridx)
            lnthdf <- rbind(lnthdf,fnlnth)
            break
          }
        }
      }
      else{}
    }
  }
  return(lnthdf)
}