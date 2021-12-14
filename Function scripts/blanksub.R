.blankSubtract <- function(IDdata, blankdata){
  IDdata$blank <- NA
  for(i in seq_along(IDdata$Compound.ID)){
    feat <- IDdata[i,]
    mz <- feat$m.z
    b <- 0
    #set tolerance range for mz matching
    tol <- (5/10^6)*mz
    mz.range <- c((mz-tol), (mz+tol))
    names(mz.range) <- c("low", "high")
    #find matching features
    mzmatch <- blankdata[c(which(blankdata$m.z>mz.range["low"] & blankdata$m.z<mz.range["high"])),]
    if(nrow(mzmatch) != 0){
     b <- 1
    }
    IDdata[i,]$blank <- b
  }
  message(paste(length(unique(subset(IDdata, blank == 1)$Compound)), "features removed."))
  IDdata <- IDdata[-which(IDdata$blank == 1),]
  rm(id, mz, tol, mz.range, mzmatch)
  IDdata
}
  