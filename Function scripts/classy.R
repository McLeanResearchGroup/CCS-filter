#eventually, hopefully can use the classyfireR package to classify in bulk.
library(webchem)
library(pbapply)

#this function will add classifications of all lipids classified by lipid blast. 
#(PE, PC, etc. abbrevations)
.blastClassy <- function(data, classfile){
  lipclass <- read.csv(classfile, header = TRUE)
  for(i in seq_along(lipclass$Abbr)){
    abbr <- as.character(lipclass$Abbr)[i]
    data$Class <- ifelse((startsWith(data$Name, abbr)), 
                           as.character(lipclass$Class)[i], data$Class)
    data$Subclass <- ifelse((startsWith(data$Name, abbr)), 
                           as.character(lipclass$Subclass)[i], data$Subclass)
  }
  data$SuperClass <- "Lipids and lipid-like molecules"
  message(paste(length(which(data$Class != "NA")), "IDs classified.", 
                length(which(is.na(data$Class))), "remaining."))
  return(data)
}

#this function will try to add inchi keys to all ccsids without an existing classification.
#it works in step sizes of 50 in order to allow easy monitoring of progress, as well as easy
#termination if needed (all progress won't be lost if terminated). 
.addInchi <- function(data){
  chemnames <- unique(subset(data, is.na(Class))$Name)
  inchied <- as.data.frame(chemnames)
  #add InChi Key
  inchied$inchis <- NA
  n = 1
  
  #works in step sizes of 50.
  for(i in 1:ceiling(length(chemnames)/50)){
    m <- i*50
    inchied$inchis[n:m] <- pbsapply(chemnames[n:m], cts_convert, from = "Chemical Name", to = "InChIKey", 
                                    choices = 1, verbose = FALSE, USE.NAMES = FALSE)
    n <- m
  }

    #match inchis to original data.
  for(i in seq_along(inchied$chemnames)){
    data$inchi <- ifelse(data$Name == as.character(inchied$chemnames[i]), 
                         inchied$inchis[i], data$inchi)
  }
  
  message(paste(nrow(subset(data, !is.na(inchi)))), " out of ", nrow(subset(data, is.na(Class))), 
                     " remaining IDs have been assigned an InChI Key.")
 return(data) 
}


