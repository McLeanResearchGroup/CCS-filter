.broadFilterbase <- function(set, model){
  for(j in 1:nrow(set)){
    #find the x values on the curve on either side of the feature
    lo <- model[which(model$`X curve` == max(subset(model, model$`X curve` < set[j,]$mz)$`X curve`)),]
    hi <- model[which(model$`X curve` == min(subset(model, model$`X curve` > set[j,]$mz)$`X curve`)),]
    range <- as.data.frame(rbind(lo, hi))
    #fit this to a line
    hifit <- lm(range$`PI hi` ~ range$`X curve`)
    lofit <- lm(range$`PI lo` ~ range$`X curve`)
    #predict the position on the PI of the feature's mass
    pihi <- hifit$coefficients[2]*set[j,]$mz + hifit$coefficients[1]
    pilo <- lofit$coefficients[2]*set[j,]$mz + lofit$coefficients[1]
    if(!is.na(pihi) & !is.na(pilo)){
      if(set[j,]$ccsz > pilo & set[j,]$ccsz < pihi){
        set[j,]$Conf <- 1
      } else {
        set[j,]$Conf <- -1
      }
    } else {
      set[j,]$Conf <- -1
    }
  }
  return(set)
}


.broadFilter <- function(ids, models){
  classes <- unique(ids$Class)
  subclasses <- unique(ids$Subclass)
  testmodels <- models[c(which(names(models) %in% classes))]
  subtestmodels <- models[c(which(names(models) %in% subclasses))]
  ids$Conf <- 0
  #class analysis
  for(i in 1:length(testmodels)){
    set <- subset(ids, Class == names(testmodels)[i])
    model <- testmodels[[i]]
    set <- .broadFilterbase(set, model)
    for(n in 1:nrow(set)){
      id <- set[n,]
      hit <- which(ids$Feature == id$Feature & ids$ID == id$ID & 
                     ids$Name == id$Name)
      #assign Confidence for class analysis
      ids[hit,]$Conf <- id$Conf
    }
  }
  #subclass analysis
  for(i in 1:length(subtestmodels)){
    set <- subset(ids, Subclass == names(subtestmodels)[i])
    model <- subtestmodels[[i]]
    set <- .broadFilterbase(set, model)
    for(n in 1:nrow(set)){
      id <- set[n,]
      hit <- which(ids$Feature == id$Feature & ids$ID == id$ID & 
                     ids$Name == id$Name)
      #assign Confidence for subclass analysis
      ids[hit,]$Conf <- id$Conf
    }
  }
  message(length(which(ids$Conf == -1)), " IDs rejected out of ", nrow(ids))
  return(ids)
}


.fineFilter <- function(ids, models, AICs){
  
  for(i in 1:length(unique(ids$Feature))){
    #grab all ids for one feature
    matches <- subset(ids, Feature == unique(ids$Feature)[i])
    #is there > 1 ID remaining from > 1 subclass (and therefore also class)?
    
    if(nrow(matches > 1) & length(unique(matches$Subclass)) > 1){
      #is there > 1 class WITH a model
      ccs <- matches[1,]$ccsz
      mz <- matches[1,]$mz
      
      if(length(which(names(models) %in% matches$Class)) > 1){
        #filter1- class
        #which models to test: 
        test <- models[c(which(names(models) %in% unique(matches$Class)))]
        pred <- as.data.frame(names(test))
        
        #predict the mean of the model at the mz of the feature
        pred$predccs <- lapply(test, predict, newdata = data.frame(x = mz))
        
        #subtract mean from actual ccs
        pred$delta <- lapply(pred$predccs, function(x, c) x - c, c = ccs)
        
        #pull the std err of the model
        pred$stderr <- AICs[which(AICs[, 1] %in% pred[, 1]), 7]
        
        #calculate distance in sigma of the ccs from the model mean
        pred$dist <- as.numeric(pred[, 3])/pred[, 4]
        likely <- pred$`names(test)`[which(abs(pred$dist) == min(abs(pred$dist)))]
        
        #if pass, conf = 2 
        matches$Conf <- ifelse(matches$Class == likely, 3, matches$Conf)
        matches <- subset(matches, Conf == 3)
        #if passing features still have more than one subclass with a model, 
        #put them back into the filter. 
        }
      #is there > 1 subclass WITH a model
      if(length(which(names(models) %in% matches$Subclass)) > 1){
        #filter2
        #which models to test: 
        test <- models[c(which(names(models) %in% unique(matches$Subclass)))]
        pred <- as.data.frame(names(test))
        
        #predict the mean of the model at the mz of the feature
        pred$predccs <- as.numeric(lapply(test, predict, newdata = data.frame(x = mz)))
        
        #subtract mean from actual ccs
        pred$delta <- as.numeric(lapply(pred$predccs, function(x, c) x - c, c = ccs))
        
        #pull the std err of the model
        pred$stderr <- AICs[which(AICs[, 1] %in% pred[, 1]), 7]
        
        #calculate distance in sigma of the ccs from the model mean
        pred$dist <- as.numeric(pred[, 3])/pred[, 4] 
        
        # #calculate pct difference of experimental and predicted CCS
        # pred$pctdiff <- as.numeric((abs(pred$delta)/pred$predccs)*100)
        # 
        # #what needs to happen: if all pctdiff values are < 2, don't make the call
        # #or.. maybe.. if all dist values are within 2 std dev of each other, 
        # #don't make the call?
        
        likely <- pred[which(abs(pred$dist) == min(abs(pred$dist))), 1]
        
        #if pass, conf = 2 
        matches$Conf <- ifelse(matches$Subclass == likely, matches$Conf + 2, matches$Conf)
        matches <- subset(matches, Conf > 1)
      }
      
      #assign confidence scores to original dataframe
      for(j in 1:nrow(matches)){
        m <- matches[j,]
        ids[which(ids$Feature == m$Feature & ids$Name == m$Name),]$Conf <- m$Conf
      }
    } else{
      #if the feature already has ids from only one subclass, increase its confidence
      #this ensures these highly confident ids do not get penalized
      ids$Conf <- ifelse(ids$Feature == matches$Feature[1], 2, ids$Conf)
    }
  }
  rm(ccs, mz, pred, matches, test)
  message(length(which(ids$Conf > 2)), " IDs passed the filter.")
  return(ids)
}

