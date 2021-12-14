#Model Generator
#Adapted from Compendium Rmd file

library(basicTrendline)

#import data
all <- read.csv(here("Function scripts/ccs-compendium.csv"))

all$Subclass <- as.factor(all$Subclass)
classes <- levels(all$Subclass)
fulldata = list()
for(i in seq_along(classes)) {
  class <- classes[i]
  dat <- as.data.frame(subset(all, Subclass == class))
  fulldata[[i]] <- dat
}
names(fulldata) <- classes
names(fulldata)[[1]] <- "none"

#initialize functions

.yTOz <- function(y){
  bottom <- as.numeric(quantile(y, .025, na.rm=TRUE))
  top <- as.numeric(quantile(y, .975, na.rm=TRUE))
  z <- (y - bottom)/(top - bottom)
  z[z<=0] <- 0.05; z[z>=1] <- 0.95
  z
}

.lmlz <- function(x, y){
  z <- .yTOz(y)
  z <- as.vector(by(z, x, mean, na.rm = TRUE))
  lz <- log(z/(1-z))
  lm(lz ~ unique(x))
}

.estimScal <- function(x, y){
  model <- .lmlz(x, y)
  return(coef(model)[2])
}
.estimMid <- function(x, y){
  model <- .lmlz(x, y)
  fit <- model$fitted.values
  predict(model, data.frame(x = median(fit, na.rm = TRUE)))  
}

.initPars <- function(x, y, npars){
  if(npars<5) s <- 1 else s <- 1
  if(npars<4) bottom <- 0 else bottom <- quantile(y, .05, na.rm=TRUE)
  if(npars<3) top <-1 else top <- quantile(y, .95, na.rm=TRUE)
  xmid <- (max(x, na.rm = TRUE) + min(x, na.rm = TRUE))/2
  scal <- .estimScal(x, y)
  return(c(bottom, top, xmid, scal, s))
}

CI <- function(stdErr, yobs, newy){
  n <- length(yobs)
  ybar <- mean(yobs, na.rm = TRUE)
  t <- qt(.99, n-2)
  ci <- t*stdErr*sqrt((1/n+(newy - ybar)^2/sum((newy - ybar)^2)))
  lo <- newy - ci
  hi <- newy + ci
  return(list(lo = lo, hi = hi))
}

PI <- function(stdErr, yobs, newy){
  n <- length(yobs)
  ybar <- mean(yobs, na.rm = TRUE)
  t <- qt(.99, n-2)
  pi <- t*stdErr*sqrt((1+1/n+(newy - ybar)^2/sum((newy - ybar)^2)))
  lo <- newy - pi
  hi <- newy + pi
  return(list(lo = lo, hi = hi))
}

#generate fits

#throw out classes where n < 10
data <- list()
for(i in seq_along(fulldata)){
  if(nrow(fulldata[[i]]) > 9){
    name <- names(fulldata)[[i]]
    data[[name]] <- fulldata[[i]]
  }
}
data[['none']] <- NULL

pf2 <- list()
for(i in seq_along(data)){
  tryCatch( {
    name <- names(data)[[i]]
    x <- data[[i]][, "mz"]
    y <- data[[i]][, "CCS.z"]
    xy <- data.frame(cbind(x, y))
    fit <- nls(y ~ SSpower2P(x, a, b))
    pf2[[name]] <- fit 
  }, error = function(e){cat("Error: ", names(data)[[i]], conditionMessage(e), "\n")})
} 

pf3 <- list()
for(i in seq_along(data)){
  tryCatch( {
    name <- names(data)[[i]]
    x <- data[[i]][, "mz"]
    y <- data[[i]][, "CCS.z"]
    xy <- data.frame(cbind(x, y))
    fit <- nls(y ~ SSpower3P(x, a, b, c))
    pf3[[name]] <- fit 
  }, error = function(e){cat("Error: ", names(data)[[i]], conditionMessage(e), "\n")})
} 


sig.4P.nls <- list()
for(i in seq_along(data)){
  tryCatch( {
    name <- names(data)[[i]]
    x <- data[[i]][, "mz"]
    y <- data[[i]][, "CCS.z"]
    xy <- as.data.frame(cbind(x, y))
    startvals <- .initPars(x, y, 4) 
    start <- startvals[1:4]
    names(start) <- c("b", "top", "xmid", "scal")
    fit <- nls(y ~ b+(top-b)/(1+10^((xmid-x)*scal)), data = xy, 
               start = start)
    sig.4P.nls[[name]] <- fit 
  }, error = function(e){cat("Error: ", names(data)[[i]], conditionMessage(e), "\n")})
}

sig.5P.nls <- list()
for(i in seq_along(data)){
  tryCatch( {
    name <- names(data)[[i]]
    x <- data[[i]][, "mz"]
    y <- data[[i]][, "CCS.z"]
    xy <- as.data.frame(cbind(x, y))
    startvals <- .initPars(x, y, 5) 
    start <- startvals[1:5]
    names(start) <- c("b", "top", "xmid", "scal", "s")
    fit <- nls(y ~ b+(top-b)/(1+10^((xmid-x)*scal))^s, data = xy, 
               start = start)
    sig.5P.nls[[name]] <- fit 
  }, error = function(e){cat("Error: ", names(data)[[i]], conditionMessage(e), "\n")})
}

#choose a model

subAICs <- data.frame(names(data))

subAICs$PF2 <- NA
for(i in seq_along(pf2)){
  AIC <- AIC(pf2[[i]])
  c <- (2*3^2 + 2*3)/(length(residuals(pf2[[i]]))-3-1)
  AICc <- AIC + c
  name <- names(pf2)[[i]]
  ind <- which(subAICs[, 1] == name)
  subAICs[ind, 2] <- AICc
}

subAICs$PF3 <- NA
for(i in seq_along(pf3)){
  AIC <- AIC(pf3[[i]])
  c <- (2*3^2 + 2*3)/(length(residuals(pf3[[i]]))-3-1)
  AICc <- AIC + c
  name <- names(pf3)[[i]]
  ind <- which(subAICs[, 1] == name)
  subAICs[ind, 3] <- AICc
}

subAICs$Sig4 <- NA
for(i in seq_along(sig.4P.nls)){
  AIC <- AIC(sig.4P.nls[[i]])
  c <- (2*3^2 + 2*3)/(length(residuals(sig.4P.nls[[i]]))-3-1)
  AICc <- AIC + c
  name <- names(sig.4P.nls)[[i]]
  ind <- which(subAICs[, 1] == name)
  subAICs[ind, 4] <- AICc
}

subAICs$Sig5 <- NA
for(i in seq_along(sig.5P.nls)){
  AIC <- AIC(sig.5P.nls[[i]])
  c <- (2*3^2 + 2*3)/(length(residuals(sig.5P.nls[[i]]))-3-1)
  AICc <- AIC + c
  name <- names(sig.5P.nls)[[i]]
  ind <- which(subAICs[, 1] == name)
  subAICs[ind, 5] <- AICc
}

for(i in seq_along(subAICs[, 1])){
  if(!(is.na(subAICs[i, 2]) & is.na(subAICs[i, 3]) & is.na(subAICs[i, 4]))){
    p2 <- subAICs[i, 2]
    p3 <- subAICs[i, 3]
    r <- subAICs[i, 4]
    v <- subAICs[i, 5]
    f <- c(p2 = p2, p3 = p3, r = r, v = v)
    fit <- names(which(f == min(f, na.rm = TRUE)))
    subAICs[i, 6] <- fit
  }
}

subfitteddata <- list()
for(i in seq_along(data)){
  if(!is.na(subAICs[i, 6])) {
    name = names(data)[[i]]
    subfitteddata[[name]] <- data[[i]]
    subfitteddata[[name]] <- subfitteddata[[name]][, c("Compound", "mz", "CCS.z", "RSD")]
  } else{}
}


#add to classes final list

subAICs$stderr <- NA

for(i in seq_along(subfitteddata)){
  name <- names(subfitteddata)[[i]]
  fit <- subAICs[which(subAICs[, 1] == name), 6]
  y <- subfitteddata[[i]][, 3]
  if(fit == 'p2'){
    models[[name]] <- pf2[[name]]
    yfit <- predict(pf2[[name]], newdata = data.frame(x = subfitteddata[[i]][, 2]))
  }
  if(fit == 'p3'){
    models[[name]] <- pf3[[name]]
    yfit <- predict(pf3[[name]], newdata = data.frame(x = subfitteddata[[i]][, 2]))
  }
  if(fit == 'r'){
    models[[name]] <- sig.4P.nls[[name]]
    yfit <- predict(sig.4P.nls[[name]], newdata = data.frame(x = subfitteddata[[i]][, 2]))
  } 
  if(fit == 'v'){
    models[[name]] <- sig.5P.nls[[name]]
    yfit <- predict(sig.5P.nls[[name]], newdata = data.frame(x = subfitteddata[[i]][, 2]))
  } 
  stdErr <- sqrt(1/(length(y)-2)*sum((yfit - y)^2))
  subAICs[which(subAICs[, 1] == name), ]$stderr <- stdErr 
}

subfitteddata <- list()
for(i in seq_along(data)){
  if(!is.na(subAICs[i, 6])) {
    name = names(data)[[i]]
    subfitteddata[[name]] <- data[[i]]
    subfitteddata[[name]] <- subfitteddata[[name]][, c("Compound", "mz", "CCS.z", "RSD")]
  } else{}
}

rm(data, fulldata, all)

#Generate tables of parameters
power2table <- data.frame()
power3table <- data.frame()
fourtable <- data.frame()
fivetable <- data.frame()

for(i in seq_along(subfitteddata)){
  name <- names(subfitteddata)[[i]]
  fit <- subAICs[which(subAICs[, 1] == name), 6]
  if(fit == 'p2'){
    pars <- coef(pf2[[name]])
    row <- as.data.frame(t(c(name, pars)))
    power2table <- rbind(power2table, row)
  }
  if(fit == 'p3'){
    pars <- coef(pf3[[name]])
    row <- as.data.frame(t(c(name, pars)))
    power3table <- rbind(power3table, row)
  }
  if(fit == 'r'){
    pars <- coef(sig.4P.nls[[name]])
    row <- as.data.frame(t(c(name, pars)))
    fourtable <- rbind(fourtable, row)
  }
  if(fit == 'v'){
    pars <- coef(sig.5P.nls[[name]])
    row <- as.data.frame(t(c(name, pars)))
    fivetable <- rbind(fivetable, row)
  }
}

colnames(power2table) <- c('class', 'a', 'k')
colnames(power3table) <- c('class', 'a', 'k', 'y0')
colnames(fourtable) <- c('class', 'y0', 'ymax', 'y50', 'H')
tryCatch({colnames(fivetable) <- c('class', 'y0', 'ymax', 'y50', 'H', 's')}, 
         error = function(e){cat("no 5P fits today")})

#Generate curves
for(i in seq_along(subfitteddata)){
  x <- subfitteddata[[i]]$mz
  name <- names(subfitteddata)[[i]]
  fit <- subAICs[which(subAICs[, 1] == name), 6]
  subfitteddata[[i]]$Xcurve <- seq(from = min(x), to = max(x), length.out = length(x))
  if(fit == 'p2'){
    subfitteddata[[i]]$Ycurve <- predict(pf2[[name]], newdata = data.frame(x = subfitteddata[[i]]$Xcurve))
  }
  if(fit == 'p3'){
    subfitteddata[[i]]$Ycurve <- predict(pf3[[name]], newdata = data.frame(x = subfitteddata[[i]]$Xcurve))
  }
  if(fit == 'r'){
    subfitteddata[[i]]$Ycurve <- predict(sig.4P.nls[[name]], newdata = data.frame(x = subfitteddata[[i]]$Xcurve))  
  }
  if(fit == 'v'){
    subfitteddata[[i]]$Ycurve <- predict(sig.5P.nls[[name]], newdata = data.frame(x = subfitteddata[[i]]$Xcurve))
  }
}

#Add CIs and PIs
for(i in seq_along(subfitteddata)){
  name <- names(subfitteddata)[[i]]
  fit <- subAICs[which(subAICs[, 1] == name), 6]
  y <- subfitteddata[[i]][, 3]
  if(fit == 'p2'){
    yfit <- predict(pf2[[name]], newdata = data.frame(x = subfitteddata[[i]][, 2]))
  }
  if(fit == 'p3'){
    yfit <- predict(pf3[[name]], newdata = data.frame(x = subfitteddata[[i]][, 2]))
  }
  if(fit == 'r'){
    yfit <- predict(sig.4P.nls[[name]], newdata = data.frame(x = subfitteddata[[i]][, 2]))
  }
  if(fit == 'v'){
    yfit <- predict(sig.5P.nls[[name]], newdata = data.frame(x = subfitteddata[[i]][, 2]))
  }
  stdErr <- sqrt(1/(length(y)-2)*sum((yfit - y)^2))
  stdDev <- sd(subfitteddata[[i]][, "CCS.z"])
  ci <- CI(stdErr, y, subfitteddata[[i]][, "Ycurve"])
  pi <- PI(stdErr, y, subfitteddata[[i]][, "Ycurve"])
  subfitteddata[[i]] <- cbind(subfitteddata[[i]], ci, pi, stdErr, stdDev)
  colnames(subfitteddata[[i]]) <- c('Compound', 'm/z', 'CCS/z', 'rsd', 'X curve', 'Y curve', 
                                 'CI lo', 'CI hi', 'PI lo', 'PI hi', 'std error', 'SD')
}

rm(pf2, pf3, sig.4P.nls, sig.5P.nls)

subAICs <- subset(subAICs, !is.na(subAICs$V6))

rm(AIC, AICc, c, ci, CI, class, classes, dat, f, fit, fivetable, fourtable, 
   i, ind, name, p2, p3, pars, pi, PI, power2table, power3table, r, row, start, 
   startvals, stdDev, stdErr, v, x, xy, y, yfit)

