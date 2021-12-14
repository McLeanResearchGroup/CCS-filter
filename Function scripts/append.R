#append CCS data- returns dataset with appended CCS
.appendCCS <- function(IDdata, CCSdata){
  #initialize ccs column
  IDdata$CCS <- as.numeric(NA)

  #de-duplicate features
  features <- IDdata %>% 
    distinct(Compound, .keep_all = TRUE)
  
  #enter the loop
  for(i in seq_along(features$Compound)){
    #initialize CCS
    ccs <- NA
    #grab feature name
    feat <- features[i,]$Compound
    #grab mz
    mz <- as.numeric(features[i,]$m.z)
    #grab RT
    rt <- as.numeric(features[i,]$RT)
    #set RT tol range
    rtmax <- rt + rttol
    rtmin <- rt - rttol
    #set mz tol range
    mzmax <- mz + mztol
    mzmin <- mz - mztol 
    #find matches
    matches <- CCSdata %>% 
      filter(RT <= rtmax & RT >= rtmin) %>% 
      filter(`m/z` <= mzmax & `m/z` >= mzmin)
    #pick minimum RSD
    if(nrow(matches) > 0){
      ccs <- matches[1,] %>% slice_min(RSD, n = 1) %>% pull(CCSz)
    }
    #assign ccs value
    IDdata$CCS <- ifelse(IDdata$Compound == feat, ccs, IDdata$CCS)
  }
  ccsdata <- drop_na(IDdata, CCS)
  message(paste(length(unique(ccsdata$Compound)), "features with a CCS value out of", 
                length(unique(IDdata$Compound))))
  return(ccsdata)
}
