makeEcomData <- function(pathMamOccs, pathEcoTraits){
  require(e1071)
  Occs <- read.csv(pathMamOccs, header=T) # read occurrence matrix
  Occs <- Occs[,2:ncol(Occs)] #trim off first column ## This has already been done in Excel
  jdat <- Occs[2:(ncol(Occs)-2)]  #"just data columns without HID label
  jdat[jdat == 0] <- NA
  
  ecomt <- read.csv(pathEcoTraits, header=T) ## trait data

  #Find species count by point
  Occs$GenCount <- rowSums(!(is.na(jdat)))  #save species count in occurrence sheet
  
  ## Following nested for-loop fills in a matrix with trait values (sums) and
  ## N of species with data
  
  EM <- matrix(nrow=nrow(Occs),ncol=65) #creating empty matrix
  
  for(k in 1:nrow(Occs)){   #for each HID
    BMASS <- rep(NA, len=Occs$GenCount[k])  #empty vectors for trait values per HID   
    HYP <- rep(NA, len=Occs$GenCount[k])  
    LOP <- rep(NA, len=Occs$GenCount[k])
    OL <- rep(NA, len=Occs$GenCount[k])
    SF <- rep(NA, len=Occs$GenCount[k])
    OT <- rep(NA, len=Occs$GenCount[k])  
    BUN <- rep(NA, len=Occs$GenCount[k])
    ALX <- rep(NA, len=Occs$GenCount[k])
    OO <- rep(NA, len=Occs$GenCount[k])
    
    indx <- which(jdat[k,]=='1', arr.ind=TRUE)  #which entries in row k have 1 for that species
    species <- colnames(jdat)[indx[,2]]     #list of species at point
    
    
    if (length(species)>=1){
      ## BMASS loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species
        logl <- ecomt[ecomt$Genus==spec, "BodyMass_mean"] ## body MASS value for a given species 
        BMASS[p] <- logl[1]
      }
      
      ## HYP loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Genus==spec,"HYP_mode"] ## hypsodonty value for a given species 
        HYP[p] <- logl[1]
      }
      ## LOP loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Genus==spec,"LOP_mode"]
        LOP[p] <- logl[1]
      }
      ## OT loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Genus==spec,"OT_mode"]
        OT[p] <- logl[1]
      }
      
      ## OL loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Genus==spec,"OL_mode"]
        OL[p] <- logl[1]
      }
      
      ## SF loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Genus==spec,"SF_mode"]
        SF[p] <- logl[1]
      }
      
      ## BUN loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Genus==spec,"BUN_mode"]
        BUN[p] <- logl[1]
      }
      ## Exclusively acute loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Genus==spec,"AL_excl_mode"]
        ALX[p] <- logl[1]
      }
      ## OO goat type loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Genus==spec,"OO_goat_mode"]
        OO[p] <- logl[1]
      }
    }
    
    ## combining data for a HID
    EM[k, ] <- c(
      Occs[k,1], Occs[k, "GenCount"],
      mean(BMASS, na.rm=TRUE),
      sd(BMASS, na.rm=TRUE),
      min(BMASS, na.rm=TRUE),
      max(BMASS, na.rm=TRUE),
      max(BMASS, na.rm=TRUE)-min(BMASS, na.rm=TRUE),
      skewness(BMASS, na.rm=TRUE),
      kurtosis(BMASS, na.rm=TRUE, type = 1),
      mean(HYP, na.rm=TRUE),
      sd(HYP, na.rm=TRUE),
      min(HYP, na.rm=TRUE),
      max(HYP, na.rm=TRUE),
      max(HYP, na.rm=TRUE)-min(HYP, na.rm=TRUE),
      skewness(HYP, na.rm=TRUE),
      kurtosis(HYP, na.rm=TRUE, type = 1),
      mean(LOP, na.rm=TRUE),
      sd(LOP, na.rm=TRUE),
      min(LOP, na.rm=TRUE),
      max(LOP, na.rm=TRUE),
      max(LOP, na.rm=TRUE)-min(LOP, na.rm=TRUE),
      skewness(LOP, na.rm=TRUE),
      kurtosis(LOP, na.rm=TRUE, type = 1),
      mean(ALX, na.rm=TRUE),
      sd(ALX, na.rm=TRUE),
      min(ALX, na.rm=TRUE),
      max(ALX, na.rm=TRUE),
      max(ALX, na.rm=TRUE)-min(ALX, na.rm=TRUE),
      skewness(ALX, na.rm=TRUE),
      kurtosis(ALX, na.rm=TRUE, type = 1),
      mean(OL, na.rm=TRUE),
      sd(OL, na.rm=TRUE),
      min(OL, na.rm=TRUE),
      max(OL, na.rm=TRUE),
      max(OL, na.rm=TRUE)-min(OL, na.rm=TRUE),
      skewness(OL, na.rm=TRUE),
      kurtosis(OL, na.rm=TRUE, type = 1),
      mean(SF, na.rm=TRUE),
      sd(SF, na.rm=TRUE),
      min(SF, na.rm=TRUE),
      max(SF, na.rm=TRUE),
      max(SF, na.rm=TRUE)-min(SF, na.rm=TRUE),
      skewness(SF, na.rm=TRUE),
      kurtosis(SF, na.rm=TRUE, type = 1),
      mean(BUN, na.rm=TRUE),
      sd(BUN, na.rm=TRUE),
      min(BUN, na.rm=TRUE),
      max(BUN, na.rm=TRUE),
      max(BUN, na.rm=TRUE)-min(BUN, na.rm=TRUE),
      skewness(BUN, na.rm=TRUE),
      kurtosis(BUN, na.rm=TRUE, type = 1),
      mean(OT, na.rm=TRUE),
      sd(OT, na.rm=TRUE),
      min(OT, na.rm=TRUE),
      max(OT, na.rm=TRUE),
      max(OT, na.rm=TRUE)-min(OT, na.rm=TRUE),
      skewness(OT, na.rm=TRUE),
      kurtosis(OT, na.rm=TRUE, type = 1),
      mean(OO, na.rm=TRUE),
      sd(OO, na.rm=TRUE),
      min(OO, na.rm=TRUE),
      max(OO, na.rm=TRUE),
      max(OO, na.rm=TRUE)-min(OO, na.rm=TRUE),
      skewness(OO, na.rm=TRUE),
      kurtosis(OO, na.rm=TRUE, type = 1)
    )
    print(k)
    Sys.sleep(0.01)
    flush.console()
  }
  colnames(EM) <- c("HID", "GenCount", "BMASS", "BMASSsd", "BMASSmin", "BMASSmax", "BMASSrange",
                    "BMASSskew", "BMASSkurt", "HYP", "HYPsd", "HYPmin", "HYPmax", "HYPrange",
                    "HYPskew", "HYPkurt", "LOP", "LOPsd", "LOPmin", "LOPmax", "LOPrange",
                    "LOPskew", "LOPkurt", "ALX", "ALXsd", "ALXmin", "ALXmax", "ALXrange",
                    "ALXskew", "ALXkurt", "OL", "OLsd", "OLmin", "OLmax", "OLrange",
                    "OLskew", "OLkurt", "SF", "SFsd", "SFmin", "SFmax", "SFrange",
                    "SFskew", "SFkurt", "BUN", "BUNsd", "BUNmin", "BUNmax", "BUNrange",
                    "BUNskew", "BUNkurt", "OT", "OTsd", "OTmin", "OTmax", "OTrange",
                    "OTskew", "OTkurt", "OO", "OOsd", "OOmin", "OOmax", "OOrange",
                    "OOskew", "OOkurt")
  EM <- as.data.frame(EM)
  row.names(EM) <- NULL
  return(EM)
}
