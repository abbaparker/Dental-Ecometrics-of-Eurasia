#Initial trait sampling code from Miikka

makeEcomData <- function(pathMamOccs, pathEcoTraits, massVar){
  require(e1071)
  Occs <- read.csv(pathMamOccs, header=T) # read occurrence matrix
  Occs <- Occs[,2:ncol(Occs)] #trim off first column ## This has already been done in Excel
  jdat <- Occs[2:ncol(Occs)]  #"just data columns without HID label
  jdat[jdat == 0] <- NA

  ecomt <- read.csv(pathEcoTraits, header=T) ## trait data
  ## enam <- names(ecomt) No need yet in this version...
  #Find species count by point
  Occs$SpCount <- rowSums(!(is.na(jdat)))  #save species count in occurrence sheet

  ## Following nested for-loop fills in a matrix with trait values (sums) and
  ## N of species with data

  EM <- matrix(nrow=nrow(Occs),ncol=58) #creating empty matrix

  for(k in 1:nrow(Occs)){   #for each HID
    BMASS <- rep(NA, len=Occs$SpCount[k])  #empty vectors for trait values per HID   
    HYP <- rep(NA, len=Occs$SpCount[k])  
    LOP <- rep(NA, len=Occs$SpCount[k])
    OL <- rep(NA, len=Occs$SpCount[k])
    SF <- rep(NA, len=Occs$SpCount[k])
    AL <- rep(NA, len=Occs$SpCount[k])  
    BUN <- rep(NA, len=Occs$SpCount[k])
    OT <- rep(NA, len=Occs$SpCount[k])
    AA <- rep(NA, len=Occs$SpCount[k])

    ## V[1] <-  #HID as first entry
    indx <- which(jdat[k,]=='1', arr.ind=TRUE)  #which entries in row k have 1 for that species
    species <- colnames(jdat)[indx[,2]]     #list of species at point
    
    
    if (length(species)>=1){
    ## BMASS loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species
        logl <- ecomt[ecomt$Binomial==spec, "massg"] ## body MASS value for a given species 
        BMASS[p] <- logl[1]
      }

    ## HYP loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species namesum(AL, na.rm=TRUE),
        logl <- ecomt[ecomt$Binomial==spec,"HYP"] ## hypsodonty value for a given species 
        HYP[p] <- logl[1]
      }
    ## LOP loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Binomial==spec,"LOP"]
        LOP[p] <- logl[1]
      }
    ## AL loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Binomial==spec,"AL"]
        AL[p] <- logl[1]
      }

    ## OL loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Binomial==spec,"OL"]
        OL[p] <- logl[1]
      }
  
    ## SF loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Binomial==spec,"SF"]
        SF[p] <- logl[1]
      }
 
    ## BUN loop
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Binomial==spec,"BUN"]
        BUN[p] <- logl[1]
      }
    ## OT
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Binomial==spec,"OT"]
        OT[p] <- logl[1]
      }      
      # Exclusively acute
      for (p in (1:length(species))) {
        spec <- as.character(species[p])  #pulls species name
        logl <- ecomt[ecomt$Binomial==spec,"Excl_AL"]
        AA[p] <- logl[1]
      }
    }

    ## combining data for a HID
    EM[k, ] <- c(
      Occs[k,1], Occs[k, "SpCount"],
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
      mean(AL, na.rm=TRUE),
      sd(AL, na.rm=TRUE),
      min(AL, na.rm=TRUE),
      max(AL, na.rm=TRUE),
      max(AL, na.rm=TRUE)-min(AL, na.rm=TRUE),
      skewness(AL, na.rm=TRUE),
      kurtosis(AL, na.rm=TRUE, type = 1),
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
      kurtosis(OT, na.rm=TRUE, type = 1)
      mean(AA, na.rm=TRUE),
      sd(AA, na.rm=TRUE),
      min(AA, na.rm=TRUE),
      max(AA, na.rm=TRUE),
      max(AA, na.rm=TRUE)-min(OT, na.rm=TRUE),
      skewness(AA, na.rm=TRUE),
      kurtosis(AA, na.rm=TRUE, type = 1)
    )
    print(k)
    Sys.sleep(0.01)
    flush.console()
  }
  colnames(EM) <- c("HID", "SpCount", "BMASS", "BMASSsd", "BMASSmin", "BMASSmax", "BMASSrange",
  "BMASSskew", "BMASSkurt", "HYP", "HYPsd", "HYPmin", "HYPmax", "HYPrange",
  "HYPskew", "HYPkurt", "LOP", "LOPsd", "LOPmin", "LOPmax", "LOPrange",
  "LOPskew", "LOPkurt", "AL", "ALsd", "ALmin", "ALmax", "ALrange",
  "ALskew", "ALkurt", "OL", "OLsd", "OLmin", "OLmax", "OLrange",
  "OLskew", "OLkurt", "SF", "SFsd", "SFmin", "SFmax", "SFrange",
  "SFskew", "SFkurt", "BUN", "BUNsd", "BUNmin", "BUNmax", "BUNrange",
  "BUNskew", "BUNkurt", "OT", "OTsd", "OTmin", "OTmax", "OTrange",
  "OTskew", "OTkurt", "AA", "AAsd", "AAmin", "AAmax", "AArange",
  "AAskew", "AAkurt")
  EM <- as.data.frame(EM)
  row.names(EM) <- NULL
  return(EM)
}


## Plotting function for relationships between climate and ecometric
## variables across continents
plotEcomRegions <- function(data, clivar, ecomvar, continents=TRUE){
  xlim <- c(floor(min(data[, clivar])), ceiling(max(data[, clivar]))) 
  ylim <- c(floor(min(data[, ecomvar], na.rm=TRUE)), ceiling(max(data[, ecomvar], na.rm=TRUE))) 

  if(continents==TRUE){
    par(mfrow=c(2,2))
    plot(data[data$Africa==TRUE, ][, clivar], data[data$Africa==TRUE, ][, ecomvar],
    col="red", xlim=xlim, ylim=ylim, xlab=clivar, ylab=ecomvar, main="Africa")
    
    plot(data[data$Eurasia==TRUE, ][, clivar], data[data$Eurasia==TRUE, ][, ecomvar],
    col="blue", xlim=xlim, ylim=ylim, xlab=clivar, ylab=ecomvar, main="Eurasia")
    
    plot(data[data$North_America==TRUE, ][, clivar], data[data$North_America==TRUE, ][, ecomvar],
    col="green", xlim=xlim, ylim=ylim, xlab=clivar, ylab=ecomvar, main="North America")
    
    plot(data[data$South_America==TRUE, ][, clivar], data[data$South_America==TRUE, ][, ecomvar],
    col="brown", xlim=xlim, ylim=ylim, xlab=clivar, ylab=ecomvar, main="South America")

  } else {
    par(mfrow=c(2,3))
    plot(data[data$Realm_Mode=="Afrotropic", ][, clivar], data[data$Realm_Mode=="Afrotropic", ][, ecomvar],
    col="red", xlim=xlim, ylim=ylim, xlab=clivar, ylab=ecomvar, main="Afrotropic")
    
    plot(data[data$Realm_Mode=="Neotropic", ][, clivar], data[data$Realm_Mode=="Neotropic", ][, ecomvar],
    col="brown", xlim=xlim, ylim=ylim, xlab=clivar, ylab=ecomvar, main="Neotropic")
    
    plot(data[data$Realm_Mode=="Indomalaya", ][, clivar], data[data$Realm_Mode=="Indomalaya", ][, ecomvar],
    col="blue", xlim=xlim, ylim=ylim, xlab=clivar, ylab=ecomvar, main="Indomalaya")
    
    plot(data[data$Realm_Mode=="Nearctic", ][, clivar], data[data$Realm_Mode=="Nearctic", ][, ecomvar],
    col="green", xlim=xlim, ylim=ylim, xlab=clivar, ylab=ecomvar, main="Nearctic")
      
    plot(data[data$Realm_Mode=="Palearctic", ][, clivar], data[data$Realm_Mode=="Palearctic", ][, ecomvar],
    col="black", xlim=xlim, ylim=ylim, xlab=clivar, ylab=ecomvar, main="Palearctic")
  }
}



## movePoints function to move points that hit ocean to the nearest land point.
## this function is used within the next function. Highlight and run.

## FUNCTION STARTS HERE
movePoints <- function(ras, locs, distance=NULL){
    require(terra)
    require(sf)
    require(class)


    locsV <- vect(as.matrix(locs))
    crs(locsV) <- crs(ras)

    if (nlyr(ras) > 1) ras <- rast(ras[[1]])
    
    ## Get NA cells
    rasvals <- terra::extract(ras, locsV)
    miss <- is.na(rasvals)[, 2]
    
    if (sum(miss) > 0) {
        coord.miss <- st_as_sf(locsV[miss, ])  # points without data
        cells.notNA <- terra::as.points(ras, values=FALSE)  # get coordinates of cells with data
        coord.ras <- st_as_sf(cells.notNA)
        cell.id <- factor(seq_len(nrow(coord.ras)))

        coord.miss2 <- st_coordinates(coord.miss)
        coord.ras2 <- st_coordinates(coord.ras)

        # find the nearest raster cell for each point with missing data
        nearest.cell <- class::knn1(coord.ras2, coord.miss2, cell.id)

        new.coords <- as.data.frame(matrix(coord.ras2[nearest.cell, ], ncol = 2))
        rownames(new.coords) <- NULL
        colnames(new.coords) <- c("x", "y")

        if (!is.null(distance)){
        new.coordsSF <- st_as_sf(new.coords, coords=c("x", "y"))
        st_crs(new.coordsSF) <- st_crs(coord.miss)
        # calculate distances between old and new coordinates
        distances <- as.vector(diag(st_distance(coord.miss, new.coordsSF)))
        
        # if distance below threshold, accept, otherwise keep old coordinates
        x <- ifelse(distances < distance, new.coords[,1], coord.miss2[,1])
        y <- ifelse(distances < distance, new.coords[,2], coord.miss2[,2])
        new.coords <- cbind(x, y)

        }
        # assign new coordinates to those points
        locs[miss, ] <- new.coords
        return(locs)

    } else {
        message("All points fall within a raster cell")
        return(locs)

    }

}
## FUNCTION ENDS HERE


## extractClimateRange
## This function extracts datasets of climate variables for each location and saves them into a list
## 

## FUNCTIONS STARTS HERE
extractClimateRange <- function(precipPaths, tempPaths, locs, minAges, maxAges, distance=NULL, bufferSize, locNames=NULL){
    require(terra)
    require(ncdf4)
    require(class)
    require(sf)
    
    ## Read in NetCDF files
    precip1 <- precipPaths[1]
    precip2 <- precipPaths[2]
    precip3 <- precipPaths[3]
    temp1 <- tempPaths[1]
    temp2 <- tempPaths[2]
    temp3 <- tempPaths[3]

    nc_precip1 <- nc_open(precip1)
    nc_precip2 <- nc_open(precip2)
    nc_precip3 <- nc_open(precip3)
    nc_temp1 <- nc_open(temp1)
    nc_temp2 <- nc_open(temp2)
    nc_temp3 <- nc_open(temp3)

    lon <- ncvar_get(nc_precip1, "lon")
    lat <- ncvar_get(nc_precip1, "lat", verbose = F)
    t1 <- rev(ncvar_get(nc_precip1, "time"))
    t2 <- rev(ncvar_get(nc_precip2, "time"))
    t3 <- rev(ncvar_get(nc_precip3, "time"))
    fillvalueP1 <- ncatt_get(nc_precip1, "precip_mm_srf", "_FillValue")
    fillvalueP2 <- ncatt_get(nc_precip2, "precip_mm_srf", "_FillValue")
    fillvalueP3 <- ncatt_get(nc_precip3, "precip_mm_srf", "_FillValue")
    fillvalueT1 <- ncatt_get(nc_temp1, "temp_mm_1_5m", "_FillValue")
    fillvalueT2 <- ncatt_get(nc_temp2, "temp_mm_1_5m", "_FillValue")
    fillvalueT3 <- ncatt_get(nc_temp3, "temp_mm_1_5m", "_FillValue")

    precip.array1 <- ncvar_get(nc_precip1, "precip_mm_srf")
    nc_close(nc_precip1)
    precip.array1[precip.array1 == fillvalueP1$value] <- NA
    precip.array1 <- precip.array1[, , dim(precip.array1)[3]:1]

    precip.array2 <- ncvar_get(nc_precip2, "precip_mm_srf")
    nc_close(nc_precip2)
    precip.array2[precip.array2 == fillvalueP2$value] <- NA
    precip.array2 <- precip.array2[, , dim(precip.array2)[3]:1]

    precip.array3 <- ncvar_get(nc_precip3, "precip_mm_srf")
    nc_close(nc_precip3)
    precip.array3[precip.array3 == fillvalueP3$value] <- NA
    precip.array3 <- precip.array3[, , dim(precip.array3)[3]:1]

    ###

    temp.array1 <- ncvar_get(nc_temp1, "temp_mm_1_5m")
    nc_close(nc_temp1)
    temp.array1[temp.array1 == fillvalueT1$value] <- NA
    temp.array1 <- temp.array1[, , dim(temp.array1)[3]:1]

    temp.array2 <- ncvar_get(nc_temp2, "temp_mm_1_5m")
    nc_close(nc_temp2)
    temp.array2[temp.array2 == fillvalueT2$value] <- NA
    temp.array2 <- temp.array2[, , dim(temp.array2)[3]:1]

    temp.array3 <- ncvar_get(nc_temp3, "temp_mm_1_5m")
    nc_close(nc_temp3)
    temp.array3[temp.array3 == fillvalueT3$value] <- NA
    temp.array3 <- temp.array3[, , dim(temp.array3)[3]:1]

    ###

    ## First we need to check, if some of the localities
    ## hit ocean points in the climate layers. If they do
    ## we will move them to the closest land cell. We will
    ## use the modern conditions as a land mask, because the 
    ## highest sea level and the smalles land area.
    ## 

    ## Take the first time slice from the precipitation array.
    ## Note that we have reversed the order of layers in the 
    ## array, so that the youngest layer is the first.
    landMask <- flip(t(rast(precip.array1[, , 1], extent=ext(c(min(lat), max(lat), min(lon), max(lon))),
                crs=crs(rast()))), direction="vertical")*360

    ## Then we use movePoints -function to check if some of the points hit
    ## ocean, and if so move them to the closest land cell.

    newCoords <- movePoints(ras=landMask, locs=locs)

    climDataList <- list()
    
    ## FOR LOOP HERE
    for(i in 1:nrow(locs)){

        LON <- newCoords[i, 1]  ## column 1 here is longitude in locs matrix
        LAT <- newCoords[i, 2] ## column 2 here is latitude in locs matrix
        max_a <- maxAges[i]*1000000  ## ages in years, not Ma
        min_a <- minAges[i]*1000000 ##
        max_a <- round(max_a, digits=-3) #find oldest 1000 years in age range
        min_a <- round(min_a, digits=-3)  #find youngest 1000 years in age range

        ## Define the area 10 x 10 deg area around the locality

        lonMAX <- LON + 5
        lonMIN <- LON - 5
        latMAX <- LAT + 5
        latMIN <- LAT - 5

        ## and find the closest (corner) grid cells in the climate data
        lonMAXindx <- which.min(abs(lon - lonMAX))
        lonMINindx <- which.min(abs(lon - lonMIN))
        latMAXindx <- which.min(abs(lat - latMAX))
        latMINindx <- which.min(abs(lat - latMIN))


        extr <- ext(c(latMIN, latMAX, lonMIN, lonMAX))

        locB <- buffer(vect(cbind(LON, LAT)), bufferSize)

        #get data, from whichever 1 Ma bin site falls in (5 combinations)
        if (max_a < 1001000){  #just 1-0Ma bin
            t_indx1 <- which(t1 <= max_a & t1 >= min_a)  #define age range to pull data from
            precipR1 <- flip(t(rast(precip.array1[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")*360
            tempR1 <- flip(t(rast(temp.array1[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")
        }
        
        if (max_a < 2001000 & min_a >= 1001000){  #just 2-1Ma bin
            t_indx1 <- which(t2 <= max_a & t2 >= min_a)  #define age range to pull data from
            precipR1 <- flip(t(rast(precip.array2[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")*360
            tempR1 <- flip(t(rast(temp.array2[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")
        }

        if (min_a >= 2001000){  #just 3-2Ma bin
            t_indx1 <- which(t3 <= max_a & t3 >= min_a)  #define age range to pull data from
            precipR1 <- flip(t(rast(precip.array3[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")*360
            tempR1 <- flip(t(rast(temp.array3[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")
        }
        
        if (max_a > 1001000 & max_a < 2001000 & min_a < 1001000){  #just 2-1 and 1-0 Ma bins
            t_indx1 <- which(t1 <= max_a & t1 >= min_a)  #define age range to pull data from
            t_indx2 <- which(t2 <= max_a & t2 >= min_a) # for second bin
            precipR1 <- flip(t(rast(precip.array1[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")*360
            tempR1 <- flip(t(rast(temp.array1[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")
            precipR2 <- flip(t(rast(precip.array2[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx2], extent=extr, crs=crs(rast()))), direction="vertical")*360
            tempR2 <- flip(t(rast(temp.array2[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx2], extent=extr, crs=crs(rast()))), direction="vertical")
            precipR1 <- c(precipR1, precipR2)
            tempR1 <- c(tempR1, tempR2)
        }
        
        if (max_a > 2001000 & min_a < 2001000 & min_a >= 1001000 ){  #just 3-2 and 2-1 Ma bins
            t_indx1 <- which(t2 <= max_a & t2 >= min_a)  #define age range to pull data from
            t_indx2 <- which(t3 <= max_a & t3 >= min_a)   # for second bin
            precipR1 <- flip(t(rast(precip.array2[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")*360
            tempR1 <- flip(t(rast(temp.array2[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")
            precipR2 <- flip(t(rast(precip.array3[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx2], extent=extr, crs=crs(rast()))), direction="vertical")*360
            tempR2 <- flip(t(rast(temp.array3[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx2], extent=extr, crs=crs(rast()))), direction="vertical")
            
            precipR1 <- c(precipR1, precipR2)
            tempR1 <- c(tempR1, tempR2)

        }

        if (max_a > 2001000 & min_a < 1001000){  #all three bins (widest age range)
            print("age range too long for anything useful")
            t_indx1 <- which(t1 <= max_a & t1 >= min_a)  #define age range to pull data from
            t_indx2 <- which(t2 <= max_a & t2 >= min_a)
            t_indx3 <- which(t3 <= max_a & t3 >= min_a)
            precipR1 <- flip(t(rast(precip.array1[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")*360
            tempR1 <- flip(t(rast(temp.array1[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx1], extent=extr, crs=crs(rast()))), direction="vertical")
            precipR2 <- flip(t(rast(precip.array2[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx2], extent=extr, crs=crs(rast()))), direction="vertical")*360
            tempR2 <- flip(t(rast(temp.array2[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx2], extent=extr, crs=crs(rast()))), direction="vertical")
            precipR3 <- flip(t(rast(precip.array3[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx3], extent=extr, crs=crs(rast()))), direction="vertical")*360
            tempR3 <- flip(t(rast(temp.array3[lonMINindx:lonMAXindx, latMINindx:latMAXindx, t_indx3], extent=extr, crs=crs(rast()))), direction="vertical")

            precipR1 <- c(precipR1, precipR2, precipR3)
            tempR1 <- c(tempR1, tempR2, tempR3) 
        }
        climDataList[[i]] <- cbind(seq(min_a, max_a, 1000), 
                                    as.vector(t(extract(precipR1, locB, fun=mean, na.rm=TRUE, ID=FALSE))),
                                    as.vector(t(extract(tempR1, locB, fun=mean, na.rm=TRUE, ID=FALSE))))
        
         print(i)

    }   ## FOR LOOP ENDS HERE

    if(!is.null(locNames)){
        names(climDataList) <- locNames
    }

    rm(precip.array1, precip.array2, precip.array3, temp.array1, temp.array2, temp.array3,
    landMask, precipR1, precipR2, precipR3, tempR1, tempR2, tempR3)
    gc()

    return(climDataList)

}
## FUNCTION ENDS HERE


## This function calculates number of climate variables from individual climate dataset
## created by the above function
## FUNCTIONS STARTS HERE
calculateClimateValues <- function(data, prob=c(0.1, 0.9)){

    climat <- data.frame(matrix(nrow=length(data), ncol=6))

    for(i in 1:length(data)){
        precip.d <- density(data[[i]][, 2], bw=20)
        temp.d <- density(data[[i]][, 3], bw=0.5)

        precip.mode <- precip.d$x[which.max(precip.d$y)]
        precip.min <- as.numeric(quantile(data[[i]][, 2], probs=0.1))
        precip.max <- as.numeric(quantile(data[[i]][, 2], probs=0.9))
        temp.mode <- temp.d$x[which.max(temp.d$y)]
        temp.min <- as.numeric(quantile(data[[i]][, 3], probs=0.1))
        temp.max <- as.numeric(quantile(data[[i]][, 3], probs=0.9))

        climat[i, ] <- c(precip.mode, precip.min, precip.max, temp.mode, temp.min, temp.max)
    }
    names(climat) <- c("precipMode","precipMin","precipMax","tempMode","tempMin","tempMax")
    return(climat)
}
## FUNCTIONS ENDS HERE