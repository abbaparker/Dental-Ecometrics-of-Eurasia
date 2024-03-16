library(tidyverse)
library(RColorBrewer)
library(sf)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
theme_set(theme_bw())

#First, genus count by sample point, Eco-ISEA3H resolution 9, global ####
#read in hexagon centroid samples
setwd("C:/Users/abigpark/OneDrive - University of Helsinki/")
setwd("C:/Users/abbap/OneDrive - University of Helsinki/GitEcometrics")
hids <- read.csv("AllTerrestrialHIDs_continents.csv")

#Read in Genus occurrences
setwd("C:/Users/abigpark/OneDrive - University of Helsinki/Phylacine data/GENUS 2024")
setwd("C:/Users/abbap/OneDrive - University of Helsinki/Phylacine data/GENUS 2024")
GenusMatrixPN <- read.csv("PPPA_PhylacineGenusOccurrences_PresentNatural_Mar24_NoHomo.csv", header=T, row.names=1)
GenusMatrixP <- read.csv("PPPA_PhylacineGenusOccurrences_Present_Jan24_NoHomo.csv", , header=T, row.names=1)

# select points not over 2500m in elevation
hid_low <- hids %>% filter(Elevation_Mean < 2500)

# Find genus counts
colnames(GenusMatrixP)[1:12]
jdatP <- GenusMatrixP[2:ncol(GenusMatrixP)]  #"just data' columns with HID label
jdatP[jdatP == 0] <- NA
jdatP <- jdatP %>% mutate(GenCount=rowSums(!(is.na(.))))
GenusMatrixP$GenCountPresent <- jdatP$GenCount  #save species count in occurrence sheet
max(GenusMatrixP$GenCount)
jdatP <- jdatP[1:(ncol(jdatP)-1)]  #reset just data sheet to not include species count
# present natural
colnames(GenusMatrixPN)[1:12]
jdatPN <- GenusMatrixPN[2:ncol(GenusMatrixPN)]  #"just data' columns with HID label
jdatPN[jdatPN == 0] <- NA
jdatPN <- jdatPN %>% mutate(GenCount=rowSums(!(is.na(.))))
GenusMatrixP$GenCountPresentNatural <- jdatPN$GenCount  #save species count in occurrence sheet
max(GenusMatrixP$GenCountPresentNatural)
jdatPN <- jdatPN[1:(ncol(jdatPN)-1)]  #reset just data sheet to not include species count

#genus counts
gencount <- GenusMatrixP[,c(1,154,155)]
gencount <- gencount %>% mutate(Difference=GenCountPresentNatural - GenCountPresent)
summary(gencount)

#3 or more Present genera
gencount_use <- gencount %>% filter(GenCountPresent>2)
summary(gencount_use)

#merge with elevation filtered Eco-ISEA3H climate data
hid_use <- merge(hid_low, gencount_use, by="HID")
hid_use %>% count(Realm_Mode) #Australasia under-represented, preexisting bias towards Palearctic 
#plot points and generic richnesses
world <- ne_countries(scale="large", continent = c('europe', 'asia',
                                                   'oceania', 'north america',
                                                   'south america','africa'), returnclass = "sf")
defaultCRS <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +units=m +no_defs")
points_sf <- hid_use %>%
  st_as_sf(
    coords = c("X", "Y"),
    crs = defaultCRS)

ggplot()+ geom_sf(data=world)+
  geom_sf(data=points_sf, aes(color=Difference), size=1)+
  scale_color_gradientn(colors=brewer.pal(7,"YlOrRd"))+
  ggtitle("Generic enrichment Phylacine present natural vs. present")

ggplot()+ geom_sf(data=world)+
  geom_sf(data=points_sf, aes(color=GenCountPresentNatural), size=1)+
  scale_color_gradientn(colors=brewer.pal(7,"YlGnBu"))+
  ggtitle("Genus count present natural")

#Read in dental traits
setwd("C:/Users/abigpark/OneDrive - University of Helsinki/From Miikka/Jan24")
setwd("C:/Users/abbap/OneDrive - University of Helsinki/From Miikka/Jan24")
# For genera (use mode of species within genera for modern taxa)
Gtraits <- read.csv("DentalTraits_Genus_PPPA_withfossils.csv", header=T)
colnames(Gtraits)
PPPAgen <- Gtraits$Genus #list of genera in these 4 orders
#For species
Straits <- read.csv("DentalTraits_Species_120324.csv", header=T)

#Define PPPA taxa from Phylacine
setwd("C:/Users/abbap/OneDrive - University of Helsinki/GitEcometrics")
phy <- read.csv("Phylacine_trait_data.csv")

phy_pppa <- phy %>% filter(Order.1.2 %in% c("Perissodactyla", "Primates",
                                              "Proboscidea", "Cetartiodactyla"))

#remove aquatic species 
marine <- phy_pppa %>% filter(Marine==1)
marine <- marine %>% filter(!(Family.1.2 =="Hippopotamidae")) #leave hippos

# add freshwater fully-aquatic animals
fwat <- phy_pppa %>% filter(Family.1.2 %in% c("Delphiniidae","Delphinidae", "Iniidae", "Platanistidae",
                                         "Phocoenidae", "Trichechidae"))
aq <- rbind(marine, fwat) #add marine, freshwater families, plus these 3 awquatic species
aqs_sp <- aq$Binomial.1.2

phy_pppa_terr <- phy_pppa %>% filter(!(Binomial.1.2 %in% aqs_sp))
pppa_species <- phy_pppa_terr$Binomial.1.2

# Species count for global Eco-ISEA3H sample points
setwd("C:/Users/abigpark/OneDrive - University of Helsinki/Phylacine data/Occurrence Matrices")
setwd("C:/Users/abbap/OneDrive - University of Helsinki/Phylacine data/Occurrence Matrices")
# read in site/occurrence matrix for all terrestrial mammals
gmat <- read.table("PhylacineRes9_PresentNatural_NoFlying.txt", sep="\t", header=T)
gmat <- read.csv("PhylacineRes9_PresentNatural_AllMammalOccurrences.csv", header=T)
colnames(gmat) <- colnames(gmat) %>% str_replace("_Centroid", "")

#Present natural species occurrences in Perissodactyla, Artiodactyla, Proboscidea, and Primates
gmatPPPA_PN <- gmat %>% dplyr::select(any_of(c("HID", Straits$Binomial)))
gmatPPPA_PN <- gmat %>% dplyr::select(any_of(c("HID", pppa_species))) #what about using 814 species?

#same for present only
gmat_P <- read.table("PhylacineRes9_Present_NoFlying.txt", sep="\t", header=T)
gmat_P <- read.table("PhylacineRes9_Present_AllMammalOccurrences.txt", sep="\t", header=T)
colnames(gmat_P) <- colnames(gmat_P) %>% str_replace("_Centroid", "")
gmatPPPA_P <- gmat_P %>% dplyr::select(any_of(c("HID", Straits$Binomial)))
gmatPPPA_P <- gmat_P %>% dplyr::select(any_of(c("HID", pppa_species)))  #phylacine species list instead?

write.csv(gmatPPPA_PN, "PPPA_SpeciesOccurrences_PresentNatural_DentalTraitTaxa.csv")
write.csv(gmatPPPA_P, "PPPA_SpeciesOccurrences_Present_DentalTraitTaxa.csv")

#map differences in species richness, as for genera above
# from same set of sites with 3+ genera, under 2500m elevation
gmatPPPA_PN_use <- gmatPPPA_PN %>% filter(HID %in% hid_use$HID)
gmatPPPA_P_use <- gmatPPPA_P %>% filter(HID %in% hid_use$HID)

#PN species counts
colnames(gmatPPPA_PN_use)[1:12]
jdatPN <- gmatPPPA_PN_use[2:ncol(gmatPPPA_PN_use)]  #"just data' columns with HID label
jdatPN[jdatPN == 0] <- NA
jdatPN <- jdatPN %>% mutate(SpCount=rowSums(!(is.na(.))))
gmatPPPA_PN_use$SpCountPresentNatural <- jdatPN$SpCount  #save species count in occurrence sheet
max(gmatPPPA_PN_use$SpCountPresentNatural)
jdatPN <- jdatPN[1:(ncol(jdatPN)-1)]

#P species counts
jdatP <- gmatPPPA_P_use[2:ncol(gmatPPPA_P_use)]  #"just data' columns with HID label
jdatP[jdatP == 0] <- NA
jdatP <- jdatP %>% mutate(SpCount=rowSums(!(is.na(.))))
gmatPPPA_PN_use$SpCountPresent <- jdatP$SpCount  #save species count in occurrence sheet
max(gmatPPPA_PN_use$SpCountPresent)
jdatP <- jdatP[1:(ncol(jdatP)-1)]

spcount <- gmatPPPA_PN_use[,c(1,815,816)]   #1, 815,816)]          #
spcount <- spcount %>% mutate(Sp_Difference=SpCountPresentNatural - SpCountPresent)
summary(spcount)

#merge with climate data and genus counts
hid_use_sp <- merge(hid_use, spcount, by="HID")

points_sf <- hid_use_sp %>%
  st_as_sf(
    coords = c("X", "Y"),
    crs = defaultCRS)

ggplot()+ geom_sf(data=world)+
  geom_sf(data=points_sf, aes(color=Sp_Difference), size=1)+
  scale_color_gradientn(colors=brewer.pal(7,"YlOrRd"))+
  ggtitle("Species enrichment Phylacine present natural vs. present")

ggplot()+ geom_sf(data=world)+
  geom_sf(data=points_sf, aes(color=SpCountPresent), size=1)+
  scale_color_gradientn(colors=brewer.pal(7,"YlGnBu"))+
  ggtitle("Species count present")

#which species are added outside dental data?
phy_addt <- phy_pppa_terr %>% filter(!(Binomial.1.2 %in% Straits$Binomial))

phy_addt_large <- phy_addt %>% filter(Mass.g>=1000)
phy_addt_largeH <- phy_addt_large %>% filter(Diet.Plant>50) #%>% filter(Diet.Vertebrate<50)
phy_addt_largeH_mod <- phy_addt_largeH %>% filter(!(IUCN.Status.1.2 %in% c("EP", "EX")))
phy_addt_largeH_ex <- phy_addt_largeH %>% filter(IUCN.Status.1.2 %in% c("EP", "EX")) #these 122 need to be added to best use PN data
# which extinct genera in PN dataset have no dental data? 61 genera
unique(phy_addt_largeH_ex$Genus.1.2[!(phy_addt_largeH_ex$Genus.1.2 %in% PPPAgen)])

