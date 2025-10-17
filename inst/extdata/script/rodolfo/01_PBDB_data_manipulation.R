##### Palaeobiology Data Base curatorial work ####

rm(list = ls());gc()

library(dplyr)
library(cowplot)

#potentially useful links of PBDB
#https://paleobiodb.org/data1.2/
#https://paleobiodb.org/data1.2/occs/geosum_doc.html
#https://paleobiodb.org/data1.2/occs/list_doc.html

# After downloading the data files and scripts, it is advised to create a new folder and paste all to the folder, here called MAIN:
dir.create("./MAIN/", recursive = T)
setwd("./MAIN/")

#import all Canidae data downloaded from PBDB
PBDB_NEW_Canidae_NA_full  <- read.csv(file="./pbdb_data_Canidae_NA_full.csv", header = TRUE, sep = ",")

colnames(PBDB_NEW_Canidae_NA_full)
unique(PBDB_NEW_Canidae_NA_full$cc)
#note that it includes occurrences  from Central America
# those are all that have "SV" "HN" "PA" "CU" "GL" "GT"  in the column "cc"
which(PBDB_NEW_Canidae_NA_full$cc == "SV")
which(PBDB_NEW_Canidae_NA_full$cc == "HN")
which(PBDB_NEW_Canidae_NA_full$cc == "PA")
which(PBDB_NEW_Canidae_NA_full$cc == "CU")
which(PBDB_NEW_Canidae_NA_full$cc == "GL")
which(PBDB_NEW_Canidae_NA_full$cc == "GT")
#"SV" is El Salvador
#"HN" is Honduras
#"PA" is Panama
#"CU" is Cuba
#"GL" is Greenland
#"GT" is Guatemala

#remove occurrences only identified to genus, family, subfamily, and tribe

unique(PBDB_NEW_Canidae_NA_full$accepted_rank)
length(which(PBDB_NEW_Canidae_NA_full$accepted_rank == "species"))
length(which(PBDB_NEW_Canidae_NA_full$accepted_rank == "genus"))
length(which(PBDB_NEW_Canidae_NA_full$accepted_rank == "family"))
length(which(PBDB_NEW_Canidae_NA_full$accepted_rank == "subfamily"))
length(which(PBDB_NEW_Canidae_NA_full$accepted_rank == "tribe"))
length(which(PBDB_NEW_Canidae_NA_full$accepted_rank == "subspecies"))

PBDB_NEW_Canidae_NA_full <- PBDB_NEW_Canidae_NA_full  %>%
  filter(! (accepted_rank %in% c("genus", "family", "subfamily", "tribe")))

#remove occurrences that either have "cf." or "?"

PBDB_NEW_Canidae_NA_full <- PBDB_NEW_Canidae_NA_full %>%
  filter(! (species_reso  %in% c("cf.","?")))




#remove the occurrences of Canis familiaris

PBDB_NEW_Canidae_NA_full <- PBDB_NEW_Canidae_NA_full %>%
  filter(PBDB_NEW_Canidae_NA_full$accepted_name != "Canis familiaris")


#replace the species name of some occurrences which have subspecies to only have the species name
#Urocyon cinereoargenteus townsendi; only two occurrences
#Vulpes vulpes macroura; only one occurrence
which(PBDB_NEW_Canidae_NA_full$accepted_name == "Urocyon cinereoargenteus townsendi")
which(PBDB_NEW_Canidae_NA_full$accepted_name == "Vulpes vulpes macroura")

PBDB_NEW_Canidae_NA_full$accepted_name[PBDB_NEW_Canidae_NA_full$accepted_name == "Urocyon cinereoargenteus townsendi"] <- "Urocyon cinereoargenteus"

PBDB_NEW_Canidae_NA_full$accepted_name[PBDB_NEW_Canidae_NA_full$accepted_name == "Vulpes vulpes macroura"] <- "Vulpes vulpes"

# remove several columns that we do not really need, to check the data
PBDB_NEW_Canidae_NA_simpler <- PBDB_NEW_Canidae_NA_full %>%
  select(occurrence_no, collection_no, accepted_name, accepted_rank, accepted_no, max_ma, min_ma, genus, cc, lng, lat, primary_name, primary_reso, species_name, species_reso)

# check the comments on the species_resolution
unique(PBDB_NEW_Canidae_NA_simpler$species_reso)


# remove more columns
PBDB_NEW_Canidae_NA_simpler <- PBDB_NEW_Canidae_NA_full %>%
  select(occurrence_no, collection_no, accepted_name, max_ma, min_ma, genus, cc, lng, lat)


## check the lat & long of occurrences to see if any fits a museum or any other locality that we should be careful
library(CoordinateCleaner)

#cc_inst
#Identify Records in the Vicinity of Biodiversity Institutions
#there aare no records that match the lat & long of Biodiversity institutions
test1= cc_inst(
  PBDB_NEW_Canidae_NA_simpler,
  lon = "lng",
  lat = "lat",
  species = "accepted_name",
  buffer = 100,
  geod = TRUE,
  ref = NULL,
  verify = FALSE,
  verify_mltpl = 10,
  value = "clean",
  verbose = TRUE
)



#######################################################################

#Identify Coordinates match the Gibif headquarters 
#there are no coordinates that match it
test2=cc_gbif(
  PBDB_NEW_Canidae_NA_simpler,
  lon = "lng",
  lat = "lat",
  species = "accepted_name",
  buffer = 1000,
  geod = TRUE,
  verify = FALSE,
  value = "clean",
  verbose = TRUE
)


#######################################################################

#Identify Coordinates in Vicinity of Country and Province Centroids
#there are no coordinates in the vivicinty of centroids
test3 <- cc_cen(
  PBDB_NEW_Canidae_NA_simpler,
  lon = "lng",
  lat = "lat",
  species = "accepted_name",
  buffer = 1000,
  geod = TRUE,
  test = "both",
  ref = NULL,
  verify = FALSE,
  value = "clean",
  verbose = TRUE
)

######
#identifies if any occurrence is at the ocean
#althugh the test flags some occurrences, many are not in fact on the ocean
#but rather on very small islands really near the continent (see details below)
#others are in the ocean but really, really close to the continent, e.g. 2 meters (see below)
# and others are placed in the ocean because of lack of resolution
#hence we decided to keep all of the 28 occurrences that were identified as being "at the ocean""
#Moreover, most of our analysis do not rely on lat & long, and the ones that do we explicitly discuss the limiations

test4 <- cc_sea(
  PBDB_NEW_Canidae_NA_simpler,
  lon = "lng",
  lat = "lat",
  ref = NULL,
  scale = 110,
  value = "flagged",
  speedup = FALSE,
  verbose = TRUE,
  buffer = NULL
)

#checking those occurrences that supposedly are on the ocean
length(PBDB_NEW_Canidae_NA_simpler[which(test4 == "FALSE"), ]$lat)

occurrences_at_the_sea <-PBDB_NEW_Canidae_NA_simpler[which(test4 == "FALSE"), ]

occurrences_at_the_sea <- occurrences_at_the_sea %>%
  arrange(lng, lat)

occurrences_at_the_sea
 
### plot those points in a map to visually inspect them
library(ggplot2)
library(maps)

world_map <- map_data("world")
na_map <- map_data("world", region = c("USA", "Canada", "Mexico", "Greenland", "Guatemala", "Honduras", "Panama", "El Salvador", "Cuba", "Belize", "Nicaragua", "Costa Rica", "Jamaica", "Dominican Republic", "Haiti", "Bahamas", "Barbados", "Trinidad and Tobago"))

ggplot() +
  geom_polygon(data = na_map, aes(x = long, y = lat, group = group), fill = "lightgray") +
  geom_point(data = occurrences_at_the_sea, aes(x = lng, y = lat), color = "red", size = 2) +
  coord_fixed(1.3) +
  theme_minimal() +
  labs(title = "Points on Map", x = "lng", y = "lat")

#detailed explanation for each flagged occurrence that were, many times mistankely identified to be at the sea
# occurrence 1429286 is at Saint Paul Island of Alaska, which is quite far from the main land; 57.17750 -170.34473
# occurrences 1614787 & 1614788 are not at the ocean but at a small islet called Herschel Island near the shore of  of Alaska; 69.57361 -138.89833
# occurrence 766675 is actually, according to googlemaps not at the ocean but at Alaska; 56.34 -133.58
#occurrence 1534974 is  actually, according to googlemaps, not at the ocean but at a small islet near Skagit bay that is easily connected to the main land. 48.21230 122.74992
#occurrence 1630722 is actually, according to googlemaps, not at the ocean but at a small islet called Santa Rosa Island near the shore of Los Angeles, USA; 34.00389 -120.19250
#occurrence 1630723 is actually, according to googlemaps, not at the ocean but at a small islet called San Nicolas Island near the shore of Los Angeles, USA;  33.27917 -119.53028
#occurrence 1622357 is actually, according to googlemaps, not at the ocean but at a small islet called Somerset Island near the shore (basically connected) to North Canada;  72.46670 -93.50460
#occurrence 187973 and 377062 are at the ocean but really close (about 2 km; but note that coordinates seem to be rounded) to the coast of Florida;   29.30000 -83.20000
#occurrence 198145 is at the ocean but really close (about 5 km; but note that coordinates seem to be rounded) to the coast of Florida; 28.10000 -80.50000
#occurrence 745057 is is actually, according to googlemaps, not at the ocean but at a small islet called Edisto Island near the shore (basically connected) of South Carolina, USA  32.47944 -80.33472
####occurrence 745056 is at the ocean but really close (about 10 km) to the coast of South Carolina   ; 33.01028 -79.08778   
#occurrence 1622316 is actually, according to googlemaps, not at the ocean but at a small islet called Ellesmere Island near the shore (basically connected) of North Canada  ; 78.75210 -75.67980
#occurrence 1374955 is is actually, according to googlemaps, at rally near the Davis Strait (about 50 meters) near the shore (basically connected) of Greenland 68.37144 -53.28046
#occurrence 1611242  is, according to googlemaps, at the sea but really close to Greenland (about 100 meters);   82.18361 -31.21944
#occurrence 1611355 is, according to googlemaps, at the sea but really close to Greenland (about 100 meters);    82.15667 -30.26667
#occurrence 1611369 is actually, according to googlemaps, not at the ocean but at Greenland; 82.15300 -30.20890
#occurrence 1611387 is actually, according to googlemaps, not at the ocean but at Greenland; 82.10520 -30.02760
#occurrence 1611843 is actually, according to googlemaps, not at the ocean but at Greenland 82.15090 -29.92880
#occurrence 1230888 is, according to googlemaps, at the sea but really close to Greenland (about 10 meters); 82.12500 -29.88833
#occurrence 1230889 is, according to googlemaps, at the sea but really close (about 10 meters) to Greenland 82.12500 -29.88833
#occurrence 1611867 is, according to googlemaps, at the sea but really close (about 100 meters) to Greenland  81.99583 -24.56472
#occurrence 1625005 is actually, according to googlemaps, not at the ocean but at the Clavering Island which is in fact connected to Greenland;  74.10917 -21.12083
#occurrence 1625019 is actually, according to googlemaps, not at the ocean but at the Clavering Island which is in fact connected to Greenland 74.12778 -20.77445
#occurrence 1624987 is actually, according to googlemaps, not at the ocean but at the Sabine Øer Island which is in fact rally close to Greenland;  74.55055 -18.79361
#occurrence 1624695 is actually, according to googlemaps, not at the ocean but at the Sabine Øer Island which is in fact rally close to Greenland; 74.54250 -18.78320
#occurrence 1611891 is according to googlemaps, at the sea but really close (about 500 meters) to on of the islands in Greenland 77.61083 -18.08861


### 
# 12 occurrences flagged for "artificial hotspots occurrence inventory.
#We will not remove those because manual check does not suggest any problem. 
#In addition most of our results do not rely on precise coordinates. 

test5 <- cc_aohi(
  PBDB_NEW_Canidae_NA_simpler,
  lon = "lng",
  lat = "lat",
  species = "Species",
  value = "flagged")

which(test5 == "FALSE")
PBDB_NEW_Canidae_NA_simpler[878, ]
#occurrence number 187875  lat & long  29.55 -82.51667 

PBDB_NEW_Canidae_NA_simpler[1159, ]
#occurrence number 193774   lat & long 28.9553 -82.6769

PBDB_NEW_Canidae_NA_simpler[1283, ]
#occurrence number 197793  lat & long    29 -82.68

PBDB_NEW_Canidae_NA_simpler[1284, ]
#occurrence number 197794  lat & long    29 -82.68

PBDB_NEW_Canidae_NA_simpler[1288, ]
#occurrence number 197990  lat & long   27.7 -82.5

PBDB_NEW_Canidae_NA_simpler[1289, ]
#occurrence number 197991  lat & long  27.7 -82.5

PBDB_NEW_Canidae_NA_simpler[1290, ]
#occurrence number 197992  lat & long  27.7 -82.5

PBDB_NEW_Canidae_NA_simpler[1291, ]
#occurrence number 198018  lat & long  27.7 -82.5

PBDB_NEW_Canidae_NA_simpler[1380, ]
#occurrence number 709947  lat & long  32.36667 -104.45 

PBDB_NEW_Canidae_NA_simpler[1381, ]
#occurrence number 709948  lat & long  32.36667 -104.45 

PBDB_NEW_Canidae_NA_simpler[1382, ]
#occurrence number 709950  lat & long  32.36667 -104.45 

PBDB_NEW_Canidae_NA_simpler[1596, ]
#occurrence number 1217425  lat & long  28.9553 -82.6769 



####test for country capital
## no flagged occurrence
test6 <- cc_cap(
  PBDB_NEW_Canidae_NA_simpler,
  lon = "lng",
  lat = "lat",
  species = "Species",
  value = "flagged")

# Testing equal lat/lon
#none showed any flagged occurrences
test7 <- cc_equ(
  PBDB_NEW_Canidae_NA_simpler,
  lon = "lng",
  lat = "lat",
  test = "absolute",
  value = "flagged",
  verbose = TRUE
)


# Testing coordinate validity
#none showed any flagged occurrences
test8 <- cc_val(
  PBDB_NEW_Canidae_NA_simpler,
  lon = "lng",
  lat = "lat",
  value = "flagged",
  verbose = TRUE
)



#####
# test for duplicate records. added the "Site" column as well. 
#It represents the PBDB "collection_no" and hence a way to get to a "local assemblage"
#Note that some species might have the same lat & long and time interval but 
#have different "collection_no". This is because some times the lat and long are not for
#the site itself but for a nearby place, for example a county. Because the "collection_no" is different, 
#those are treated differently.
#that said there were indeed some few replicated record where one species had the same 
#lat, long, and "collection_no"". Those represent duplicates and were removed, in this case
#those were only 11 records.

test9 <- cc_dupl(
  PBDB_NEW_Canidae_NA_simpler,
  lon = "lng",
  lat = "lat",
  species = "accepted_name",
  additions = "collection_no",
  value = "flagged",
  verbose = TRUE
)

# check which ones are duplicates
which(test9 == "FALSE")
PBDB_NEW_Canidae_NA_simpler[c(852, 1420, 1454, 1459, 1461, 1538, 1570, 1661, 1666, 1669, 1678), ]

PBDB_NEW_Canidae_NA_simpler[c(851, 852), ]

PBDB_NEW_Canidae_NA_simpler[c(1419, 1420), ]

PBDB_NEW_Canidae_NA_simpler[c(1454, 1262), ]

PBDB_NEW_Canidae_NA_simpler[c(1460, 1461), ]

PBDB_NEW_Canidae_NA_simpler[c(1458, 1459), ]

PBDB_NEW_Canidae_NA_simpler[c(1538, 1009), ]

PBDB_NEW_Canidae_NA_simpler[c(1570, 983), ]

PBDB_NEW_Canidae_NA_simpler[c(1661, 613), ]

PBDB_NEW_Canidae_NA_simpler[c(1666, 898), ]
 
PBDB_NEW_Canidae_NA_simpler[c(1669, 899), ]

PBDB_NEW_Canidae_NA_simpler[c(1678, 892), ] 


#remove those duplicates
PBDB_NEW_Canidae_NA_simpler <- cc_dupl(
  PBDB_NEW_Canidae_NA_simpler,
  lon = "lng",
  lat = "lat",
  species = "accepted_name",
  additions = "collection_no",
  value = "clean",
  verbose = TRUE
)

###########################################
# remove occurrences of particular species

#remove the single occurrence of Cubacyon transversidens because it is a dubious
#identification and likely synonym of Canis lupus familiaris (see Morgan and Woods 1986)
#Morgan and Woods (1986). Extinction and the zoogeography of West Indian land mammals. Hiolo,cical Journal of the Linnean Socieg (1986),28: 167-203. 
which(PBDB_NEW_Canidae_NA_simpler$accepted_name == "Cubacyon transversidens")

PBDB_NEW_Canidae_NA_simpler <- PBDB_NEW_Canidae_NA_simpler %>%
  filter(PBDB_NEW_Canidae_NA_simpler$accepted_name != "Cubacyon transversidens")

#remove the single occurrence of Prohesperocyon wilsoni.
#There is considerable doubt if this species belongs to Canidae or not.
#Some authors consider it to be a Miacidae
#Even when it is considered to be a Canidae, it is considered to the sister to all other Canids (Wang & Tedford 2008)

which(PBDB_NEW_Canidae_NA_simpler$accepted_name == "Prohesperocyon wilsoni")

PBDB_NEW_Canidae_NA_simpler <- PBDB_NEW_Canidae_NA_simpler %>%
  filter(PBDB_NEW_Canidae_NA_simpler$accepted_name != "Prohesperocyon wilsoni")


# remove two species (Canis feneus & Cynarctoides roii) for which 
#we do not have morphological data.
#they represent only 10 occurrences

which(PBDB_NEW_Canidae_NA_simpler$accepted_name == "Canis feneus")
which(PBDB_NEW_Canidae_NA_simpler$accepted_name == "Cynarctoides roii")

PBDB_NEW_Canidae_NA_simpler <- PBDB_NEW_Canidae_NA_simpler %>%
  filter(! (accepted_name  %in% c("Canis feneus","Cynarctoides roii")))



###########################################
#prepare dataset to use in analysis
#Include an underscore between the Genus and  species name 
#use the following column labels: Species	MaxT	MinT	lng	lat	Site

PBDB_NEW_Canidae_NA_simpler$accepted_name <-gsub(" ", "_", PBDB_NEW_Canidae_NA_simpler$accepted_name)


PBDB_NEW_Canidae_NA_simpler_PyRate <- PBDB_NEW_Canidae_NA_simpler %>%
  select(accepted_name, max_ma, min_ma, lng, lat, collection_no)

colnames(PBDB_NEW_Canidae_NA_simpler_PyRate) <- c("Species", "MaxT", "MinT", "lng", "lat", "Site")

# Figure with the ranges to remove - ESM 
pdf("./fig_s1.pdf", width = 8, height = 6)
ggplot(data = PBDB_NEW_Canidae_NA_simpler_PyRate, 
       aes(x = MaxT - MinT, fill = MaxT - MinT > 6)) +
  geom_histogram(bins = 17, color = "black",  
                 show.legend = F, linewidth = 1) +
  scale_fill_manual(values = c("white", "red")) +
  geom_vline(xintercept = 6, linetype = 2, linewidth = 2) + 
  labs(title = "Cut point of occurrences range", 
       x = "Range (Myr)", 
       y = "Count") +
  theme_cowplot() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 30),
        plot.title = element_text(size = 30))
dev.off()

#remove occurrences that have stratigraphic interval longer than 6MY
PBDB_NEW_Canidae_NA_simpler_PyRate_high_resol <- PBDB_NEW_Canidae_NA_simpler_PyRate %>%
  filter((PBDB_NEW_Canidae_NA_simpler_PyRate$MaxT-PBDB_NEW_Canidae_NA_simpler_PyRate$MinT) <= 6)


# create directory for curated data
dir.create("./PBDB/", recursive = T)

# write datasets to merge
write.table(PBDB_NEW_Canidae_NA_simpler_PyRate_high_resol, file="./PBDB/PBDB_NEW_Canidae_NA_simpler_PyRate_high_resol.csv", row.names = F, sep = ",")


#write csv to check the list of species

PBDB_NEW_Canidae_NA_simpler_PyRate_high_resol_Species_List <- sort(unique(PBDB_NEW_Canidae_NA_simpler_PyRate_high_resol$Species))
write.table(PBDB_NEW_Canidae_NA_simpler_PyRate_high_resol_Species_List, file="./PBDB/PBDB_NEW_Canidae_NA_simpler_PyRate_high_resol_Species_List.csv", row.names = F, sep = ",")

# end