library(rgbif)
library(spocc)
library(mapr)
library(usethis)
library(geodata)
library(raster)
library(sp)
library(ggplot2)
#Make sure to run the functions provided in the GBIFWC.R prior to using this. 
#This script is based on the A. triloba tutorial that can be found at:
#Modify GBIF credentials
usethis::edit_r_environ()
#Set working directory for Worldclim/GBIF downloads
setwd()
#Define Working Directory
direc<-getwd()
#Import GBIF Data
#Include Genus and specific epithet for the species
#Get Taxon key corresponding to species in question
taxonKey <- name_backbone("Asimina triloba")$usageKey
#import species data using occ_download with taxon key
GBIF_spdata<-occ_download(pred("taxonKey", taxonKey), user = 'user', pwd = 'password')

#Get download code from output to monitor download
GBIF_spdata
#Check download
occ_download_wait('downloadkey')
#Get data once download is ready
gbif_data <- occ_download_get('downloadkey') %>%
  occ_download_import()

#source worldclim bioclimatee variables data
bioclim<-worldclim_global('bio', 2.5, direc, version="2.1")
#prepare bioclim for stacking with elevation data
bioclim<- bioclim[[c(1:19)]] 
#get elevation data
bioelev<-worldclim_global('elev', 2.5, direc, version="2.1")
bioclimd<-stack(c(bioclim, bioelev))
stack

#Filter GBIF data
GBIFilt<-GBIFilt_country(gbif_data, inat = F, 'US')
#remove samples from westcoast and africa by getting samples above -97 longitude
GBIFilt<-dplyr::filter(GBIFilt, decimalLongitude > -97)
GBIFilt<-dplyr::filter(GBIFilt, decimalLongitude < 0)
#Checking maximum longitude
max(GBIFilt$decimalLongitude)
#Checking Unique country codes to ensure only samples from the US
unique(GBIFilt$countryCode)

#Generate extent with buffer to crop raster
extent <- calculate_extent(GBIFilt, lon_column = "decimalLongitude",
                           lat_column = "decimalLatitude",
                           buffer = 4)
print(extent)
bioclimd<-crop(bioclimd, extent)
#check extent by plotting
plot(bioclimd$wc2.1_2.5m_bio_1)

#project sample points
coords<-(cbind(GBIFilt$decimalLongitude,GBIFilt$decimalLatitude))
samplepoints<-SpatialPoints(coords, proj4string = crs(bioclimd))

#plot presence points
plot(samplepoints, add=T)

#make longitude groups
lonlabel<-add_geo_labels(biovals,group_column = 'decimalLongitude',
                         interval = 6.5,
                         geolabel = 'long')
#plot longitude groups on map
plotgeolabels(lonlabel, bioclimd, geolabel = 'long',layer = 'wc2.1_2.5m_bio_1',
              plot_title = "A. triloba Distribution based on Longitude Groups",
              x_label = 'Longitude',
              y_label = 'Latitude')
#plot density distributions of bioclim variables
plot_bioclim2(lonlabel,value_column = 'wc2.1_2.5m_bio_1',
              color_column = 'long',
              plot_type = 'density')
