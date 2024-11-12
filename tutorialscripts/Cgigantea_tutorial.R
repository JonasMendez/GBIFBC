#GBIFBioClim R Toolkit - Saguaro Tutorial
#visit www.endemicbio.info/GBIFBioClimR for documentation and tutorials
library(rgbif)
library(spocc)
library(mapr)
library(usethis)
library(geodata)
library(raster)
library(sp)
library(dplyr)
library(ggplot2)
library(dplyr)
library(lme4)
library(rstatix)

#Modify GBIF credentials
usethis::edit_r_environ()
#Set working directory for Worldclim/GBIF downloads
setwd()
#Define Working Directory
direc<-getwd()
#Import GBIF Data
#Below is where we source GBIF data for the species in question,
#Include Genus and specific epithet for the species
#Get Taxon key corresponding to species in question
taxonKey <- name_backbone("Carnegiea gigantea")$usageKey
#import species data using occ_download with taxon key
GBIF_spdata<-occ_download(pred("taxonKey", taxonKey), user = 'jonasmr', pwd = 'Mexicano12?')

#Get download code from output to monitor download
GBIF_spdata
#Check download
occ_download_wait('0003554-241107131044228')
#Get data once download is ready
gbif_data <- occ_download_get('0003554-241107131044228') %>%
  occ_download_import()

#source worldclim bioclimate variables data
bioclim<-worldclim_global('bio', 2.5, direc, version="2.1")
#prepare bioclim for stacking with elevation data
bioclim<- bioclim[[c(1:19)]] 
#get elevation data
bioelev<-worldclim_global('elev', 2.5, direc, version="2.1")
bioclimd<-stack(c(bioclim, bioelev))

#Filter to US samples, with Inaturalist data
GBIFilt<-GBIFilter_world(gbif_data, inat = T)

#filter old world samples
GBIFilt<-dplyr::filter(GBIFilt, decimalLongitude < 0)

#calculate extent
extent <- calculate_extent(GBIFilt, lon_column = "decimalLongitude", 
                           lat_column = "decimalLatitude", buffer = 1)
print(extent)

#crop raster with extent
bioclimd<-crop(bioclimd, extent)
#check extent by plotting raster
plot(bioclimd$wc2.1_2.5m_bio_1)

#project sample points
coords<-(cbind(GBIFilt$decimalLongitude,GBIFilt$decimalLatitude))
samplepoints<-SpatialPoints(coords, proj4string = crs(bioclimd))

#plot presence points on top of raster
plot(samplepoints, add=T)

#Extract bioclim data onto new dataframe based on lon/lat
biovals<-raster_extract(GBIFilt, bioclimd)
#remove sample with NA data and > 0m for elevation (we can't use negative values for lme4 models)
biovals<-dplyr::filter(biovals, !is.na(biovals$wc2.1_2.5m_elev))
biovals<-dplyr::filter(biovals, wc2.1_2.5m_elev > 0)

#make latitude and elevation groups
#check latitude range
max(biovals$decimalLatitude)
min(biovals$decimalLatitude)
#Latitude between 26-36, lets use 2.5 latitude intervals
biovals<-add_geo_labels(biovals,group_column = 'decimalLatitude',
                          interval = 2.5,
                          geolabel = 'lat')
#Check elevation range:
max(biovals$wc2.1_2.5m_elev)
min(biovals$wc2.1_2.5m_elev)
#We have an elevation range between 0 and 2500ft, so lets use 500ft intervals
biovals<-add_geo_labels(biovals, group_column = 'wc2.1_2.5m_elev',
                        interval = 500,
                        geolabel = 'elev')

#Plot latitude groups 
plotgeolabels(biovals, bioclimd, geolabel = 'lat',layer = 'wc2.1_2.5m_bio_1',
              plot_title = "C. gigantea Distribution based on Latitude Groups",
              x_label = 'Longitude',
              y_label = 'Latitude')

#Plot Elevation groups
plotgeolabels(biovals, bioclimd, geolabel = 'elev',layer = 'wc2.1_2.5m_bio_1',
              plot_title = "C. gigantea Distribution based on Elevation Groups",
              x_label = 'Longitude',
              y_label = 'Latitude')

#Assess relationship between latitude and elevation
# we will use the plotting function made by the Johnston lab to visualize
# the correlation between continuous variables: 
#(https://sejohnston.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

#plot Linear model to assess how elevation changes across latitude
ggplotRegression(lm(biovals$wc2.1_2.5m_elev ~ biovals$decimalLatitude))

#Positive correlation between latitude and elevation,
#i.e. higher latitudes also host higher elevation populations

#Elevation vs Latitude on BIO12 and BIO5
#Bio 12- Latitude
plot_bioclim2(biovals,value_column = 'wc2.1_2.5m_bio_12',
             color_column = 'lat',
             plot_type = 'density')
#Plot linear model
ggplotRegression(lm(biovals$wc2.1_2.5m_bio_12 ~ biovals$decimalLatitude))

#Bio12 - Elevation
plot_bioclim(biovals,value_column = 'wc2.1_2.5m_bio_12',
             color_column = 'elev',
             plot_type = 'density')
#Plot linear model
ggplotRegression(lm(biovals$wc2.1_2.5m_bio_12 ~ biovals$wc2.1_2.5m_elev))

#Bio5 - Latitude
plot_bioclim2(biovals,value_column = 'wc2.1_2.5m_bio_5',
             color_column = 'lat',
             plot_type = 'density')
#Plot linear model
ggplotRegression(lm(biovals$wc2.1_2.5m_bio_5 ~ biovals$decimalLatitude))

#Bio5 - Elevation
plot_bioclim(biovals,value_column = 'wc2.1_2.5m_bio_5',
             color_column = 'elev',
             plot_type = 'density')
#Plot linear model
ggplotRegression(lm(biovals$wc2.1_2.5m_bio_5 ~ biovals$wc2.1_2.5m_elev))

#BIO12~Elevation with latitdue as random effect
lmer1<-lmer(wc2.1_2.5m_bio_12 ~ wc2.1_2.5m_elev + (1 | lat), data = biovals)
summary(lmer1)
Anova(lmer1)

#BIO12~Elevation with latitude as a random effect
ip<- ggplot(data=biovals, aes(x = wc2.1_2.5m_elev , y= wc2.1_2.5m_bio_12, color = lat))+
  geom_point()+
  stat_smooth(method='lm')+
  #facet_grid(.~lat)+
  scale_x_continuous(name = 'Elevation')+
  theme(axis.title.x=element_text(face = 'bold', color='darkred'))
ip

#BIO5~Elevation with latitude as random effect
lmer2<-lmer(wc2.1_2.5m_bio_5 ~ wc2.1_2.5m_elev + (1 | lat), data = biovals)
summary(lmer2)
Anova(lmer2)

#BIO5~Elevation with latitude as a random effect
ip<- ggplot(data=biovals, aes(x = wc2.1_2.5m_elev , y= wc2.1_2.5m_bio_5, color = lat))+
  geom_point()+
  stat_smooth(method='lm')+
  #facet_grid(.~lat)+
  scale_x_continuous(name = 'Elevation')+
  theme(axis.title.x=element_text(face = 'bold', color='darkred'))
ip
