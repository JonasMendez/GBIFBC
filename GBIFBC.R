#GBIFBioClim R Toolkit
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

#Set up to get GBIF data
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
GBIF_spdata<-occ_download(pred("taxonKey", taxonKey), user = 'user', pwd = 'password')

#Get download code from output to monitor download
GBIF_spdata
#Check download
occ_download_wait('GBIF Download ID')
#Get data once download is ready
gbif_data <- occ_download_get('GBIF Download ID') %>%
  occ_download_import()

#Getting bioclim data
#Change the resolution if needed (default 2.5)

#Source worldclim bioclimate variables data
bioclim<-worldclim_global('bio', 2.5, direc, version="2.1")
#prepare bioclim for stacking with elevation data
bioclim<- bioclim[[c(1:19)]] 
#get elevation datam ensure the same resolution for stacking
bioelev<-worldclim_global('elev', 2.5, direc, version="2.1")
#Stack bioclim and elevation variables
bioclimd<-stack(c(bioclim, bioelev))

#Functions
#Filter GBIG data by lon/lat, include/exclude inaturalist, and or filter by country code
#This version of the function allows automatic filtering by country code
GBIFilter_country <- function(gbifdf, inat, ccode){
  if(inat == T){#include inaturalist collections
    sp<-dplyr::filter(gbifdf, !is.na(decimalLongitude), !is.na(decimalLatitude), countryCode == ccode)
    #Remove Dupes
    dups<- duplicated(sp[, c('decimalLongitude', 'decimalLatitude')])
    sp <- sp[!dups, ] #remove dupes
    return(sp) }else{#exclude inaturalist collections
      sp<-dplyr::filter(gbifdf, !is.na(decimalLongitude), !is.na(decimalLatitude), institutionCode != "iNaturalist", countryCode == ccode)
      #Remove Dupes
      dups<- duplicated(sp[, c('decimalLongitude', 'decimalLatitude')])
      sp <- sp[!dups, ] #remove dupes
      return(sp)
    }
}

#This version samples worldwide distributions
GBIFilter_world <- function(gbifdf, inat){
  if(inat == T){#include inaturalist collections
    sp<-dplyr::filter(gbifdf, !is.na(decimalLongitude))
    #Remove Dupes
    dups<- duplicated(sp[, c('decimalLongitude', 'decimalLatitude')])
    sp <- sp[!dups, ] #remove dupes
    return(sp) }else{#exclude inaturalist collections
      sp<-dplyr::filter(gbifdf, !is.na(decimalLongitude), !is.na(decimalLatitude), institutionCode != "iNaturalist")
      #Remove Dupes
      dups<- duplicated(sp[, c('decimalLongitude', 'decimalLatitude')])
      sp <- sp[!dups, ] #remove dupes
      return(sp)
    }
}

#Function to get extent + buffer of GBIF sampling
calculate_extent <- function(data, 
                             lon_column = "longitude", 
                             lat_column = "latitude", 
                             buffer = 0) {
  # Find the minimum and maximum values for longitude and latitude
  lon_min <- min(data[[lon_column]], na.rm = TRUE)
  lon_max <- max(data[[lon_column]], na.rm = TRUE)
  lat_min <- min(data[[lat_column]], na.rm = TRUE)
  lat_max <- max(data[[lat_column]], na.rm = TRUE)
  
  # Apply the buffer to the extent
  lon_min <- lon_min - buffer
  lon_max <- lon_max + buffer
  lat_min <- lat_min - buffer
  lat_max <- lat_max + buffer
  
  # Create the extent list
  extent <- c(lon_min, lon_max, lat_min, lat_max)
  
  return(extent)
}

#Function to extract bioclim data from GBIF lon/lat data
raster_extract<- function(gbifdist, bioclimd) {
  spoints<-sp::SpatialPoints(spdist[,c('decimalLongitude', 'decimalLatitude')], proj4string = crs(bioclimd))
  wclimvals <- raster::extract(bioclimd, spoints)
  spvals <- data.frame(cbind(spdist, wclimvals))
  #filtvals<-dplyr::filter(spvals, !is.na(bio1), !is.na(bio2), !is.na(bio3), !is.na(bio4), !is.na(bio5), !is.na(bio6), !is.na(bio7), !is.na(bio8), !is.na(bio9), !is.na(bio10), !is.na(bio11), !is.na(bio12), !is.na(bio13), !is.na(bio14), !is.na(bio15), !is.na(bio16), !is.na(bio17), !is.na(bio18), !is.na(bio19), !is.na(alt))
  return(spvals)
}

#Function to make lon/lat/elevation or other numeric value labels based on intervals
add_geo_labels <- function(data, group_column = "group_column", interval = 10, geolabel = 'geolabel') {
  # Remove rows with NA or non-numeric values in the group_column
  data <- data[!is.na(data[[group_column]]) & is.numeric(data[[group_column]]), ]
  
  # Find the minimum and maximum group values
  group_min <- floor(min(data[[group_column]], na.rm = TRUE))
  group_max <- max(data[[group_column]], na.rm = TRUE)  # Do not ceiling here to avoid exceeding max value
  
  # Create breaks for value intervals
  breaks <- seq(group_min, group_max, by = interval)
  
  # If the last break is less than the max group value, add an additional interval that ends at group_max
  if (tail(breaks, 1) < group_max) {
    breaks <- c(breaks, group_max)
  }
  
  # Assign each sample to a value range
  data[[geolabel]] <- cut(data[[group_column]], breaks = breaks, include.lowest = TRUE, right = FALSE,
                          labels = paste0(breaks[-length(breaks)], "_", breaks[-1]))
  
  # Adjust the last label to reflect the correct maximum value
  last_label <- paste0(breaks[length(breaks)-1], "_", group_max)
  data[[geolabel]] <- as.character(data[[geolabel]])
  data[[geolabel]][data[[geolabel]] == paste0(breaks[length(breaks)-1], "_", breaks[length(breaks)])] <- last_label
  
  return(data)
}

#Define function to group samples by grid, the group label has the header 'gridlabel' instead of default geolabel
add_grid_labels <- function(data, lon_column = "longitude", lat_column = "latitude", grid_size = 1) {
  # Find the minimum and maximum longitude and latitude
  lon_min <- min(data[[lon_column]], na.rm = TRUE)
  lon_max <- max(data[[lon_column]], na.rm = TRUE)
  lat_min <- min(data[[lat_column]], na.rm = TRUE)
  lat_max <- max(data[[lat_column]], na.rm = TRUE)
  
  # Generate sequences for grid lines based on grid size
  lon_seq <- seq(lon_min, lon_max, by = grid_size)
  lat_seq <- seq(lat_max, lat_min, by = -grid_size)
  
  # Initialize grid labels
  grid_labels <- expand.grid(lon_seq, lat_seq)
  grid_labels <- grid_labels[order(-grid_labels$Var2, grid_labels$Var1), ]
  grid_labels$label <- paste0("G", 1:nrow(grid_labels))
  
  # Assign each sample to a grid square
  data$gridlabel <- apply(data[, c(lon_column, lat_column)], 1, function(coord) {
    lon_index <- findInterval(coord[1], lon_seq)
    lat_index <- findInterval(-coord[2], -lat_seq)
    if (lon_index == 0 | lat_index == 0) {
      return(NA)  # Outside the grid
    }
    return(grid_labels$label[grid_labels$Var1 == lon_seq[lon_index] & grid_labels$Var2 == lat_seq[lat_index]])
  })
  
  return(data)
}

#Plot points on raster map based on grouping
plotgeolabels <- function(df, bioclimd, layer, geolabel = 'geolabel', plot_title = "Plot Title", 
                          x_label = "Longitude", y_label = "Latitude") {
  
  # Generate colors for each unique geolabel
  colors <- rainbow(length(unique(df[[geolabel]])))
  names(colors) <- unique(df[[geolabel]])
  
  # Access the specified layer dynamically using the provided layer name
  selected_layer <- bioclimd[[layer]]
  
  # Get the extent of the selected raster layer
  raster_extent <- extent(selected_layer)
  
  # Set up the plot without going beyond the raster extent
  par(mar = c(5, 4, 4, 2) + 0.1)  # Adjust margins
  plot(selected_layer, main = plot_title, xlab = x_label, ylab = y_label,
       xlim = c(raster_extent@xmin, raster_extent@xmax),
       ylim = c(raster_extent@ymin, raster_extent@ymax))
  
  # Plot the spatial points with colors based on the geolabel column
  points(df$decimalLongitude, df$decimalLatitude, 
         col = colors[df[[geolabel]]], 
         pch = 18,   # pch is the plotting symbol
         cex = 0.5)  # cex controls the size of the points 
}

#Functions to plot density distribution of bioclim variables based on groups
#with Mean and SD on distribution peaks, good for differentiated groups
plot_bioclim <- function(data, value_column, color_column, plot_type = c("density", "histogram"), 
                         title = value_column, x_label = "Values", 
                         y_label = "Density", legend_title = "Group") {
  plot_type <- match.arg(plot_type)
  
  # Convert the color_column to a factor to ensure proper grouping and coloring
  data[[color_column]] <- as.factor(data[[color_column]])
  
  # Calculate mean and standard deviation for each group
  stats <- data %>%
    group_by(!!sym(color_column)) %>%
    summarise(mean = mean(!!sym(value_column), na.rm = TRUE),
              sd = sd(!!sym(value_column), na.rm = TRUE),
              max_density = ifelse(plot_type == "density",
                                   max(density(!!sym(value_column), na.rm = TRUE)$y),
                                   NA))
  
  # Create the plot based on the selected plot_type
  p <- ggplot(data, aes_string(x = value_column, fill = color_column)) +
    labs(title = title, x = x_label, y = y_label, fill = legend_title) +
    theme_minimal()
  
  if (plot_type == "density") {
    p <- p + geom_density(alpha = 0.5)
    
    # Add centered text annotations within a white box with 75% opacity for mean and standard deviation
    p <- p + geom_label(data = stats, aes(x = mean, y = max_density * 0.95, 
                                          label = paste0("Mean: ", round(mean, 2), "\nSD: ", round(sd, 2))),
                        vjust = 0.5, hjust = 0.5, size = 3, color = "black", fill = "white", alpha = 0.75)
  } else if (plot_type == "histogram") {
    p <- p + geom_histogram(position = "identity", alpha = 0.5, bins = 30)
    
    # Add centered text annotations within a white box with 75% opacity for mean and standard deviation
    stats <- stats %>%
      mutate(max_density = max(data[[value_column]], na.rm = TRUE) * 1.05)
    
    p <- p + geom_label(data = stats, aes(x = mean, y = max_density * 0.95, 
                                          label = paste0("Mean: ", round(mean, 2), "\nSD: ", round(sd, 2))),
                        vjust = 0.5, hjust = 0.5, size = 3, color = "black", fill = "white", alpha = 0.75)
  }
  
  print(p)
}

#with Mean and SD next to group legend, good when most groups overlap
plot_bioclim2 <- function(data, value_column, color_column, plot_type = c("density", "histogram"), 
                                           title = value_column, x_label = "Values", y_label = "Density", 
                                           legend_title = "Group") {
  plot_type <- match.arg(plot_type)
  
  # Convert the color_column to a factor to ensure proper grouping and coloring
  data[[color_column]] <- as.factor(data[[color_column]])
  
  # Calculate mean and standard deviation for each group
  stats <- data %>%
    group_by(!!sym(color_column)) %>%
    summarise(mean = mean(!!sym(value_column), na.rm = TRUE),
              sd = sd(!!sym(value_column), na.rm = TRUE))
  
  # Create a new legend label with mean and SD included
  stats <- stats %>%
    mutate(label = paste0(!!sym(color_column), " (Mean: ", round(mean, 2), ", SD: ", round(sd, 2), ")"))
  
  # Replace the color column with the new label
  data[[color_column]] <- factor(data[[color_column]], labels = stats$label)
  
  # Create the plot based on the selected plot_type
  p <- ggplot(data, aes_string(x = value_column, fill = color_column)) +
    labs(title = title, x = x_label, y = y_label, fill = legend_title) +
    theme_minimal() +
    theme(legend.position = c(0.95, 0.95),              # Position legend inside plot area (top right)
          legend.justification = c("right", "top"),     # Align legend to the top right of its position
          legend.box.just = "right",                    # Justify legend box to the right
          legend.background = element_rect(fill = alpha("white", 0.75)),  # Semi-transparent legend background
          legend.direction = "vertical")                # Arrange legend items vertically
  
  if (plot_type == "density") {
    p <- p + geom_density(alpha = 0.5)
  } else if (plot_type == "histogram") {
    p <- p + geom_histogram(position = "identity", alpha = 0.5, bins = 30)
  }
  
  print(p)
}


