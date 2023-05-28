# Obtaining and processing spatial observation data of the species. 

occurrences.get.data <- function(occurrences = "",
                                 species = "",
                                 range.shapefile = "",
                                 from = c('gbif', 'vertnet', 'inat', 'idigbio'),
                                 limit = 1,
                                 has_coords = TRUE,
                                 return.occurrences.raw = FALSE,
                                 return.occurrences.filtered = FALSE) {
  
  if(species == ""){
    stop("*Error: species is not defined.")
  }
  
  if(range.shapefile == ""){
    stop("*Error: range.shapefile is not defined.")
  }
  
  if(return.occurrences.raw == TRUE & return.occurrences.filtered == TRUE) {
    stop("*Error: please select only one of return.occurrences.raw or return.occurrences.filtered.")
  }
  
  #species <- gsub(" ", "_", species)

  if(occurrences == ""){
    occurrences.raw <- spocc::occ(query = species, from = from, limit = limit, has_coords = has_coords)
    
    occurrences.raw <- data.frame(database = c(rep('GBIF', nrow(occurrences.raw$gbif$data[[1]])),
                                               rep('iNat', nrow(occurrences.raw$inat$data[[1]])),
                                               rep('VertNet', nrow(occurrences.raw$vertnet$data[[1]])),
                                               rep('iDigBio', nrow(occurrences.raw$idigbio$data[[1]]))),
                                  
                                  species = c(occurrences.raw$gbif$data[[1]]$scientificName, 
                                              occurrences.raw$inat$data[[1]]$name,
                                              occurrences.raw$vertnet$data[[1]]$name,
                                              occurrences.raw$idigbio$data[[1]]$dwc.scientificName),
                                  
                                  year = c(occurrences.raw$gbif$data[[1]]$year, 
                                           occurrences.raw$inat$data[[1]]$observed_on_details.year,
                                           occurrences.raw$vertnet$data[[1]]$year,
                                           occurrences.raw$idigbio$data[[1]]$dwc.year),
                                  
                                  longitude = c(occurrences.raw$gbif$data[[1]]$longitude, 
                                                occurrences.raw$inat$data[[1]]$longitude,
                                                occurrences.raw$vertnet$data[[1]]$longitude,
                                                occurrences.raw$idigbio$data[[1]]$longitude),
                                  
                                  latitude = c(occurrences.raw$gbif$data[[1]]$latitude, 
                                               occurrences.raw$inat$data[[1]]$latitude,
                                               occurrences.raw$vertnet$data[[1]]$latitude,
                                               occurrences.raw$idigbio$data[[1]]$latitude))
  } else {
    occurrences.raw <- readRDS(occurrences)
  }

  
  occurrences.raw$longitude <- stringr::str_replace_all(occurrences.raw$longitude, ",", ".")
  occurrences.raw$latitude <- stringr::str_replace_all(occurrences.raw$latitude, ",", ".")
  
  cat(paste("Number of raw occurrences:", nrow(occurrences.raw), "\n"))
  
  occurrences.filtered <- occurrences.raw[!duplicated(occurrences.raw),]
  cat(paste("Number of occurrences after elimination of duplicates:", nrow(occurrences.filtered), "\n"))
  
  occurrences.filtered <- occurrences.filtered[occurrences.filtered$longitude != 0 & occurrences.filtered$latitude != 0, ]
  cat(paste("Number of occurrences after excluding occurrences with zero values of longitude and latitude:", nrow(occurrences.filtered), "\n"))
  
  occurrences.filtered$year = substr(occurrences.filtered$year, 1, 4)
  occurrences.filtered <- occurrences.filtered[occurrences.filtered$year>=1960 & occurrences.filtered$year<=2023,]
  cat(paste("Number of occurrences after exclusion of occurrences collected earlier 1960-2023:", nrow(occurrences.filtered), "\n"))
  
  decimalplaces <- function(x) {
    if (abs(x - round(x)) > .Machine$double.eps^0.5) {
      nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }
  
  occurrences.filtered <- na.omit(occurrences.filtered)
  occurrences.filtered$longitude <- as.numeric(gsub("\\,", "", as.character(occurrences.filtered$longitude)))
  occurrences.filtered$latitude <- as.numeric(gsub("\\,", "", as.character(occurrences.filtered$latitude)))
  occurrences.filtered <- occurrences.filtered[sapply(occurrences.filtered$longitude, decimalplaces) > 2 & sapply(occurrences.filtered$latitude, decimalplaces) > 2, ]
  cat(paste("Number of occurrences after excluding occurrences with low level of accuracy:", nrow(occurrences.filtered), "\n"))
  
  occurrences.filtered$year <- as.numeric(gsub("\\,", "", as.character(occurrences.filtered$year)))
  
  range.map <- rgdal::readOGR(range.shapefile, 'eurasia')
  occurrences.inout <- data.frame(occurrences.filtered, inside = raster::extract(range.map, occurrences.filtered[, 4:5]))
  occurrences.inout$inside <- ifelse(!is.na(occurrences.inout$inside.ID_NO), TRUE, FALSE)
  occurrences.filtered <- occurrences.filtered[occurrences.inout$inside == TRUE, ]
  cat(paste("After removal of observations outside the range:", nrow(occurrences.filtered), "\n"))
  
  
  if(return.occurrences.raw == TRUE){
    return(occurrences.raw)
  }
  
  if(return.occurrences.filtered == TRUE){
    return(occurrences.filtered)
  }

}






stack.create.fun <- function(path, prefix, years, output, variable) {
  for(y in years){
    for(pr in prefix){
      prefixes <- as.factor(paste0(pr, y, '-'))
      for(p in prefixes){
        s <- stack(list.files(path = path, pattern = rep(p, 12), full.names=TRUE))
        saveRDS(s, paste0(output, y, '.Rds'))
        #writeRaster(s, filename = paste0(output, variable, '.', y, '.'), format = 'GTiff', overwrite=TRUE, progress = 'text')
      } 
    }
  }
}

biovars.dismo <- function(tmin.path, tmax.path, prec.path, years, output) {
  for(y in years){
    tmin <- readRDS(paste0(tmin.path, y, '.Rds'))
    tmax <- readRDS(paste0(tmax.path, y, '.Rds'))
    prec <- readRDS(paste0(prec.path, y, '.Rds'))
    variables <- as.factor(paste0('bio_', 1:19))
    bioclims <- biovars(prec = prec, tmin = tmin, tmax = tmax)
    saveRDS(bioclims, paste0(output, y, '.Rds'))
    for (i in 1:length(variables)) {
      writeRaster(bioclims[[i]], filename = paste0(output, y, '.', variables[i], '.'), format = 'GTiff', overwrite=TRUE, progress = 'text')
      writeRaster(bioclims[[i]], filename = paste0(output, y, '.', variables[i], '.'), format = 'ascii', overwrite=TRUE, progress = 'text')
    }
  }
}


# Усреднение
biostacks <- function(inputpath, outputpath, bioprefix, bionumb) {
  cat('\nInput path with 19 bioclimatic data for each specific year: ', inputpath, '\n')
  cat('Output path, where files with averaged bioclimatic data for all years will be create: ', outputpath, '\n')
  cat('Bioclimatic data prefix in inputpath: ', bioprefix, '\n')
  cat('Number of bioclimatic data: ', length(bionumb), '\n\n')
  for(n in bionumb){
      clim <- as.factor(paste0('.', bioprefix, n, '.tif'))
      for(c in clim){
        s <- stack(list.files(path = inputpath, pattern = c, full.names=TRUE))
        cat('Bioclimatic data averaging: ', '\n')
        print(s)
        s <- mean(s)
        cat('\nOutput file preparation: ', bioprefix, n, '.tif', '\n')
        writeRaster(s, filename = paste0(outputpath, bioprefix, n, '.'), format = 'GTiff', overwrite=TRUE, progress = 'text')
      }
  }
}


# Clim data preparation

clim.preparation <- function(path.tifs = '', prefix.tifs = '', path.range = '', prefix.range = '', path.output = '', return.tif_crop = FALSE, return.asc = FALSE, return.asc_crop = FALSE){

  if(return.tif_crop == TRUE & return.asc_crop == TRUE || return.tif_crop == TRUE & return.asc == TRUE || return.asc_crop == TRUE & return.asc == TRUE) {
    stop("*Error: please select only one of return.tif_crop, return.asc_crop or return.asc")
  }
  clim <- lapply(list(paste0(path.tifs, prefix.tifs, 1:19, '.tif')), function(file) stack(lapply(file, raster)))
  variables <- as.factor(paste0(prefix.tifs, 1:19))
  if (return.asc == TRUE || return.tif_crop == TRUE || return.asc_crop == TRUE) {
    if (return.tif_crop == TRUE || return.asc_crop == TRUE) {
      range.map <- readOGR(path.range, prefix.range)
      clim.crop <- stack(crop(clim[[1]], range.map))
      if (return.tif_crop == TRUE){
        for (i in 1:length(variables)) {
          writeRaster(clim.crop[[i]], filename = paste0(path.output, "crop_", variables[i], "."), format = "GTiff", overwrite=TRUE)
        }
      }
      if (return.asc_crop == TRUE) {
        for (i in 1:length(variables)) {
          writeRaster(clim.crop[[i]], filename = paste0(path.output, "crop_", variables[i], "."), format = "ascii", overwrite=TRUE)
        }
      }
    }
    if (return.asc == TRUE) {
      clim <- stack(clim[[1]])
      for (i in 1:length(variables)) {
        writeRaster(clim[[i]], filename = paste0(path.output, variables[i], "."), format = "ascii", overwrite=TRUE)
      }
    } 
  } else {
    stop("*Error: please select one of return.tif_crop, return.asc_crop or return.asc")
  }
}













# Plots

occurrences.stats.plot <- function(occurrences, plot.save = FALSE, path.output = '', output.prefix = ''){

  occ.points <- ggplot() +
    geom_point(data = occurrences, aes(x = longitude, y = latitude), pch='.', alpha = 0.5) +
    theme_classic() +
    coord_equal() +
    ggtitle(paste0('Number of occurrences: ', nrow(occurrences), '\n')) +
    theme(axis.line = element_line(colour = "grey90"), plot.title = element_text(size = 12, hjust = 0.0))
  
  occ.hist <- ggplot(occurrences, aes(x = year)) +
    stat_count() +
    geom_histogram(bins = length(unique(occurrences$year))) +
    theme_classic() +
    ggtitle(paste0('Collection year (', min(occurrences$year), '-', max(occurrences$year), ') ', nrow(occurrences), ' occurrences\n')) +
    labs(x = "Collection year", y = "Number of occurrences") +
    theme(axis.line = element_line(colour = "grey90"), plot.title = element_text(size = 12, hjust = 0.0))
  
  if(plot.save == TRUE){
    if(path.output == ''){
      stop("*Error: path.output is not defined")
    } else{
          output <- paste0(path.output, "occurrences.", output.prefix, min(occurrences$year), 'to', max(occurrences$year), '.png')
          ggsave(file = output, plot = grid.arrange(occ.points, occ.hist, nrow = 2), dpi = 600)
    }
  } else{
    grid.arrange(occ.points, occ.hist, nrow = 2)
  }
}














# Пример получения биоварс, который нашла Диана
setwd('C:/Users/User/Downloads/Downloads/phd_project/worldclim')
one_year_tmax <- list.files('tmax_2000_2009/output',pattern='2000.*\\.asc', full=TRUE)
one_year_tmin <- list.files('tmin_2000_2009/output',pattern='2000.*\\.asc', full=TRUE)
one_year_prec <- list.files('prec_2000_2009/output',pattern='2000.*\\.asc', full=TRUE)


ras_tmax <- lapply(one_year_tmax,raster) 
ras_tmin <- lapply(one_year_tmin,raster)
ras_prec <- lapply(one_year_prec,raster) 

tmax_stack <- stack(ras_tmax)
tmin_stack <- stack(ras_tmin)
prec_stack <- stack(ras_prec)

b <- biovars(prec_stack, tmin_stack, tmax_stack)

writeRaster(b[[1]], filename = "C:/Users/User/Downloads/Downloads/phd_project/worldclim/bioclimatic/bio1_2000.asc", format='ascii')





