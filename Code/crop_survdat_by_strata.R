library(sf); library(data.table)


all_data <- fread('data derived/survdat_names_sed.csv')

## Drop length
all_data[, c('length', 'numlen') := NULL]

## Drop redundant
all_data <- unique(all_data, by = names(all_data)[1:21])

## Turn into spatial and crop
all_data <- st_as_sf(all_data,
                     coords = c('lon', 'lat'),
                     remove = F,
                     crs = 4326) %>%
  st_transform(crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')



strata_info <- fread('data/svdbs_supporttables/svdbs_svmstrata.csv')

strata <- read_sf('data/mapping/strata/strata.shp')
setDT(strata)

strata <- strata[strata_info, on = c('STRATA' = 'stratum')]
strata <- strata[grepl('MAB|SNE', stratum_name)]
strata[, geometry :=
         st_transform(geometry, crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')]

strata <- strata[, .(geometry = st_union(geometry))]

# Buffer by 5k
strata[, geometry := st_buffer(geometry, 5000)]

setDT(all_data)
all_data <- all_data[as.vector(st_intersects(geometry,
                                             strata$geometry,
                                             sparse = F))]


fwrite(all_data, 'data derived/survat_mabsne_only.csv')
