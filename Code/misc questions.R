# Task: how far apart are stations in a year?

library(data.table); library(sf)


all_data <- fread("data derived/survdat_mabsne_only.csv")

station_key <- unique(all_data,
                      by = c('cruise6', 'station', 'stratum', 'tow',
                             'year', 'season'))

station_key <- st_as_sf(station_key,
                         coords = c('lon', 'lat'),
                         remove = F,
                         crs = 4326) %>%
  st_transform(crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')


setDT(station_key)

k <- station_key[, st_distance(.SD$geometry),
                 by = c('year', 'season')]

View(k[V1 > units::set_units(0, 'm'),
  min(V1),
  by = c('year', 'season')])
