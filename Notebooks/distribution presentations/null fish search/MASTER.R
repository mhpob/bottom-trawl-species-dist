library(rmarkdown)
library(concaveman); library(sf); library(data.table)
library(ggplot2); library(patchwork)

all_data <- fread('data derived/survdat_names_sed.csv')

## Drop length
all_data[, c('length', 'numlen') := NULL]

## Drop redundant
all_data <- unique(all_data, by = names(all_data)[1:21])

station_key <- unique(all_data,
                      by = c('cruise6', 'station', 'stratum', 'tow',
                             'year', 'season'))
station_key[, 18:23 := NULL]

all_data <- all_data[station_key, on = names(station_key)]

all_data[, abundance := fifelse(is.na(abundance), 0, abundance)]

all_data <- st_as_sf(all_data,
                     coords = c('lon', 'lat'),
                     remove = F,
                     crs = 4326) %>%
  st_transform(crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')

all_data <- setDT(all_data)
# [
#   complete.cases(species_data[, .(season, bottemp, botsalin,
#                                   grpsed, year, depth,
#                                   lat, lon)])]

all_data <- all_data[abundance > 0 & year >= 1967]
all_data[, ':='(X = st_coordinates(geometry)[, 1],
                Y = st_coordinates(geometry)[, 2],
                season = factor(season,
                                ordered = T,
                                levels = c('SPRING', 'FALL')))]

# target_species <- c('horseshoe crab',
#                     'smallmouth flounder',
#                     'shortnose greeneye',
#                     'planehead filefish',
#                     'moustache sculpin',
#                     'lumpfish',
#                     'longnose greeneye',
#                     'fawn cusk-eel',
#                     'conger eel')

target_species <- c('blue runner',
                    'chain dogfish',
                    'gulf stream flounder')



all_data <- all_data[comname %in% target_species]




coast <- read_sf('data derived/mapping/coast_crop.shp') %>%
  st_transform(st_crs('+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74'))


strata_info <- fread('data/svdbs_supporttables/svdbs_svmstrata.csv')

strata <- read_sf('data/mapping/strata/strata.shp')
setDT(strata)

strata <- strata[strata_info, on = c('STRATA' = 'stratum')]
# strata <- strata[grepl('MAB|SNE', stratum_name)]
strata[, geometry :=
         st_transform(geometry, crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')]

strata <- strata[, .(geometry = st_union(geometry))]







grid <- st_make_grid(strata$geometry, cellsize = c(50*1000, 50*1000))

grid <- grid[st_covered_by(grid, strata$geometry, sparse = F) |
               st_overlaps(grid, strata$geometry, sparse = F)]

grid_dt <- data.table(geometry = grid,
                      id = seq_along(grid))

con_hull <- all_data[, concaveman(st_as_sf(geometry), 3, 1000),
                     by = c('comname', 'year', 'season')]

com <- all_data[, .(x_wt = sum(abundance * X)/sum(abundance),
               y_wt = sum(abundance * Y)/sum(abundance),
               dep_wt = sum(abundance * depth)/sum(abundance)),
           by = c('comname', 'year', 'season')]

com <- com[con_hull, on = c('comname', 'year', 'season')]

com <- st_as_sf(com, coords = c('x_wt', 'y_wt'),
                crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74',
                remove = F)

setDT(com)

for(i in 1:length(target_species)){
  spec <- all_data[comname == target_species[i]]

  spec_grid <- st_join(st_as_sf(grid_dt), st_as_sf(spec))


  abun_grid <- setDT(spec_grid)[, .(med_abun = median(abundance, na.rm = T),
                                    geometry = geometry),
                                by = c('id', 'year', 'season')]

  com_spec <-com[comname == target_species[i]]

  render(input = file.path("notebooks",
                           "distribution presentations",
                           "null fish search",
                           "template.Rmd"),
         params = list(species = target_species[i]),
         output_file = paste0(target_species[i], ".html"))
}
