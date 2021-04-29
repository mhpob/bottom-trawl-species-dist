library(rmarkdown)
library(concaveman); library(sf); library(data.table)
library(ggplot2); library(patchwork)

# Import station key
station_key <- fread('data derived/survdat_names_sed.csv')
station_key <- unique(station_key,
                      by = c('cruise6', 'station', 'stratum', 'tow',
                             'year', 'season'))
station_key[, 18:25 := NULL]

# IMport demersal species selected via inclusion rule in "com comparison table.R"
demersal_spp <- fread('data derived/demersal.csv')[
  comname %in%
    fread('data derived/COM slopes_dem75pctincl.csv',
          select = 'comname')$comname]


# Some species (scallops) have negative abundance
#   Assume this is a code and the abundance should be positive
demersal_spp[, abundance := abs(abundance)]

# Sum across sex/age classes
demersal_spp[, ':='(abundance = sum(abundance, na.rm = T),
                    biomass = sum(biomass, na.rm = T)),
             by = c('cruise6', 'station', 'stratum', 'tow',
                    'year', 'season', 'comname')]

# Pick unique tows
demersal_spp <- unique(demersal_spp,
                       by = c('cruise6', 'station', 'stratum', 'tow',
                              'year', 'season', 'comname'))


# Add in all tows where catch was 0 for each species
demersal_spp <- demersal_spp[, .SD[station_key, on = names(station_key)], by = 'comname']
demersal_spp[, ':='(abundance = fifelse(is.na(abundance), 0, abundance),
                    biomass = fifelse(is.na(biomass), 0, biomass))]


demersal_spp <- st_as_sf(demersal_spp,
                     coords = c('lon', 'lat'),
                     remove = F,
                     crs = 4326) %>%
  st_transform(crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')

demersal_spp <- setDT(demersal_spp)


demersal_spp <- demersal_spp[abundance > 0 & year >= 1967]
demersal_spp[, ':='(X = st_coordinates(geometry)[, 1],
                Y = st_coordinates(geometry)[, 2],
                season = factor(season,
                                ordered = T,
                                levels = c('SPRING', 'FALL')))]


unique_species <- unique(demersal_spp$comname)

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



demersal_spp[, dist := st_distance(geometry, st_union(coast))]



grid <- st_make_grid(strata$geometry, cellsize = c(50*1000, 50*1000))

grid <- grid[st_covered_by(grid, strata$geometry, sparse = F) |
               st_overlaps(grid, strata$geometry, sparse = F)]

grid_dt <- data.table(geometry = grid,
                      id = seq_along(grid))

con_hull <- demersal_spp[, concaveman(st_as_sf(geometry), 3, 1000),
                         by = c('comname', 'year', 'season')]

com <- demersal_spp[, .(x_wt = sum(abundance * X, na.rm = T)/sum(abundance, na.rm = T),
                        y_wt = sum(abundance * Y, na.rm = T)/sum(abundance, na.rm = T),
                        dist_wt = sum(abundance * as.numeric(dist), na.rm = T)/sum(abundance, na.rm = T),
                        dep_wt = sum(abundance * depth, na.rm = T)/sum(abundance, na.rm = T),
                        bt_wt = sum(abundance * bottemp, na.rm = T)/sum(abundance, na.rm = T)),
                    by = c('comname', 'year', 'season')]



com <- st_as_sf(com, coords = c('x_wt', 'y_wt'),
                crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74',
                remove = F)

setDT(com)

com[, com_dist_line := st_nearest_points(geometry, st_union(coast))]
com[, com_dist := as.numeric(st_length(com_dist_line))]

com <- com[con_hull, on = c('comname', 'year', 'season')]

for(i in 1:length(unique_species)){
  spec <- demersal_spp[comname == unique_species[i]]

  spec_grid <- st_join(st_as_sf(grid_dt), st_as_sf(spec))


  abun_grid <- setDT(spec_grid)[, .(med_abun = median(abundance, na.rm = T),
                                    geometry = geometry),
                                by = c('id', 'year', 'season')]

  com_spec <-com[comname == unique_species[i]]

  render(input = file.path("notebooks",
                           "distribution presentations",
                           "com-selected species",
                           "template.Rmd"),
         params = list(species = unique_species[i]),
         output_file = paste0(unique_species[i], ".html"))
}
