library(sdm)

library(sf); library(data.table)


strata_info <- fread('data/svdbs_supporttables/svdbs_svmstrata.csv')


# subset strata to midatl/sne
strata <- read_sf('data/mapping/strata/strata.shp')
setDT(strata)

strata <- strata[strata_info, on = c('STRATA' = 'stratum')]
strata <- strata[grepl('MAB|SNE', stratum_name)]
strata[, geometry := st_transform(geometry, crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')]

strata <- strata[, .(geometry = st_union(geometry))]



survdat <- fread('data derived/survdat_names_sed.csv')
subs <- survdat[comname == 'atlantic croaker']
subs <- unique(subs, by = c('cruise6', 'station', 'stratum', 'tow'))
subs[, yr_fac := as.factor(year)]
subs <- st_as_sf(subs, coords = c('lon', 'lat'), crs = 4326)
subs <- setDT(subs)[, geometry := st_transform(geometry, crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')]
subs <- cbind(subs, st_coordinates(subs$geometry))


subs <- subs[as.vector(st_intersects(geometry, strata$geometry, sparse = F))]

train <- subs[year < 2014]


k <- sdmData(
  formula = abundance ~ season + bottemp + coords(X+Y) + g(yr_fac),
  train = subs[year < 2014]
)


m1 <- sdm(abundance ~ f(season) + bottemp + coords(X+Y) + g(yr_fac),
          data = k,
          method = 'glm')
?sf::as_Spatial
