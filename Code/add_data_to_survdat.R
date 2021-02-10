library(data.table)

# Load and manipulate NOAA EDAB's compiled survey data ----
load('data/noaa_edab_survdat.rdata')

setnames(survdat, tolower)


# Add species names ----
species_list <- fread(
  input = 'data/svdbs_supporttables/svdbs_svspecies_list.csv',
  col.names = tolower,
  fill = T
)

## Quick data cleaning
species_list[!grepl('^\\d', svspp), ':='(
  comname = paste(comname, svspp, v1),
  svspp = v2)]
species_list[, c('v1', 'v2') := NULL]
species_list[, svspp := as.numeric(svspp)]
species_list[, c('sciname', 'comname') :=
               lapply(.SD, tolower), .SDcols = c('sciname', 'comname')]

## Join species names
### Note that 726 observations don't have an associated species name
survdat <- species_list[survdat, on = 'svspp']




# Add sediment data ----
## Unique info by trawl
key <- unique(survdat, by = c('cruise6', 'station', 'stratum', 'tow'))

## Pick needed columns
key <- key[, .(cruise6, station, stratum, tow, season, lat, lon)]


## Load bottom sediment type
library(sf)

## Other shapefiles -- not used at the moment
# bot_mab <- st_read('data/mapping/benthic_habitats/benhab_mab.shp',
#                    query = 'SELECT
#                               SEDIMENT as sediment,
#                               ECOLOGICAL as ecological
#                             FROM benhab_mab')
#
# bot_sne <- st_read('data/mapping/benthic_habitats/benhab_sne.shp',
#                    query = 'SELECT
#                               SEDIMENT as sediment,
#                               ECOLOGICAL as ecological
#                             FROM benhab_sne')
#
# bottom_type <- rbind(bot_mab, bot_sne)
# bottom_type <- st_make_valid(bottom_type)
#
# st_write(bottom_type, 'data derived/compiled_bottom_type.shp')

bottom_type <- st_read('data/mapping/tncbenthic/benthic.gdb',
                       layer = 'benthic_sediment')

key <- st_as_sf(key,
                coords = c('lon', 'lat'),
                crs = 4326)
key <- st_transform(key,
                    st_crs(bottom_type))



key <- st_join(key, bottom_type, join = st_intersects)
key <- st_transform(key, 4326)
key$lat <- st_coordinates(key)[, 2]
key$lon <- st_coordinates(key)[, 1]

key <- data.table(key)
key <- key[, !c('ID', 'GRIDCODE', 'GLOBALID', 'SHAPE_Length', 'SHAPE_Area',
                'geometry', 'lat', 'lon')]



survdat <- survdat[key, on = c('cruise6', 'station', 'stratum', 'tow', 'season')]



# Get ready to export ----
setnames(survdat, tolower)

survdat <- survdat[, .(cruise6, station, stratum, tow, svvessel, year, season, lat, lon,
                       est_towdate, depth, sediment, grpsed, surftemp, surfsalin, bottemp,
                       botsalin, sciname, comname, svspp, catchsex, abundance, biomass,
                       length, numlen)]

fwrite(survdat, 'data derived/survdat_names_sed.gz')
