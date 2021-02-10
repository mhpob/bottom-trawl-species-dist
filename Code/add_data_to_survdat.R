library(sf)
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

library(data.table)

load('data/noaa_edab_survdat.rdata')

setnames(survdat, tolower)

# Drop length
survdat[, c('length', 'numlen') := NULL]

# Unique info by trawl
survdat <- unique(survdat, by = c('cruise6', 'station', 'stratum', 'tow'))


sp_survdat <- st_as_sf(survdat,
                       coords = c('lon', 'lat'),
                       crs = 4326)
sp_survdat <- st_transform(sp_survdat,
                           st_crs(bottom_type))



sp_survdat <- st_join(sp_survdat, bottom_type, join = st_intersects)
sp_survdat <- st_transform(sp_survdat, 4326)
sp_survdat$lat <- st_coordinates(sp_survdat)[, 2]
sp_survdat$lon <- st_coordinates(sp_survdat)[, 1]
sp_survdat <- data.table(sp_survdat)
sp_survdat <- sp_survdat[, c(1:4, 20:21)]
load('data/noaa_edab_survdat.rdata')

setnames(survdat, tolower)


sp_survdat <- survdat[sp_survdat, on = c('cruise6', 'station', 'stratum', 'tow')]

fwrite(sp_survdat, 'data derived/survdat_w_sediment.gz')
