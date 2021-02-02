library(data.table)

list.files('data/22560_FSCSTables')

svcat <- fread(
  grep('SVCAT', list.files('data/22560_FSCSTables', full.names = T), value = T),
  col.names = tolower, fill = T
  )


svcat[!grepl('^\\d', svspp), ':='(logged_species_name = paste(logged_species_name, svspp, catchsex),
                                  svspp = expcatchnum,
                                  catchsex = expcatchwt,
                                  expcatchnum = v11,
                                  expcatchwt = v12)]
svcat[, c('v11', 'v12') := NULL]


svsta <- fread(
  grep('SVSTA', list.files('data/22560_FSCSTables', full.names = T), value = T),
  col.names = tolower)



subs <- svcat[grepl('alewife', logged_species_name, ignore.case = T)]
subs <- svsta[subs, on = c('cruise6', 'stratum', 'tow', 'station', 'id')]




library(sf)
# stations <- svsta[!is.na(decdeg_endlon)][
#   , {
#     geometry <- sf::st_linestring(x = matrix(c(decdeg_beglon, decdeg_endlon, decdeg_beglat, decdeg_endlat), ncol = 2))
#     geometry <- sf::st_sfc(geometry)
#     geometry <- sf::st_sf(geometry = geometry)
#   },
#   by = id
# ]
#
# stations <- st_as_sf(stations,
#                      crs = 4326)

stations <- st_as_sf(svsta,
                     coords = c('decdeg_beglon', 'decdeg_beglat'),
                     remove = F,
                     crs = 4326)


subs <- st_as_sf(subs,
                     coords = c('decdeg_beglon', 'decdeg_beglat'),
                     remove = F,
                     crs = 4326)


coast <- read_sf('data/mapping/natural earth',
                 wkt_filter = st_as_text(st_as_sfc(st_bbox(mack))))

library(ggplot2)

ggplot() +
  geom_sf(data = coast) +
  geom_sf(data = stations) +
  coord_sf(xlim = c(-82, -63), ylim = c(34, 48)) +
  facet_wrap(~est_year, nrow = 5)
