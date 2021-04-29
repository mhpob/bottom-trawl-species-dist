library(sf); library(data.table)


demersal <- fread('data derived/demersal.csv')

# Some abundances (scallops) are negative.
#   Assuming that was an internal note and should be positive.
demersal[, abundance := abs(abundance)]

# Some sexes/stages are listed independently. Sum all of these
demersal[, ':='(abundance = sum(abundance, na.rm = T),
                biomass = sum(biomass, na.rm = T)),
         by = c('cruise6', 'station', 'stratum', 'tow',
                'year', 'season', 'comname')]

demersal <- unique(demersal, by = c('cruise6', 'station', 'stratum', 'tow',
                                    'year', 'season', 'comname'))

## Drop length
demersal[, c('length', 'numlen') := NULL]


station_key <- unique(demersal,
                      by = c('cruise6', 'station', 'stratum', 'tow',
                             'year', 'season'))
station_key[, 18:23 := NULL]

demersal <- demersal[station_key, on = names(station_key)]

demersal[, abundance := fifelse(is.na(abundance), 0, abundance)]

demersal <- st_as_sf(demersal,
                     coords = c('lon', 'lat'),
                     remove = F,
                     crs = 4326) %>%
  st_transform(crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')

demersal <- setDT(demersal)


demersal <- demersal[abundance > 0]
demersal[, ':='(X = st_coordinates(geometry)[, 1],
                Y = st_coordinates(geometry)[, 2],
                season = factor(season,
                                ordered = T,
                                levels = c('SPRING', 'FALL')))]


# Select trawls after 1967
com <- demersal[year >= 1967]

com <- com[, .(x = sum(abundance * X)/sum(abundance),
               y = sum(abundance * Y)/sum(abundance),
               depth = sum(abundance * depth)/sum(abundance)),
           by = c('year', 'season', 'comname')]

# Select only species where 75% of COM is in MAB (Hatteras to Nantucket)
target_spp <- st_as_sf(com,
                       coords = c('x', 'y'),
                       crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74',
                       remove = F) %>%
  st_transform(4326)
setDT(target_spp)[, ':='(lon = st_coordinates(geometry)[, 1],
                  lat = st_coordinates(geometry)[, 2])]


target_spp <- target_spp[lat %between% c(35.2, 41.2)]

target_spp <- target_spp[, .N, by = .(comname, season)]
target_spp <- unique(target_spp[N > (2019-1967) * 0.75]$comname)

com <- com[comname %in% target_spp]

com <- melt(com, measure.vars = c('x', 'y', 'depth'))


# Get slope and slope p-value
## {} allows creation of an intermediate object
com <- com[, {
  tmp = coef(summary(lm(value ~ year, data = .SD)))
  list(slope = tmp[2],
       slope_p = tmp[8],
       n = .N)
}, by = c('comname', 'variable', 'season')]


# Only select those with >10 observations per yr x season combo
com <- com[n > 10]
com <- dcast(com, comname + season + n ~ variable, value.var = c('slope', 'slope_p'))

com <- com[order(comname, season)]
# com <- com[order(-abs(slope_y_wt)), .SD, by = comname]

com <- com[, .(comname, season, n, slope_y, slope_p_y, slope_x,
               slope_p_x, slope_depth, slope_p_depth)]

fwrite(com, 'data derived/COM slopes_dem75pctincl.csv')

# com <- st_as_sf(com, coords = c('x_wt', 'y_wt'),
                # crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74',
                # remove = F)