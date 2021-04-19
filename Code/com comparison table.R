library(sf); library(data.table)

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


all_data <- all_data[abundance > 0 & !grepl('(^| )un[ck]', comname)]
all_data[, ':='(X = st_coordinates(geometry)[, 1],
                Y = st_coordinates(geometry)[, 2],
                season = factor(season,
                                ordered = T,
                                levels = c('SPRING', 'FALL')))]


# target <- all_data[comname %in% c('atlantic croaker',
#                                   'summer flounder',
#                                   'black sea bass',
#                                   'scup',
#                                   'little skate',
#                                   'spotted hake',
#                                   'northern searobin',
#                                   'silver hake',
#                                   'weakfish',
#                                   'smooth dogfish',
#                                   'red hake',
#                                   'winter skate',
#                                   'bluntnose stingray',
#                                   'american lobster')]

# Select trawls after 1967 that have a species listed.
com <- all_data[year >= 1967 & comname != ""]

# Select only species that have been caught off Delmarva at least once
com <- com[comname %in% unique(com[between(lat, 37, 39), comname])]

com <- com[, .(x = sum(abundance * X)/sum(abundance),
               y = sum(abundance * Y)/sum(abundance),
               depth = sum(abundance * depth)/sum(abundance)),
           by = c('year', 'season', 'comname')]


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

fwrite(com, 'data derived/COM slopes.csv')

# com <- st_as_sf(com, coords = c('x_wt', 'y_wt'),
                # crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74',
                # remove = F)