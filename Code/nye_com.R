library(data.table)
library(sf)
nye <- fread('data derived/nye_svdata.gz')


nye_sp <- st_as_sf(nye,
                   coords = c('lon', 'lat'),
                   crs = 4326) %>%
  # transform to Lambert projection to get distances
  st_transform(5072)

nye_trans <- data.table(nye_sp)
coords <- st_coordinates(nye_trans$geometry)
nye_trans[, ':='(x = coords[, 1],
                 y = coords[, 2])]


nye_trans <- rbind(
  nye_trans[, lapply(.SD, function (.) sum(abundance * ., na.rm = T) /
                     sum(abundance, na.rm = T)),
    by = .(svspp, species, year, season), .SDcols = c('x', 'y', 'depth')],
  nye_trans[, lapply(.SD, function (.) sum(biomass * ., na.rm = T) /
                       sum(biomass, na.rm = T)),
    by = .(svspp, species, year, season), .SDcols = c('x', 'y', 'depth')]
)

nye_trans <- data.table(
  nye_trans,
  com_type = rep(c('abundance', 'biomass'), each  = nrow(nye_trans) / 2)
)

nye_trans <- st_as_sf(nye_trans[complete.cases(nye_trans)],
                      coords = c('x', 'y'),
                      remove = F,
                      crs = 5072) %>%
  st_transform(crs = 4326)

coords <- st_coordinates(nye_trans)
nye_trans$lon <- coords[, 1]
nye_trans$lat <- coords[, 2]

nye_out <- data.table(nye_trans)[, c(2, 1, 3:4, 8, 5:7, 10:11)]

fwrite(nye_out, 'data derived/nye_COM.csv')


coast <- read_sf('data/mapping/natural earth')
%>%
  st_transform(5072)




library(ggplot2)

ggplot() +
  geom_sf(data = coast) +
  geom_sf(data = nye_trans[nye_trans$svspp == 15,],
          aes(color = year, shape = season)) +
  geom_path(data = nye_trans[nye_trans$svspp == 15,],
            aes(x = x, y = y, color = year, group = season)) +
  # coord_sf(xlim = c(1.6e6, 2.4e6), ylim = c(1848755, 2527557)) +
  coord_sf(xlim = c(-76.5, -69), ylim = c(35, 45)) +
  labs(x = NULL, y = NULL, color = 'Year', shape = 'Season') +
  scale_color_viridis_c() +
  facet_wrap(~ type) +
  theme_bw()

k <- melt(data.table(nye_trans)[, c(1:3, 6:7, 9:10)],
          id.vars = c('svspp', 'year', 'season', 'type'))
k[variable == 'depth', value := -value]


ggplot(data = k[svspp == 15]) +
  geom_line(aes(x = year, y = value, color = season)) +
  facet_grid(variable ~ type, scales ='free_y')
