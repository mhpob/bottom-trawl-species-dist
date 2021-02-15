library(sf)

coast <- read_sf('data/mapping/natural earth')

coast <- coast %>%
  st_crop(st_bbox(c(xmin = -83.5, ymin = 24, xmax = -59, ymax = 50)))

st_write(coast, 'data derived/mapping/coast_crop.shp')
