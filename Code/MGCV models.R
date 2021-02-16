library(data.table)

survdat <- fread('data derived/survdat_names_sed.gz')

cruise_key <- unique(survdat[, cruise6:botsalin])
cruise_key <- cruise_key[lat > 35]
setkeyv(cruise_key, names(cruise_key))

crk <- survdat[comname == 'atlantic croaker' &
                 lat > 35]
crk <- unique(crk, by = c('cruise6', 'station', 'stratum', 'tow'))
setkeyv(crk, names(cruise_key))


crk <- crk[cruise_key]
crk[, c('abundance', 'biomass') := lapply(.SD, function(.) fifelse(is.na(.), 0, .)),
    .SDcols = c('abundance', 'biomass')]

crk <- crk[, .(year, season, lat, lon, depth, grpsed, bottemp, abundance, biomass)]
crk <- crk[complete.cases(crk)]

crk[, ':='(season = as.factor(season),
           yr_fac = as.factor(year))]

library(mgcv)


m <- bam(abundance ~
           # variation in abundance due to sampling season
           season +

           # long-term trend in abundance
           # s(year, bs = 'cs', k = 45) +

           #spatial variation in abundance
           s(lon, lat, bs = 'ds', k = 100, m = c(1, 0.5)) +

           # variation in abundance due to bottom temperature
           s(bottemp, k = 15) +

           # how spatial trend varies over time
           # ti(lon, lat, yr_fac, bs = c('ds', 're'), d = c(2, 1),
           #    m = list(c(1, 0.5), NA),
           #    k = c(20, 10)) +

           s(yr_fac, bs = 're')
         ,
         data = crk[year <= 2008],
         discrete = T,
         samfrac = 0.1,
         family = 'tw')


new_data <- crk[year > 2008]
new_data[, yr_fac := '2000']

pred <- predict(m, new_data, exclude = 's(yr_fac)', type = 'link')
pred <- as.data.table(pred)

new_data <- data.table(
  new_data,
  pred
)

new_data[, ':='(pred = exp(pred))]
new_data[, err := abundance - pred]

library(ggplot2)

ggplot(data = new_data) +
  geom_point(aes(x = lon, y = lat, size = fit)) +
  scale_color_viridis_c() +
  facet_grid(season ~ year)


# RMSE of zeros
new_data[abundance == 0,
         round(sqrt( mean( (abundance - pred) ^ 2)), 0), by = c('year', 'season')]

#RMSE of >0
new_data[abundance > 0,
         round(sqrt( mean( (abundance - pred) ^ 2)), 0), by = c('year', 'season')]









k <- gam(
  list(
    abundance ~ s(bottemp) + s(yr_fac, bs = 're'),
    ~ season + s(lon, lat, bs = 'ds', k = 100, m = c(1, 0.5))
  ),
  family = 'ziplss',
  data = crk,
  method = 'REML',
  control = gam.control(nthreads = 3)
)


new_data <- crk[year > 2008]

pred <- predict(k, new_data, exclude = 's(yr_fac)', type = 'link')
pred <- as.data.table(pred)

new_data <- data.table(
  new_data,
  pred
)

new_data[, ':='(fit =exp(V1) * binomial(link = "cloglog")$linkinv(V2))]
new_data[, err := abundance - fit]

library(ggplot2)

ggplot(data = new_data) +
  geom_point(aes(x = lon, y = lat, size = fit)) +
  scale_color_viridis_c() +
  facet_grid(season ~ year)

# RMSE of zeros
new_data[abundance == 0,
         round(sqrt( mean( (abundance - fit) ^ 2)), 0), by = c('year', 'season')]

#RMSE of >0
new_data[abundance > 0,
         round(sqrt( mean( (abundance - fit) ^ 2)), 0), by = c('year', 'season')]
