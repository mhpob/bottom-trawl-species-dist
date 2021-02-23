library(data.table); library(mgcv)
source("miller et al. supp/mgcv_spde_smooth.R") # note that this loads INLA



# load strata shp
strata_info <- fread('data/svdbs_supporttables/svdbs_svmstrata.csv')

library(sf)



# subset strata to midatl/sne
strata <- read_sf('data/mapping/strata/strata.shp')
setDT(strata)

strata <- strata[strata_info, on = c('STRATA' = 'stratum')]
strata <- strata[grepl('MAB|SNE', stratum_name)]
strata[, geometry := st_transform(geometry, crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')]

strata <- strata[, .(geometry = st_union(geometry))]

strata[, geometry := st_simplify(geometry, preserveTopology = F, dTolerance = 5000)]



boundary <- inla.sp2segment(as_Spatial(strata$geometry))
# make mesh using meshbuilder 'cause it's finnicky
# meshbuilder()

mesh <- inla.mesh.2d(boundary=boundary,
                     max.edge=c(22000, 108000),
                     min.angle=c(30, 21),
                     max.n=c(48000, 16000), ## Safeguard against large meshes.
                     max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                     cutoff=8000, ## Filter away adjacent points.
                     offset=c(1000, 55000))

plot(mesh, asp = 1)

#load species data

survdat <- fread('data derived/survdat_names_sed.csv')

station_key <- unique(survdat,
                              by = c('cruise6', 'station', 'stratum', 'tow'))
station_key[, 18:25 := NULL]




subs <- survdat[comname == 'atlantic croaker']
subs <- unique(subs, by = c('cruise6', 'station', 'stratum', 'tow'))

subs <- subs[station_key, on = names(station_key)]

subs <- subs[complete.cases(subs[, .(season, bottemp, year, lat, lon)])]
subs[, ':='(yr_fac = as.factor(year),
            abundance = fifelse(is.na(abundance), 0, abundance),
            biomass = fifelse(is.na(biomass), 0, biomass))]
subs <- st_as_sf(subs, coords = c('lon', 'lat'), crs = 4326)
subs <- setDT(subs)[, geometry := st_transform(geometry, crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')]
subs <- cbind(subs, st_coordinates(subs$geometry))


# select trawls within the boundary
subs <- subs[as.vector(st_intersects(geometry, strata$geometry, sparse = F))]


# Biomass
# mod_tw <- bam(biomass ~ season + s(bottemp, bs = 'cr) +
#              s(yr_fac, bs = 're') +
#              s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
#            data = subs[year < 2014], family = tw(),
#            control =  gam.control(scalePenalty = FALSE),
#            discrete = T, nthreads = 11)
#
# mod_nb <- bam(biomass ~ season + s(bottemp, bs = 'cr') +
#                 s(yr_fac, bs = 're') +
#                 s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
#               data = subs[year < 2014], family = nb(),
#               control =  gam.control(scalePenalty = FALSE),
#               discrete = T, nthreads = 11)

# Abundance
mod_tw_abun <- bam(abundance ~ season + grpsed +
                     s(bottemp, bs = 'cr', k = 10) +
                     s(year, bs = 'cr', k = 10) +
                     s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
                   data = subs[between(year, 1967, 2014)], family = tw(),
                   control =  gam.control(scalePenalty = FALSE),
                   discrete = T, nthreads = 11)

mod_poi_abun <- bam(abundance ~ season + grpsed +
                      s(bottemp, bs = 'cr', k = 10) +
                      s(year, bs = 'cr', k = 10) +
                      s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
                    data = subs[between(year, 1967, 2014)], family = poisson(),
                    control =  gam.control(scalePenalty = FALSE),
                    discrete = T, nthreads = 11)

mod_nb_abun <- bam(abundance ~ season + grpsed +
                     s(bottemp, bs = 'cr', k = 10) +
                     s(year, bs = 'cr', k = 10) +
                     s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
                   data = subs[between(year, 1967, 2014)], family = nb(),
                   control =  gam.control(scalePenalty = FALSE),
                   discrete = T, nthreads = 11)

mod_zip_abun <- bam(abundance ~ season + grpsed +
                      s(bottemp, bs = 'cr', k = 10) +
                      s(year, bs = 'cr', k = 10) +
                      s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
                    data = subs[between(year, 1967, 2014)], family = ziP(),
                    control =  gam.control(scalePenalty = FALSE),
                    discrete = T, nthreads = 11)

abun <- list(nb = mod_nb_abun,
             poi = mod_poi_abun,
             tw = mod_tw_abun,
             zip = mod_zip_abun)
#
saveRDS(abun, 'data derived/model output/bam_abun_bottemp_xy_yrran_with0_cryr.RDS')


abun <- readRDS('data derived/model output/bam_abun_bottemp_xy_yrran.RDS')
lapply(abun, AIC)

k <- predict(abun$nb, newdata = abun$nb$model, type = 'terms', se = T,
             exclude = 'season')

kk <- setDT(cbind(
  abun$nb$model,
  pred = exp(k$fit),
  lci = exp(k$fit - 1.96 * k$se.fit),
  uci = exp(k$fit + 1.96 * k$se.fit)
))


library(ggplot2)

ggplot(data = kk, aes(x = bottemp)) +
  geom_ribbon(aes(ymin = `lci.s(bottemp)`, ymax = `uci.s(bottemp)`),
              fill = 'lightgray') +
  geom_line(aes(y = `pred.s(bottemp)`)) +
  geom_rug() +
  labs(x = 'Bottom temperature (°C)', y = 'Individuals',
       title = 'Partial effect of bottom temperature') +
  theme_bw()


ggplot(data = kk, aes(x = yr_fac)) +
  geom_pointrange(aes(y = `pred.s(yr_fac)`,
                      ymin = `lci.s(yr_fac)`, ymax = `uci.s(yr_fac)`)) +
  geom_rug() +
  labs(x = NULL, y = 'Individuals',
       title = 'Partial effect of year; I.I.D. random intercept') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


ggplot(data = melt(kk,
                   id.vars = c('abundance','X','Y'),
                   measure.vars = c('pred.s(X,Y)', 'lci.s(X,Y)', 'uci.s(X,Y)')),
       aes(x = X, y = Y)) +
  geom_point(aes(color = value)) +
  scale_color_viridis_c() +
  facet_wrap(~variable)

coast <- st_read('data derived/mapping/coast_crop.shp') %>%
  st_transform(crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')

spat_dat <- melt(kk,
                 id.vars = c('abundance','X','Y'),
                 measure.vars = c('pred.s(X,Y)', 'lci.s(X,Y)', 'uci.s(X,Y)'))
spat_dat[, variable := fcase(variable %like% 'pred', 'Predicted',
                             variable %like% 'lci', 'Lower confidence interval',
                             variable %like% 'uci', 'Upper confidence interval')]


ggplot() +

  geom_point(data = setorder(spat_dat, value),
             aes(x = X, y = Y, color = value)) +
  geom_sf(data = coast) +
  coord_sf(xlim = c(-2e5, 4.2e5), ylim = c(3.6e6, 4.38e6)) +
  scale_color_viridis_c(trans = 'log10') +
  labs(x = NULL, y = NULL, color = NULL,
       title = 'Atlantic croaker modeled abundance',
       subtitle = 'Partial effect of X-Y space') +
  facet_wrap(~variable) +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.4))



grid <- st_make_grid(strata$geometry, cellsize = 10*1000)

grid <- grid[st_covered_by(grid, strata$geometry, sparse = F) |
                st_overlaps(grid, strata$geometry, sparse = F)]

grid <- st_centroid(grid)
grid <- data.table(X = st_coordinates(grid)[,1],
                   Y = st_coordinates(grid)[,2],
                   season = 'fall')

spat_pred <- predict(mod_nb_abun, cbind(grid, season = c('FALL', 'SPRING'),
                                        grpsed = c('Sand', 'Gravel'),
                                        bottemp = 20, year = 1963),
                     type = 'terms', terms = 's(X,Y)', se = T,
                     newdata.guaranteed = T)
spat_pred <- data.table(
  data.table(grid),
  pred = exp(spat_pred$fit),
  lci = exp(spat_pred$fit - 1.96 * spat_pred$se.fit),
  uci = exp(spat_pred$fit + 1.96 * spat_pred$se.fit)
)

bla <- cbind(spat_pred, grid)

ggplot() +
  # geom_sf(data = bla, aes(geometry = grid, fill  = `pred.s(X,Y)`)) +
  geom_contour_filled(data = bla, aes(x = X, y = Y, z  = `pred.s(X,Y)`)) +
  scale_fill_viridis_c(trans = 'log10')



## Cross validation ----
### Model fits

cv_5y <- lapply(2015:2019,
                function(.){
                  tryCatch({
                    bam(formula = abundance ~ season + grpsed +
                          s(bottemp, bs = "cr", k = 10) +
                          s(year, bs = "cr", k = 10) +
                          s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
                        family = nb(), data = croaker[between(year, 1967, .)],
                        control = gam.control(scalePenalty = FALSE),
                        discrete = T, nthreads = 11)},
                    error = function(e) NULL
                  )
                })

# cv_5y <- c(abun$nb, cv_5y)

names(cv_5y) <- c('m15', 'm16', 'm17', 'm18', 'm19')




saveRDS(cv_5y,  file = 'data derived/model output/bam_abun_bottemp_xy_yr_seas_sed_cvyrly.RDS')


pred15 <- predict(abun$nb, newdata = croaker[year == 2015],
                  type = 'link',
                  exclude = c('s(year)', 's(X,Y)'))

pred15 <- setDT(
  cbind(
    croaker[year == 2015],
    pred = exp(pred15)
  )
)


cv <- pred15[, .(rmse = sqrt(mean((abundance - pred)^2)),
                 prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
                 year = 'yr2015'),
       by = abundance > 0]


pred16 <- predict(cv_5y$m15, newdata = croaker[year == 2016],
                  type = 'link',
                  exclude = c('s(year)', 's(X,Y)'))

pred16 <- setDT(
  cbind(
    croaker[year == 2016],
    pred = exp(pred16)
  )
)

cv <- rbind(cv,
            pred16[, .(rmse = sqrt(mean((abundance - pred)^2)),
                       prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
                       year = 'yr2016'),
                   by = abundance > 0]
)

pred17 <- predict(cv_5y$m16, newdata = croaker[year == 2017],
                  type = 'link',
                  exclude = c('s(year)', 's(X,Y)'))

pred17 <- setDT(
  cbind(
    croaker[year == 2017],
    pred = exp(pred17)
  )
)

cv <- rbind(cv,
            pred17[, .(rmse = sqrt(mean((abundance - pred)^2)),
                       prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
                       year = 'yr2017'),
                   by = abundance > 0]
)



pred18 <- predict(cv_5y$m17, newdata = croaker[year == 2018],
                  type = 'link',
                  exclude = c('s(year)', 's(X,Y)'))

pred18 <- setDT(
  cbind(
    croaker[year == 2018],
    pred = exp(pred18)
  )
)

cv <- rbind(cv,
            pred18[, .(rmse = sqrt(mean((abundance - pred)^2)),
                       prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
                       year = 'yr2018'),
                   by = abundance > 0]
)


fake <- croaker[year == 2019]
fake <- rbind(fake, fake[1])
fake[nrow(fake), season := 'FALL']


pred19 <- predict(cv_5y$m18, newdata = fake,
                  type = 'link',
                  exclude = c('s(year)', 's(X,Y)', 'seasonSPRING'))

pred19 <- setDT(
  cbind(
    croaker[year == 2019],
    pred = exp(pred19[-length(pred19)])
  )
)

cv <- rbind(cv,
            pred19[, .(rmse = sqrt(mean((abundance - pred)^2)),
                       prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
                       year = 'yr2019'),
                   by = abundance > 0]
)



## cv using 2014 and predicting forward
pred15 <- predict(abun$nb, newdata = croaker[year == 2015],
                  type = 'link',
                  exclude = c('s(year)', 's(X,Y)'))

pred15 <- setDT(
  cbind(
    croaker[year == 2015],
    pred = exp(pred15)
  )
)


cv <- pred15[, .(rmse = sqrt(mean((abundance - pred)^2)),
                 prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
                 year = 'yr2015'),
             by = abundance > 0]




pred16 <- predict(abun$nb, newdata = croaker[year == 2016],
                  type = 'link',
                  exclude = c('s(year)', 's(X,Y)'))

pred16 <- setDT(
  cbind(
    croaker[year == 2016],
    pred = exp(pred16)
  )
)


cv <- rbind(cv,
            pred16[, .(rmse = sqrt(mean((abundance - pred)^2)),
                 prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
                 year = 'yr2016'),
             by = abundance > 0]
)


pred17 <- predict(abun$nb, newdata = croaker[year == 2017],
                  type = 'link',
                  exclude = c('s(year)', 's(X,Y)'))

pred17 <- setDT(
  cbind(
    croaker[year == 2017],
    pred = exp(pred17)
  )
)


cv <- rbind(cv,
            pred17[, .(rmse = sqrt(mean((abundance - pred)^2)),
                 prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
                 year = 'yr2017'),
             by = abundance > 0])



pred18<- predict(abun$nb, newdata = croaker[year == 2018],
                  type = 'link',
                  exclude = c('s(year)', 's(X,Y)'))

pred18 <- setDT(
  cbind(
    croaker[year == 2018],
    pred = exp(pred18)
  )
)


cv <- rbind(cv,
            pred18[, .(rmse = sqrt(mean((abundance - pred)^2)),
                 prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
                 year = 'yr2018'),
             by = abundance > 0])



fake <- croaker[year == 2019]
fake <- rbind(fake, fake[1])
fake[nrow(fake), season := 'FALL']


pred19 <- predict(abun$nb, newdata = fake,
                  type = 'link',
                  exclude = c('s(year)', 's(X,Y)', 'seasonSPRING'))

pred19 <- setDT(
  cbind(
    croaker[year == 2019],
    pred = exp(pred19[-length(pred19)])
  )
)

cv <- rbind(cv,
            pred19[, .(rmse = sqrt(mean((abundance - pred)^2)),
                       prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
                       year = 'yr2019'),
                   by = abundance > 0]
)






mods <- readRDS('data derived/model output/bam_abun_bottemp_xy_yrran_yrly.RDS')

mod_2014 <- mods$m14
mod_2019 <- mods$m19

### Yearly forecast RMSE
forecast_rmse <- function(yr){
  df <- copy(subs)[year == yr]
  # need to have a factor level that was in the model fitting
  df[, yr_fac := '1967']


  # If the subset has only 1 of "FALL" or "SPRING", predict() won't work.
  if(length(unique(df$season)) == 1){
    df <- rbind(df, df[1])
    df[nrow(df), season :=
         c('SPRING', 'FALL')[!c('SPRING', 'FALL') %in% unique(df$season)]]
    df[, dummy := T]
  }

  preds <- predict(mod_2014,
                   newdata = df, type = 'response',
                   exclude = 's(yr_fac)')

  df <- setDT(cbind(
    df,
    pred = preds
  ))

  if('dummy' %in% names(df)){
    df[-nrow(df)]
  }

  df[, .(rmse = sqrt(mean((abundance - pred)^2)),
         prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
         year = paste0('yr', yr - 2014)),
     by = abundance > 0]

}

rbind(
  forecast_rmse(2015),
  forecast_rmse(2016),
  forecast_rmse(2017),
  forecast_rmse(2018))

,
  forecast_rmse(2019),
)




#P/A-----
subs[, present := fifelse(abundance > 0, T, F)]


mod_bin <- bam(present ~ season + grpsed +
                     s(bottemp, bs = 'cr') +
                     s(yr_fac, bs = 're') +
                     s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
                   data = subs[between(year, 1967, 2014)], family = binomial(),
                   control =  gam.control(scalePenalty = FALSE),
                   discrete = T, nthreads = 11)

