library(data.table);library(mgcv)
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
subs <- survdat[comname == 'atlantic croaker']
subs <- unique(subs, by = c('cruise6', 'station', 'stratum', 'tow'))
subs <- subs[complete.cases(subs[, .(season, bottemp, year, lat, lon)])]
subs[, yr_fac := as.factor(year)]
subs <- st_as_sf(subs, coords = c('lon', 'lat'), crs = 4326)
subs <- setDT(subs)[, geometry := st_transform(geometry, crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')]
subs <- cbind(subs, st_coordinates(subs$geometry))


# select trawls within the boundary
subs <- subs[as.vector(st_intersects(geometry, strata$geometry, sparse = F))]


# Biomass
# mod_tw <- bam(biomass ~ season + s(bottemp) +
#              s(yr_fac, bs = 're') +
#              s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
#            data = subs[year < 2014], family = tw(),
#            control =  gam.control(scalePenalty = FALSE),
#            discrete = T, nthreads = 11)
#
# mod_nb <- bam(biomass ~ season + s(bottemp) +
#                 s(yr_fac, bs = 're') +
#                 s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
#               data = subs[year < 2014], family = nb(),
#               control =  gam.control(scalePenalty = FALSE),
#               discrete = T, nthreads = 11)

# Abundance
mod_tw_abun_no921 <- bam(abundance ~ season + s(bottemp) +
                 s(yr_fac, bs = 're') +
                 s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
               data = subs[year <= 2014][-921], family = tw(),
               control =  gam.control(scalePenalty = FALSE),
               discrete = T, nthreads = 11)

mod_poi_abun <- bam(abundance ~ season + s(bottemp) +
                     s(yr_fac, bs = 're') +
                     s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
                   data = subs[year <= 2014], family = poisson(),
                   control =  gam.control(scalePenalty = FALSE),
                   discrete = T, nthreads = 11)

mod_nb_abun <- bam(abundance ~ season + s(bottemp) +
                      s(yr_fac, bs = 're') +
                      s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
                    data = subs[year <= 2014], family = nb(),
                    control =  gam.control(scalePenalty = FALSE),
                    discrete = T, nthreads = 11)

mod_zip_abun <- bam(abundance ~ season + s(bottemp) +
                     s(yr_fac, bs = 're') +
                     s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
                   data = subs[year <= 2014], family = ziP(),
                   control =  gam.control(scalePenalty = FALSE),
                   discrete = T, nthreads = 11)
#
# abun <- list(nb = mod_nb_abun,
#              poi = mod_poi_abun,
#              tw = mod_tw_abun,
#              zip = mod_zip_abun)
# #
# saveRDS(abun, 'data derived/model output/bam_abun_bottemp_xy_yrran.RDS')


abun <- readRDS('data derived/model output/bam_abun_bottemp_xy_yrran.RDS')
lapply(abun, AIC)

k <- predict(abun$nb, newdata = abun$nb$model[1:15,], type = 'response', se = T,
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
  geom_sf(data = coast) +
  geom_point(data = spat_dat,
             aes(x = X, y = Y, color = value)) +
  coord_sf(xlim = c(-2e5, 2.1e5), ylim = c(3.6e6, 4.35e6)) +
  scale_color_viridis_c() +
  labs(x = NULL, y = NULL, color = NULL,
       title = 'Atlantic croaker modeled abundance',
       subtitle = 'Partial effect of X-Y space') +
  facet_wrap(~variable) +
  theme_minimal() +
  theme(legend.position = c(0.9, 0.4))



## Cross validation ----
mod <- mod_nb_abun

subs <- subs[complete.cases(subs[, .(abundance, season, bottemp, X, Y)])]

run1 <- subs[year == 2015
             ]

# need to have a factor level that was in the model fitting
run1[, yr_fac := '1967']

# If the subset has only 1 of "FALL" or "SPRING", predict() won't work.
if(length(unique(run1$season)) == 1){
  run1 <- rbind(run1, run1[1])
  run1[nrow(run1), season :=
         c('SPRING', 'FALL')[!c('SPRING', 'FALL') %in% unique(run1$season)]]
  run1[, dummy := T]
}


preds <- predict(mod, newdata = run1, type = 'link',
                 exclude = 's(yr_fac)')

run1 <- setDT(cbind(
  run1,
  pred = preds
))

if('dummy' %in% names(run1)){
  run1[-nrow(run1)]
}

run1[, .(rmse = sqrt(mean((abundance - pred)^2)),
         prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
         fold = 'fold1'),
     by = abundance > 0]



mod_nb_cv <- bam(abundance ~ season + s(bottemp) +
                   s(yr_fac, bs = 're') +
                   s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
                 data = subs[year <= 2015], family = nb(),
                 control =  gam.control(scalePenalty = FALSE),
                 discrete = T, nthreads = 11)



