library(parallel); library(pbapply); library(sf)
library(data.table); library(mgcv)
source("miller et al. supp/mgcv_spde_smooth.R") #note: loads INLA

load_species_data <- function(species){

  all_data <- fread("data derived/survdat_mabsne_only.csv")

  station_key <- unique(all_data,
                        by = c('cruise6', 'station', 'stratum', 'tow',
                               'year', 'season'))
  station_key[, 18:23 := NULL]

  species_data <- all_data[comname == species]

  species_data <- species_data[station_key, on = names(station_key)]

  species_data[, abundance := fifelse(is.na(abundance), 0, abundance)]

  species_data <- st_as_sf(species_data,
                           coords = c('lon', 'lat'),
                           remove = F,
                           crs = 4326) %>%
    st_transform(crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')

  species_data <- setDT(species_data)[
    complete.cases(species_data[, .(season, bottemp, botsalin,
                                    grpsed, year, depth,
                                    lat, lon)])]

  species_data[, ':='(X = st_coordinates(geometry)[, 1],
                      Y = st_coordinates(geometry)[, 2])]
}




# Load mesh ----
mesh <- readRDS('data derived/mab_sne_mesh.RDS')


# Load species ----
spec_abun <- load_species_data('atlantic croaker')
spec_abun <- spec_abun[abundance > 0]


# Model error structure selection ----
cl <- makeCluster(detectCores(logical = F) - 2)
clusterEvalQ(cl = cl, c(library(mgcv), library(INLA), library(data.table)))
clusterExport(cl, c('spec_abun', 'mesh', 'Predict.matrix.spde.smooth',
                    'smooth.construct.spde.smooth.spec'))

models <- pbsapply(cl = cl,
                   X = list(tw = 'tw', pois = 'poisson',
                            negbin = 'nb'),
                   FUN = function(.){
                     bam(abundance ~ season + grpsed +
                           s(depth, bs = 'cs', k = 10) +
                           s(botsalin, bs = 'cs', k = 10) +
                           s(bottemp, bs = 'cs', k = 10) +
                           s(year, bs = 'cs', k = 10) +
                           s(X, Y, bs = 'ds', m = c(1, 0.5)),
                         data = spec_abun[year <= 2014], family = .,
                         control =  gam.control(scalePenalty = FALSE),
                         discrete = T)
                   },
                   simplify = F,
                   USE.NAMES = T)

stopCluster(cl)


saveRDS(models, 'data derived/model output/ac_PO_familyselection.RDS')

sapply(models, AIC)[order(sapply(models, AIC))]


## Model term selection ----
forms <- list(ac_f = "abundance ~ season + grpsed +
              s(depth, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_0 = "abundance ~ s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_1 = "abundance ~ s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_2 = "abundance ~
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_3 = "abundance ~ season +
              s(depth, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_4 = "abundance ~ grpsed +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_5 = "abundance ~ season +
              s(depth, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_6 = "abundance ~ season + grpsed +
              s(depth,  bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_7 = "abundance ~ grpsed +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_8 = "abundance ~
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_9 = "abundance ~ season +
              s(depth,  bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_10 = "abundance ~ season + grpsed +
              s(depth, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_11 = "abundance ~ grpsed +
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_12 = "abundance ~ season +
              s(depth, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_13 = "abundance ~ season + grpsed +
              s(depth, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_14 = "abundance ~ grpsed +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_15 = "abundance ~
              s(depth, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_16 = "abundance ~ grpsed +
              s(depth, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_17 = "abundance ~ season +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              ac_18 = "abundance ~ season + grpsed +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))")




forms <- sapply(forms, function(.) gsub('\n|\\s+', ' ', .),
                simplify = F, USE.NAMES = T)

model_dat <- spec_abun[year <= 2014]
model_dat[, season := as.factor(season)]



cl <- makeCluster(11)
clusterEvalQ(cl, c(library(mgcv), library(INLA)))
clusterExport(cl, c('model_dat', 'mesh', 'Predict.matrix.spde.smooth',
                    'smooth.construct.spde.smooth.spec'))

st <- Sys.time()

models <- pbsapply(cl = cl,
                   forms,
                   function(.){
                     bam(formula = as.formula(.),
                         data = model_dat,
                         family = nb(),
                         control = gam.control(scalePenalty = F),
                         discrete = T)
                   },
                   simplify = F,
                   USE.NAMES = T)

end <- Sys.time()

stopCluster(cl)

end - st #22.25 min

saveRDS(models, 'data derived/model output/ac_PO_modelselection.rds')

sapply(models, AIC)[order(sapply(models, AIC))]

