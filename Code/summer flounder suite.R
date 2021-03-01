library(ggplot2); library(sf); library(data.table)

all_data <- fread("data derived/survdat_mabsne_only.csv")

station_key <- unique(all_data,
                      by = c('cruise6', 'station', 'stratum', 'tow',
                             'year', 'season'))
station_key[, 18:23 := NULL]

sumfl <- all_data[comname == 'summer flounder']

sumfl <- sumfl[station_key, on = names(station_key)]

sumfl[, abundance := fifelse(is.na(abundance), 0, abundance)]

sumfl <- st_as_sf(sumfl,
                  coords = c('lon', 'lat'),
                  remove = F,
                  crs = 4326) %>%
  st_transform(crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')

sumfl <- setDT(sumfl)[complete.cases(sumfl[, .(season, bottemp, botsalin,
                                               grpsed, year, depth,
                                               lat, lon)])]

sumfl[, ':='(X = st_coordinates(geometry)[, 1],
               Y = st_coordinates(geometry)[, 2])]


library(mgcv)
source("miller et al. supp/mgcv_spde_smooth.R") #note: loads INLA


mesh <- readRDS('data derived/mab_sne_mesh.RDS')



# Model error structure selection ----
library(parallel)
library(pbapply)

cl <- makeCluster(detectCores(logical = F) - 2)
clusterEvalQ(cl = cl, c(library(mgcv), library(INLA), library(data.table)))
clusterExport(cl, c('sumfl', 'mesh', 'Predict.matrix.spde.smooth',
                    'smooth.construct.spde.smooth.spec'))

models <- pbsapply(cl = cl,
                   X = list(tw = 'tw', pois = 'poisson', negbin = 'nb', zipois = 'ziP'),
                   FUN = function(.){
                     bam(abundance ~ season + grpsed +
                           s(depth, bs = 'cr', k = 10) + s(botsalin, bs = 'cr', k = 10) +
                           s(bottemp, bs = 'cr', k = 10) +
                           s(year, bs = 'cr', k = 10) +
                           s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
                         data = sumfl[year <= 2014], family = .,
                         control =  gam.control(scalePenalty = FALSE),
                         discrete = T)
                   },
                   simplify = F,
                   USE.NAMES = T)

stopCluster(cl)


saveRDS(models, 'data derived/model output/sumfl_familyselection.RDS')


## Model term selection ----
forms <- list(sfm_f = "abundance ~ season + grpsed +
              s(depth, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_0 = "abundance ~ s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_1 = "abundance ~ s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_2 = "abundance ~
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_3 = "abundance ~ season +
              s(depth, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_4 = "abundance ~ grpsed +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_5 = "abundance ~ season +
              s(depth, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_6 = "abundance ~ season + grpsed +
              s(depth,  bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_7 = "abundance ~ grpsed +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_8 = "abundance ~
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_9 = "abundance ~ season +
              s(depth,  bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_10 = "abundance ~ season + grpsed +
              s(depth, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_11 = "abundance ~ grpsed +
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_12 = "abundance ~ season +
              s(depth, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_13 = "abundance ~ season + grpsed +
              s(depth, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
              sfm_14 = "abundance ~ grpsed +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))")




forms <- sapply(forms, function(.) gsub('\n|\\s+', ' ', .),
                simplify = F, USE.NAMES = T)

model_dat <- sumfl[year <= 2014]
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

saveRDS(models, 'data derived/model output/summfl_modelselection.rds')


