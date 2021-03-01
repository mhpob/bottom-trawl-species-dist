library(sf); library(data.table); library(mgcv)

all_data <- fread("data derived/survdat_mabsne_only.csv")

station_key <- unique(all_data,
                      by = c('cruise6', 'station', 'stratum', 'tow',
                             'year', 'season'))
station_key[, 18:23 := NULL]

croaker <- all_data[comname == 'atlantic croaker']

croaker <- croaker[station_key, on = names(station_key)]

croaker[, abundance := fifelse(is.na(abundance), 0, abundance)]

croaker <- st_as_sf(croaker,
                    coords = c('lon', 'lat'),
                    remove = F,
                    crs = 4326) %>%
  st_transform(crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')

croaker <- setDT(croaker)[complete.cases(croaker[, .(season, bottemp, botsalin,
                                                     depth, grpsed, year,
                                                     lat, lon)])]

croaker[, ':='(X = st_coordinates(geometry)[, 1],
               Y = st_coordinates(geometry)[, 2])]


library(mgcv)
source("miller et al. supp/mgcv_spde_smooth.R") #note: loads INLA

mesh <- readRDS('data derived/mab_sne_mesh.RDS')


forms <- list(cm_f = "abundance ~ season + grpsed +
              s(depth, by = season, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_0 = "abundance ~ s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_1 = "abundance ~ s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_2 = "abundance ~
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_3 = "abundance ~ season +
              s(depth, by = season, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_4 = "abundance ~ grpsed +
              s(botsalin, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_5 = "abundance ~ season +
              s(depth, by = season, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_6 = "abundance ~ season + grpsed +
              s(depth, by = season, bs = 'cs', k = 10) +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_7 = "abundance ~ grpsed +
              s(bottemp, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_8 = "abundance ~
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_9 = "abundance ~ season +
              s(depth, by = season, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_10 = "abundance ~ season + grpsed +
              s(depth, by = season, bs = 'cs', k = 10) +
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_11 = "abundance ~ grpsed +
              s(botsalin, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_12 = "abundance ~ season +
              s(depth, by = season, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_13 = "abundance ~ season + grpsed +
              s(depth, by = season, bs = 'cs', k = 10) +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))",
           cm_14 = "abundance ~ grpsed +
              s(year, bs = 'cs', k = 10) +
              s(X, Y, bs = 'spde', k = mesh$n, xt = list(mesh = mesh))")




forms <- gsub('\n', ' ', forms)

model_dat <- croaker[year <= 2014]
model_dat[, season := as.factor(season)]



cl <- makeCluster(11)
clusterEvalQ(cl, library(mgcv))
clusterEvalQ(cl, library(INLA))
clusterExport(cl, 'model_dat')
clusterExport(cl, 'mesh')
clusterExport(cl, 'Predict.matrix.spde.smooth')
clusterExport(cl, 'smooth.construct.spde.smooth.spec')

st <- Sys.time()

models <- parLapply(cl = cl,
                    forms,
                    function(.){
                      bam(formula = as.formula(.),
                          data = model_dat,
                          family = nb(),
                          control = gam.control(scalePenalty = F),
                          discrete = T)
                    })

end <- Sys.time()

stopCluster(cl)

saveRDS(models, 'data derived/model output/croaker_nb_suite_shrinkage.rds')

#yr as cr = 38.28167 min