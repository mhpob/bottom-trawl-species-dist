library(data.table);library(mgcv)
source("miller et al. supp/mgcv_spde_smooth.R") # note that this loads INLA



# load strata shp
strata_info <- fread('data/svdbs_supporttables/svdbs_svmstrata.csv')

library(sf); library(dplyr)



# subset strata to midatl/sne
strata <- read_sf('data/mapping/strata/strata.shp')
setDT(strata)

strata <- strata[strata_info, on = c('STRATA' = 'stratum')]
strata <- strata[grepl('MAB|SNE', stratum_name)]
strata[, geometry := st_transform(geometry, crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')]

strata <- strata[, .(geometry = st_union(geometry))]

strata[, geometry := st_simplify(geometry, preserveTopology = F, dTolerance = 5000)]



boundary <- inla.sp2segment(as_Spatial(strata$geometry))
# make mesh
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
subs <- st_as_sf(subs, coords = c('lon', 'lat'), crs = 4326)
subs <- setDT(subs)[, geometry := st_transform(geometry, crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')]
subs <- cbind(subs, st_coordinates(subs$geometry))

# select trawls within the boundary

mod <- bam(abundance ~ s(bottemp) + s(X, Y, bs = "spde", k = mesh$n, xt = list(mesh = mesh)),
           data = subs, family = poisson(),
           control =  gam.control(scalePenalty = FALSE),
           discrete = T)
