library(INLA); library(data.table); library(sf)

strata_info <- fread('data/svdbs_supporttables/svdbs_svmstrata.csv')

strata <- read_sf('data/mapping/strata/strata.shp')
setDT(strata)

strata <- strata[strata_info, on = c('STRATA' = 'stratum')]
strata <- strata[grepl('MAB|SNE', stratum_name)]
strata[, geometry :=
         st_transform(geometry, crs = '+proj=aea +lat_1=33 +lat_2=45 +lon_0=-74')]

strata <- strata[, .(geometry = st_union(geometry))]

strata[, geometry := st_simplify(geometry,
                                 preserveTopology = F,
                                 dTolerance = 5000)]

boundary <- inla.sp2segment(as_Spatial(strata$geometry))
# make mesh using meshbuilder 'cause it's finnicky
# meshbuilder()

mesh <- inla.mesh.2d(boundary=boundary,
                     max.edge=c(22000, 108000),
                     min.angle=c(30, 21),

                     ## Safeguard against large meshes.
                     max.n=c(48000, 16000),

                     ## Don't build a huge mesh!
                     max.n.strict=c(128000, 128000),

                     ## Filter away adjacent points.
                     cutoff=8000,
                     offset=c(1000, 55000))


saveRDS(mesh, 'data derived/mab_sne_mesh.RDS')
