library(data.table)
nye_fish <- fread('data derived/nye_svdata.gz')




spp15 <- nye_species[svspp == 15]






library(mgcv)

m1 <- bam(abundance ~ te(lat,lon,year),
          data = nye_species,
          subset = svspp == 15,
          family = 'tw',
          discrete = T)

m2 <- bam(abundance ~ te(lat, lon, year) + s(bottemp),
          data = nye_species,
          subset = svspp == 15,
          family = 'tw',
          discrete = T)


m3 <- bam(abundance ~ s(lon, lat) + s(year) + ti(lon, lat, year, d = c(2,1)),
          data = spp15,
          family = 'tw',
          discrete = T)

# https://github.com/eric-pedersen/mgcv-esa-workshop/blob/master/example-spatio-temporal%20data.Rmd


m4 <- bam(abundance ~ s(lon, lat) + s(year) + ti(lon, lat, year, d = c(2,1)) +
            s(bottemp),
          data = nye_species,
          subset = svspp == 15,
          family = 'tw',
          discrete = T)


m5 <- bam(abundance ~ s(lon, lat) + s(year) + ti(lon, lat, year, d = c(2,1)) +
            s(bottemp) + s(depth),
          data = spp15,
          family = 'tw',
          discrete = T)


m5_pois <- bam(abundance ~ s(lon, lat) + s(year) + ti(lon, lat, year, d = c(2,1)) +
            s(bottemp) + s(depth),
          data = spp15,
          family = 'ziP',
          discrete = T)


newdat <-





