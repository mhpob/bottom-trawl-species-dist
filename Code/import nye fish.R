library(data.table)

survdat <- fread('data derived/survdat_w_sediment.gz')

setnames(survdat, tolower)

# Drop length and sex
survdat[, c('length', 'numlen', 'catchsex') := NULL]

# Unique info by trawl
survdat <- unique(survdat, by = c('cruise6', 'station', 'stratum', 'tow', 'svspp'))


# Read in species used by Nye et al. 2009
nye_species <- fread('data derived/nye_species.csv',
                     col.names = function(.) tolower(gsub(' ', '_', .)))
nye_species[, 1:6 := lapply(.SD, tolower), .SDcol = 1:6]

# Names by which to join
nn <- names(survdat)[c(1:4, 6:16, 19:20)]

# Make sure each species is listed in each trawl, even if catch was 0
nye_species <- survdat[, .SD[nye_species, on = 'svspp'],
                       by = nn]

# Set NA abundance and biomass to 0
nye_species[, c('abundance', 'biomass') := lapply(.SD, nafill, type = 'const', fill = 0),
            .SDcols = c('abundance', 'biomass')]


# Export
fwrite(nye_species, 'data derived/nye_svdata.gz')
