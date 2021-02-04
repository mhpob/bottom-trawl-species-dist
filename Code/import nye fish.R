library(data.table)

load('data/noaa_edab_survdat.rdata')

setnames(survdat, tolower)

# Drop length
survdat[, c('length', 'numlen') := NULL]

# Unique info by trawl
survdat <- unique(survdat, by = c('cruise6', 'station', 'stratum', 'tow', 'svspp'))


# Read in species used by Nye et al. 2009
nye_species <- fread('data derived/nye_species.csv',
                     col.names = function(.) tolower(gsub(' ', '_', .)))
nye_species[, 1:6 := lapply(.SD, tolower), .SDcol = 1:6]


# Make sure each species is listed in each trawl, even if catch was 0
nye_species <- survdat[, .SD[nye_species, on = 'svspp'],
                       by = c('cruise6', 'station', 'stratum', 'tow')]


# Fill in necessary NA values with tow information
nye_species[, 7:17 := lapply(.SD, function(.) fifelse(is.na(.), unique(.)[1], .)),
            by = c('cruise6', 'station', 'stratum', 'tow'),
            .SDcols = 7:17]

# Set NA abundance and biomass to 0
nye_species[, 18:19 := lapply(.SD, nafill, type = 'const', fill = 0),
            .SDcols = 18:19]


# Export
fwrite(nye_species, 'data derived/nye_svdata.gz')