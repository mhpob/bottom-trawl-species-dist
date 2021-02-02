# Refer to "code/download trawl data" for data access.

# Packages ----
library(data.table)





# Species table ----
species_list <- fread(input = 'data/svdbs_supporttables/svdbs_svspecies_list.csv',
                      col.names = tolower,
                      fill = T
)

# A few species names that should be in "logged_species_name" span multiple columns,
#   forcing the columns to bump further to the right. The following lines fix that.
species_list[!grepl('^\\d', svspp), ':='(
  comname = paste(comname, svspp, v1),
  svspp = v2)]
species_list[, c('v1', 'v2') := NULL]
species_list[, svspp := as.numeric(svspp)]
species_list[, c('sciname', 'comname') := lapply(.SD, tolower), .SDcols = c('sciname', 'comname')]





# Spring botom trawl survey ----
## Data import ----
spring_catch <- fread(
  grep('SVCAT', list.files('data/22561_FSCSTables', full.names = T), value = T),
  col.names = tolower, fill = T
)

## Data cleaning ----
##  A few species names that should be in "logged_species_name" span multiple columns,
##    forcing the columns to bump further to the right
spring_catch[!grepl('^\\d', svspp), ':='(
  logged_species_name = paste(logged_species_name, svspp, catchsex),
  svspp = expcatchnum,
  catchsex = expcatchwt,
  expcatchnum = v1,
  expcatchwt = v2)]
spring_catch[, c('v1', 'v2') := NULL]
spring_catch[, ':='(logged_species_name = tolower(logged_species_name),
                    svspp = as.numeric(svspp),
                    # ID can wind up being a 17-digit integer, which can cause issues
                    #   on 32bit systems convert to character to avoid this
                    id = as.character(id),
                    season = 'spring')]


## Join species names ----
##  Join catch data and species list to get consistent species naming
spring_catch <- species_list[spring_catch, on = 'svspp']


## Join station data ----
spring_station <- fread(
  grep('SVSTA', list.files('data/22561_FSCSTables', full.names = T), value = T),
  col.names = tolower
)

# IDs are saved as long numeric code -- causes issues on 32bit systems. Convert to
#   character since keeping it as a number isn't necessary.
spring_station[, id := as.character(id)]


spring <- spring_station[spring_catch, on = c('cruise6', 'stratum', 'tow', 'station', 'id')]





# Fall bottom trawl survey ----
## Data import ----
fall_catch <- fread(
  grep('SVCAT', list.files('data/22560_FSCSTables', full.names = T), value = T),
  col.names = tolower, fill = T
)


## Data cleaning ----
##  A few species names that should be in "logged_species_name" span multiple columns,
##    forcing the columns to bump further to the right
fall_catch[!grepl('^\\d', svspp), ':='(
  logged_species_name = paste(logged_species_name, svspp, catchsex),
  svspp = expcatchnum,
  catchsex = expcatchwt,
  expcatchnum = v11,
  expcatchwt = v12)]
fall_catch[, c('v11', 'v12') := NULL]
fall_catch[, ':='(logged_species_name = tolower(logged_species_name),
                  svspp = as.numeric(svspp),
                  # ID can wind up being a 17-digit integer, which can cause issues
                  #   on 32bit systems convert to character to avoid this
                  id = as.character(id),
                  season = 'fall')]


## Join species names ----
##  Join catch data and species list to get consistent species naming
fall_catch <- species_list[fall_catch, on = 'svspp']


## Join station data ----
fall_station <- fread(
  grep('SVSTA', list.files('data/22560_FSCSTables', full.names = T), value = T),
  col.names = tolower
)

# IDs are saved as long numeric code -- causes issues on 32bit systems. Convert to
#   character since keeping it as a number isn't necessary.
fall_station[, id := as.character(id)]

fall <- fall_station[fall_catch, on = c('cruise6', 'stratum', 'tow', 'station', 'id')]





# Export aggregated data ----
all <- rbind(
  spring,
  fall
)

fwrite(all, 'data derived/trawl_data.gz')