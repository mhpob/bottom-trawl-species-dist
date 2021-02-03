# Option 1: Download fully-aggregated file ----
#    Source is NOAA Ecosystem Dynamics and Assessment Branch. Data is nested
#     within a non-exported section of their {ecodata} package
#     https://github.com/NOAA-EDAB/ecodata

download.file('https://github.com/NOAA-EDAB/ecodata/raw/master/data-raw/Survdat.RData',
              destfile = 'data/noaa_edab_survdat.rdata')









# Option 2: Do it yourself ----
## Support tables ----
#   These will be used the get consistent species names, etc., and can be used for
#   both the spring and fall bottom trawl surveys.

### Download ----
download.file(paste0('ftp://ftp.nefsc.noaa.gov/pub/dropoff/',
                     'PARR/PEMAD/ESB/SVDBS/SVDBS_SupportTables.zip'),
              destfile = 'data/SVDBS_SupportTables.zip')

### Unzip ----
unzip('data/SVDBS_SupportTables.zip',
      exdir = 'data/SVDBS_SupportTables')
file.remove('data/SVDBS_SupportTables.zip')



## Spring bottom trawl data ----

### Download ----
download.file(paste0('ftp://ftp.nefsc.noaa.gov/pub/dropoff/',
                     'PARR/PEMAD/ESB/22561/22561_FSCSTables.zip'),
              destfile = 'data/22561_FSCSTables.zip')

### Unzip ----
unzip('data/22561_FSCSTables.zip',
      exdir = 'data/22561_FSCSTables')
file.remove('data/22561_FSCSTables.zip')





## Fall bottom trawl data ----

### Download ----
download.file(paste0('ftp://ftp.nefsc.noaa.gov/pub/dropoff/',
                     'PARR/PEMAD/ESB/22560/22560_FSCSTables.zip'),
              destfile = 'data/22560_FSCSTables.zip')


### Unzip ----
unzip('data/22560_FSCSTables.zip',
      exdir = 'data/22560_FSCSTables')
file.remove('data/22561_FSCSTables.zip')