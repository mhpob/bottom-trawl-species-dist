
library(data.table)
all_data <- fread('data derived/survdat_names_sed.csv')

all_data <- all_data[comname != '' | sciname != '']

capfun <- Vectorize(
  function(x) paste(toupper(substring(x, 1, 1)), substring(x, 2),
                    sep = "", collapse = " "),
  USE.NAMES = FALSE
)

spp_sci <- capfun(unique(all_data$sciname))
spp_com <- capfun(unique(all_data$comname))

library(rfishbase)

fulltable_fb <- species(fields = c('Species', 'FBname', 'DemersPelag'))
fulltable_sl <- species(fields = c('Species', 'FBname', 'DemersPelag'),
                        server = 'sealifebase')
setDT(fulltable_fb)
setDT(fulltable_sl)

spp_sci_fb <- fulltable_fb[Species %in% spp_sci]
spp_com_fb <- fulltable_fb[FBname %in% spp_com]

spp_sci_sl <- fulltable_sl[Species %in% spp_sci]
spp_com_sl <- fulltable_sl[FBname %in% spp_com]

spp <- rbind(spp_sci_fb, spp_com_fb, spp_sci_sl, spp_com_sl)

unique(spp$DemersPelag)
# [1] "reef-associated" "bathypelagic"    "pelagic-oceanic" "pelagic-neritic" "demersal"        "benthopelagic"
# [7] "bathydemersal"   "pelagic"         "benthic"

# Not selecting "bathy-" (likely not shelf), "pelagic-" (want bottom-oriented spp).

spp <- unique(spp)

spp <- spp[DemersPelag %in% c('reef-associated', 'demersal', 'benthic')]


demersal <- all_data[sciname %in% tolower(spp$Species) |
                       comname %in% tolower(spp$FBname)]

fwrite(demersal, 'data derived/demersal.csv')
