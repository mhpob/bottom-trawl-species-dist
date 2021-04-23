
library(data.table)
all_data <- fread('data derived/survdat_names_sed.csv')

capfun <- Vectorize(
  function(x) paste(toupper(substring(x, 1, 1)), substring(x, 2),
                    sep = "", collapse = " "),
  USE.NAMES = FALSE
)

spp <- capfun(unique(all_data$sciname))
spp_com <- capfun(unique(all_data$comname))

library(rfishbase)

full_table <- species(fields = c('Species', 'FBname', 'DemersPelag'))
setDT(full_table)

spp <- full_table[Species %in% spp]
spp_com <- full_table[FBname %in% spp_com]

spp <- rbind(spp, spp_com)

spp <- unique(spp, by = 'Species')
spp <- unique(spp, by = 'FBname')

# unique(spp$DemersPelag)
# [1] "reef-associated" "bathypelagic"    "pelagic-oceanic" "pelagic-neritic" "demersal"        "benthopelagic"
# [7] "bathydemersal"   "pelagic"

# Not selecting "bathy-" (likely not shelf), "pelagic-" (want bottom-oriented spp).


p <- spp[DemersPelag %in% c('reef-associated', 'demersal')]


demersal <- all_data[sciname %in% tolower(p$Species) |
                       comname %in% tolower(p$FBname)]

fwrite(demersal, 'data derived/demersal.csv')
