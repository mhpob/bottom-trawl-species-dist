library(data.table)


all_data <- fread('data derived/survdat_names_sed.csv')

station_key <- unique(all_data,
                      by = c('cruise6', 'station', 'stratum', 'tow',
                             'year', 'season'))
station_key[, 18:25 := NULL]

croaker <- all_data[comname == 'atlantic croaker']
croaker <- croaker[station_key, on = names(station_key)]
croaker <- unique(croaker, by = c('cruise6', 'station', 'stratum', 'tow',
                                  'year', 'season'))


croaker[, abundance := fifelse(is.na(abundance), 0, abundance)]
croaker[, period := as.factor(fifelse(year <= 1998, 'start', 'end'))]




wq <- croaker[period == 'start']$bottemp
det <- croaker[period == 'start']$abundance

# Create breaks.
minval <- min(wq, na.rm = T)
maxval <- max(wq, na.rm = T)

brks <- seq(minval, maxval, bin_width)
brks <- if(maxval > max(brks)) c(brks, max(brks) + bin_width) else brks

# Create grouping bins.
bins <- cut(wq, brks,
            include.lowest = T,
            ordered_result = T)

# Aggregate by environmental bins.
agg.func <- function(x){
  if(pres_abs == T){
    # Only count number over 0 if presence/absence
    x <- x > 0
  }
  hold <- as.data.frame(xtabs(x ~ bins), responseName = 'det')
  hold$bins <- as.ordered(hold$bins)
  hold
}

boot_func <- function(x, index){
  x <- x[index]
  boot_det <- agg.func(x)
  # Bootstrapped pMe
  as.vector(boot_det$det)/sum(boot_det$det)
}

strap <- boot::boot(det, boot_func, R)

#Confidence Interval
ci <- matrix(nrow = length(strap$t0), ncol = 2)
for(i in 1:length(strap$t0)){
  ci[i,] <- boot::boot.ci(strap, type = 'perc', index = i)$percent[4:5]
}

fish <- agg.func(det)

station <- as.data.frame(table(bins), responseName = 'wq')
station$bins <- as.ordered(station$bins)


# Merge data and correctly order bins.
q_an <- merge(fish, station, by = 'bins', sort = F)

# Quotient analysis
q_an$pme <- q_an$det / sum(q_an$det)
q_an$pse <- q_an$wq / sum(q_an$wq)

q_an$qe <- q_an$pme / q_an$pse

q_an$ci.025 <- ci[, 1] / q_an$pse
q_an$ci.975 <- ci[, 2] / q_an$pse

q_an$result.95 <- ifelse(q_an$qe < q_an$ci.025, 'avoidace',
                         ifelse(q_an$qe > q_an$ci.975, 'attraction',
                                'ns'))

names(q_an) <- c('bin', 'detections', 'wq.var', 'pMe', 'pSe',
                 'Qe', 'CI_0.025', 'CI_0.975', 'sig_0.95')
q_an


library(mgcv)
croaker[, yr_fac := as.factor(year)]
croaker[, season := as.factor(season)]
jj <- bam(abundance ~ te(year, bottemp) + te(lon, lat, by =season, bs = 'cr'),
          discrete = T,
          data = croaker[abundance >0],
          family = nb())
draw(jj)

k <- difference_smooths(jj, 's(bottemp)')
draw(k)



jk <- bam(abundance ~ s(bottemp) + s(bottemp, by = period) + te(lon, lat, bs = 'cr'),
          discrete = T,
          data = croaker[abundance >0 & season == 'FALL'],
          family = nb())
