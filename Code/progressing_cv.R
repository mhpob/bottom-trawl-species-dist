progressing_cv <- function(model, full_data, yr_seq){

  library(data.table); library(mgcv)

  data_prep <- function(subset_data){

    # Need to have a factor level for year that was in the original model fitting
    # Since year was a random effect, this will be removed by "exclude" in predict()
    subset_data[, yr_fac := model$model$yr_fac[1]]


    # If the subset has only 1 of "FALL" or "SPRING", predict() won't work.
    # Check to see if all seasons are accounted for. If not, then tack a dummy
    #   observation onto the end
    unique_seasons <- unique(subset_data$season)

    if(length(unique_seasons) == 1){
      subset_data <- rbind(subset_data, subset_data[1])

      subset_data[nrow(subset_data), season :=
                    c('SPRING', 'FALL')[!c('SPRING', 'FALL') %in% unique_seasons]]

      subset_data[, dummy := T]
    }

    subset_data
  }


  pb <- txtProgressBar(max = length(yr_seq))

  cv_results <- lapply(yr_seq, function(.){

    setTxtProgressBar(pb, which(yr_seq == .))

    cv_fold <- full_data[year == .]

    bam.call <- model$call

    mod <- bam(formula = bam.call$formula,
               data = full_data[year < .],
               family = eval(bam.call$family),
               control = eval(bam.call$control),
               discrete = bam.call$discrete,
               nthreads = bam.call$nthreads)


    cv_fold <- data_prep(cv_fold)


    preds <- predict(mod, newdata = cv_fold,
                     type = 'response',
                     exclude = 's(yr_fac)')

    cv_fold <- setDT(cbind(
      cv_fold,
      pred = preds
    ))

    if('dummy' %in% names(cv_fold)){
      cv_fold <- cv_fold[-nrow(cvfold)]
    }

    cv_fold <- cv_fold[, .(rmse = sqrt(mean((abundance - pred)^2)),
                           prop_max = sqrt(mean((abundance - pred)^2)) / max(abundance),
                           fold = paste0('fold', .)),
                       by = abundance > 0]



  }
  )

  cv_results <- rbindlist(cv_results)
  cv_results

}


j <- progressing_cv(mod, subs, 2015:2019)

saveRDS(j, 'data derived/model output/nb_bam_abun_progcv.rds')
