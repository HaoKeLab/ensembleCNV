
output_gatk_result_LRR <- function(dt_LRRBAF) {

  dt_LRRBAF$LRRBAF0 <- -10*log(dt_LRRBAF$LRR0)/log(10)
  dt_LRRBAF$LRRBAF1 <- -10*log(dt_LRRBAF$LRR1)/log(10)  ## save V1
  dt_LRRBAF$LRRBAF2 <- -10*log(dt_LRRBAF$LRR2)/log(10)
  dt_LRRBAF$LRRBAF3 <- -10*log(dt_LRRBAF$LRR3)/log(10)
  
  ## must be deal with
  
  dt_sub <- dt_LRRBAF[, c("LRRBAF0", "LRRBAF1", "LRRBAF2", "LRRBAF3")]
  
  value_GQs <- unlist(lapply(1:nrow(dt_sub), FUN = function(k) {
    v1 <- unlist(dt_sub[k, ])
    v1 <- sort(v1)
    gq1 <- v1[2] - v1[1]
    gq1
  }))
  
  # mean(GQs)
  dt_LRRBAF$value_GQ <- value_GQs
  
  CN_gatk_preds <- unlist(lapply(1:nrow(dt_sub), FUN = function(k) {
    v1 <- unlist(dt_sub[k, ])
    idx1 <- which.min(v1)
    return(idx1 - 1)
  }))
  
  dt_LRRBAF$CN_gatk_pred <- CN_gatk_preds
  
  ## add BAF0 1 2 3 column
  dt_LRRBAF$BAF0 = NA
  dt_LRRBAF$BAF1 = NA
  dt_LRRBAF$BAF2 = NA
  dt_LRRBAF$BAF3 = NA
  
  dt_LRRBAF # return result
}

# calculate gatk result
output_gatk_result <- function(dt_LRRBAF) {
  
  dt_LRRBAF$LRRBAF0 <- -10*log(dt_LRRBAF$LRR0*dt_LRRBAF$BAF0)/log(10)
  dt_LRRBAF$LRRBAF1 <- -10*log(dt_LRRBAF$LRR1*dt_LRRBAF$BAF1)/log(10)  ## save V1
  dt_LRRBAF$LRRBAF2 <- -10*log(dt_LRRBAF$LRR2*dt_LRRBAF$BAF2)/log(10)
  dt_LRRBAF$LRRBAF3 <- -10*log(dt_LRRBAF$LRR3*dt_LRRBAF$BAF3)/log(10)
  
  dt_sub <- dt_LRRBAF[, c("LRRBAF0", "LRRBAF1", "LRRBAF2", "LRRBAF3")]
  
  value_GQs <- unlist(lapply(1:nrow(dt_sub), FUN = function(k) {
    v1 <- unlist(dt_sub[k, ])
    v1 <- sort(v1)
    gq1 <- v1[2] - v1[1]
    gq1
  }))
  
  # mean(GQs)
  dt_LRRBAF$value_GQ <- value_GQs
  
  CN_gatk_preds <- unlist(lapply(1:nrow(dt_sub), FUN = function(k) {
    v1 <- unlist(dt_sub[k, ])
    idx1 <- which.min(v1)
    return(idx1 - 1)
  }))
  
  dt_LRRBAF$CN_gatk_pred <- CN_gatk_preds
  
  dt_LRRBAF # return result
}
