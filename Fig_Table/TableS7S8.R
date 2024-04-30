




#
#--- LOAD DATA ---
#




cv.f = list()
for (f.period in c(4, 13)) {
  # f.period = 4
  cv.table = read.csv(list.files("program/ILI/arima/", sprintf("arima_fore_%02dw_CV.csv", f.period), full.names = T))
  cv.table = cv.table[, c("model", "v.WIS", "v.RMSE", "v.LNQ", "v.MAE")]

  ILIpS.lag = sapply(
    strsplit(cv.table[, "model"], "log_ILIpS_lag_|-1"), 
    function(x) { suppressWarnings(max(as.numeric(x), na.rm = T)) } 
  )
  
  cv.table[, "model"] = stringr::str_replace(cv.table[, "model"], "log_ILIpS.*-1", sprintf("ILI(%s)", ILIpS.lag))
  cv.table[, "model"] = with(cv.table, gsub("t_lag_1_M_bs.*12", "bs(12)", model) )
  
  cv.table = dplyr::arrange(cv.table, v.WIS)
  for (k in c("v.WIS", "v.RMSE", "v.LNQ", "v.MAE")) {
    cv.table[,k] = with(cv.table, sprintf("%.3f", get(k)/min(get(k))))
  }
  
  cv.f[[paste(f.period)]] = rbind(f.period, cv.table)
  
}

cv.f = do.call(rbind, cv.f)

write.csv(cv.f, file = "TableS7_cv.csv", row.names = F)





cv.f = list()
for (f.period in c(4, 13)) {
  # f.period = 4
  cv.table = read.csv(list.files("program/hosp/arima/", sprintf("arima_fore_%02dw_CV.csv", f.period), full.names = T))
  cv.table = cv.table[, c("model", "v.WIS", "v.RMSE", "v.LNQ", "v.MAE")]
  
  ILIpS.lag = sapply(
    strsplit(cv.table[, "model"], "log_ILIpS_lag_|-1"), 
    function(x) { suppressWarnings(max(as.numeric(x), na.rm = T)) } 
  )
  
  cv.table[, "model"] = stringr::str_replace(cv.table[, "model"], "log_ILIpS.*-1", sprintf("ILI(%s)", ILIpS.lag))
  cv.table[, "model"] = with(cv.table, gsub("t_lag_1_M_bs.*12", "bs(12)", model) )
  
  cv.table = dplyr::arrange(cv.table, v.WIS)
  for (k in c("v.WIS", "v.RMSE", "v.LNQ", "v.MAE")) {
    cv.table[,k] = with(cv.table, sprintf("%.3f", get(k)/min(get(k))))
  }
  
  cv.f[[paste(f.period)]] = rbind(f.period, cv.table)
  
}

cv.f = do.call(rbind, cv.f)

write.csv(cv.f, file = "TableS8_cv.csv", row.names = F)





#
#
#--- END ---
#
#




