




#
#--- SETUP ---
#




os = Sys.info()[["sysname"]]
if (os == "Linux") { setwd("DanielLau/regression/") }



load(list.files(pattern = "HK_flu_reg", full.names = T))
load(list.files(pattern = "hosp_data", full.names = T))
setwd("20240415/")
source("program/20240415_reg_sep_fore_Xw_fx.R")





#
#--- DATA INPUT ---
#





f.yr = 2020
f.period = 4

data = merge(
  x     = data,
  y     = hosp.data,
  by.x  = c("WeekEnd", "Year", "Week"),
  by.y  = c("weekend", "year", "week"),
  all.y = T
)

data[, "ILIpS"] = round(data[, "admission.1e6"])

data.fx(data, f.yr, f.period, 26)
model.list.fx()

data.2 = rbind(train.d, forec.d)
data.2 = arrange(data.2, WeekEnd)





#
#--- REGRESSION ---
#





#-- tsCV --
# 13: 2014-04-01; 26: 2014-07-01
train.ed = seq(as.Date("2018-01-01"), as.Date("2019-12-31")-7*f.period, by = "week")
train.st = seq(as.Date("2014-07-01"), by = "week", length.out = length(train.ed))
CV.w     = data.frame("t.st" = train.st, "t.ed" = train.ed)





start.t = Sys.time()

cl = makeCluster(n.cores, type = "SOCK")
registerDoSNOW(cl)



alpha    = c(0.02, 0.05, seq(0.1, 0.9, 0.1))
iter     = length(model.list)
pb       = txtProgressBar(max = iter * nrow(CV.w), style = 3)
progress = function(n) setTxtProgressBar(pb, n)


cv.table = foreach(
  i             = 1:iter,
  .packages     = c("MASS"),
  .combine      = "rbind",
  .options.snow = list(progress = progress)
) %:% foreach (
  k             = 1:nrow(CV.w),
  .combine      = "rbind"
) %dopar% {
  # i = iter; k = 6
  
  model.i    = model.list[i]
  
  t.st.k = CV.w[k, "t.st"]
  t.ed.k = CV.w[k, "t.ed"]
  
  train.k = subset(data.2, WeekEnd  > t.st.k & WeekEnd <= t.ed.k)
  valid.k = subset(data.2, WeekEnd  > t.ed.k & WeekEnd <= t.ed.k + 7)
  
  pred.valid.k = pred.fx(model.i, train.k, valid.k, alpha, keep.mc = F)
  valid.cost.k = with(pred.valid.k, cost.fx(pred.y.sum, alpha, "v") )
  
  cv.ik = data.frame("model" = model.i, "tsCV.ind" = k)
  cv.ik = cbind(cv.ik, rbind(valid.cost.k))
  
  return(cv.ik)
  
}



stopCluster(cl)

end.t = Sys.time()
end.t - start.t

save(cv.table, file = file.fx( sprintf("reg_sep_fore_%02dw_%d_CV_result.Rdata", f.period, f.yr) ) )





#
#--- PREDICTION ---
#





( linux.file = list.files("program/Xw_linux", sprintf("%02dw.*CV_result.*Rdata", f.period), full.names = T) )
cv.table = lapply(linux.file, function(x) { load(x); return(cv.table) } )
cv.table = do.call(rbind, cv.table)
save(cv.table, file = file.fx( sprintf("reg_sep_fore_%02dw_%d_CV_result.Rdata", f.period, f.yr) ) )



load(list.files(, sprintf("reg_sep_fore_%02dw_%d_CV_result.Rdata", f.period, f.yr), full.names = T))

cv.table = aggregate(. ~ model, data = cv.table, FUN = mean)
cv.table = subset(cv.table, select = -tsCV.ind)
cv.table = criteria.rank.fx(cv.table, rank.w.value = F)

( RMSE.model = with(cv.table, model[which.min(v.RMSE     )]) )
( LNQ.model  = with(cv.table, model[which.min(v.LNQ      )]) )
( MAE.model  = with(cv.table, model[which.min(v.MAE      )]) ) 
( WIS.model  = with(cv.table, model[which.min(v.WIS      )]) ) 





#
#--- FORECAST ---
#




# 10000 mc
# f.period = 4:   <1  mins
# f.period = 13:   2  mins 



start.t = Sys.time()


cl = makeCluster(n.cores, type = "SOCK")
registerDoSNOW(cl)

start.d = seq(as.Date(sprintf("%s-01-01", f.yr)), by = "week", length.out = 8)

f.list = foreach(
  t = 1:length(start.d),
  .packages = c("MASS", "Matrix"),
  .export  = sprintf("%s.model", c("RMSE", "MAE", "LNQ", "WIS")),
  .combine = c
) %dopar% {

  # t = 1
  f.start = start.d[[t]]

  data.2  = rbind(train.d, forec.d)
  train.t = train.d
  forec.t = head(subset(data.2, WeekEnd > f.start), 1)
  
  alpha   = c(0.02, 0.05, seq(0.1, 0.9, 0.1))
  

    
  f.list.t = list()
  print(sprintf("Forecasting progress: %s / %s", t, length(start.d)) )
  
  for (j in c("RMSE", "MAE", "LNQ", "WIS")) {
    # j = "RMSE"
    print(sprintf(" - Working on model: %s", j) )
    model.j  = get(sprintf("%s.model", j))
    
    pred.t.j = pred.fx(model.j, train.t, train.t, alpha, keep.mc = F) 

    fore.t.j = pred.fx(model.j, train.t, forec.t, alpha, keep.mc = T)
    fore.t.j = within(
      fore.t.j, {
        pred.y.sum = do.call(rbind, pred.y.sum)
        pred.y.mc  = do.call(cbind, pred.y.mc )
      }
    )
    
    Sum.j = list(
      "type"  = j,
      "model" = model.j, 
      "pred"  = pred.t.j,
      "forec" = fore.t.j
    )
    
    f.list.t[[j]] = Sum.j
  }
  
  return(
    list( list("start.d" = f.start, "result"  = f.list.t) )
  )
  
}

stopCluster(cl)

end.t = Sys.time()
end.t - start.t

save(f.list, file = file.fx( sprintf("reg_sep_fore_%02dw_%d_forecast.Rdata", f.period, f.yr) ) )





#
#--- RESULT ---
#





load(list.files(, sprintf("reg_sep_fore_%02dw_%d_forecast.Rdata", f.period, f.yr), full.names = T))



pdf(
  width = 16, height = 9, 
  file = file.fx(sprintf("reg_sep_fore_%02dw_%s_forecast_plot.pdf", f.period, f.yr))
)

# windows(width = 16, height = 12)
# quartz( width = 16, height = 12)

par(mfcol = c(4, 2), mar = c(3,4,3,1), las = 1)

for (j in c("RMSE", "MAE", "LNQ", "WIS")) { 
  
  for (t in 1:length(f.list)) {

    # t=1; j="RMSE"
    f.t       = f.list[[t]]
    forec.sum = f.t[["result"]][[j]][["forec"]][["pred.y.sum"]]
    start.d   = f.t[["start.d"]]
    
    title.j.t = paste(
      sprintf("Forecast for %d weeks ahead (Fitted model: %s criteria)", f.period, j)
    )
    
    fore.plot.fx(
      forec.sum = forec.sum, 
      title     = title.j.t, 
      xlim.yr   = c(2019, 2021)
    )
    
  }
}

dev.off()

  



#
#--- END ---
#




