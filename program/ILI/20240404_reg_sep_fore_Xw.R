




#
#--- SETUP ---
#




os = Sys.info()[["sysname"]]



load(list.files("data", "HK_flu_reg", full.names = T))
source("program/ILI/20240404_reg_sep_fore_Xw_fx.R")





#
#--- DATA INPUT ---
#





f.yr = 2020
f.period = 4
data.fx(data, f.yr, f.period, 26)
model.list.fx()

data.2 = rbind(train.d, forec.d)
data.2 = arrange(data.2, WeekEnd)






#
#--- REGRESSION ---
#





#-- tsCV --

train.ed = seq(as.Date("2017-01-14"), as.Date("2019-12-31")-7*f.period, by = "week")
train.st = seq(as.Date("2010-01-14"), by = "week", length.out = length(train.ed))
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





load(list.files("program/ILI", sprintf("reg_sep_fore_%02dw_%d_CV_result.Rdata", f.period, f.yr), full.names = T))

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





forec.st = seq(as.Date(sprintf("%s-01-01", f.yr)), by = "week", length.out = 8)
train.st = zoo::as.Date( sapply(forec.st, function(d) { seq(d, by = "-7 years", length.out = 2)[2] } ) )
forec.w  = data.frame("t.st" = train.st, "t.ed" = forec.st)





start.t = Sys.time()

cl = makeCluster(n.cores, type = "SOCK")
registerDoSNOW(cl)



alpha   = c(0.02, 0.05, seq(0.1, 0.9, 0.1))

f.list = foreach(
  t = 1:nrow(forec.w),
  .packages = "MASS",
  .export  = sprintf("%s.model", c("RMSE", "MAE", "LNQ", "WIS")),
  .combine = c
) %dopar% {

  # t = 1
  t.st.t = forec.w[t, "t.st"]
  t.ed.t = forec.w[t, "t.ed"]
  
  train.t = subset(data.2, WeekEnd  > t.st.t & WeekEnd <= t.ed.t)
  forec.t = subset(data.2, WeekEnd  > t.ed.t & WeekEnd <= t.ed.t + 7)

  
    
  f.list.t = list()

  for (j in c("RMSE", "MAE", "LNQ", "WIS")) {
    # j = "RMSE"
    # print(sprintf(" - Working on model: %s", j) )
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
    list( list("start.d" = t.ed.t, "result"  = f.list.t) )
  )
  
}

stopCluster(cl)

end.t = Sys.time()
end.t - start.t

save(f.list, file = file.fx( sprintf("reg_sep_fore_%02dw_%d_forecast.Rdata", f.period, f.yr) ) )

  



#
#--- END ---
#




