




#
#--- SETUP ---
#





load(list.files("data", "HK_flu_reg", full.names = T))
source("program/ILI/20240404_reg_sep_fore_Xw_fx.R")
library("forecast")
setwd("program/ILI/arima")



f.yr = 2020
f.period = 4
data.fx(data, f.yr, f.period, 26)
model.list.fx()

data.2 = rbind(train.d, forec.d)
data.2 = arrange(data.2, WeekEnd)



train.ed = seq(as.Date("2017-01-14"), as.Date("2019-12-31")-7*f.period, by = "week")
train.st = seq(as.Date("2010-01-14"), by = "week", length.out = length(train.ed))
valid.ed = train.ed + 7*f.period

CV.w = data.frame("t.st" = train.st, "t.ed" = train.ed, "v.ed" = valid.ed)

alpha = c(0.02, 0.05, seq(0.1, 0.9, 0.1))





#
#--- RESULT FROM GLM ---
#





load(list.files("..", sprintf("reg_sep_fore_%02dw_%d_CV_result.Rdata", f.period, f.yr), full.names = T))
cv.table = aggregate(. ~ model, data = cv.table, FUN = mean)
cv.table = subset(cv.table, select = -tsCV.ind)
cv.table = criteria.rank.fx(cv.table, rank.w.value = F)

( RMSE.model = with(cv.table, model[which.min(v.RMSE)]) )
( LNQ.model  = with(cv.table, model[which.min(v.LNQ )]) )
( MAE.model  = with(cv.table, model[which.min(v.MAE )]) ) 
( WIS.model  = with(cv.table, model[which.min(v.WIS )]) ) 

reg.x = unlist( strsplit(c(RMSE.model, LNQ.model, MAE.model, WIS.model), " [+] ") )
reg.x = reg.x[!grepl("ILIpS", reg.x)]
reg.x = unique(reg.x)



tCV.WIS.fx = function(y, ts.f) {
  
  Q50 = t(rbind(ts.f[["mean"]]))
  LB  = ts.f[["lower"]]; LB = LB[, c(ncol(LB):1)]
  UB  = ts.f[["upper"]]
  
  q.sum = cbind(LB, Q50, UB); colnames(q.sum) = sprintf("PI.%d%%", 1:ncol(q.sum))
  q.sum = as.data.frame(q.sum)
  q.sum = cbind(data.frame("f.obs" = y), q.sum)
  
  WIS.t = WIS.fx(q.sum, alpha)
  return(WIS.t)
}





#
#--- ARIMA ---
#





start.t = Sys.time()

cl = makeCluster(n.cores, type = "SOCK")
registerDoSNOW(cl)

arima.cv = foreach(
  i         = 1:nrow(CV.w),
  .packages = "forecast",
  .combine  = "rbind"
) %dopar% {
  
  # i = 1
  # i = nrow(CV.w)
  t.st.i = CV.w[i, "t.st"]
  t.ed.i = CV.w[i, "t.ed"]
  v.ed.i = CV.w[i, "v.ed"]
  
  train.i = subset(data.2, WeekEnd  > t.st.i & WeekEnd <= t.ed.i)
  valid.i = subset(data.2, WeekEnd  > t.ed.i & WeekEnd <= v.ed.i)

  
  
  t.reg.i = as.matrix( train.i[, reg.x] )
  v.reg.i = as.matrix( valid.i[, reg.x] )

  t.reg.i = t.reg.i[, !grepl("M_bs", colnames(t.reg.i))]
  v.reg.i = v.reg.i[, !grepl("M_bs", colnames(v.reg.i))]

  CV.fit.i = auto.arima(
    y         = ts(train.i[, "ILIpS"], frequency = 52), 
    max.p     = 26,
    max.q     = 26,
    lambda    = 0, 
    xreg      = t.reg.i
  )

  AR.i  = CV.fit.i[["arma"]][1]
  MA.i  = CV.fit.i[["arma"]][2]
  sAR.i = CV.fit.i[["arma"]][3]
  sMA.i = CV.fit.i[["arma"]][4]
  D.i   = CV.fit.i[["arma"]][6]
  sD.i  = CV.fit.i[["arma"]][7]

  CV.f.i   = forecast(object = CV.fit.i, h = f.period, xreg = v.reg.i, level = rev(1-alpha) * 100)
  
  
  
  cv.i = data.frame(
    "i"    = i, "f" = f.period,
    "AR"   = AR.i,  "D"  = D.i,  "MA"  = MA.i,
    "sAR"  = sAR.i, "sD" = sD.i, "sMA" = sMA.i,
    "bs"   = 0,
    "WIS"  = mean( tCV.WIS.fx(valid.i[, "ILIpS"], CV.f.i) ),
    "RMSE" = sqrt(mean( (    valid.i[, "ILIpS"]  -     CV.f.i[["mean"]] )^2 )),
    "LNQ"  = sqrt(mean( (log(valid.i[, "ILIpS"]) - log(CV.f.i[["mean"]]))^2 )),
    "MAE"  = mean( abs(      valid.i[, "ILIpS"]  -     CV.f.i[["mean"]]     ))
  )
  
  
  
  if (D.i == 0) {
    
    t.reg.i = as.matrix( train.i[, reg.x] )
    v.reg.i = as.matrix( valid.i[, reg.x] )
    n.x     = ncol(t.reg.i)
    
    order.i    = c(AR.i, D.i, MA.i)
    seasonal.i = c(sAR.i, sD.i, sMA.i)
    n.coef     = sum(order.i[c(1,3)], seasonal.i[c(1,3)]) + n.x
    
    ndeps   = 1e-6
    fit.ind = T
    
    while (fit.ind) {
      CV.fit.i2 = try(Arima(
        y         = ts(train.i[, "ILIpS"], frequency = 52), 
        order     = order.i,
        seasonal  = seasonal.i,
        lambda    = 0, 
        xreg      = t.reg.i, 
        include.mean = F,
        optim.control = list("ndeps" = rep(ndeps, n.coef))
      ), silent = T
      )
      fit.ind = (length(CV.fit.i2) == 1) & (ndeps < 1)
      ndeps   = ndeps * 10^0.1
    }

    CV.f.i2   = forecast(object = CV.fit.i2, h = f.period, xreg = v.reg.i, level = rev(1-alpha) * 100)

    cv.i = rbind(
      cv.i,
      data.frame(
        "i"    = i, "f" = f.period,
        "AR"   = AR.i,  "D"  = D.i,  "MA"  = MA.i,
        "sAR"  = sAR.i, "sD" = sD.i, "sMA" = sMA.i,
        "bs"   = 1,
        "WIS"  = mean( tCV.WIS.fx(valid.i[, "ILIpS"], CV.f.i2) ),
        "RMSE" = sqrt(mean( (    valid.i[, "ILIpS"]  -     CV.f.i2[["mean"]] )^2 )),
        "LNQ"  = sqrt(mean( (log(valid.i[, "ILIpS"]) - log(CV.f.i2[["mean"]]))^2 )),
        "MAE"  = mean( abs(      valid.i[, "ILIpS"]  -     CV.f.i2[["mean"]]     ))
      )
    )
    
  }

  return(cv.i)
}

stopCluster(cl)

end.t = Sys.time()
end.t - start.t


save(
  arima.cv, 
  file = file.fx(sprintf("arima_fore_%02dw_CV.Rdata", f.period))
)





#
#--- SELECT ARIMA MODEL ---
#





file = list.files(, sprintf("arima_fore_%02dw_CV.*Rdata", f.period), full.names = T)
load(file)

arima.cv = split(arima.cv, ~ i)
arima.cv = lapply(arima.cv, function(x) { y = arrange(x, WIS)[1, ]; return(y) } )
arima.cv = do.call(rbind, arima.cv)

arima.model = aggregate(len ~ AR+D+MA+sAR+sD+sMA+bs, data = cbind(arima.cv, len = 1), FUN = sum)
arima.model = lapply(1:5, function(i) { as.numeric(arrange(arima.model, -len)[i, 1:7]) } )
(arima.model)





start.t = Sys.time()

cl = makeCluster(n.cores, type = "SOCK")
registerDoSNOW(cl)

arima.cv = foreach(
  k         = 1:length(arima.model),
  .packages = "forecast",
  .combine  = "rbind"
) %:% foreach(
  i         = 1:nrow(CV.w),
  .combine  = "rbind"
) %dopar% {
  
  # k = 1; i = 1
  t.st.i = CV.w[i, "t.st"]
  t.ed.i = CV.w[i, "t.ed"]
  v.ed.i = CV.w[i, "v.ed"]
  
  train.i = subset(data.2, WeekEnd  > t.st.i & WeekEnd <= t.ed.i)
  valid.i = subset(data.2, WeekEnd  > t.ed.i & WeekEnd <= v.ed.i)

  t.reg.i = as.matrix( train.i[, reg.x] )
  v.reg.i = as.matrix( valid.i[, reg.x] )
  
  
  
  order.k    = arima.model[[k]][1:3]
  seasonal.k = arima.model[[k]][4:6]
  bs.k       = arima.model[[k]][7]
  
  if (bs.k == 0) {
    t.reg.i = t.reg.i[, !grepl("M_bs", colnames(t.reg.i))]
    v.reg.i = v.reg.i[, !grepl("M_bs", colnames(v.reg.i))]
  }
  n.x = ncol(t.reg.i)
  
  
  
  mu.ind = ifelse(bs.k == 1 | order.k[2] > 0, F, T)
  n.coef = sum(order.k[c(1,3)], seasonal.k[c(1,3)]) + n.x + 1*mu.ind
  
  ndeps   = 1e-6
  fit.ind = T
  
  while (fit.ind) {
    CV.fit.ki = try(Arima(
      y             = ts(train.i[, "ILIpS"], frequency = 52), 
      order         = order.k,
      seasonal      = seasonal.k,
      lambda        = 0, 
      xreg          = t.reg.i, 
      include.mean  = mu.ind,
      optim.control = list("ndeps" = rep(ndeps, n.coef))
    ), silent = T
    )
    fit.ind = (length(CV.fit.ki) == 1) & (ndeps < 1)
    ndeps   = ndeps * 10^0.1
  }
  
  AR.k  = CV.fit.ki[["arma"]][1]
  MA.k  = CV.fit.ki[["arma"]][2]
  sAR.k = CV.fit.ki[["arma"]][3]
  sMA.k = CV.fit.ki[["arma"]][4]
  D.k   = CV.fit.ki[["arma"]][6]
  sD.k  = CV.fit.ki[["arma"]][7]
  
  CV.f.ki   = forecast(object = CV.fit.ki, h = f.period, xreg = v.reg.i, level = rev(1-alpha) * 100)
  
  
  
  cv.ki = data.frame(
    "i"    = i, "f" = f.period,
    "AR"   = AR.k,  "D"  = D.k,  "MA"  = MA.k,
    "sAR"  = sAR.k, "sD" = sD.k, "sMA" = sMA.k,
    "bs"   = bs.k,
    "WIS"  = mean( tCV.WIS.fx(valid.i[, "ILIpS"], CV.f.ki) ),
    "RMSE" = sqrt(mean( (    valid.i[, "ILIpS"]  -     CV.f.ki[["mean"]] )^2 )),
    "LNQ"  = sqrt(mean( (log(valid.i[, "ILIpS"]) - log(CV.f.ki[["mean"]]))^2 )),
    "MAE"  = mean( abs(      valid.i[, "ILIpS"]  -     CV.f.ki[["mean"]]     ))
  )
  return(cv.ki)
}

stopCluster(cl)

end.t = Sys.time()
end.t - start.t


save(
  arima.cv, 
  file = file.fx(sprintf("arima_fore_%02dw_CV.Rdata", f.period))
)





#
#--- COMPARE ---
#





file = list.files(, sprintf("arima_fore_%02dw_CV.*Rdata", f.period), full.names = T)
load(file)

arima.cv    = aggregate(cbind(WIS, RMSE, LNQ, MAE) ~ AR+D+MA+sAR+sD+sMA+bs, data = arima.cv, FUN = mean)

arima.cv = with(
  arima.cv,
  data.frame(
    "model"  = sprintf("ARIMA (%d,%d,%d)(%d,%d,%d)+%s", AR,D,MA,sAR,sD,sMA, paste(reg.x, collapse = "+")), 
    "v.RMSE" = RMSE, 
    "v.LNQ"  = LNQ, 
    "v.MAE"  = MAE, 
    "v.WIS"  = WIS)
)
cv.table = rbind(cv.table, arima.cv)
rownames(cv.table) = NULL

write.csv(
  cv.table, 
  file = file.fx(sprintf("arima_fore_%02dw_CV.csv", f.period))
)





#
#--- END ---
#




