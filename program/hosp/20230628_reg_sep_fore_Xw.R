




#
#--- SETUP ---
#




os = Sys.info()[["sysname"]]



load(list.files("data", "HK_flu_reg", full.names = T))
load(list.files("data", "HK_hosp_data", full.names = T))
source("program/20230628_reg_sep_fore_Xw_fx.R")





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

data.fx(data, f.yr, f.period, 10)
model.list.fx()





#
#--- REGRESSION ---
#





#-- 10-fold CV --

set.seed(555)
K = 10
CV.gp = sample(1:K, size = nrow(train.d), replace = T)
train.list = split(train.d, CV.gp)

# 10000 mc (valid)
# f.period = 4 ;  64 cores;  54m
# f.period = 13:  64 cores; 130m


start.t = Sys.time()

cl = makeCluster(n.cores, type = "SOCK")
registerDoSNOW(cl)



iter     = length(model.list)
pb       = txtProgressBar(max = iter, style = 3)
progress = function(n) setTxtProgressBar(pb, n)

cv.list = foreach(
  i             = 1:iter,
  .packages     = c("MASS", "Matrix"),
  .combine      = c,
  .options.snow = list(progress = progress)
) %dopar% {
  # i = 1
  
  model.i    = model.list[i]
  train.cost = valid.cost = list()
  alpha      = c(0.02, 0.05, seq(0.1, 0.9, 0.1))
  
  # start.t = Sys.time()
  
  for (k in 1:K) {
    
    # k=1
    train.k = do.call(rbind, train.list[-k])
    valid.k = train.list[[k]]
    
    # pred.train.k = pred.fx(model.i, train.k, train.k, alpha, keep.mc = F)
    pred.valid.k = pred.fx(model.i, train.k, valid.k, alpha, keep.mc = F)
    
    # train.cost.k = with(pred.train.k, cost.fx(pred.y.sum, alpha, "t") )
    valid.cost.k = with(pred.valid.k, cost.fx(pred.y.sum, alpha, "v") )
    
    # train.cost[[k]] = train.cost.k
    valid.cost[[k]] = valid.cost.k
    
  }
  
  # train.cost = do.call(rbind, train.cost)
  # train.cost = apply(train.cost, 2, mean)
  
  valid.cost = do.call(rbind, valid.cost)
  valid.cost = apply(valid.cost, 2, mean)

  cv.i = data.frame("model" = model.i)
  # cv.i = cbind(cv.i, rbind(train.cost), rbind(valid.cost))
  cv.i = cbind(cv.i, rbind(valid.cost))
  row.names(cv.i) = NULL
  
  return(list(cv.i))
  
}


stopCluster(cl)

end.t = Sys.time()
end.t - start.t

cv.table = do.call("rbind", cv.list)


save(cv.table, file = file.fx( sprintf("reg_sep_fore_%02dw_%d_CV_result.Rdata", f.period, f.yr) ) )





#
#--- PREDICTION ---
#





load(list.files("program", sprintf("reg_sep_fore_%02dw_%d_CV_result.Rdata", f.period, f.yr), full.names = T))
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
#--- END ---
#




