




#
#--- SETUP ---
#




os = Sys.info()[["sysname"]]



load(list.files("data", "HK_flu_reg", full.names = T))
source("program/ILI/20240404_reg_sep_fore_Xw_fx.R")
source("program/ILI/20240404_Rt_fx.R")





#
#--- RESULT ---
#





f.yr     = 2020
f.period = 13
load(list.files("program/ILI", sprintf("reg_sep_fore_%02dw_%d_forecast.Rdata", f.period, f.yr), full.names = T))




data.2020 = subset(data, WeekEnd >= as.Date("2019-12-01") & WeekEnd <= as.Date("2020-03-31"))
pred.sum  = f.list[[3]][["result"]][["WIS"]][["forec"]][["pred.y.sum"]]

data.2020[, "dec.date"] = decimal_date(data.2020[, "WeekEnd"])
pred.sum = transform(
  pred.sum,
  dec.date = decimal_date(f.WeekEnd),
  fit.LB   = fit.mean + qnorm(0.025) * fit.mean.se,
  fit.UB   = fit.mean + qnorm(0.975) * fit.mean.se,
  PI.LB    = `PI.2.5%`,
  PI.UB    = `PI.97.5%`
)



current.w = subset(data.2020, WeekEnd == (min(pred.sum[, "f.WeekEnd"]) - 7))
current.w = with(
  current.w, data.frame(
    "WeekEnd"   = WeekEnd,
    "dec.date"  = decimal_date(WeekEnd),
    "fore"      = 0,
    "f.WeekEnd" = WeekEnd,
    "fit.mean"  = ILIpS,
    "fit.LB"    = ILIpS,
    "fit.UB"    = ILIpS,
    "PI.mean"   = ILIpS,
    "PI.LB"     = ILIpS,
    "PI.UB"     = ILIpS
  )
)


pred.sum = plyr::rbind.fill(current.w, pred.sum)







inc.data = merge(
  x     = data.2020[, c("WeekEnd", "ILIpS")],
  y     = pred.sum[, c("f.WeekEnd", "fit.mean")],
  by.x  = "WeekEnd",
  by.y  = "f.WeekEnd",
  all.x = T
)

inc.data = transform(
  inc.data, 
  fit.mean = ifelse(is.na(fit.mean), ILIpS, fit.mean)
)





#
#--- RESULT ---
#





Rt.compare.fx = function(pred.Rt) {
  
  for (i in c("Cori.obs.Rt", "Cori.fit.Rt")) { # }, "WT.obs.Rt", "WT.fit.Rt")) {
    
    Rt.i = lapply(pred.Rt, function(x) { x[[i]] } )
    Rt.i = do.call(rbind, Rt.i)
    Rt.i = split(Rt.i, Rt.i[, "w"])
    
    Rt.mc  = list()
    Rt.sum = list()
    
    
    for (w in 1:length(Rt.i)) {
      
      Rt.i.w = Rt.i[[w]]
      
      Rt.mc.w = sapply(
        1:nrow(Rt.i.w), 
        function(j) { set.seed(j); with(Rt.i.w, rlnorm(1, log.mu[j], log.sd[j]) ) } 
      )
      
      Rt.sum.w = data.frame(
        "weekend"  = unique(Rt.i.w[, "weekend"]),
        "dec.date" = unique(Rt.i.w[, "dec.date"]),
        "mu"       = mean(Rt.mc.w),
        "sd"       = sd(Rt.mc.w),
        "LB"       = quantile(Rt.mc.w, 0.025),
        "UB"       = quantile(Rt.mc.w, 0.975)
      )
      
      Rt.mc[[w]]  = Rt.mc.w
      Rt.sum[[w]] = Rt.sum.w
      
    }
    
    Rt.sum = do.call(rbind, Rt.sum)
    rownames(Rt.sum) = NULL
    
    assign(sprintf("%s.mc",  i), Rt.mc,  envir = .GlobalEnv)
    assign(sprintf("%s.sum", i), Rt.sum, envir = .GlobalEnv)
    
  }

  Rt.compare = lapply(
    c("Cori"), # , "WT"),
    function(k) {
      
      obs.Rt.sum = get(sprintf("%s.obs.Rt.sum", k))
      obs.Rt.k   = get(sprintf("%s.obs.Rt.mc", k))
      fit.Rt.sum = get(sprintf("%s.fit.Rt.sum", k))
      fit.Rt.k   = get(sprintf("%s.fit.Rt.mc", k))
      
      obs.Rt.sum = obs.Rt.sum[1:13, ]
      fit.Rt.sum = fit.Rt.sum[1:13, ]
      obs.Rt.k   = do.call(rbind, obs.Rt.k[1:13])
      fit.Rt.k   = do.call(rbind, fit.Rt.k[1:13])
      
      Rt.redu = (1 - obs.Rt.k / fit.Rt.k) * 100
      Rt.redu = apply(Rt.redu, 1, CI95.fx)
      Rt.redu = do.call(rbind, Rt.redu)
      
      Rt.redu = data.frame(
        "obs.Rt"  = with(obs.Rt.sum, sprintf("%.2f (%.2f, %.2f)",   mu,     LB,      UB) ),
        "fit.Rt"  = with(fit.Rt.sum, sprintf("%.2f (%.2f, %.2f)",   mu,     LB,      UB) ),
        "Rt.redu" = with(Rt.redu,    sprintf("%.2f%% (%.2f%%, %.2f%%)", mean, `2.5%`, `97.5%`) )
      )
      colnames(Rt.redu) = paste0(k, ".", colnames(Rt.redu))
      
      
      Rt.redu = cbind(
        obs.Rt.sum[1:13, c("weekend", "dec.date")],
        Rt.redu
      )
      
      return(Rt.redu)
    }
  )
  
  Rt.compare = Rt.compare[[1]]
  return(Rt.compare)
}









Rt.compare = lapply(
  1:4,
  function(f.start) {
    # f.start = 3
    load(list.files("program/ILI", "forecast_RtPI.Rdata", full.names = T))
    RtPI.compare.f = Rt.compare.fx(pred.Rt[[f.start]])
    load(list.files("program/ILI", "forecast_RtCI.Rdata", full.names = T))
    RtCI.compare.f = Rt.compare.fx(pred.Rt[[f.start]])

    Rt.compare.f = merge(
      x   = RtPI.compare.f[, 1:4],
      y   = RtCI.compare.f[, c(1,2,5)],
      by  = c("weekend", "dec.date"),
      all = T
    )
    colnames(Rt.compare.f) = c("WeekEnd", "dec.date", "Obs.Rt", sprintf("%s.Rt.%d", c("CounterPI", "ReducedCI"), f.start))
    Rt.compare.f = subset(Rt.compare.f, dec.date > 2020)
    return(Rt.compare.f)  
  }
)

Rt.compare = Reduce(function(x, y) merge(x, y, by = c("WeekEnd", "dec.date", "Obs.Rt"), all = T), Rt.compare)  





write.csv( Rt.compare, file = file.fx("TableS3.csv"), row.names = F )





#
#--- END ---
#




