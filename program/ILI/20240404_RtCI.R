




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
#--- RT ---
#





start.t = Sys.time()

cl = makeCluster(n.cores, type = "SOCK")
registerDoSNOW(cl)



N = nrow( f.list[[1]][["result"]][["WIS"]][["forec"]][["pred.y.sum"]] )

pred.Rt = foreach(
  f.start   = 1:4, 
  .packages = c("EpiEstim", "lubridate", "pracma"),
  .options.snow = pb.fx(4 * N)
) %:% foreach(
  i         = 1:N
) %dopar% {
  
  
  pred.mu.mc = f.list[[f.start]][["result"]][["WIS"]][["forec"]][["pred.y.sum"]]
  pred.mu.mc = sapply(
    1:nrow(pred.mu.mc),
    function(i) { with(pred.mu.mc, rnorm(1e4, fit.mean[i], fit.mean.se[i])) }
  )
  
  inc.data.i = inc.data
  ind = (4 + f.start):nrow(inc.data.i)
  inc.data.i[ind, "fit.mean"] = pred.mu.mc[i, 1:length(ind)]
  
  daily.inc.i = merge(
    x  = daily.incid.fx(inc.data.i, "ILIpS"),
    y  = daily.incid.fx(inc.data.i, "fit.mean")
  )
  
  
  
  # Rt setup
  t.start = 2:(nrow(daily.inc.i)-6)
  t.end   = t.start + 6
  # https://academic.oup.com/aje/article/180/9/865/2739204
  SI.mu   = 2.8
  SI.sd   = 1.2
  
  Rt.config = list(
    t_start    = t.start,
    t_end      = t.end,
    mean_si    = SI.mu, 
    std_si     = SI.sd, 
    mean_prior = 1, 
    std_prior  = 1000
  )
  
  
  
  Cori.obs.Rt.i = w.Rt.fx(daily.inc.i, "ILIpS",    "Cori", Rt.config)
  Cori.fit.Rt.i = w.Rt.fx(daily.inc.i, "fit.mean", "Cori", Rt.config)
  # WT.obs.Rt.i   = w.Rt.fx(daily.inc.i, "ILIpS",    "WT", Rt.config)
  # WT.fit.Rt.i   = w.Rt.fx(daily.inc.i, "fit.mean", "WT", Rt.config)
  
  
  output = list(
    "Cori.obs.Rt" = Cori.obs.Rt.i,
    "Cori.fit.Rt" = Cori.fit.Rt.i
    # "WT.obs.Rt"   = WT.obs.Rt.i,
    # "WT.fit.Rt"   = WT.fit.Rt.i
  )
  
  return(output)
}



stopCluster(cl)

end.t = Sys.time()

cat("\n")
print(end.t - start.t)



save(
  pred.Rt,
  file = file.fx(sprintf("reg_sep_fore_%02dw_%d_forecast_RtCI.Rdata", f.period, f.yr))
)





#
#--- RESULT ---
#





load(list.files("program/ILI", "forecast_RtCI.Rdata", full.names = T))

f.start = 3
pred.Rt = pred.Rt[[f.start]]



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
# Reduce(function(x, y) merge(x, y, all = T), Rt.compare)

write.csv(
  Rt.compare, 
  file = file.fx(sprintf("reg_sep_fore_13w_2020_forecast_%dw_RtCI.csv", f.start)), 
  row.names = F
)





#
#--- END ---
#




