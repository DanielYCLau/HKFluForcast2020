




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





f.start = 3

load(list.files("program/ILI", "forecast_RtPI.Rdata", full.names = T))
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
  
  assign(sprintf("%s.PI.mc",  i), Rt.mc,  envir = .GlobalEnv)
  assign(sprintf("%s.PI.sum", i), Rt.sum, envir = .GlobalEnv)
  
}





load(list.files("program/ILI", "forecast_RtCI.Rdata", full.names = T))
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
  
  assign(sprintf("%s.CI.mc",  i), Rt.mc,  envir = .GlobalEnv)
  assign(sprintf("%s.CI.sum", i), Rt.sum, envir = .GlobalEnv)
  
}



#
#--- PLOT ---
#






pdf(width = 7, height = 5, file = file.fx("FigS6_Rt.pdf"))

# quartz(width = 7, height = 5)
# windows(width = 7, height = 5)

par(mar = c(3,6,1,6), las = 1)

x.lim = decimal_date(c(as.Date("2019-12-01"), as.Date("2020-04-01")))



for (k in c("Cori")) { # "WT", "Cori")) {
  plot(NULL, xlim = x.lim, ylim = c(0, 5e4), yaxt = "n", xaxt = "n", ylab = "", xlab = "")
  
  axis.date = seq(as.Date("2019-12-01"), as.Date("2020-04-01"), by = "month")
  axis.t    = decimal_date(axis.date)
  axis(1, at = axis.t, labels = NA)
  axis(1, at = head(axis.t, -1) + 0.5 * diff(axis.t), labels = format(head(axis.date, -1), "%b %Y"), tick = F)
  
  axis(2, at = seq(0, 3e4, 1e4), labels = seq(0, 0.03, 0.01))
  
  
  lines(ILIpS    ~ dec.date, data = data.2020, lwd = 5)
  lines(fit.mean ~ dec.date, data = pred.sum,  col = col.fx(1), lty = 2, lwd = 2)
  with(pred.sum, polygon.fx(t = dec.date, LB = fit.LB, UB = fit.UB, col = col.fx(1, 0.1)))
  with(pred.sum, polygon.fx(t = dec.date, LB = PI.LB,  UB = PI.UB,  col = col.fx(1, 0.1)))
  
  abline(v = axis.t, col = alpha(1, 0.1))
  
  
  
  mtext(
    text = "Influenza Activty (ILI+ proxy)", 
    font = 2,
    side = 2, 
    line = 4, 
    las  = 0
  )
  
  
  
  
  par(new = T)
  
  # title  = sprintf("%s Rt comparison\n(Forecast based on WIS model)", ifelse(k == "WT", "Case", "Instantaneous") )
  # title  = sprintf("%s Rt", ifelse(k == "WT", "Case", "Instantaneous") )
  obs.Rt = get(sprintf("%s.obs.Rt.PI.sum", k))[1:11, ]
  fit.Rt = head(get(sprintf("%s.fit.Rt.PI.sum", k)), -1)
  
  plot(NULL, xlim = x.lim, ylim = c(-2, 2), axes = F, ylab = "", xlab = "")
  axis(4, at = seq(0,2,0.5))
  lines(mu ~ dec.date, data = obs.Rt, col = 1, lwd = 2)
  with(obs.Rt, polygon.fx(t = dec.date, LB = LB, UB = UB, col = col.fx(10, 0.1)))
  
  lines(mu ~ dec.date, data = fit.Rt, col = col.fx(10), lty = 2, lwd = 2)
  with(fit.Rt, polygon.fx(t = dec.date, LB = LB, UB = UB, col = col.fx(10, 0.1)))
  
  
  
  obs.Rt = get(sprintf("%s.obs.Rt.CI.sum", k))[1:11, ]
  fit.Rt = head(get(sprintf("%s.fit.Rt.CI.sum", k)), -1)
  with(fit.Rt, polygon.fx(t = dec.date, LB = LB, UB = UB, col = col.fx(10, 0.2)))
  
  
  
  
  
  abline(h = 1, col = alpha(1, 0.2), lty = 3)
  
  text(
    x      = 2020.32, 
    y      = -0.5, 
    labels = bquote(bold(.("Effective Reproduction Number"))), 
    srt    = 270, 
    xpd    = T
  )
  
  abline(v = decimal_date(as.Date("2020-01-22")), col = 2)
  
}

dev.off()





#
#--- END ---
#




