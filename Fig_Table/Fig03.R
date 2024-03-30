




#
#--- SETUP ---
#





load(list.files("data", "HK_flu_reg", full.names = T))
load(list.files("data", "hosp_data", full.names = T))
source("program/hosp/20230628_reg_sep_fore_Xw_fx.R")





f.yr = 2020
models = c("WIS", "RMSE", "LNQ", "MAE")





data = merge(
  x     = data,
  y     = hosp.data,
  by.x  = c("WeekEnd", "Year", "Week"),
  by.y  = c("weekend", "year", "week"),
  all.y = T
)
data[, "ILIpS"] = round(data[, "admission.1e6"])





#
#--- PLOT SETTNG ---
#





t.lim = c(decimal_date(as.Date("2019-12-01")), f.yr + 0.325)

data[, "t"] = decimal_date(data[, "WeekEnd"])
data = subset(data, t >= decimal_date(as.Date("2019-12-01")-7) & t <= t.lim[2]) 

axis.m = seq(as.Date(sprintf("%d-01-02", f.yr)), by = "month", length.out = 24)
axis.t = c(t.lim[1], decimal_date(axis.m), t.lim[2])

axis.t = axis.t[order(axis.t)]
axis.t = axis.t[axis.t <= t.lim[2]]

axis.t.label = head(month.abb[month(date_decimal(axis.t))], -1)
# axis.t.label = sprintf("%s %d", axis.t.label, head(year(date_decimal(axis.t)), -1) )



event.fx = function(event, event.t, event.y, ...) {
  text(x = event.t, y = event.y, labels = event, pos = 4, ...)
  segments(x0 = event.t, y0 = -100, y1 = event.y)
  shape::Arrows(
    x0 = event.t, 
    x1 = event.t + 0.02, 
    y0 = event.y, 
    y1 = event.y, 
    arr.type = "triangle", 
    arr.length = 0.15,
    arr.width = 0.15
  )
}





#
#--- PLOT ---
#





pdf(
  height = 8, 
  width  = 6, 
  file   = file.fx("Fig3_short_long_forec.pdf")
)


# windows(width = 12, height = 16)
# quartz(width = 8, height = 8)
par(mfrow = c(2, 1), mar = c(4,6,1,1), las = 1, mgp = c(4,1,0))



for (f.period in c(4, 13)) {
  
  # f.period = 13
  load(list.files("program/hosp/", sprintf("reg_sep_fore_%02dw_%d_forecast.Rdata", f.period, f.yr), full.names = T))
  f.list = f.list[1:4]
  
  
  
  k = 3
  # same for every plot
  plot(
    ILIpS/1e4 ~ t, data = data,
    ylim = c(0, 0.03), xlim = t.lim,
    ylab = "", #, bquote("      Influenza-associated\nhospital admission ( per"~10^6~")"),
    xlab = "",
    type = "l", xaxt = "n", yaxt = "n" #, ann = F
  )
  axis(2, at = seq(0, 0.03, 0.005))
  
  if (k == 1) {
    title = sprintf("Forecast 1-%s weeks ahead in %s", f.period, f.yr)
    mtext(
      text = bquote(bold(.(title))), 
      side = 3, 
      line = 1.5, 
      outer = TRUE
    )
  }
  
  axis(1, at = axis.t, labels = NA)
  
  axis(1, at = axis.t[c(1,2,length(axis.t))], labels = NA, tck = -0.15)
  axis(1, at = head(axis.t, -1) + diff(axis.t)/2, labels = axis.t.label, tick = F)
  axis(1, at = quantile(head(axis.t, -1) + diff(axis.t)/2, c(0, 0.625)), labels = 2019:2020, line = 1.25, tick = F)
  
  
  polygon.fx(t = c(f.yr-1, f.yr), LB = -1e4, UB = 1e4, col = alpha(1, 0.1))
  axis.w = decimal_date(seq(min(data[1, "WeekEnd"]), by = "week", length.out = 100))
  abline(v = axis.w, col = alpha(1, 0.2), lty = 3)
  abline(h = seq(0, 0.03, 0.005), col = alpha(1, 0.2), lty = 3)
  
  
  
  # start loop 
  i = 1
  m.i = models[[i]]
  
  f.list.k = f.list[[k]]
  
  forec.y.ki = f.list.k[["result"]][[m.i]][["forec"]][["pred.y.sum"]]
  forec.y.ki[, "PI.LB"] = forec.y.ki[, "PI.2.5%"]
  forec.y.ki[, "PI.UB"] = forec.y.ki[, "PI.97.5%"]
  
  forec.y.ki[, "fit.mean.LB"] = with(forec.y.ki, fit.mean + qnorm(0.025) * fit.mean.se)
  forec.y.ki[, "fit.mean.UB"] = with(forec.y.ki, fit.mean + qnorm(0.975) * fit.mean.se)
  
  
  
  current.w = subset(data, WeekEnd == unique(forec.y.ki[, "WeekEnd"]) - 7)
  current.w = with(
    current.w, {
      data.frame(
        "WeekEnd"     = WeekEnd,
        "fore"        = 0,
        "f.WeekEnd"   = WeekEnd,
        "f.obs"       = ILIpS,
        "fit.mean"    = ILIpS,
        "fit.mean.LB" = ILIpS,
        "fit.mean.UB" = ILIpS,
        "PI.LB"       = ILIpS,
        "PI.UB"       = ILIpS
      )
    }
  )
  
  forec.y.ki = bind_rows(current.w, forec.y.ki)
  forec.y.ki[, "t"] = decimal_date(forec.y.ki[, "f.WeekEnd"])
  # for ( x in grep("fit|PI", colnames(forec.y.ki)) ) {
  #   forec.y.ki[, x] = forec.y.ki[, x] / 1e4
  # }
  points(fit.mean / 1e4 ~ t, data = forec.y.ki, col = col.fx(i), pch = i)
  lines(fit.mean / 1e4 ~ t, data = forec.y.ki, col = col.fx(i), lty = i)
  with(forec.y.ki, polygon.fx(t, PI.LB / 1e4, PI.UB / 1e4, col = col.fx(i, 0.15)))
  points(fit.mean / 1e4 ~ t, data = head(forec.y.ki, 1), col = 1, lwd = 2, pch = 19)
  
  with(forec.y.ki, polygon.fx(t, fit.mean.LB / 1e4, fit.mean.UB / 1e4, col = col.fx(i, 0.25)))
  
  
  if (f.period == 4) {
    mtext("A", side = 3, line = -1.5, adj = -0.275, font = 2, cex = 3)
  } else {
    mtext("B", side = 3, line = -1.5, adj = -0.275, font = 2, cex = 3)
  }
  
  
  
  event.fx(
    event   = "      First COVID-19 case in Wuhan\n      on 1st Dec 2019",
    event.t = decimal_date(as.Date("2019-12-01")),
    event.y = 0.028,
    cex     = 0.65
  )
  
  
  event.fx(
    event   = "      Screening at airports and train stations with\n      connects to Wuhan on 6th Jan 2020",
    event.t = decimal_date(as.Date("2020-01-06")),
    event.y = 0.024,
    cex     = 0.65
  )
  
  event.fx(
    event   = "      First COVID-19 case in Hong Kong\n      on 22nd Jan 2020",
    event.t = decimal_date(as.Date("2020-01-22")),
    event.y = 0.020,
    cex     = 0.65
  )
  
  mtext(bquote("Influenza-associated"), side = 2, line = 4.5, outer = F, las = 0)
  mtext(bquote("hospital admission rate"), side = 2, line = 3.5, outer = F, las = 0)
  
}

dev.off()





#
#--- END ---
#




