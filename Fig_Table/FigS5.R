




#
#--- SETUP ---
#





load(list.files("data", "HK_flu_reg", full.names = T))
source("program/ILI/20240404_reg_sep_fore_Xw_fx.R")

substrRight = function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }




f.yr = 2020
f.period = 13
load(list.files("program/ILI", sprintf("reg_sep_fore_%02dw_%d_forecast.Rdata", f.period, f.yr), full.names = T))
f.list = f.list[1:4]
models = c("WIS", "RMSE", "LNQ", "MAE")





#
#--- PLOT SETTNG ---
#





t.lim = c(f.yr - 0.05, f.yr + 0.325)

data[, "t"] = decimal_date(data[, "WeekEnd"])
data = subset(data, t >= t.lim[1] & t <= t.lim[2]) 

axis.m = seq(as.Date(sprintf("%d-01-02", f.yr)), by = "month", length.out = 24)
axis.t = c(t.lim[1], decimal_date(axis.m), t.lim[2])

axis.t = axis.t[order(axis.t)]
axis.t = axis.t[axis.t <= t.lim[2]]

axis.t.label = head(month.abb[month(date_decimal(axis.t))], -1)
# axis.t.label = sprintf("%s %d", axis.t.label, head(year(date_decimal(axis.t)), -1) )




#
#--- PLOT ---
#





pdf(height = 6, width = 5, file = file.fx("FigS5_long_forec.pdf"))

# windows(width = 8, height = 8)
# quartz(width = 8, height = 8)
par(mfrow = c(4,1), mar = c(0,5,0,1), las = 1, oma = c(6,2,4,1))
for (k in 1:length(f.list)) {
  
  # k = 1
  # same for every plot
  plot(
    ILIpS/1e6 ~ t, data = data,
    ylim = c(0, max(data[, "ILIpS"]/1e6) * 6), xlim = t.lim,
    ylab = "", type = "l", xaxt = "n", yaxt = "n", ann = F
  )
  axis(2, at = seq(0, 0.03, 0.01))
  
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
  if (k == length(f.list)) {
    axis(1, at = axis.t[c(1,2,length(axis.t))], labels = NA, tck = -0.15)
    axis(1, at = head(axis.t, -1) + diff(axis.t)/2, labels = axis.t.label, tick = F)
    axis(1, at = quantile(head(axis.t, -1) + diff(axis.t)/2, c(0, 0.625)), labels = 2019:2020, line = 1.25, tick = F)
  }
  
  polygon.fx(t = c(f.yr-1, f.yr), LB = -1e6, UB = 1e6, col = alpha(1, 0.1))
  axis.w = decimal_date(seq(min(data[1, "WeekEnd"]), by = "week", length.out = 100))
  abline(v = axis.w, col = alpha(1, 0.2), lty = 3)
  abline(h = seq(0, 1, 0.01), col = alpha(1, 0.2), lty = 3)
  
  abline(v = decimal_date(as.Date("2020-01-22")), col = 2)
  # legend(
  #   "topright", 
  #   legend = "First COVID case\nin Hong Kong\non 22nd Jan 2020", 
  #   lty = 1, 
  #   col = 2,
  #   cex = 0.75,
  #   bty = "n"
  # )
  
  # start loop 
  for ( i in 1:length(models) ) {
    
    m.i = models[[i]]
    
    f.list.k = f.list[[k]]
    
    forec.y.ki = f.list.k[["result"]][[m.i]][["forec"]][["pred.y.sum"]]
    forec.y.ki[, "PI.LB"] = forec.y.ki[, "PI.2.5%"]
    forec.y.ki[, "PI.UB"] = forec.y.ki[, "PI.97.5%"]
    
    
    current.w = subset(data, WeekEnd == unique(forec.y.ki[, "WeekEnd"]) - 7)
    current.w = with(
      current.w, {
        data.frame(
          "WeekEnd"   = WeekEnd,
          "fore"      = 0,
          "f.WeekEnd" = WeekEnd,
          "f.obs"     = ILIpS,
          "fit.mean"  = ILIpS,
          "PI.LB"     = ILIpS,
          "PI.UB"     = ILIpS
        )
      }
    )
    
    forec.y.ki = bind_rows(current.w, forec.y.ki)
    forec.y.ki[, "t"] = decimal_date(forec.y.ki[, "f.WeekEnd"])
    for ( x in grep("fit|PI", colnames(forec.y.ki)) ) {
      forec.y.ki[, x] = forec.y.ki[, x] / 1e6
    }
    points(fit.mean ~ t, data = forec.y.ki, col = col.fx(i), pch = i)
    lines(fit.mean ~ t, data = forec.y.ki, col = col.fx(i), lty = i)
    with(forec.y.ki, polygon.fx(t, PI.LB, PI.UB, col = col.fx(i, 0.15)))
    points(fit.mean ~ t, data = head(forec.y.ki, 1), col = 1, lwd = 2, pch = 19)
    
    
    if (i == 1) {
      f.d.min = forec.y.ki[2, "f.WeekEnd"]
      f.d.max = max( forec.y.ki[, "f.WeekEnd"] )
      
      legend(
        "topleft", 
        legend = bquote(
          paste(
            .(week(f.d.min))^.(substrRight(ordinal(week(f.d.min)), 2)), "week of ", .(year(f.d.min)), " - ",
            .(week(f.d.max))^.(substrRight(ordinal(week(f.d.max)), 2)), "week of ", .(year(f.d.max))
          )
        ),
        bty = "n"
      )
      
    }
  }
  
}

legend(
  "bottom", 
  xpd    = NA, 
  inset  = c(0, -0.65), 
  legend = sprintf("%s       ", gsub("LNQ", "RMSLE", models)),
  lty    = 1:length(models), 
  pch    = 1:length(models), 
  col    = col.fx(1:length(models)),
  horiz  = T, 
  bty    = "n"
)

mtext("Influenza Activity (ILI+ proxy)", side = 2, line = 0, outer = T, las = 0)


dev.off()





#
#--- END ---
#




