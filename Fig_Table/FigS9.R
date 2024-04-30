




#
#--- SETUP ---
#




os = Sys.info()[["sysname"]]
# if (os == "Linux") { setwd("DanielLau/regression/") }



load(list.files("data", "HK_flu_reg", full.names = T))
load(list.files("data", "hosp_data", full.names = T))
source("program/ILI/20240404_reg_sep_fore_Xw_fx.R")
source("program/hosp/20230818_hosp_fx.R")

f.yr     = 2020
f.period = 13
load(list.files("program/ILI", sprintf("reg_sep_fore_%02dw_%d_forecast.Rdata", f.period, f.yr), full.names = T))
# load("../../../Data/hospitalization/20230619_hosp_data.Rdata")





#
#--- DATA INPUT ---
#




hosp.C       = hosp.C.fx(data, hosp.data) # scaling constant by week no.
hosp.list    = hosp.fx(hosp.data, f.list, hosp.C) # first 4 forecast by week
f.hosp.list  = f.hosp.list.fx()
hosp.compare = hosp.compare.fx(hosp.C)




#
#--- PLOT ---
#





pdf(
  width = 5, height = 10,
  file = file.fx(sprintf("FigS9.pdf", f.period, f.yr))
)

# quartz(width = 6, height = 12)
# windows(width = 6, height = 12)
par(mfrow = c(4, 1), mar = c(0,3,0,1), oma = c(3,5,1,0), las = 1)

for (k in 1:4) {
  
  plot(
    NULL,
    xlim = c(0, 120),
    ylim = c(0, 300),
    main = "", # Admission rates in public hospitals with\nprincipal diagnosis of influenza per 1e6 population",
    xaxt = "n",
    yaxt = "n",
    ylab = ""
  )
  
  axis.date = seq(as.Date("2019-12-01"), as.Date("2020-04-01"), by = "month")
  axis.t    = axis.date - min(axis.date)
  
  if (k == 4) {
    axis(1, at = axis.t, labels = NA)
    axis(1, at = head(axis.t, -1) + 0.5 * diff(axis.t), labels = format(head(axis.date, -1), "%b %Y"), tick = F)
  } 

  
  mtext(
    "Influenza-associated hospital admission rate",
    side = 2, line = 2, outer = T, las = 0, font = 2
  )
  axis(2, at = seq(0, 300, 50), labels = sprintf("%.3f", seq(0, 300, 50) / 1e4))
  
  abline(v = axis.t, col = alpha(1, 0.1))
  
  
  # data
  lines(
    admission.1e6 ~ I(weekend - as.Date("2019-12-01") + 1),
    data = subset(hosp.data, weekend >= as.Date("2019-12-01")),
    col = 1, lwd = 3
  )
  
  
  
  # scaled ILI proxy
  hosp.k = hosp.list[[k]]
  hosp.k[, "t"] = with(hosp.k, f.weekend - as.Date("2019-12-01") + 1)
  lines( fit.mean ~ t, data = hosp.k, lwd = 3, col = col.fx(1) , lty = 3)
  # with( hosp.k, polygon.fx( t   = t, LB = fit.LB, UB = fit.UB, col = col.fx(1, 0.1) ) )
  # with( hosp.k, polygon.fx( t   = t, LB = PI.LB, UB = PI.UB, col = col.fx(1, 0.1) ) )
  with( hosp.k, polygon.fx( t   = t, LB = fit.LB, UB = fit.UB, col = col.fx(1, 0.1) ) )

  # abline(v = min(hosp.k[, "t"]), col = alpha(1, 0.5), lty = 3, lwd = 2)
  
  
  
  # forecasted hosp
  f.hosp.k = f.hosp.list[[k]]
  f.hosp.k[, "t"] = with(f.hosp.k, f.weekend - as.Date("2019-12-01") + 1)
  lines( fit.mean ~ t, data = f.hosp.k, lwd = 3, col = col.fx(7) , lty = 2)
  # with( f.hosp.k, polygon.fx( t   = t, LB = fit.LB, UB = fit.UB, col = col.fx(7, 0.1) ) )
  # with( f.hosp.k, polygon.fx( t   = t, LB = PI.LB, UB = PI.UB, col = col.fx(7, 0.2) ) )
  with( f.hosp.k, polygon.fx( t   = t, LB = fit.LB, UB = fit.UB, col = col.fx(7, 0.2) ) )
  
  points(fit.mean ~ t, data = hosp.k[1,], pch = 19, cex = 2)
  
  abline(v = as.Date("2020-01-22") - min(axis.date), col = 2)
  
}


dev.off()





write.csv(
  hosp.compare, 
  file = file.fx(sprintf("TableX_%02dw_hosp_compare.csv", f.period)), 
  na = "", 
  row.names = F
)





#
#--- END ---
#




