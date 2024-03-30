




#
#--- SETUP ---
#





load(list.files("data", "HK_flu_reg", full.names = T))
load(list.files("data", "hosp_data", full.names = T))
source("program/20230628_reg_sep_fore_Xw_fx.R")





f.yr = 2020
f.period = 4
data.fx(data, f.yr, f.period, 10)
data[, "t"] = decimal_date(data[, "WeekEnd"])


data = merge(
  x = data,
  y = hosp.data,
  by.x  = c("WeekEnd", "Year", "Week"),
  by.y  = c("weekend", "year", "week"),
  all.x = T
)






# plot function
plot.fx = function(x, y.lab, panel, n = 1) {

  with(
    data, {
      plot(
        NULL,
        ylim = range(get(x), na.rm = T)/n,
        xlim = range(t), 
        type = "l",
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = y.lab
      )
    }
  )
  
  axis(2)
  axis(1, at = 2010:2030, labels = F)
  axis(1, at = 2010:2030 + 0.5, labels = 2010:2030, tick = F)
  abline(v = 2010:2030, col = alpha(1, 0.1))
  
  with(data, lines(x = t, y = get(x)/n))
  
  mtext(panel, side = 3, adj = -0.1, cex = 2, font = 2)
  
}



# school closure
ScClos.dt = as.Date(
  c(
    "15/03/2008", "22/03/2008", "29/03/2008", "13/06/2009", "20/06/2009", "27/06/2009", 
    "04/07/2009", "10/02/2018",  "17/02/2018", "26/01/2019", "02/02/2019","22/01/2020", "23/01/2020"
  ),
  format = "%d/%m/%Y"
)





# output

pdf(height = 10, width = 12, file = file.fx("Fig1.pdf"))

# windows(height = 10, width = 12)
# quartz(height = 10, width = 12)
par(mfrow = c(5,1), mar = c(3,6,2,1), mgp = c(4, 1, 0), oma = c(0, 5, 1, 0), las = 1)

plot.fx("ILIpS", "ILI+ proxy\n(influenza activity)", "A", 1e6)
with(data, polygon.fx(t  = t, UB = ifelse(ScHD > 2, 1, -1), LB = -1, col = alpha("lightblue", 0.5) ) )
abline(v = decimal_date(ScClos.dt), col = alpha(2, 0.5), lty = 1 , lwd = 2)
with(data, lines(x = t, y = ILIpS/1e6))

plot.fx("admission.1e6", "Influenza-associated\nhospital admission rate", "B", 1e4)
with(data, polygon.fx(t  = t, UB = ifelse(ScHD > 2, 1, -1), LB = -1, col = alpha("lightblue", 0.5) ) )
abline(v = decimal_date(ScClos.dt), col = alpha(2, 0.5), lty = 1 , lwd = 2)
with(data, lines(x = t, y = admission.1e6/1e4))

plot.fx("Mean_Temp", bquote("Mean temperature (" *degree*"C)"), "C")
plot.fx(       "AH", bquote("Absolute humidiy (mg/m"^"3"*")"),  "D")
plot.fx(       "o3", bquote("Ozone ("*mu*"g/m"^"3"*")"),        "E")

dev.off()





#
#--- END ---
#




