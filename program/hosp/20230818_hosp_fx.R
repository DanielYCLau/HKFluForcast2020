




library("splines")
library("pomp")

hosp.C.fx = function(ILIpS, hosp) {
  
  # scaling from ILIpS to hospitailization
  # ILIpS = data; hosp = hosp.data; yr.from = 2014; yr.to = 2019
  hosp.C = merge(
    x    = ILIpS[, c("WeekEnd", "Year", "Week", "ILIpS")],
    y    = hosp, # subset(hosp, weekend >= as.Date("2019-12-01") & weekend <= (as.Date("2019-12-31") + (week.no-1) * 7)),
    # y    = subset(hosp, weekend >= as.Date("2019-12-01") & weekend <= as.Date("2019-12-31")),
    by.x = c("WeekEnd", "Year", "Week"),
    by.y = c("weekend", "year", "week"),
    all.y = T
  )
  hosp.C = subset(hosp.C, Year %in% c(2014:2016, 2018:2019))
  
  hosp.C[, "Week"]    = with(hosp.C, ifelse(Week == 53, 52, Week))
  hosp.C[, "scaling"] = with(hosp.C, admission.1e6 / ILIpS)
  # hosp.C = aggregate(scaling ~ Week, data = hosp.C, FUN = function(x) { mean(x[is.finite(x)]) } )

  # tmp = aggregate(scaling ~ Week, data = hosp.C, FUN = function(x) { mean(x[is.finite(x)]) } )
  # plot(scaling ~ Week, data = hosp.C)
  # lines(scaling ~ Week, data = tmp, col = 2, lwd = 2)
  
  # pbsp.X    = with(hosp.C, bs(x = Week, knots = c(13, 26, 39)))
  # pbsp.fit  = lm(hosp.C$scaling ~ pbsp.X)
  # pbsp.pred = cbind(1, bs(x = 1:52, knots = c(13, 26, 39))) %*% coef(pbsp.fit)
  
  nb = 6
  pbsp.X    = with(hosp.C, periodic_bspline_basis(x = Week, nbasis = nb, period = 52))
  pbsp.fit  = lm(hosp.C$scaling ~ -1 + pbsp.X)
  pbsp.pred = periodic_bspline_basis(x = 1:52, nbasis = nb, period = 52) %*% coef(pbsp.fit)

  hosp.C = data.frame("Week" = 1:52, "scaling" = pbsp.pred)
  # lines(x = 1:52, y = pbsp.pred, col = 3, lwd = 2)
  
  return(hosp.C)
  
}





hosp.fx = function(hosp.data, f.list, hosp.C) {
  
  
  
  # forecast from scaling of ILIpS's forecast
  hosp.list = list()
  for (k in 1:4) {
    
    f.k = f.list[[k]][["result"]][["WIS"]][["forec"]][["pred.y.sum"]]
    f.k[, "Week"] = week(f.k[, "f.WeekEnd"])
    f.k = merge(x = f.k, y = hosp.C)
    t0 = min(f.k[, "f.WeekEnd"]) - 7

    
    
    current.hosp = subset(hosp.data, weekend == t0)
    current.hosp = with(
      current.hosp, data.frame(
        "f.weekend"  = weekend,
        "f.dec.date" = dec.date,
        "fit.mean"   = admission.1e6,
        "fit.LB"     = admission.1e6,
        "fit.UB"     = admission.1e6,
        "PI.mean"    = admission.1e6,
        "PI.LB"      = admission.1e6,
        "PI.UB"      = admission.1e6
      )
    )
    
    
    
    
    
    scale.t0 = merge(
      x = subset(data,      WeekEnd %in% (t0 - 7*0:3))[, c("WeekEnd", "ILIpS")],
      y = subset(hosp.data, weekend %in% (t0 - 7*0:3))[, c("weekend", "admission.1e6")],
      by.x = "WeekEnd",
      by.y = "weekend"
    )
    scale.t0 = with(scale.t0, admission.1e6 / ILIpS)
    f.k = transform(f.k, scaling = scaling / scaling[1] * mean(scale.t0))
    
    
    # f.adm.k = with(
    #   f.k, data.frame(
    #     "f.weekend"  = f.WeekEnd,
    #     "f.dec.date" = decimal_date(f.WeekEnd),
    #     "fit.mean"   = fit.mean   * scaling,
    #     "fit.LB"     = (fit.mean + qnorm(0.025) * fit.mean.se) * scaling,
    #     "fit.UB"     = (fit.mean + qnorm(0.975) * fit.mean.se) * scaling,
    #     "PI.mean"    = PI.mean    * scaling,
    #     "PI.LB"      = `PI.2.5%`  * scaling,
    #     "PI.UB"      = `PI.97.5%` * scaling
    #   )
    # )
    
    f.adm.k = with(
      f.k, data.frame(
        "f.weekend"  = f.WeekEnd,
        "f.dec.date" = decimal_date(f.WeekEnd),
        "fit.mean"   = fit.mean   * scaling,
        "fit.LB"     = (fit.mean + qnorm(0.025) * fit.mean.se) * scaling,
        "fit.UB"     = (fit.mean + qnorm(0.975) * fit.mean.se) * scaling,
        "PI.mean"    = PI.mean    * scaling,
        "PI.LB"      = `PI.2.5.`  * scaling,
        "PI.UB"      = `PI.97.5.` * scaling
      )
    )
    
    
    
    f.adm.k = rbind(current.hosp, f.adm.k)

    hosp.list[[k]] = f.adm.k
  }
  
  
  
  # output
  return(hosp.list)

}




f.hosp.list.fx = function() {
  
  load(list.files("program/hosp", sprintf("reg_sep_fore_%02dw_%d_forecast.Rdata", f.period, f.yr), full.names = T))
  
  f.hosp.list = lapply(
    1:4,
    function(k) {
      f.k = f.list[[k]][["result"]][["WIS"]][["forec"]][["pred.y.sum"]]
      f.k = with(
        f.k, data.frame(
          "f.weekend"  = f.WeekEnd,
          "f.dec.date" = decimal_date(f.WeekEnd),
          "fit.mean"   = fit.mean,
          "fit.LB"     = (fit.mean + qnorm(0.025) * fit.mean.se),
          "fit.UB"     = (fit.mean + qnorm(0.975) * fit.mean.se),
          "PI.mean"    = PI.mean,
          "PI.LB"      = `PI.2.5%`,
          "PI.UB"      = `PI.97.5%`
        )
      )
      
      current.hosp = subset(hosp.data, weekend == (min(f.k[, "f.weekend"]) - 7))
      current.hosp = with(
        current.hosp, data.frame(
          "f.weekend"  = weekend,
          "f.dec.date" = dec.date,
          "fit.mean"   = admission.1e6,
          "fit.LB"     = admission.1e6,
          "fit.UB"     = admission.1e6,
          "PI.mean"    = admission.1e6,
          "PI.LB"      = admission.1e6,
          "PI.UB"      = admission.1e6
        )
      )
      
      f.k = rbind(current.hosp, f.k)
      
    }
  )
  
  return(f.hosp.list)
}






hosp.compare.fx = function(hosp.C) {
  
  
  
  # forecast from scaling of ILIpS's forecast
  hosp.compare.list = list()
  for (k in 1:4) {
    
    load(list.files("program/hosp", sprintf("reg_sep_fore_%02dw_%d_forecast.Rdata", f.period, f.yr), full.names = T))
    # f.hosp.k = f.list[[k]][["result"]][["WIS"]][["forec"]][["pred.y.mc"]]
    f.hosp.k = f.list[[k]][["result"]][["WIS"]][["forec"]][["pred.y.sum"]]
    
    load(list.files("program/ILI", sprintf("reg_sep_fore_%02dw_%d_forecast.Rdata", f.period, f.yr), full.names = T))
    # f.sILI.k = f.list[[k]][["result"]][["WIS"]][["forec"]][["pred.y.mc"]]
    f.sILI.k = f.list[[k]][["result"]][["WIS"]][["forec"]][["pred.y.sum"]]
    
    tmp.k = f.list[[k]][["result"]][["WIS"]][["forec"]][["pred.y.sum"]]
    t0    = min(tmp.k[, "f.WeekEnd"]) - 7
    scaling = merge(
      x = transform(tmp.k, Week = week(tmp.k[, "f.WeekEnd"])),
      y = hosp.C
    )
    scaling = scaling[, "scaling"]
    
    scale.t0 = merge(
      x = subset(data,      WeekEnd %in% (t0 - 7*0:3))[, c("WeekEnd", "ILIpS")],
      y = subset(hosp.data, weekend %in% (t0 - 7*0:3))[, c("weekend", "admission.1e6")],
      by.x = "WeekEnd",
      by.y = "weekend"
    )
    scale.t0 = with(scale.t0, admission.1e6 / ILIpS)
    
    
    
    
    # peak.sILI.k = t(t(f.sILI.k) * scaling / scaling[1] * mean(scale.t0))
    # peak.sILI.k = peak.sILI.k[, which.max(apply(peak.sILI.k, 2, mean))] 
    # peak.hosp.k = f.hosp.k[, which.max(apply(f.hosp.k, 2, mean))] 
    # peak.k      = peak.sILI.k / peak.hosp.k 
# 
#     output.k = rbind(
#       data.frame(t0, "par" = "sILI",  "Q50" = median(peak.sILI.k), CI95.fx(peak.sILI.k) ),
#       data.frame(t0, "par" = "hosp",  "Q50" = median(peak.hosp.k), CI95.fx(peak.hosp.k) ),
#       data.frame(t0, "par" = "ratio", "Q50" = median(peak.k), CI95.fx(peak.k) )
#     )
    
    
    
    peak.sILI.k = transform(
      f.sILI.k, 
      fit.mean    = fit.mean    * scaling / scaling[1] * mean(scale.t0),
      fit.mean.se = fit.mean.se * scaling / scaling[1] * mean(scale.t0)
    )
    peak.sILI.k = with(
      subset(peak.sILI.k, fit.mean == max(fit.mean)),
      rnorm(1e6, fit.mean, fit.mean.se)
    )
    peak.hosp.k = with(
      subset(f.hosp.k, fit.mean == max(fit.mean)),
      rnorm(1e6, fit.mean, fit.mean.se)
    )
    peak.k      = peak.sILI.k / peak.hosp.k 
    
    output.k = rbind(
      data.frame(t0, "par" = "ratio", CI95.fx(peak.k) )
    )
    
    for (i in grep("t0|par", colnames(output.k), invert = T)) {
      output.k[,i] = sprintf("%.2f", output.k[,i])
    }
    
    output.k = rbind(output.k, NA)
    
    
    hosp.compare.list[[k]] = output.k
  }
  hosp.compare = do.call(rbind, hosp.compare.list)
  
  
  

  
  # output
  return(hosp.compare)
  
}
