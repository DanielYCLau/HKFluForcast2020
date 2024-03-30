




#
#--- SETUP ---
#




# also for S6

library("lubridate")
load(list.files("data", "HK_flu_reg", full.names = T))
load(list.files("data", "hosp_data", full.names = T))

data = merge(
  x     = data,
  y     = hosp.data,
  by.x  = c("WeekEnd", "Year", "Week"),
  by.y  = c("weekend", "year", "week"),
  all.y = T
)
data[, "ILIpS"] = round(data[, "admission.1e6"])

ILIpS  = data[, c("WeekEnd", "ILIpS")]
season = with(ILIpS, ifelse(month(WeekEnd) >= 10, year(WeekEnd), year(WeekEnd) - 1))
ILIpS  = split(ILIpS, season)
ILIpS  = ILIpS[["2019"]]
ILIpS[, "f.WeekEnd"] = ILIpS[, "WeekEnd"]

source("program/hosp/20230628_reg_sep_fore_Xw_fx.R")





#
#--- DATA INPUT ---
#





f.yr = 2020
f.period = 13
data.fx(data, f.yr, f.period, 10)
model.list.fx()





#
#--- RESULT ---
#





load(list.files("program/hosp", sprintf("reg_sep_fore_%02dw_%d_forecast.Rdata", f.period, f.yr), full.names = T))


att.fm = list()

for (k in 1:4) {
  
  # k = 1
  stdt.k   = as.Date(f.list[[k]][["start.d"]])
  f.list.k = f.list[[k]][["result"]] # with highest incidence
  
  
  att.fm.k = list()
  
  for (i in c("WIS", "RMSE", "LNQ", "MAE")) { 
    
    # i = "WIS"
    forec.y  = f.list.k[[i]][["forec"]][["pred.y.sum"]]
    forec.mc = f.list.k[[i]][["forec"]][["pred.y.mc"]]
    
    pred.d = subset(data, WeekEnd < unique(forec.y[, "WeekEnd"]))
    expected.y = with(pred.d, data.frame("f.WeekEnd" = WeekEnd, "fit.mean" = ILIpS, "fit.mean.se" = 0))
    expected.y = plyr::rbind.fill(expected.y, forec.y)
    expected.y[, "PI.LB"] = expected.y[, "PI.2.5%"]
    expected.y[, "PI.UB"] = expected.y[, "PI.97.5%"]
    
    # expected.y = expected.y[colSums(is.na(expected.y)) == 0]
    
    
    # define season
    season     = with(expected.y, ifelse(month(f.WeekEnd) >= 10, year(f.WeekEnd), year(f.WeekEnd) - 1))
    expected.y = split(expected.y, season)
    
    
    
    # calculate attack rate
    a.r = sapply(
      expected.y[as.character(2013:2018)], 
      function(x) { 
        x = subset(x, month(f.WeekEnd) %in% c(12, 1:3))
        y = sum(x[, "fit.mean"]) # / 1e6
        return(y)
      }
    )
    median(a.r); range(a.r)
    
    a.r.1920 = expected.y[["2019"]]
    
    d.tmp = subset(a.r.1920, fit.mean.se == 0)
    f.tmp = subset(a.r.1920, fit.mean.se  > 0)
    
    cat(sprintf("\n[ --- Model: %s --- ]:", i))
    
    period = c("winter")
    for (p in period) {
      
      # p = "winter"
      subset.p = switch(
        p,
        "summer" = 'f.WeekEnd >= as.Date("2020-05-01") & f.WeekEnd < as.Date("2020-10-01")',
        "winter" = 'f.WeekEnd >= as.Date("2019-12-01") & f.WeekEnd < as.Date("2020-04-01")',
        "all"    = 'f.WeekEnd >= as.Date("2019-10-01") & f.WeekEnd < as.Date("2020-10-01")'
      )
      
      f.ind = with(f.tmp, eval(parse(text = subset.p)))
      
      f.Inc   = f.tmp[f.ind, ]
      f.Inc.t = f.Inc[which.max(f.Inc[, "fit.mean"]), "f.WeekEnd"]
      f.Inc   = f.Inc[which.max(f.Inc[, "fit.mean"]), c("fit.mean", "PI.LB", "PI.UB")]
      f.Inc   = f.Inc / 1e4 * 100

      
      
      f.Inc.pk = f.tmp[f.ind, ]
      f.Inc.pk = f.Inc.pk[which.max(f.Inc.pk[, "fit.mean"]), c("fit.mean", "fit.mean.se")]

      f.Inc.pk = with(
        f.Inc.pk, c(
          fit.mean,
          fit.mean + qnorm(0.025) * fit.mean.se,
          fit.mean + qnorm(0.975) * fit.mean.se
        )
      )
      f.Inc.pk = f.Inc.pk / 1e4 * 100

      obs.peak = subset(data, WeekEnd >= as.Date("2019-12-01") & WeekEnd <= as.Date("2020-04-01"))
      obs.peak = max(obs.peak[, "ILIpS"] / 1e4 * 100)
      peak.reduce = -(obs.peak / f.Inc.pk - 1) * 100

      
      f.att = f.tmp[f.ind, ]
      f.att = sum(f.att[, "fit.mean"])
      
      f.mc  = forec.mc[, c(f.ind, rep(F, ncol(forec.mc) - length(f.ind)))]
      f.mc  = rowSums(f.mc)
      f.mc  = quantile(f.mc, c(0.025, 0.975))
      
      f.se  = with(f.tmp[f.ind, ], sqrt(sum((fit.mean.se)^2)))
      f.CI  = f.att + qnorm(c(0.025, 0.975)) * f.se
      
      
      
      if (p %in% c("winter", "all")) {
        
        d.ind = with(d.tmp, eval(parse(text = subset.p)))
        d.att = d.tmp[d.ind, ]
        d.att = sum(d.att[, "fit.mean"])
        
        f.att.PI = (d.att + c(f.att, f.mc)) / 1e4 * 100
        f.att.CI = (d.att + c(f.att, f.CI)) / 1e4 * 100
      } else {
        f.att.PI = c(f.att, f.mc) / 1e4 * 100
        f.att.CI = c(f.att, f.CI) / 1e4 * 100
      }
      
      # reduction in attack rate
      ILIpS.p = subset(ILIpS, eval(parse(text = subset.p)))
      d.att = sum(ILIpS.p[, "ILIpS"]) / 1e4 * 100
      f.att.reduce.PI = (1 - d.att/f.att.PI) * 100
      f.att.reduce.CI = (1 - d.att/f.att.CI) * 100
      
      
      # summarise
      att.i.p = data.frame(
        "stdt"  = stdt.k,
        "model" = i,
        "p"     = p,
        "Inc.t" = paste(format(c(f.Inc.t - 6, f.Inc.t), "%d%B"), collapse = "-"),
        "Inc"   = sprintf("%.1f (%.1f - %.1f)", f.Inc[1], f.Inc[2], f.Inc[3]),
        "R.pk"  = sprintf("%.1f (%.1f - %.1f)", peak.reduce[1], peak.reduce[2], peak.reduce[3]),
        "Att"   = sprintf("%.1f (%.1f - %.1f)", f.att.PI[1], f.att.PI[2], f.att.PI[3]),
        # "R.Att" = sprintf("%.1f (%.1f - %.1f)", f.att.reduce.PI[1], f.att.reduce.PI[2], f.att.reduce.PI[3]),
        "R.Att" = sprintf("%.1f (%.1f - %.1f)", f.att.reduce.CI[1], f.att.reduce.CI[2], f.att.reduce.CI[3]),
        "NA"    = NA
      )
      
      att.fm.k[[i]][[p]] = t(att.i.p)
    }
    
  }
  
  att.fm.k = lapply(att.fm.k, function(x) { do.call(rbind, x) } )
  att.fm.k = do.call(rbind, att.fm.k)
  att.fm[[k]] = att.fm.k
  
}

att.fm = do.call(cbind, att.fm)



write.csv(
  att.fm, 
  file = file.fx(sprintf("Table2_%02dw_att.csv", f.period)),
  na = ""
)





#
#--- END ---
#




