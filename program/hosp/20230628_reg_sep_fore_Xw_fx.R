


#
#--- PACKAGES ---
#



library("stringr")
library("zoo")
library("lubridate")
library("scales")
library("dplyr")
library("MASS")
library("pomp")
library("doSNOW")
library("parallel")
n.cores = detectCores() - 1

os = Sys.info()[['sysname']]





#
#--- FUNCTIONS ---
#



data.fx = function(data, f.yr, f.period, lag, start = 2010) {
  data = within(
    data, {
      
      ILIpS_fore_01 = ILIpS
      
      AH_lag_1      = lag(AH,        1 )
      ozone_lag_1   = lag(o3,        1 )
      ScHD_lag_1    = lag(ScHD,      1 )
      temp_lag_1    = lag(Mean_Temp, 1 )
      
    }
  )
  
  for (i in 1:lag) {
    data[, sprintf("ILIpS_lag_%02d", i)] = lag(data[, "ILIpS"], i)
  }
  
  for (x in grep("^ILIpS_lag", colnames(data)) ) {
    col.x = colnames(data)[x]
    data[, paste0("log_", col.x)] = log( ifelse( data[, x] == 0, 0.5, data[, x] ) )
  }

  # forecast
  for (i in 2:f.period) {
    data[, sprintf("ILIpS_fore_%02d", i)] = 
      with(data, c( tail( ILIpS, -(i-1) ), rep(NA, (i-1)) ) )
    if (i == f.period) {
      data[, "f.WeekEnd"] = data[, "WeekEnd"] + f.period * 7
    }
  }

  dec.t.lag.1 = lag( decimal_date(data[, "WeekEnd"]) )
  
  
  t.lag.1.M.bs  = periodic_bspline_basis(dec.t.lag.1, nbasis = 12, names = "t_lag_1_M_bs_%02d")
  
  M.bs.coln = colnames(t.lag.1.M.bs)
  data      = cbind(data, t.lag.1.M.bs)
  
  data    = subset(data, Year >= start)
  train.d = subset(data, Year >= start & year(f.WeekEnd) <  f.yr)
  forec.d = subset(data, year(f.WeekEnd) >= f.yr)

  na.id = apply(train.d, 1, function(x) { any(is.na(x)) } )
  if (any(na.id)) {
    forec.d = rbind(train.d[na.id, ], forec.d)
    train.d = train.d[!na.id, ]
  }

  
  
  for (x in colnames(data)[grep("temp|ScHD|ozone|AH", colnames(data))] ) {
    train.d[, paste0(x, "_sq")]  = train.d[, x]^2
    train.d[, paste0(x, "_Pow")] = log(train.d[, x])
    forec.d[, paste0(x, "_sq")]  = forec.d[, x]^2
    forec.d[, paste0(x, "_Pow")] = log(forec.d[, x])
  }

  assign("data",      data,      envir = .GlobalEnv)
  assign("train.d",   train.d,   envir = .GlobalEnv)
  assign("forec.d",   forec.d,   envir = .GlobalEnv)
  assign("M.bs.coln", M.bs.coln, envir = .GlobalEnv)
  
}



model.list.fx = function() {

  lag.n = length( grep("log_ILIpS_lag_[0-9][0-9]$", colnames(data)) )
  ILI.lag.coln = sprintf("log_ILIpS_lag_%02d", 1:lag.n)
  
  ILI.lag.coln = lapply(
    1:lag.n, function(i) {
      x.i = ILI.lag.coln[1:i]
      list(
        paste0(x.i, collapse = " + ")
      )
    } 
  )
  
  model.list = expand.grid(
    "lag"    = sprintf("%s", unlist(ILI.lag.coln)),
    "month"  = paste( c("-1", M.bs.coln), collapse = " + " ),
    "AH"     = c("", " + AH_lag_1"   , " + AH_lag_1    + AH_lag_1_sq"    , " + AH_lag_1_Pow"    ),
    "ozone"  = c("", " + ozone_lag_1", " + ozone_lag_1 + ozone_lag_1_sq" , " + ozone_lag_1_Pow" ),
    "temp"   = c("", " + temp_lag_1" , " + temp_lag_1  + temp_lag_1_sq"  , " + temp_lag_1_Pow"  ),
    "ScHD"   = c("", " + ScHD_lag_1" , " + ScHD_lag_1  + ScHD_lag_1_sq"  , " + ScHD_lag_1_Pow"  )
  )
  
  model.list = apply( model.list, 1, function(x) { paste0(x, collapse = "") } )
  model.list = gsub("    ", " ", model.list)
  model.list = gsub("   ",  " ", model.list)
  model.list = gsub("  ",   " ", model.list)
  
  assign( "model.list", model.list, envir = .GlobalEnv )
}





col.fx = function(ind, alpha = 1) {
  col = c(
    "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd", 
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
  )
  return(alpha(col[ind], alpha))
}

file.fx = function(filename) {
  sprintf( "%s_%s", format(Sys.Date(), "%Y%m%d"), filename )
}

polygon.fx = function(t, LB, UB, col, ...) {
  if (length(LB) == 1) { LB = rep(LB, length(t)) }
  if (length(UB) == 1) { UB = rep(UB, length(t)) }
  polygon(x = c(t, rev(t)), y = c(LB, rev(UB)), border = NA, col = col, ...)
}

CI95.fx = function(x, p = c(0.025, 0.5, 0.975)) {
  
  x.q = quantile(x, p)

  output = data.frame(
    "mean" = mean(x),
    "sd"   = sd(x)
  )
  output = cbind(output, t(x.q))
  
  return(output)
}





WIS.fx = function(y.sum, alpha) {
  
  K = length(alpha)
  
  y = y.sum[, "f.obs"]
  
  h.alpha = 0.5 * alpha
  alpha.Q = y.sum[, grep("[%]", colnames(y.sum))]
  
  m  = alpha.Q[, K+1]
  LB = alpha.Q[, 1:K]
  UB = alpha.Q[, ncol(alpha.Q):(K+2)]  
  
  
  WIS = sapply(
    1:nrow(y.sum),
    function(i) {
      # i=1
      
      y.i = y[i]
      m.i = m[i]

      tmp.i = data.frame(
        "h.alpha" = h.alpha,
        "y"       = y.i,
        "LB"      = as.numeric(LB[i,]),
        "UB"      = as.numeric(UB[i,])
      )
      
      wis.alpha.i = with(tmp.i, h.alpha*(UB - LB) + (LB - y) * (y < LB) + (y - UB) * (y > UB))
      
      WIS.i = 1/(K + 0.5) * ( 0.5 * abs(y.i - m.i) + sum( wis.alpha.i ) )
      return(WIS.i)
    }
  )
  
  return( WIS )
  
}





cost.fx = function(y.sum, alpha, prefix) {
  
  # y.sum = pred.train.k[["pred.y.sum"]]
  
  y.sum = do.call(rbind, y.sum)
  
  y.sum = within(
    y.sum, {
      e2     = (f.obs - fit.mean)^2
      abs.e  = abs(f.obs - fit.mean)
      e2.log = ( log(ifelse(f.obs == 0, 0.5, f.obs)) - log(fit.mean) )^2
    }
  )
  
  y.sum[, "wis"] = WIS.fx(y.sum, alpha)
  
  
  
  RMSE = aggregate(e2     ~ WeekEnd, data = y.sum, FUN = function(x) { sqrt(mean(x)) } )
  MAE  = aggregate(abs.e  ~ WeekEnd, data = y.sum, FUN = mean)
  LNQ  = aggregate(e2.log ~ WeekEnd, data = y.sum, FUN = function(x) { sqrt(mean(x)) } )
  WIS  = aggregate(wis    ~ WeekEnd, data = y.sum, FUN = mean)
  
  
  
  cost = data.frame(
    "RMSE" = RMSE[, 2], 
    "LNQ"  =  LNQ[, 2], 
    "MAE"  =  MAE[, 2],
    "WIS"  =  WIS[, 2]
  )
  
  colnames(cost) = sprintf("%s.%s", prefix, colnames(cost))
  
  return(cost)
}





nb.llh.fx = function(pars, y, X) {
  
  beta = head(pars, -1)
  disp = exp( tail(pars,  1) )
  
  mu = exp( beta %*% t(X) )

  llh = dnbinom(x = y, mu = mu, size = disp, log = T)
  llh = ifelse(is.infinite(llh), -1000, llh)
  
  sum.llh = sum(llh)
  return(-sum.llh)
  
}





pred.fx = function(model, data, pred.X, alpha, keep.mc) {
  
  
  
  #-- model fitting --
  
  #  model = model.i; data = train.k; pred.X = train.k
  #  model = model.i; data = train.k; pred.X = valid.k; keep.mc = F
  #  model = model.j; data = train.t; pred.X = forec.t
  pred.y.sum = list()
  pred.y.mc  = list()
  
  len   = length(grep("ILIpS_fore", colnames(data)))

  
  pred.y.list = lapply(
    1:len,
    function(t) {
      
      # for (t in 1:len) {
      
      #
      #--- fitting ---
      #
      
      # t=1
      model.t  = as.formula( sprintf("ILIpS_fore_%02d ~ %s", t, model) )
      
      
      fit.t     = glm(model.t, data = data, x = T, family="poisson", maxit = 1000)
      # fit.t.est = coef(fit.t)
      # fit.t.cov = vcov(fit.t)

      # fit.t    = suppressWarnings( glm.nb(model.t, data = data, x = T, maxit = 1000) )

      
      
      
      # fit.t.2 = optim(
      #   par     = c(coef(fit.t), log(0.1)), # log( fit.t[["theta"]] ) ),
      #   fn      = nb.llh.fx,
      #   method  = "L-BFGS-B",
      #   hessian = T,
      #   y       = fit.t[["y"]],
      #   X       = fit.t[["x"]],
      #   upper   = c(rep(Inf, length(coef(fit.t))), log(1e8)), 
      #   control = list("ndeps" = rep(1e-6, length(coef(fit.t))+1 ) )
      # )
      # 
      # fit.t.est = fit.t.2[["par"]]
      # fit.t.cov = solve(fit.t.2[["hessian"]])
      # 
      # if (any(eigen(fit.t.cov)[["values"]] < 0)) {
      #   fit.t.cov = as.matrix(Matrix::nearPD(fit.t.cov)[["mat"]])
      # }

      
      
      fit.t.2 = try(
        optim(
          par     = c(coef(fit.t), log(0.1)), # log( fit.t[["theta"]] ) ),
          fn      = nb.llh.fx,
          method  = "L-BFGS-B",
          hessian = T,
          y       = fit.t[["y"]],
          X       = fit.t[["x"]],
          upper   = c(rep(Inf, length(coef(fit.t))), log(1e9)), 
          control = list("ndeps" = rep(1e-6, length(coef(fit.t))+1 ) )
        ), silent = T 
      )
      
      if (is.character(fit.t.2)) {
        fit.t.est = coef(fit.t)
        fit.t.cov = vcov(fit.t)
      } else {
        fit.t.est = fit.t.2[["par"]]
        fit.t.cov = solve(fit.t.2[["hessian"]])
        
        if (any(eigen(fit.t.cov)[["values"]] < 0)) {
          fit.t.cov = as.matrix(Matrix::nearPD(fit.t.cov)[["mat"]])
        }
      }
      

      
      #
      #--- generate prediction sample ---
      #
      
      set.seed(555 * t)
      N          = 10000
      est.t.mc   = mvrnorm(N, mu = fit.t.est, Sigma = fit.t.cov)
      
      if (is.character(fit.t.2)) {
        beta.t.mc  = est.t.mc
      } else {
        beta.t.mc  = est.t.mc[, -ncol(est.t.mc)]
        theta.t.mc = exp( est.t.mc[,  ncol(est.t.mc)] )
      }
        
      
      if ( any( grepl("Intercept", colnames(beta.t.mc)) ) ) {
        pred.t.X = cbind(1, pred.X[, tail(colnames(beta.t.mc), -1)])
      } else {
        pred.t.X = pred.X[, colnames(beta.t.mc)]
      }
      mu.t.mc = exp( beta.t.mc %*% t(pred.t.X) )
      
      
      
      if ( nrow(pred.X) == 1 ) {
        
        mu.t    = mean(mu.t.mc)
        mu.t.se = sd(mu.t.mc)
        
        
        set.seed(555 * t)
        if (is.character(fit.t.2)) {
          pred.y.t.mc = rpois(n = N, lambda = mu.t.mc)
        } else {
          pred.y.t.mc = rnbinom(n = N, mu = mu.t.mc, size = theta.t.mc)
        }

        pred.y.t.sum = CI95.fx(pred.y.t.mc, p = c(0.5 * alpha, 0.5, 1-rev(0.5*alpha)) )
        colnames(pred.y.t.sum) = sprintf("PI.%s", colnames(pred.y.t.sum))
        
      } else {
        
        mu.t.mc = lapply(1:ncol(mu.t.mc), function(t) { mu.t.mc[, t] } )
        mu.t    = sapply(mu.t.mc, mean)
        mu.t.se = sapply(mu.t.mc, sd)
        
        # start.t = Sys.time()
        
        pred.y.t.mc = lapply(
          1:length(mu.t.mc),
          function(t) {
            
            set.seed(555 * t)
            if (is.character(fit.t.2)) {
              pred.y.t = rpois(n = N, lambda = mu.t.mc[[t]])
            } else {
              pred.y.t = rnbinom(n = N, mu = mu.t.mc[[t]], size = theta.t.mc)
            }
            return( pred.y.t )
          }
        )
        
        # end.t = Sys.time()
        # end.t - start.t
        
        
        
        pred.y.t.sum = lapply(pred.y.t.mc, CI95.fx, p = c(0.5 * alpha, 0.5, 1-rev(0.5*alpha)) )
        pred.y.t.sum = do.call(rbind, pred.y.t.sum)
        colnames(pred.y.t.sum) = sprintf("PI.%s", colnames(pred.y.t.sum))

      }
      
      pred.y.t.sum = cbind(
        "WeekEnd"     = pred.X[, "WeekEnd"], 
        "fore"        = t,
        "f.WeekEnd"   = pred.X[, "WeekEnd"] + 7*(t-1), 
        "f.obs"       = pred.X[, sprintf("ILIpS_fore_%02d", t)],
        "fit.mean"    = mu.t,
        "fit.mean.se" = mu.t.se, 
        pred.y.t.sum
      )
      rownames(pred.y.t.sum) = NULL
      
      
      
      # pred.y.sum[[t]] = pred.y.t.sum
      # pred.y.mc[[t]]  = pred.y.t.mc
      
      
      if (keep.mc) {
        return( list(pred.y.t.sum, pred.y.t.mc) )
      } else {
        return( list(pred.y.t.sum) )
      }
      
    }
  )
  
  
  if (keep.mc) {
    pred.y.sum = lapply(pred.y.list, function(x) { x[[1]] } )
    pred.y.mc  = lapply(pred.y.list, function(x) { x[[2]] } )
    return( list("pred.y.sum" = pred.y.sum, "pred.y.mc" = pred.y.mc) )
  } else {
    pred.y.sum = lapply(pred.y.list, function(x) { x[[1]] } )
    return( list("pred.y.sum" = pred.y.sum) )
  }

}






fore.plot.fx = function(forec.sum, title, xlim.yr) {
  
  # preparation
  forec.sum[, "PI.Q50"] = forec.sum[, "PI.50%"]
  forec.sum[, "PI.LB"]  = forec.sum[, "PI.2.5%"]
  forec.sum[, "PI.UB"]  = forec.sum[, "PI.97.5%"]
  forec.sum = forec.sum[, !grepl("[%]", colnames(forec.sum))]
  
  
  
  current.w = subset(data, WeekEnd == unique(forec.sum[, "WeekEnd"]) - 7)
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
  
  forec.sum = bind_rows(current.w, forec.sum)
  
  # plot
  plot(
    ILIpS ~ decimal_date(WeekEnd), 
    data = data, 
    type = "l", 
    xlim = xlim.yr,
    main = title
  )
  lines(fit.mean ~ decimal_date(f.WeekEnd), data = forec.sum, col = "orange")
  points(fit.mean ~ decimal_date(f.WeekEnd), data = forec.sum, col = "orange", pch = 16)
  points(f.obs ~ decimal_date(f.WeekEnd), data = current.w, pch = 16)
  
  
  with(
    forec.sum, 
    polygon.fx(
      t   = decimal_date(f.WeekEnd),
      LB  = PI.LB,
      UB  = PI.UB,
      col = alpha("orange", 0.1) 
    )
  )
  
  legend(
    "topleft",
    legend = c("When forecast is conducted", "Forecast mean"),
    pch    = 16,
    col    = c(1, "orange"),
    bty    = "n"  
  )
  
}





criteria.rank.fx = function(cv, rank.w.value) {
  # cv = cv.table
  # cv = cv[, grep("LNQ", colnames(cv), invert = T)]
  
  criteria.col.1 = grep("v[.]", colnames(cv))
  
  rank   = sapply(criteria.col.1, function(i) { rank(cv[, i]) } )
  rank.c = lapply(1:nrow(cv), function(i) { sprintf("%.f (%.01f)", rank[i,], cv[i, criteria.col.1]) } )
  rank.c = do.call(rbind, rank.c)
  
  # cv[, criteria.col.1]    = rank
  # cv[, "v.median.rank"] = apply(rank, 1, median)
  # cv[, "v.mean.rank"]   = round(apply(rank, 1, mean), 1)
  
  criteria.col.2 = grep("v[.]", colnames(cv))
  
  min.cv.ind = sapply(criteria.col.2, function(i) { which.min(cv[, i]) } )
  min.cv.ind = unique(min.cv.ind)
  cv   = cv[min.cv.ind, ]
  if (rank.w.value) {
    cv[, criteria.col.1] = rank.c[min.cv.ind, ]
    
  } 
  
  return(cv)
}





model.abb.fx = function(models, x) {
  
  # x = "ozone"
  covariate = gsub("\\s+", "", models)
  covariate = strsplit(covariate, "[+]|[-]1")
  tmp       = lapply(covariate, function(z) { z[grep(x, z)] } )
  
  if (x == "ILI") {
    tmp.sq  = sapply(tmp, function(z) { ifelse(length(grep("sq", z)) > 0, "sq", "") } )
    tmp.lag = sapply(tmp, function(z) { length(grep("[0-9][0-9]$", z)) } )
    tmp     = sprintf("%s(%s)(%s)", x, tmp.lag, tmp.sq)
  } else {
    tmp = lapply(tmp, function(z) { z[grep(x, z)] } )
    
    tmp1 = sapply(tmp, length)
    tmp1 = ifelse(tmp1 == 0, "", ifelse(tmp1 == 2, "sq", tmp1))
    tmp2 = sapply(tmp, function(x) { any(grepl("Pow", x)) } )
    
    tmp = sprintf("%s(%s)", x, ifelse(tmp2, "Pow", tmp1))
  }
  return(tmp)
}

model.abb.v2.fx = function(models) {
  
  # models = cv.table[, "model"]
  model.abb = lapply(
    c("ILI", "bs", "AH", "ozone", "temp", "ScHD"), 
    function(x) { model.abb.fx(models, x) }
  )
  
  model.abb = do.call(cbind, model.abb)
  model.abb = t(apply(model.abb, 1, function(x) { gsub("[(][])]", "", x) } ))
  model.abb = t(apply(model.abb, 1, function(x) { ifelse(grepl("[(]", x), x, "") } ))
  model.abb = apply(model.abb, 1, function(x) { paste(x, collapse = " + ") } )
  model.abb = gsub("[+]\\s+[+]",   "+", model.abb)
  model.abb = gsub("\\s+[+]\\s+$", "",  model.abb)
  return(model.abb)  
}




# hosp.C.fx = function(ILIpS, hosp, yr.range) {
#   
#   # scaling from ILIpS to hospitailization
#   # ILIpS = data; hosp = hosp.data; yr.from = 2014; yr.to = 2019
#   hosp.C = merge(
#     x    = ILIpS[, c("WeekEnd", "Year", "Week", "ILIpS")],
#     y    = subset(hosp, year %in% yr.range),
#     # y    = subset(hosp, weekend >= as.Date("2019-12-01") & weekend <= as.Date("2019-12-31")),
#     by.x = c("WeekEnd", "Year", "Week"),
#     by.y = c("weekend", "year", "week")
#   )
#   
#   # hosp.C = with(hosp.C, admission.1e6 / ILIpS)
#   # hosp.C = mean(hosp.C)
#   
#   hosp.C[, "hosp.C"] = with(hosp.C, admission.1e6 / ILIpS)
#   hosp.C = aggregate( hosp.C ~ Week, data = hosp.C, FUN = mean )
#   
#   # plot(hosp.C ~ dec.date, data = hosp.C, type = "l")
#   # plot(hosp.C ~ Week, data = hosp.C, type = "l")
#   return(hosp.C)
#   
# }



hosp.C.fx = function(ILIpS, hosp, week.no) {
  
  # scaling from ILIpS to hospitailization
  # ILIpS = data; hosp = hosp.data; yr.from = 2014; yr.to = 2019
  hosp.C = merge(
    x    = ILIpS[, c("WeekEnd", "Year", "Week", "ILIpS")],
    y    = subset(hosp, weekend >= as.Date("2019-12-01") & weekend <= (as.Date("2019-12-31") + (week.no-1) * 7)),
    # y    = subset(hosp, weekend >= as.Date("2019-12-01") & weekend <= as.Date("2019-12-31")),
    by.x = c("WeekEnd", "Year", "Week"),
    by.y = c("weekend", "year", "week")
  )
  
  hosp.C = with(hosp.C, admission.1e6 / ILIpS)
  hosp.C = mean(hosp.C)
  return(hosp.C)
  
}




hosp.fx = function(hosp.data, f.list, hosp.C) {
  
  
  
  # forecast based on week average over previous years
  hosp.week.mean = aggregate(
    admission.1e6 ~ week,
    data = subset(hosp.data, year %in% c(2014:2016, 2018:2019)),
    FUN  = mean
  )
  colnames(hosp.week.mean) = c("week", "adm.week.mean")
  
  hosp.week.mean = suppressWarnings(
    within(
      hosp.week.mean, {
        weekend       = as.Date(sprintf("2020/%d/6", week), "%Y/%W/%w") - 7
        dec.date      = decimal_date(weekend)
      }
    )
  )

  hosp.week.mean = merge(
    x  = subset(hosp.data, weekend >= as.Date("2019-12-01") & weekend <= as.Date("2020-04-01")),
    y  = hosp.week.mean,
    by = c("weekend", "dec.date", "week"),
    all.x = T
  )
  
  
  
  # forecast from scaling of ILIpS's forecast
  hosp.list = list()
  for (k in 1:4) {
    
    f.adm.k = f.list[[k]][["result"]][["WIS"]][["forec"]][["pred.y.sum"]]
    f.adm.k = with(
      f.adm.k, data.frame(
        "f.weekend"  = f.WeekEnd,
        "f.dec.date" = decimal_date(f.WeekEnd),
        "fit.mean"   = fit.mean   * hosp.C,
        "fit.LB"     = (fit.mean + qnorm(0.025) * fit.mean.se) * hosp.C,
        "fit.UB"     = (fit.mean + qnorm(0.975) * fit.mean.se) * hosp.C,
        "PI.mean"    = PI.mean    * hosp.C,
        "PI.LB"      = `PI.2.5%`  * hosp.C,
        "PI.UB"      = `PI.97.5%` * hosp.C
      )
    )
    
    
    current.hosp = subset(hosp.data, weekend == (min(f.adm.k[, "f.weekend"]) - 7))
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
    
    
    f.adm.k = rbind(current.hosp, f.adm.k)

    
    
    # f.adm.k = with(
    #   f.adm.k, data.frame(
    #     "f.weekend"  = f.WeekEnd,
    #     "f.dec.date" = decimal_date(f.WeekEnd),
    #     "fit.mean"   = fit.mean   ,
    #     "fit.LB"     = (fit.mean + qnorm(0.025) * fit.mean.se) ,
    #     "fit.UB"     = (fit.mean + qnorm(0.975) * fit.mean.se) ,
    #     "PI.mean"    = PI.mean    ,
    #     "PI.LB"      = `PI.2.5%`  ,
    #     "PI.UB"      = `PI.97.5%` 
    #   )
    # )
    
    f.adm.k = merge(
      x    = f.adm.k,
      y    = hosp.week.mean,
      by.x = c("f.weekend", "f.dec.date"),
      by.y = c("weekend", "dec.date")
    )
    
    # f.adm.k = merge(
    #   x = f.adm.k,
    #   y = hosp.C,
    #   by.x = "week",
    #   by.y = "Week",
    #   all.x = T
    # )
    # 
    # 
    # f.adm.k = transform(
    #   f.adm.k, 
    #   fit.mean = fit.mean * hosp.C,
    #   fit.LB   = fit.LB   * hosp.C,
    #   fit.UB   = fit.UB   * hosp.C,
    #   PI.mean  = PI.mean  * hosp.C,
    #   PI.LB    = PI.LB    * hosp.C,
    #   PI.UB    = PI.UB    * hosp.C
    # )

    hosp.list[[k]] = f.adm.k
  }
  
  
  
  # output
  assign("hosp.week.mean", hosp.week.mean, envir = .GlobalEnv)
  assign("hosp.list", hosp.list, envir = .GlobalEnv)
  
}

#
#--- END ---
#




