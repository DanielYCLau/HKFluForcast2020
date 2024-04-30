




#
#--- ADDITIONAL FUNCTIONS ---
#





library("EpiEstim")
library("pracma")




n.cores = parallel::detectCores()-1





pb.fx = function(iter) {
  pb = txtProgressBar(max = iter, style = 3, width = 50)
  progress = function(i) { setTxtProgressBar(pb, i) }
  opts = list(progress = progress)
  return(opts)
}





daily.incid.fx = function(data, x) {
  
  data = rbind(
    data.frame(
      "WeekEnd"  = min(data[, "WeekEnd"]) - 7,
      "ILIpS"    = 0,
      "fit.mean" = 0
    ),
    data
  )
  
  mono.sp.fx = function(t, y, pred.t) {
    # mono.sp = splinefun(x = t, y = cumsum(y), method = "hyman")
    # pred.y  = mono.sp(pred.t)
    pred.y = pchip(t, cumsum(y), pred.t)
    return(pred.y)
  }

  pred.dec.date = with(data, seq(min(WeekEnd), max(WeekEnd), by = "day"))

  pred.cumsum = with(
    data, mono.sp.fx(t = decimal_date(WeekEnd), y = get(x), pred.t = decimal_date(pred.dec.date))
  )
  # plot(cumsum(ILIpS) ~ decimal_date(WeekEnd), data = data, type = "l")
  # lines(pred.cumsum ~ pred.dec.date, col = 2, lwd = 2)
  
  output = data.frame(
    "date"  = pred.dec.date[-1],
    "incid" = diff(pred.cumsum)
  )
  colnames(output)[2] = x
  
  return(output)
  
}





log.mu.sd.fx = function(mu, sigma) {
  log.mu = log(mu^2 / sqrt(mu^2 + sigma^2))
  log.sd = sqrt(log(1 + (sigma/mu)^2))
  return(c("log.mu" = log.mu, "log.sd" = log.sd))
}





normal.mu.sd.fx = function(log.mu, log.sd) {
  log.v = log.sd^2
  mu    = exp(log.mu + 0.5 * log.v)
  sd    = sqrt( (exp(log.v)-1)*exp(2*log.mu + log.v) )
  return(c("mu" = mu, "sd" = sd))
}




w.Rt.fx = function(data, x, method, Rt.config) {
  
  # data = daily.inc; x = "ILIpS"; method = "Cori"
  
  if (method == "Cori") {
    d.Rt = estimate_R(
      incid  = data[, x],
      method = "parametric_si",
      config = make_config(Rt.config)
    )
  } else if (method == "WT") {
    d.Rt = wallinga_teunis(
      incid  = round(data[, x]),
      method = "parametric_si",
      config = c(Rt.config, list(n_sim = 100))
    )
  }
  d.Rt = d.Rt[["R"]]
  
  
  # take geometric mean of daily Rt as weekly Rt
  # mean = exp( mean( log(x) ) )
  
  log.Rt = sapply(
    1:nrow(d.Rt),
    function(i) {
      with(d.Rt[i, ], log.mu.sd.fx(mu = `Mean(R)`, sigma = `Std(R)`) ) 
    }
  )
  
  d.Rt = cbind(d.Rt, t(log.Rt))
  d.Rt[, "w"] = d.Rt[, "t_end"] %/% 7
  
  
  
  w.Rt = with(
    d.Rt,
    merge(
      x = aggregate(log.mu ~ w, FUN = mean),
      y = aggregate(log.sd ~ w, FUN = function(x) { sqrt(sum(x^2))/length(x)  } )
    )
  )
  
  gm.Rt = sapply(
    1:nrow(w.Rt),
    function(i) {
      with(w.Rt[i, ], normal.mu.sd.fx(log.mu, log.sd) )  
    }
  )
  
  w.Rt = cbind(w.Rt, t(gm.Rt))
  w.Rt = transform(
    w.Rt,
    LB       = qlnorm(0.025, log.mu, log.sd),
    UB       = qlnorm(0.975, log.mu, log.sd),
    weekend  = data[w * 7, "date"],
    dec.date = decimal_date(data[w * 7, "date"])
  )
  
  
  return(w.Rt)
  
}












#
#--- END ---
#




