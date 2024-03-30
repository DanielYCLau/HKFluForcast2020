




#
#--- SETUP ---
#





load(list.files("data", "HK_flu_reg", full.names = T))
source("program/ILI/20230619_reg_sep_fore_Xw_fx.R")

substrRight = function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }




#
#--- DATA INPUT ---
#





f.yr = 2020
f.period = 4
data.fx(data, f.yr, f.period, 10)
model.list.fx()






#
#--- PREDICTION ---
#





load(list.files("program/ILI", sprintf("reg_sep_fore_%02dw_%d_CV_result.Rdata", f.period, f.yr), full.names = T))
cv.table = criteria.rank.fx(cv.table, rank.w.value = F)

( RMSE.model = with(cv.table, model[which.min(v.RMSE     )]) )
( LNQ.model  = with(cv.table, model[which.min(v.LNQ      )]) )
( MAE.model  = with(cv.table, model[which.min(v.MAE      )]) ) 
( WIS.model  = with(cv.table, model[which.min(v.WIS      )]) ) 






#
#--- EFFECT ---
#




# 2.5 mins
# start.t = Sys.time()

cl = makeCluster(n.cores, type = "SOCK")
registerDoSNOW(cl)


nth.week = c(1,4)
effect.list = list()

model.list = sprintf("%s.model", c("WIS", "RMSE", "LNQ", "MAE"))

for ( model in model.list ) {
  
  # model = "RMSE.model"
  model.effect = list()
  par.list = c()
  for (par.i in c(sprintf("log_ILIpS_lag_%02d", 1:10), "bs", "temp", "AH", "ozone", "ScHD")) {
    if ( any( grepl(par.i, get(model)) ) ) {
      par.list = c(par.list, par.i)
    }
  }
  
  for (par in par.list) {
    
    effect.par = list()
    effect.par[["par"]] = par
    
    coln.par = switch(
      par, 
      "bs"               = "Week", 
      "log_ILIpS_lag_01" = "ILIpS",
      "log_ILIpS_lag_02" = "ILIpS",
      "log_ILIpS_lag_03" = "ILIpS",
      "log_ILIpS_lag_04" = "ILIpS",
      "log_ILIpS_lag_05" = "ILIpS",
      "log_ILIpS_lag_06" = "ILIpS",
      "log_ILIpS_lag_07" = "ILIpS",
      "log_ILIpS_lag_08" = "ILIpS",
      "log_ILIpS_lag_09" = "ILIpS",
      "log_ILIpS_lag_10" = "ILIpS",
      "AH"               = "AH", 
      "temp"             = "Mean_Temp",
      "ozone"            = "o3",
      "ScHD"             = "ScHD"
    )
    
    if (par == "bs") {
      par.x = seq(0, 53, length.out = 1e3)
    } else if ( par == "ILIpS" ) {
      par.x = log( range(train.d[, coln.par]) )
      par.x = exp( seq(par.x[1], par.x[2], length.out = 1e3) )
    } else {
      par.x = range(train.d[, coln.par])
      par.x = seq(par.x[1], par.x[2], length.out = 1e3)
    }
    
    effect.par[["x"]] = par.x
    
    effect.par.t = foreach(
      i = nth.week,
      .packages = c("MASS", "pomp"),
      .export  = c(model.list),
      .combine = c
    ) %dopar% {
      
      model.i = as.formula( sprintf("ILIpS_fore_%02d ~ %s", i, get(model)) )
      fit.i   = glm.nb(model.i, data = train.d, x = T, maxit = 1000)
      
      coln.x = names(coef(fit.i))[ grep(par, names(coef(fit.i)) ) ]
      coef.i = coef(fit.i)[coln.x]
      vcov.i = vcov(fit.i)[coln.x, coln.x]
      
      if ( any( grepl("bs", coln.x) ) ) { 
        x = seq(0, 53, length.out = 1e3)
        x = x[c(which(x >= 40), which(x < 40))]
        x = x/max(x)
        if ( any(grepl("_M_", coln.x)) ) {
          x  = periodic_bspline_basis(x, nbasis = 12, names = "t_lag_1_M_bs_%02d")
        } 
      } else {
        x = range(train.d[, coln.x[1]])
        x = seq(x[1], x[2], length.out = 1e3)
        if ( any( grepl("sq", coln.x) ) ) { x = cbind(x, x^2) }
      }
      
      set.seed(555)
      beta.i = mvrnorm(1e4, mu = coef.i, Sigma = vcov.i)
      
      y.i = as.matrix(x) %*% t(beta.i)
      y.i = apply(y.i, 1, CI95.fx)
      y.i = do.call(rbind, y.i)
      rownames(y.i) = NULL
      
      return( list(y.i) )
    }
    effect.par[["effect"]] = effect.par.t
    
    model.effect[[par]] = effect.par
  }
  
  effect.list[[model]] = model.effect
}

stopCluster(cl)

# end.t = Sys.time()
# end.t - start.t





#
#--- EFFECT ---
#



pdf(
  width  = 12, 
  height = 7, 
  file   = file.fx(sprintf("FigS1.pdf", f.period))
)

# windows()
# quartz()
par(mfrow = c(2,4), mar = c(5,4,3,1), oma = c(4,4,0,0), las = 1)

for (j in 1:length(nth.week)) { # c(1,4,8,13,26,39,52)) {
  
  par.model.list = list()
  
  for (par in c("AH", "temp", "ozone", "ScHD")) {
    
    par.list   = lapply(effect.list, function(x) { x[[par]] } )
    
    par.x      = unique( unlist( lapply(par.list, function(x) { x[["x"]] } ) ) )
    
    par.effect = lapply(par.list, function(x) { x[["effect"]][[j]] } )
    # par.effect = par.effect[!( sapply(par.effect, is.null) )]
    
    par.model.list[[par]] = list("par" = par, "x" = par.x, "effect" = par.effect)
  }
  
  par.model.list = par.model.list[sapply(par.model.list, function(x) { !is.null(x[["x"]]) } )]
  
  for (k in 1:length(par.model.list)) {
    
    par.k = par.model.list[[k]]
    
    par.name = switch(
      par.k[["par"]], 
      "AH"    = "Absolute humidity",
      "temp"  = "Mean temperature",
      "ozone" = "Ozone concentration",
      "ScHD"  = "School holiday"
    )
    
    par.lab = switch(
      par.k[["par"]], 
      "AH"    = "g/m^3",
      "temp"  = "Celsius~Degree~(C~degree)",
      "ozone" = "µg/m^3",
      "ScHD"  = "No.~of~Holiday"
    )
    
    par.x = par.k[["x"]]
    x.lim = range( par.k[["x"]] )
    y.lim = range( unlist( par.k[["effect"]] ) )
    
    
    plot(
      NULL, 
      xlim = x.lim,
      ylim = exp(y.lim),
      main = par.name,
      log  = "y",
      ylab = "",
      xlab = eval(parse(text = sprintf("expression(%s)", par.lab)))
    )
    
    legend(
      "topleft", 
      bty = "n", 
      legend = bquote( paste(.(nth.week[j])^.(substrRight(ordinal(nth.week[j]), 2)), "week") )
        # sprintf("%s week", ordinal(nth.week[j]))
    )
    abline(h = 1, col = "grey")
    
    for (m in 1:length(par.k[["effect"]])) {
      
      par.k.m = par.k[["effect"]][[m]]
      par.k.m[, "LB"] = par.k.m[, "2.5%"]
      par.k.m[, "UB"] = par.k.m[, "97.5%"]
      
      if ( !is.null(par.k.m) ) {
        with(
          par.k.m, {
            lines(x = par.x, y = exp(mean), lwd = 2, lty = m, col = col.fx(m))
            polygon.fx(t = par.x, LB = exp(LB), UB = exp(UB), col = alpha(col.fx(m), 0.1))
          }
        )
      }
      
    }
  }
}



mtext(
  text  = "Effect on predicted ILI+ proxy of n-th week ahead\n(Short-term prediction)",
  side  = 2,
  line  = 0,
  cex   = 1.25,
  outer = T,
  las   = 0
)

legend(
  "bottomright", xpd = NA, inset = c(1.5, -0.45),
  legend = sprintf("%s    ", c("WIS", "RMSE", "RMSLE", "MAE")),
  lty = 1:5, lwd = 2, col = col.fx(1:5),
  horiz = T, bty = "n"
)

dev.off()





#
#--- END ---
#




