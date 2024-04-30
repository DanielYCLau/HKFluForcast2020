




#
#--- SETUP ---
#





load(list.files("data", "HK_flu_reg", full.names = T))
source("program/ILI/20240404_reg_sep_fore_Xw_fx.R")





#
#--- DATA INPUT ---
#





# <2  mins
start.t = Sys.time()

cl = makeCluster(n.cores, type = "SOCK")
registerDoSNOW(cl)

f.r2 = list()

for (f.period in c(4, 13)) {
  
  f.yr = 2020
  # f.period = 4
  data.fx(data, f.yr, f.period, 26)
  

  
  load(list.files("program/ILI", sprintf("reg_sep_fore_%02dw_%d_CV_result.Rdata", f.period, f.yr), full.names = T))
  
  cv.table = aggregate(. ~ model, data = cv.table, FUN = mean)
  cv.table = subset(cv.table, select = -tsCV.ind)
  cv.table = criteria.rank.fx(cv.table, rank.w.value = F)
  
  ( RMSE.model = with(cv.table, model[which.min(v.RMSE     )]) )
  ( LNQ.model  = with(cv.table, model[which.min(v.LNQ      )]) )
  ( MAE.model  = with(cv.table, model[which.min(v.MAE      )]) ) 
  ( WIS.model  = with(cv.table, model[which.min(v.WIS      )]) ) 
  
  
  
  # calculate r2
  model.list = sprintf("%s.model", c("WIS", "RMSE", "LNQ", "MAE"))
  
  r2.list = foreach(
    i         = model.list,
    .packages = c("MASS"),
    .export   = model.list,
    .combine  = c
  ) %dopar% {
    
    # model.i    = get("RMSE.model")
    model.i    = get(i)
    pred.train = pred.fx(model.i, train.d, train.d, c(0.05), keep.mc = F)
    r2         = sapply(pred.train[["pred.y.sum"]], function(x) { (cor(x[, c("f.obs", "fit.mean")])[2,1])^2 })
    
    output = list(r2)
    names(output) = i
    
    return(output)
  }
  
  f.r2[[sprintf("%02dw", f.period)]] = r2.list
  
}

stopCluster(cl)

end.t = Sys.time()
end.t - start.t

save(f.r2, file = file.fx("r2.Rdata"))





#
#--- plot ---
#




# range( sapply(f.r2[["13w"]], function(i) { i[1:4]  } ) )
# range( sapply(f.r2[["13w"]], function(i) { i[5:13] } ) )


pdf(
  width  = 12,
  height = 6,
  file   = file.fx("FigS2_r2.pdf")
)



# windows(width = 12, height = 6)
# quartz(width = 12, height = 6)
par(mfrow = c(1,2), mar = c(2,5,3,1), oma = c(3,2,0,0), las = 1)

for (f.period in c(4, 13)) {
  
  plot(
    NULL,
    xlim = range(1:f.period),
    ylim = c(0, 1),
    xaxt = "n",
    ylab = ""
  )
  abline(h = 0.5, col = "grey")
  
  if (f.period == 13) {
    axis(1, at = 1:f.period, labels = NA)
    axis(1, at = seq(1, 13, 3), tick = F)
  } else {
    axis(1, at = 1:f.period)
  }
  
  
  r2.list = f.r2[[sprintf("%02dw", f.period)]]
  N       = length(r2.list)
  for (i in 1:N) {
    
    r2.i = r2.list[[i]]
    points(r2.i, col = col.fx(i), pch = i)
    lines(r2.i, col = col.fx(i), lty = i)
  }
  
  if (f.period == 4) {
    mtext("A", side = 3, line = -1.5, adj = -0.225, font = 2, cex = 3)
  } else {
    mtext("B", side = 3, line = -1.5, adj = -0.225, font = 2, cex = 3)
  }
  
  
  
  legend(
    "topright",
    legend = ifelse(f.period < 10, "Short-term prediction", "Medium-term prediction"),
    bty = "n"
  )
  
  legend(
    "bottomleft", 
    legend = sprintf("%s       ", gsub("[.]model", " ", gsub("LNQ", "RMSLE", names(r2.list)) ) ),
    lty = 1:N,
    pch = 1:N,
    col = col.fx(1:N),
    bty = "n", 
    ncol = 3
  )
  
}


mtext(
  text  = expression(R^2),
  side  = 2, 
  line  = 0, 
  cex   = 1.25, 
  outer = T, 
  las   = 0
)

mtext(
  text  = "Forecast at specific week ahead",
  side  = 1, 
  line  = 1, 
  cex   = 1.25, 
  outer = T
)

dev.off()





#
#--- END ---
#




