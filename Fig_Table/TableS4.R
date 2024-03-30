




#
#--- SETUP ---
#





source("program/hosp/20230628_reg_sep_fore_Xw_fx.R")

f.yr = 2020





#
#--- LOAD DATA ---
#




cv.f = list()
for (f.period in c(4, 13)) {
  
  load(list.files("program/hosp", sprintf("reg_sep_fore_%02dw_%d_CV_result.Rdata", f.period, f.yr), full.names = T))
  cv.table = cv.table[, grep("t[.]", colnames(cv.table), invert = T)]
  cv.table = cv.table[, c("model", "v.WIS", "v.RMSE", "v.LNQ", "v.MAE")]
  
  cv.table = criteria.rank.fx(cv.table, rank.w.value = T)
  cv.table[, "model"] = model.abb.v2.fx(cv.table[, "model"])
  
  cv.f[[paste(f.period)]] = rbind(f.period, cv.table)
  
}

cv.f = do.call(rbind, cv.f)

write.csv(
  cv.f, 
  file = file.fx("TableS4_cv.csv"),
  row.names = F
)





#
#
#--- END ---
#
#




