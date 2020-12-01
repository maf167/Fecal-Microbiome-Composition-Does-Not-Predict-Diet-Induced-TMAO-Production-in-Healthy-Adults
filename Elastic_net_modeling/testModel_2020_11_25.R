### Prediction analysis for the Approach metagenomics paper
# Adding cutC as a predictor
## Marc Ferrell 11/25/2020

rm(list=ls())

setwd("//smb-isi1.lerner.ccf.org/CCG/Ferrell/Figures/Model-Building-Analysis-2/cutC")

## Method Imports

source("loadData_2020_11_25.R")
source("preprocessing_2020_11_25.R")
source("VarExplained_2020_11_25.R")
source("fitModels_2020_11_25.R")


filepaths  = c(
  "clin"   = "CLIN.txt",
  "PID"   = "clinical_data.txt",
  "lib"    = "SeqDepth.txt",
  "cutC"   = "cutC-len2544-74-tot-count-2020-06-22.txt",
  "fdata"  = "FECAL.APP.Rdata",
  "pdata"  = "PLASMA.APP.Rdata",
  "smap"   = "Smap.txt"
)

t <- makeDataTables(filepaths)

dt <- mkDataTableDelta(filepaths["PID"], t$dat, t$cols)

cutC.perform <- structure(
  list("Plasma TMAO" = structure(list(), class="VEList"),
       "Urine TMAO" = structure(list(), class="VEList"),
       "Plasma TMAO - Abs Change" = structure(list(), class="VEList"),
       "Urine TMAO - Abs Change" = structure(list(), class="VEList"),
       "Plasma TMAO - Percent Change" = structure(list(), class="VEList"),
       "Urine TMAO - Percent Change" = structure(list(), class="VEList"),
       "Fecal d6 Choline-18hr" = structure(list(), class="VEList"),
       "Fecal d6 TMA-18hr" = structure(list(), class="VEList")
), class="VEList")

for(i in 1:100){
  print(i)
  # Plasma TMAO
  cutC.perform[["Plasma TMAO"]][[i]] <- fitModels(
                        dat = as.matrix(t$dat),
                        pred = list(
                          "Clinical" = which(t$cols == "Base"),
                          "cutC" = which(t$cols == "cutC"),
                           "Clinical + cutC" = which(t$cols %in% c("Base", "cutC"))
                        ),
                    y = which(colnames(t$dat) == "Plasma.TMAO...然.")
                    )
  # Urine TMAO
  cutC.perform[["Urine TMAO"]][[i]] <-  fitModels(
                        dat = as.matrix(t$dat),
                        pred = list(
                          "Clinical" = which(t$cols == "Base"),
                          "cutC" = which(t$cols == "cutC"),
                          "Clinical + cutC" = which(t$cols %in% c("Base", "cutC"))
                        ),
                        y = which(colnames(t$dat) == "Urine TMAO")
                        )
  # Plasma TMAO - Absolute Change
  cutC.perform[["Plasma TMAO - Abs Change"]][[i]] <-  fitModels(
                        dat = as.matrix(dt$dtab),
                        pred = list(
                          "Clinical" = which(dt$cols == "Base"),
                          "cutC" = which(dt$cols == "cutC"),
                          "Clinical + cutC" = which(dt$cols %in% c("Base", "cutC"))
                        ),
                        y = which(colnames(dt$dtab) == "Plasma.TMAO...然.")
                        )
  # Urine TMAO - Absolute Change
  cutC.perform[["Urine TMAO - Abs Change"]][[i]] <-  fitModels(
                        dat = as.matrix(dt$dtab),
                        pred = list(
                          "Clinical" = which(dt$cols == "Base"),
                          "cutC" = which(dt$cols == "cutC"),
                          "Clinical + cutC" = which(dt$cols %in% c("Base", "cutC"))
                        ),
                        y = which(colnames(dt$dtab) == "Urine TMAO")
                        )
  # Plasma TMAO - Percent Change
  cutC.perform[["Plasma TMAO - Percent Change"]][[i]] <-  fitModels(
                        dat = as.matrix(dt$ddtab),
                        pred = list(
                          "Clinical" = which(dt$cols == "Base"),
                          "cutC" = which(dt$cols == "cutC"),
                          "Clinical + cutC" = which(dt$cols %in% c("Base", "cutC"))
                        ),
                        y = which(colnames(dt$ddtab) == "Plasma.TMAO...然.")
                        )
  # Urine TMAO - Percent Change
  cutC.perform[["Urine TMAO - Percent Change"]][[i]] <-  fitModels(
                        dat = as.matrix(dt$ddtab),
                        pred = list(
                          "Clinical" = which(dt$cols == "Base"),
                          "cutC" = which(dt$cols == "cutC"),
                          "Clinical + cutC" = which(dt$cols %in% c("Base", "cutC"))
                        ),
                        y = which(colnames(dt$ddtab) == "Urine TMAO")
                        )
  
  # Fecal d6 Choline - 18hr
  cutC.perform[["Fecal d6 Choline-18hr"]][[i]] <-  fitModels(
                        dat = as.matrix(t$dat),
                        pred = list(
                          "Clinical"        = which(t$cols == "Base"),
                          "cutC"            = which(t$cols == "cutC"),
                          "Clinical + cutC" = which(t$cols %in% c("Base", "cutC"))
                        ),
                        y = which(colnames(t$dat) == "d6.Choline..18.Hr..然.")
                        )
  # Fecal d6 TMA - 18hr
  cutC.perform[["Fecal d6 TMA-18hr"]][[i]] <-  fitModels(
                        dat = as.matrix(t$dat),
                        pred = list(
                          "Clinical"        = which(t$cols == "Base"),
                          "cutC"            = which(t$cols == "cutC"),
                          "Clinical + cutC" = which(t$cols %in% c("Base", "cutC"))
                        ),
                        y = which(colnames(t$dat) == "d6.TMA..18.Hr..然.")
                        )
  
}
save(cutC.perform, file="cutC.perform-112620.Rdata")

vs.cutC <- lapply(cutC.perform, function(x) VESum(x))
vs.cutC <- do.call("rbind", vs.cutC)

vs.cutC$Outcome[1201:1800] <- "Plasma TMAO - Abs Change"
vs.cutC$Outcome[1801:2400] <- "Urine TMAO - Abs Change"
vs.cutC$Outcome[2401:3000] <- "Plasma TMAO - Percent Change"
vs.cutC$Outcome[3001:3600] <- "Urine TMAO - Percent Change"

write.table(vs.cutC, file="vs.cutC-112620.txt", sep="\t", quote=FALSE,
            row.names = FALSE)

## cutC lowers RMSE for fecal d6 choline

Ot <- "d6.Choline..18.Hr..然."
et <- "Rsq"

ggplot(vs.cutC[vs.cutC$Err_Type == et & vs.cutC$Outcome == Ot,], aes(x=factor(Pred), y=Err)) + 
  geom_violin() + geom_boxplot(width=0.1) +ylim(c(-.01,.01)) 



ggplot(vs.cutC[vs.cutC$Err_Type == et & vs.cutC$Outcome == Ot & vs.cutC$Pred == "Clinical",],
       aes(x = Err)) + 
  geom_density()

qqnorm(vs.cutC[vs.cutC$Err_Type == et & vs.cutC$Outcome == Ot & vs.cutC$Pred == "cutC", "Err"])
qqline(vs.cutC[vs.cutC$Err_Type == et & vs.cutC$Outcome == Ot & vs.cutC$Pred == "cutC", "Err"])

wilcox.test(Err ~ Pred, data = vs.cutC, 
            subset = Err_Type == et & Outcome == Ot & Pred %in% c("Clinical","Clinical + cutC"))