library(lme4)
library(ggplot2)
library(GLMMadaptive)
library(dplyr)
library(statmod)
library(stringr)
celltype_count <- read.csv('All_celltypes_counts_LMM.csv', row.names = 1)
celltype_count_orig = celltype_count
celltype_count = celltype_count[celltype_count$timepoint %in% c("A", "B"), ]
colnames(celltype_count)[2] = 'size'

cc1 = celltype_count[celltype_count$dist %in% c("0um",  "0-5um"), ] %>% group_by(leuk_ID, celltype, patient, FOV, timepoint, Response) %>% summarise(count = sum(count))
cc1$dist = '0-5um'
cc2 = celltype_count[celltype_count$dist %in% c("5-10um",  "10-15um"), ] %>% group_by(leuk_ID, celltype, patient, FOV, timepoint, Response) %>% summarise(count = sum(count))
cc2$dist = "5-15um"
cc3 = celltype_count[celltype_count$dist %in% c("15-20um", "20-25um"), ] %>% group_by(leuk_ID, celltype, patient, FOV, timepoint, Response) %>% summarise(count = sum(count))
cc3$dist = "15-25um"
cc4 = celltype_count[celltype_count$dist %in% c("25-30um", "30-35um"), ] %>% group_by(leuk_ID, celltype, patient, FOV, timepoint, Response) %>% summarise(count = sum(count))
cc4$dist = "25-35um"
cc5 = celltype_count[celltype_count$dist %in% c("35-40um", "40-45um"), ] %>% group_by(leuk_ID, celltype, patient, FOV, timepoint, Response) %>% summarise(count = sum(count))
cc5$dist = "35-45um"

cc = as.data.frame(rbind(cc1, rbind(cc2, rbind(cc3, rbind(cc4, cc5)))))
cellcount = as.data.frame(cc %>% group_by(leuk_ID, dist) %>% summarise(total = sum(count)))
rownames(cellcount) = paste0(cellcount$leuk_ID, '_', cellcount$dist)
cc$dist_totalcount = cellcount[paste0(cc$leuk_ID, '_', cc$dist), 'total']
cc$dist = factor(cc$dist, levels = c("0-5um", "5-15um",  "15-25um", "25-35um", "35-45um"))
model = list()
for (ct in unique(cc$celltype)){
  print(ct)
  a = cc[cc$celltype == ct, ]
  a$count
  fm = glmer(count ~ Response * dist * timepoint + (1|FOV), data = cc[cc$celltype == ct, ], family = 'poisson', offset = log(dist_totalcount + 0.5))
  model[[ct]] = fm
}

coefout = data.frame()
for (i in 1:length(model)){
  out = as.data.frame(summary(model[[i]])$coefficients)
  out$celltype = names(model)[[i]]
  coefout = rbind(coefout, out)
}
coefout$adj.p = coefout$`Pr(>|z|)` * (nrow(coefout)/length(unique(coefout$celltype))-1) * length(unique(coefout$celltype))
coefout$adj.p[coefout$adj.p > 1] = 1
quantile(coefout$adj.p)
coefout = coefout[-grep('Intercept', rownames(coefout)), ]
write.csv(coefout, file = 'poisson_offset_distanceInteraction.csv')

# recalculate p-values
library(epimisc)
load('poisson_offsetdist_distanceInteraction.RData')
coefout = data.frame()
for (i in 1:length(model)){
  out = as.data.frame(summary(model[[i]])$coefficients)
  # responder vs non-responder at time A
  aout = data.frame()
  aout_rowname = c()
  for (j in 8:11){
    lb = rep(0, 20)
    lb[j] = 1
    lb[2] = 1
    a = lincomR(model[[i]], lb, conflev = 0.95, digits = 6)
    aout = rbind(aout, c(a$est, a$se.est, a$wald.z, a$pvalue))
    aout_rowname = c(aout_rowname, paste(rownames(summary(model[[i]])$coefficients)[lb == 1], collapse = ' + '))
  }
  rownames(aout) = aout_rowname
  colnames(aout) = colnames(out)
  out = rbind(out, aout)
  # responder time B vs time A
  aout = data.frame()
  aout_rowname = c()
  lb = rep(0, 20)
  lb[7] = 1
  lb[12] = 1
  a = lincomR(model[[i]], lb, conflev = 0.95, digits = 6)
  aout = rbind(aout, c(a$est, a$se.est, a$wald.z, a$pvalue))
  aout_rowname = c(aout_rowname, paste(rownames(summary(model[[i]])$coefficients)[lb == 1], collapse = ' + '))
  for (j in 13:16){
    lb = rep(0, 20)
    lb[7] = 1
    lb[12] = 1
    lb[j] = 1
    lb[j+4] = 1
    a = lincomR(model[[i]], lb, conflev = 0.95, digits = 6)
    aout = rbind(aout, c(a$est, a$se.est, a$wald.z, a$pvalue))
    aout_rowname = c(aout_rowname, paste(rownames(summary(model[[i]])$coefficients)[lb == 1], collapse = ' + '))
  }
  rownames(aout) = aout_rowname
  colnames(aout) = colnames(out)
  out = rbind(out, aout)
  out$celltype = names(model)[[i]]
  coefout = rbind(coefout, out)
}
write.csv(coefout, file = 'poisson_offset_DistanceInteraction3.csv')
