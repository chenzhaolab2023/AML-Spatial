library(dplyr)
library(stringr)
setwd("/data/guig2/2203_ZhaoCollab/round3/LigandReceptor_analysis")
lr = read.csv('Ligand_receptor_input.csv', row.names = 1)
lr = lr[lr$group != '', ]
lrpair = read.table('/data/guig2/2203_ZhaoCollab/manuscript/LigandReceptor/lr_pair.txt')
# lrpair = lrpair[, c("V2", 'V1')]
pall = c(51:53, 56:58)

meta = lr[, c('patients', 'FOV', 'Timepoint')]
meta = meta[!duplicated(meta), ]
rownames(meta) = paste0('P', pall[as.numeric(gsub('p', '', meta$patients))], '_FOV', ifelse(meta$FOV < 10, paste0('0', meta$FOV), meta$FOV))

outall_list = data.frame()

for (p in 1:6){
  load(paste0('/data/guig2/2203_ZhaoCollab/round3/P', pall[p], '_raw.RData'))
  print(paste0('P', pall[p]))
  mat = mat[, -grep('Neg', colnames(mat))]
  lrp = lr[lr$patients == paste0('p', p), ]
  rownames(lrp) = paste0('P', pall[p], "_FOV", ifelse(lrp$FOV < 10, paste0('0', lrp$FOV), lrp$FOV), '_cell_', lrp$cell_id)
  mat = cbind(mat[rownames(lrp), unique(c(lrpair[,1], lrpair[,2]))], lrp[, c('Timepoint', 'celltype_detail', 'group', 'FOV')])
  outall = data.frame()
  for (gene in 1:(ncol(mat)-4)){
    out = eval(parse(text = paste0('mat %>% group_by(Timepoint, celltype_detail, group, FOV) %>% summarise(mean = mean(', colnames(mat)[gene], '))'))) %>% as.data.frame
    out$gene = colnames(mat)[gene]
    outall = rbind(outall, out)
  }
  result = data.frame()
  for (gene in 1:(ncol(mat)-4)){
    out = outall[outall$gene == unique(outall$gene)[gene], ]
    for (Timepoint in LETTERS[1:3]){
      for (celltype in unique(outall$celltype_detail)){
        a1 = out[out$Timepoint == Timepoint & out$celltype_detail == celltype & out$group == '0-5 um', 'mean']
        a2 = out[out$Timepoint == Timepoint & out$celltype_detail == celltype & out$group == '30+ um', 'mean']
        if (length(a1) > 3 & length(a2) > 3){
          if (sum(a1 > 0) > 3 & sum(a2 > 0) > 3){
            a = wilcox.test(a1, a2)
            result = rbind(result, data.frame(gene = unique(outall$gene)[gene], Timepoint = Timepoint, celltype = celltype, p = a$p.value, group0 = median(a1), group1 = median(a2)))
          }
        }
      }
    }
  }
  resultshort = result[result$p < 0.01, ]
  mat = countall[, -grep('Neg', colnames(countall))]
  all_sub = data.frame(Timepoint = meta[str_split(rownames(mat), pattern = fixed('_cell'), simplify = T)[,1], 'Timepoint'], FOV = as.numeric(gsub('FOV', '', str_split(rownames(mat), pattern = fixed('_'), simplify = T)[,2])), celltype_detail = celltypeall)
  rownames(all_sub) = rownames(mat)
  all_sub = all_sub[all_sub$celltype_detail %in% setdiff(all_sub$celltype_detail, c('Unknown', 'SmallCell', 'RBC', 'NoCellAssigned')), ]
  mat = cbind(mat[rownames(all_sub), unique(c(lrpair[,1], lrpair[,2]))], all_sub)
  outall = data.frame()
  for (gene in 1:(ncol(mat)-3)){
    out = eval(parse(text = paste0('mat %>% group_by(Timepoint, celltype_detail, FOV) %>% summarise(mean = mean(', colnames(mat)[gene], '))'))) %>% as.data.frame
    out$gene = colnames(mat)[gene]
    outall = rbind(outall, out)
  }
  medgene = outall %>% group_by(Timepoint, gene, celltype_detail) %>% summarise(med = median(mean)) %>% as.data.frame
  a = medgene %>% group_by(gene, Timepoint) %>% summarise(sd = sd(med), mean = median(med)) %>% as.data.frame
  a$malignant = medgene[medgene$celltype_detail == 'LeukemiaCell', 'med']
  lrmalig = a[a$malignant > a$mean + 2*a$sd, ]
  lrpairnew = lrpair[lrpair[,1] %in% lrmalig$gene | lrpair[,2] %in% lrmalig$gene, ]
  outall = data.frame()
  for (li in 1:nrow(lrpairnew)){
    if (lrpair[li, 1] %in% lrmalig$gene){
      genetime = lrmalig[lrmalig$gene == lrpair[li, 1], ]
      for (ti in 1:length(unique(genetime$Timepoint))){
        idall = resultshort$Timepoint == unique(genetime$Timepoint)[ti] & resultshort$gene == lrpair[li, 2]
        if (sum(idall) > 0){
          out = resultshort[idall, ]
          out$ligand_gene_malignant = lrpair[li, 1]
          outall = rbind(outall, out)
        }
      }
    }
  }
  outall$patient = p
  outall_list = rbind(outall_list, outall)
}
write.csv(outall_list, file = 'LR_result1.csv')

# switch the order of reference list and generate the result for the other direction.