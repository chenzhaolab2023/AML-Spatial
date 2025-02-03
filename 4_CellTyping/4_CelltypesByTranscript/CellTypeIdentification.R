# Cell type idenfication of spatial single-cell data for six patients
# data is a list of length 6 containing the expression matrix and metadata for six patients 
setwd("")
for (pa in 1:6){
  print(pa)
  mat = t(data[[pa]][[1]]) # get the expression matrix for patient pa
  mat = mat[, str_split(colnames(mat), pattern = fixed('_'), simplify = T)[,2] != 0]
  cd34 = data[[pa]][[2]]
  # cd34_all contains all annotated malignant cells
  # rbc_all contaiins all annotated RBC
  cd34$malignant = ifelse(rownames(cd34) %in% cd34_all[cd34_all$pa == pa, 'cell'], 'malignant', ifelse(rownames(cd34) %in% rbc_all[rbc_all$pa == pa, 'cell'], "RBC" , 'nonmalignant'))
  cd34$timepoints = meta[paste0(pa,  str_split(rownames(cd34), pattern = fixed("_"), simplify = T)[,1]), 'timepoints']
  cd34$fov = str_split(cd34$cell_ID, pattern = fixed('_'), simplify = T)[,1]
  neg = t(mat[grep("Neg", rownames(mat)), ])
  neg = neg[rowSums(mat) > 0, ]
  mat = mat[rowSums(mat) > 0, ]
  negmean <- rowMeans(neg)
  cd34 = cd34[colnames(mat),]
  ref_low = quantile(cd34[cd34$malignant == 'nonmalignant', 'CD34_median'], prob = seq(0, 1, by =0.1))
  ref_low[1] = ref_low[1]-1
  cd34$malignant2 = as.character(cut(cd34$CD34_median, ref_low, labels = paste0('nonmalig', 1:10)))
  cd34$malignant2[cd34$malignant == 'malignant'] = 'malignant'
  cd34$celltype_rmMRBC = cd34$malignant # remove morphology identified malignant and RBC
  cd34$celltype_rmRBC_addMRBC = cd34$malignant # remove morphology identified RBC, use morphology identified malignant cells as anchor cells for additional malignant cells, and add additional cells that might be RBC
  cd34$celltype_rmMRBC_detail = cd34$malignant
  cd34$celltype_rmRBC_addMRBC_detail = cd34$malignant
  cd34$anchor_rmMRBC = ''
  cd34$anchor_rmRBC_addMRBC = ''
  cd34$anchor_rmMRBC_detail = ''
  cd34$anchor_rmRBC_addMRBC_detail = ''
  # get the counts for all cells except for RBC
  counts = t(mat[-grep("Neg", rownames(mat)), cd34$malignant != 'RBC'])
  neg0 = neg[cd34$malignant != 'RBC', ]
  cohort = cd34[cd34$malignant != 'RBC', 'malignant2']
  neg1 = neg0[rowSums(counts) > 0, ]
  cohort1 = cohort[rowSums(counts) > 0]
  cd34_1 = cd34[cd34$malignant != 'RBC', ][rowSums(counts) > 0,]
  counts1 = counts[rowSums(counts) > 0, ]
  for (f in 1:length(unique(cd34_1$fov))){
    print(f)
    neg2 = neg1[cd34_1$fov == unique(cd34_1$fov)[f], ]
    counts2 = counts1[cd34_1$fov == unique(cd34_1$fov)[f], ]
    cohort2 = cohort1[cd34_1$fov == unique(cd34_1$fov)[f]]
    names(cohort2) = rownames(counts2)
    ctref = c("malignant_bygene", "RBC_bygene", 'Mega', 'CD8T', 'CD4T', 'NK', 'B', 'MonocytesCD14', 'MonocytesCD16', 'MonocytesProgenitor', 'DC')
    anchor = data.frame(matrix(rep('', nrow(counts2) * length(ctref)), nrow = nrow(counts2)))  #  * 2
    rownames(anchor) = rownames(neg2)
    colnames(anchor) = c(ctref) # , paste0(ctref, '_count')
    anchor$malignant_bygene = rownames(anchor) %in% rownames(cd34)[cd34$malignant == 'malignant']
    anchor$RBC_bygene = rowSums(counts2[,c('HBB', 'HBA1')])/rowSums(counts2) > 0.2
    anchor$Mega = rowSums(counts2[,c('PF4', 'PPBP')] > 1) == 2
    anchor$CD8T = (rowSums(counts2[,c('CD3E', 'CD3D', 'CD3G')]) > 0) + (rowSums(counts2[,c('CD8A', 'CD8B')]) > 0) + (counts2[,'CD4'] == 0) == 3
    anchor$CD4T = (rowSums(counts2[,c('CD3E', 'CD3D', 'CD3G')]) > 0) + (rowSums(counts2[,c('CD8A', 'CD8B')]) == 0) + (counts2[,'CD4'] > 0) == 3
    anchor$NK = (rowSums(counts2[,c('GNLY', 'NKG7')]) > 0) + (counts2[,'FCGR3A'] > 0) == 2
    anchor$B = (rowSums(counts2[,c('CD19', 'CD79A')]) > 1)
    anchor$MonocytesCD14 = (counts2[,'CD14'] > 0 & rowSums(counts2[,c("CSF3R", "S100A9", "CD33", "LYZ")]) > 0)
    anchor$MonocytesCD16 = (counts2[,'FCGR3A'] > 0 & rowSums(counts2[,c('GNLY', 'NKG7')]) == 0 & rowSums(counts2[,c("CD33", "LYZ")]) > 0)
    anchor$MonocytesProgenitor = (rowSums(counts2[,c('MPO', 'ELANE')] > 0) == 2)
    anchor$DC = (rowSums(counts2[,c('IL3RA', 'CD33', "LYZ")] > 0) > 1)
    anchor_final = rep('', nrow(counts2))
    names(anchor_final) = rownames(anchor)
    anchor_final[rowSums(anchor) == 1] = apply(anchor[rowSums(anchor) == 1, ], 1, function(x){return(ctref[x])})
    anchor_final[rowSums(counts2)[names(anchor_final)] < 50] = ''
    mean_mat = aggregate(counts2[anchor_final != '', ], list(anchor_final[anchor_final != '']), mean)
    profiles = t(mean_mat[,2:ncol(mean_mat)])
    colnames(profiles) = mean_mat[,1]
    profiles = profiles[, colnames(profiles) %in% names(table(anchor_final))[table(anchor_final)>1]]
    if (length(profiles) < 1000){
      next
    }
    sup1 <- insitutypeML(x = counts2,neg = rowMeans(neg2),cohort = cohort2,reference_profiles = profiles)
    cd34[names(anchor_final), 'anchor_rmRBC_addMRBC'] = anchor_final
    cd34[names(sup1$clust), 'celltype_rmRBC_addMRBC'] = sup1$clust
    cd34[names(sup1$clust), 'celltype_rmRBC_addMRBC_detail'] = sup1$clust
    # immune cell subtype
    for (subct in c('CD4T', 'CD8T', 'B')){
      print(subct)
      counts3 = counts2[names(sup1$clust)[sup1$clust == subct], ]
      cohort3 = cohort2[names(sup1$clust)[sup1$clust == subct]]
      neg3 = counts2[names(sup1$clust)[sup1$clust == subct], ]
      anchor = rep('', nrow(counts3))
      names(anchor) = rownames(neg3)
      if (subct == 'CD4T'){
        anchor[counts3[,'IL7R'] > 0 & counts3[,'CCR7'] == 0] = "CD4Tmemory"
        anchor[counts3[,'CCR7'] > 0 & counts3[,'IL7R'] == 0] = "CD4Tnaive"
      }else if (subct == 'CD8T'){
        # anchor[counts2[,'CCL5'] > 0] = "CD8Teffector"
        anchor[counts3[,'CCR7'] > 0 & counts3[,'CCL5'] == 0] = "CD8Tnaive"
        anchor[counts3[,'CCL5'] > 0 & counts3[,'GZMH'] == 0 & counts3[,'GZMK'] > 0] = "CD8TeffectorGZMK"
        anchor[counts3[,'CCL5'] > 0 & counts3[,'GZMH'] > 0 & counts3[,'GZMK'] == 0] = "CD8TeffectorGZMH"
      }else{
        anchor[counts3[,'MS4A1'] > 0 & counts3[,'TCL1A'] == 0 & counts3[,'JCHAIN'] == 0] = "MatureB"
        anchor[counts3[,'MS4A1'] == 0 & counts3[,'TCL1A'] > 0 & counts3[,'JCHAIN'] == 0] = "ProgenitorB"
        anchor[counts3[,'MS4A1'] == 0 & counts3[,'TCL1A'] == 0 & counts3[,'JCHAIN'] > 0] = "Plasma"
      }
      anchor[rowSums(counts3)[names(anchor)] < 50] = ''
      if (length(anchor) == 0 | length(setdiff(unique(anchor), '')) < 2){
        next
      }
      mean_mat = aggregate(counts3[anchor != '', ], list(anchor[anchor != '']), mean)
      profiles = t(mean_mat[,2:ncol(mean_mat)])
      colnames(profiles) = mean_mat[,1]
      if (length(profiles) < 1000){
        next
      }
      sup2 <- insitutypeML(x = counts3,neg = rowMeans(neg3),cohort = cohort3,reference_profiles = profiles)
      cd34[names(sup2$clust), 'celltype_rmRBC_addMRBC_detail'] = sup2$clust
      cd34[names(anchor), 'anchor_rmRBC_addMRBC_detail'] = anchor
    }
  }
  # get the counts for nonmalignant population
  counts = t(mat[-grep("Neg", rownames(mat)), cd34$malignant == 'nonmalignant'])
  neg0 = neg[cd34$malignant == 'nonmalignant', ]
  cohort = cd34[cd34$malignant == 'nonmalignant', 'malignant2']
  neg1 = neg0[rowSums(counts) > 0, ]
  cohort1 = cohort[rowSums(counts) > 0]
  cd34_1 = cd34[cd34$malignant == 'nonmalignant', ][rowSums(counts) > 0,]
  counts1 = counts[rowSums(counts) > 0, ]
  for (f in 1:length(unique(cd34_1$fov))){
    print(f)
    neg2 = neg1[cd34_1$fov == unique(cd34_1$fov)[f], ]
    counts2 = counts1[cd34_1$fov == unique(cd34_1$fov)[f], ]
    cohort2 = cohort1[cd34_1$fov == unique(cd34_1$fov)[f]]
    names(cohort2) = rownames(counts2)
    ctref = c('Mega', 'CD8T', 'CD4T', 'NK', 'B', 'MonocytesCD14', 'MonocytesCD16', 'MonocytesProgenitor', 'DC')
    anchor = data.frame(matrix(rep('', nrow(counts2) * length(ctref)), nrow = nrow(counts2)))  #  * 2
    rownames(anchor) = rownames(neg2)
    colnames(anchor) = c(ctref)
    anchor$Mega = rowSums(counts2[,c('PF4', 'PPBP')] > 1) == 2
    anchor$CD8T = (rowSums(counts2[,c('CD3E', 'CD3D', 'CD3G')]) > 0) + (rowSums(counts2[,c('CD8A', 'CD8B')]) > 0) + (counts2[,'CD4'] == 0) == 3
    anchor$CD4T = (rowSums(counts2[,c('CD3E', 'CD3D', 'CD3G')]) > 0) + (rowSums(counts2[,c('CD8A', 'CD8B')]) == 0) + (counts2[,'CD4'] > 0) == 3
    anchor$NK = (rowSums(counts2[,c('GNLY', 'NKG7')]) > 0) + (counts2[,'FCGR3A'] > 0) == 2
    anchor$B = (rowSums(counts2[,c('CD19', 'CD79A')]) > 1)
    anchor$MonocytesCD14 = (counts2[,'CD14'] > 0 & rowSums(counts2[,c("CSF3R", "S100A9", "CD33", "LYZ")]) > 0)
    anchor$MonocytesCD16 = (counts2[,'FCGR3A'] > 0 & rowSums(counts2[,c('GNLY', 'NKG7')]) == 0 & rowSums(counts2[,c("CD33", "LYZ")]) > 0)
    anchor$MonocytesProgenitor = (rowSums(counts2[,c('MPO', 'ELANE')] > 0) == 2)
    anchor$DC = (rowSums(counts2[,c('IL3RA', 'CD33', "LYZ")] > 0) > 1)
    anchor_final = rep('', nrow(counts2))
    names(anchor_final) = rownames(anchor)
    anchor_final[rowSums(anchor) == 1] = apply(anchor[rowSums(anchor) == 1, ], 1, function(x){return(ctref[x])})
    anchor_final[rowSums(counts2)[names(anchor_final)] < 50] = ''
    mean_mat = aggregate(counts2[anchor_final != '', ], list(anchor_final[anchor_final != '']), mean)
    profiles = t(mean_mat[,2:ncol(mean_mat)])
    colnames(profiles) = mean_mat[,1]
    profiles = profiles[, colnames(profiles) %in% names(table(anchor_final))[table(anchor_final)>1]]
    if (length(profiles) < 1000){
      next
    }
    sup1 <- insitutypeML(x = counts2,neg = rowMeans(neg2),cohort = cohort2,reference_profiles = profiles)
    cd34[names(sup1$clust), 'celltype_rmMRBC'] = sup1$clust
    cd34[names(anchor_final), 'anchor_rmMRBC'] = anchor_final
    cd34[names(sup1$clust), 'celltype_rmMRBC_detail'] = sup1$clust
    # immune cell subtype
    for (subct in c('CD4T', 'CD8T', 'B')){
      print(subct)
      counts3 = counts2[names(sup1$clust)[sup1$clust == subct], ]
      cohort3 = cohort2[names(sup1$clust)[sup1$clust == subct]]
      neg3 = counts2[names(sup1$clust)[sup1$clust == subct], ]
      anchor = rep('', nrow(counts3))
      names(anchor) = rownames(neg3)
      if (subct == 'CD4T'){
        anchor[counts3[,'IL7R'] > 0 & counts3[,'CCR7'] == 0] = "CD4Tmemory"
        anchor[counts3[,'CCR7'] > 0 & counts3[,'IL7R'] == 0] = "CD4Tnaive"
      }else if (subct == 'CD8T'){
        anchor[counts3[,'CCR7'] > 0 & counts3[,'CCL5'] == 0] = "CD8Tnaive"
        anchor[counts3[,'CCL5'] > 0 & counts3[,'GZMH'] == 0 & counts3[,'GZMK'] > 0] = "CD8TeffectorGZMK"
        anchor[counts3[,'CCL5'] > 0 & counts3[,'GZMH'] > 0 & counts3[,'GZMK'] == 0] = "CD8TeffectorGZMH"
      }else{
        anchor[counts3[,'MS4A1'] > 0 & counts3[,'TCL1A'] == 0 & counts3[,'JCHAIN'] == 0] = "MatureB"
        anchor[counts3[,'MS4A1'] == 0 & counts3[,'TCL1A'] > 0 & counts3[,'JCHAIN'] == 0] = "ProgenitorB"
        anchor[counts3[,'MS4A1'] == 0 & counts3[,'TCL1A'] == 0 & counts3[,'JCHAIN'] > 0] = "Plasma"
      }
      anchor[rowSums(counts3)[names(anchor)] < 50] = ''
      if (length(anchor) == 0 | length(setdiff(unique(anchor), '')) < 2){
        next
      }
      mean_mat = aggregate(counts3[anchor != '', ], list(anchor[anchor != '']), mean)
      profiles = t(mean_mat[,2:ncol(mean_mat)])
      colnames(profiles) = mean_mat[,1]
      if (length(profiles) < 1000){
        next
      }
      sup2 <- insitutypeML(x = counts3,neg = rowMeans(neg3),cohort = cohort3,reference_profiles = profiles)
      cd34[names(sup2$clust), 'celltype_rmMRBC_detail'] = sup2$clust
      cd34[names(anchor), 'anchor_rmMRBC_detail'] = anchor
    }
  }
  cd34$total_transcript = colSums(mat[-grep("Neg", rownames(mat)), rownames(cd34)])
  write.csv(cd34[cd34$celltype_rmMRBC_detail != 'nonmalignant', ], file= paste0('p', pa, '_celltype_v5.csv'))
}



## T cell subtype
all = list()
for (pa in 1:6){
  a = read.csv(paste0('p', pa, '_celltype_v5.csv'), row.names = 1)
  a$patient = pa
  all[[pa]] = a
}
for (pa in 1:6){
  print(pa)
  ct = all[[pa]]
  mat = t(data[[pa]][[1]])
  mat = mat[, str_split(colnames(mat), pattern = fixed('_'), simplify = T)[,2] != 0]
  cd34 = data[[pa]][[2]]
  cd34$malignant = ifelse(rownames(cd34) %in% cd34_all[cd34_all$pa == pa, 'cell'], 'malignant', ifelse(rownames(cd34) %in% rbc_all[rbc_all$pa == pa, 'cell'], "RBC" , 'nonmalignant'))
  cd34$timepoints = meta[paste0(pa,  str_split(rownames(cd34), pattern = fixed("_"), simplify = T)[,1]), 'timepoints']
  cd34$fov = str_split(cd34$cell_ID, pattern = fixed('_'), simplify = T)[,1]
  neg = t(mat[grep("Neg", rownames(mat)), ])
  neg = neg[rowSums(mat) > 0, ]
  mat = mat[rowSums(mat) > 0, ]
  negmean <- rowMeans(neg)
  cd34 = cd34[colnames(mat),]
  ref_low = quantile(cd34[cd34$malignant == 'nonmalignant', 'CD34_median'], prob = seq(0, 1, by =0.1))
  ref_low[1] = ref_low[1]-1
  cd34$malignant2 = as.character(cut(cd34$CD34_median, ref_low, labels = paste0('nonmalig', 1:10)))
  cd34$malignant2[cd34$malignant == 'malignant'] = 'malignant'
  ### get the counts for nonmalignant population ##
  counts = t(mat[-grep("Neg", rownames(mat)), cd34$malignant == 'nonmalignant'])
  neg0 = neg[cd34$malignant == 'nonmalignant', ]
  cohort = cd34[cd34$malignant == 'nonmalignant', 'malignant2']
  neg1 = neg0[rowSums(counts) > 0, ]
  cohort1 = cohort[rowSums(counts) > 0]
  cd34_1 = cd34[cd34$malignant == 'nonmalignant', ][rowSums(counts) > 0,]
  counts1 = counts[rowSums(counts) > 0, ]
  cd34$celltype = ct[rownames(cd34), "celltype"]
  cd34$anchor = ''
  for (f in 1:length(unique(cd34_1$timepoints))){
    print(f)
    counts3 = counts1[cd34_1$timepoints == unique(cd34_1$timepoints)[f], ]
    cohort3 = cohort1[cd34_1$timepoints == unique(cd34_1$timepoints)[f]]
    neg3 = neg1[cd34_1$timepoints == unique(cd34_1$timepoints)[f], ]
    for (subct in c('CD4T', 'CD8T', 'B')){
      print(subct)
      counts2 = counts3[ct[rownames(counts3), "celltype"] == subct, ]
      cohort2 = cohort3[ct[rownames(counts3), "celltype"] == subct]
      neg2 = neg3[ct[rownames(counts3), "celltype"] == subct, ]
      anchor = rep('', nrow(counts2))
      names(anchor) = rownames(neg2)
      if (subct == 'CD4T'){
        anchor[counts2[,'IL7R'] > 0 & counts2[,'CCR7'] == 0] = "CD4Tmemory"
        anchor[counts2[,'CCR7'] > 0 & counts2[,'IL7R'] == 0] = "CD4Tnaive"
      }else if (subct == 'CD8T'){
        # anchor[counts2[,'CCL5'] > 0] = "CD8Teffector"
        anchor[counts2[,'CCR7'] > 0 & counts2[,'CCL5'] == 0] = "CD8Tnaive"
        anchor[counts2[,'CCL5'] > 0 & counts2[,'GZMH'] == 0 & counts2[,'GZMK'] > 0] = "CD8TeffectorGZMK"
        anchor[counts2[,'CCL5'] > 0 & counts2[,'GZMH'] > 0 & counts2[,'GZMK'] == 0] = "CD8TeffectorGZMH"
      }else{
        anchor[counts2[,'MS4A1'] > 0 & counts2[,'TCL1A'] == 0 & counts2[,'JCHAIN'] == 0] = "MatureB"
        anchor[counts2[,'MS4A1'] == 0 & counts2[,'TCL1A'] > 0 & counts2[,'JCHAIN'] == 0] = "ProgenitorB"
        anchor[counts2[,'MS4A1'] == 0 & counts2[,'TCL1A'] == 0 & counts2[,'JCHAIN'] > 0] = "Plasma"
      }
      mean_mat = aggregate(counts2[anchor != '', ], list(anchor[anchor != '']), mean)
      profiles = t(mean_mat[,2:ncol(mean_mat)])
      colnames(profiles) = mean_mat[,1]
      sup1 <- insitutypeML(x = counts2,neg = rowMeans(neg2),cohort = cohort2,reference_profiles = profiles)
      cd34[names(sup1$clust), 'celltype'] = sup1$clust
      cd34[names(anchor), 'anchor'] = anchor
    }
  }
  cd34$celltype_detail = cd34$celltype
  cd34$anchor_detail = cd34$anchor
  cd34$celltype = ct[rownames(cd34), "celltype"]
  cd34$anchor = ct[rownames(cd34), "anchor"]
  write.csv(cd34[cd34$celltype != 'nonmalignant', ], file= paste0('p', pa, '_celltype_detail_v5.csv'))
}
