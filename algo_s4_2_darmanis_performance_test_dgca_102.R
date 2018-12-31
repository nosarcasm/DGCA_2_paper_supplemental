rm(list = ls(all.names = TRUE))

.libPaths(new="~/.Rlib")
library(profvis)
library(htmlwidgets)

.libPaths(new="~/.RlibDGCAold")
library(DGCA)

set.seed(12345)

data(darmanis)
data(design_mat)

filtmat = filterGenes(darmanis, filterTypes =  c("central", "dispersion"), keepRows = NULL,
                        filterCentralType = "median", filterDispersionType = "dispersion_index",
                        filterCentralPercentile = 0.25, filterDispersionPercentile = 0.25,
                        sequential = FALSE, allGroups = FALSE, design = design_mat)
datExpr = filtmat
nPairs = (nrow(datExpr)*nrow(datExpr)-nrow(datExpr))/2
perms = c(14,rep(15,10),rep(50,10))
minPs = c()
genesUnder = c()
cutoff = 0.05
prof = profvis({
  for(p in perms){
    print(Sys.time())
    cat(paste0("old, cores=1, perm number: ",p,"\n"))
    set.seed(12345)
    ddcor_res = ddcorAll(nPerms = p, nPairs = nPairs, inputMat = filtmat, design = design_mat,
                       compare = c("oligodendrocyte", "neuron"),corrType = "spearman")
    value = min(ddcor_res$pValDiff_adj)
    gunder = length(ddcor_res[ddcor_res$pValDiff_adj<cutoff,"pValDiff_adj"])
    minPs = c(minPs,value)
    genesUnder = c(genesUnder,gunder)
    cat(paste0("min p: ",value,"\n"))
    cat(paste0("genes under ",cutoff," cutoff: ",gunder,"\n"))
  }
})
saveWidget(prof,paste0("dgca_perm_split_profile_all_perms_nopar_linux_old.html"))

#rm(list=ls())
