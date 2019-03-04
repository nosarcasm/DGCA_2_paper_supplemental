rm(list=ls())
set.seed(12345)
library(DGCA)
library(magrittr)
library(eclust)
library(factoextra)
library(pROC)

# simulation parameters
samples = 200
total_genes = 1000

maxCor_true = 1
minCor_true = 0.3
maxCor_null = 0.1

## classes:
class1 = "+/+"
class2 = "+/0"
class3 = "+/-"
class4 = "0/+"
class5 = "0/0"
class6 = "0/-"
class7 = "-/+"
class8 = "-/0"
class9 = "-/-"

int <- function(n){as.integer(n)}

pos_cor1 = simulateModule(stats::rnorm(samples),int(total_genes/9),nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=T,geneMeans=NULL) #module is positively correlated
pos_cor2 = simulateModule(stats::rnorm(samples),int(total_genes/9),nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=T,geneMeans=NULL) #module is positively correlated
pos_cor3 = simulateModule(stats::rnorm(samples),int(total_genes/9),nNearGenes=0,minCor=minCor_true,maxCor=maxCor_true,corPower=1,signed=T,geneMeans=NULL) #module is positively correlated
neg_cor1 = simulateModule(stats::rnorm(samples),int(total_genes/9),nNearGenes=0,minCor=-minCor_true,maxCor=-maxCor_true,corPower=1,signed=T,geneMeans=NULL) #module is negatively correlated
neg_cor2 = simulateModule(stats::rnorm(samples),int(total_genes/9),nNearGenes=0,minCor=-minCor_true,maxCor=-maxCor_true,corPower=1,signed=T,geneMeans=NULL) #module is negatively correlated
neg_cor3 = simulateModule(stats::rnorm(samples),int(total_genes/9),nNearGenes=0,minCor=-minCor_true,maxCor=-maxCor_true,corPower=1,signed=T,geneMeans=NULL) #module is negatively correlated
null_cor1 = simulateModule(stats::rnorm(samples),int(total_genes/9),nNearGenes=0,minCor=-maxCor_null,maxCor=maxCor_null,corPower=1,signed=T,geneMeans=NULL) #no correlation
null_cor2 = simulateModule(stats::rnorm(samples),int(total_genes/9),nNearGenes=0,minCor=-maxCor_null,maxCor=maxCor_null,corPower=1,signed=T,geneMeans=NULL) #no correlation
null_cor3 = simulateModule(stats::rnorm(samples),int(total_genes/9),nNearGenes=0,minCor=-maxCor_null,maxCor=maxCor_null,corPower=1,signed=T,geneMeans=NULL) #no correlation

true_pos_1 = attributes(pos_cor1)$trueKME
true_pos_2 = attributes(pos_cor2)$trueKME
true_pos_3 = attributes(pos_cor3)$trueKME
true_neg_1 = attributes(neg_cor1)$trueKME
true_neg_2 = attributes(neg_cor2)$trueKME
true_neg_2 = attributes(neg_cor3)$trueKME
true_null_1 = attributes(null_cor1)$trueKME
true_null_2 = attributes(null_cor2)$trueKME
true_null_3 = attributes(null_cor3)$trueKME

order_control = c(pos_cor1,pos_cor2,pos_cor3,
                  neg_cor1,neg_cor2,neg_cor3,
                  null_cor1,null_cor2, null_cor3)

order_disease = c(pos_cor1,null_cor1,neg_cor1,
                  pos_cor2,null_cor2,neg_cor2,
                  pos_cor3,null_cor3,neg_cor3)

#======
#TODO: concat the matrices for control and disease, make matrices with information about the correlation before and after, etc.
#======

conditions = c(group1,group2)
design_mat = model.matrix(~ conditions + 0)
labels = c("cond_a", "cond_b")
colnames(design_mat) = labels

for(corrType in c("pearson","spearman","bicor","MI")){
  print(corrType)
  ddcor_start = Sys.time()
  
  ddcor_res = ddcorAll(inputMat = sim_data_merge, design = design_mat, compare = labels, 
                       corr_cutoff = 0.4, adjust = "perm", corrType = "spearman", nPerm = 15, 
                       nPairs = "all")
  
  print(paste0("time taken to run the DGCA pipeline with ", n_samples, " samples..."))
  ddcor_time = difftime(Sys.time(), ddcor_start, units = "secs")
  print(ddcor_time)
  
  ddcor_res$combos = paste(ddcor_res$Gene1, ddcor_res$Gene2, sep = " ")
  #ddcor_res$true = ddcor_res$combos %in% true_cases
}

#======
#TODO: make new ROC plotting code for these results
#======

print("finish the TODOs, Ryan.")

# print("Calculating ROC for bicor")
# observed = ddcor_res[ddcor_res["pValDiff_adj"] < 0.05,c("Gene1","Gene2")]
# observed = observed[order(observed[,"Gene1"],observed[,"Gene2"]),]
# observed_p = ddcor_res[ddcor_res["empPVals"] < 0.05,c("pValDiff_adj")]
# expected = as.data.frame(matrix(unlist(true_pairs), nrow=length(true_pairs),byrow = TRUE))
# colnames(expected) = colnames(observed)
# 
# pairs_categories = ddcor_res[,c("Gene1","Gene2","pValDiff_adj")]
# pairs_categories["obs_category"] = FALSE
# pairs_categories[pairs_categories["pValDiff_adj"] < 0.05,"obs_category"] = TRUE
# pairs_categories["exp_category"] = FALSE
# for (row in 1:nrow(pairs_categories)){
#   gene1 = pairs_categories[row,"Gene1"]
#   gene2 = pairs_categories[row,"Gene2"]
#   sorted = sort(c(gene1,gene2))
#   if (sum(expected[,"Gene1"]==sorted[1]&expected[,"Gene2"]==sorted[2])==1){
#     pairs_categories[row,"exp_category"] = TRUE
#   }
# }
# labels1 = unlist(pairs_categories[,"exp_category"])
# scores1 = as.numeric(unlist(pairs_categories[,"pValDiff_adj"]))
# print("ROC for bicor")
# print(roc(labels1, scores1))
# roc_bicor = roc(labels1, scores1)
# #========
# 
# #=========
# print("Using MI...")
# 
# ddcor_start = Sys.time()
# cl = startCluster(8)
# ddcor_res = ddcorAll(inputMat = t(X), design = design_mat, compare = labels, 
#                      corr_cutoff = 0.4, adjust = "perm", corrType = "mutualinformation", nPerm = 15, 
#                      nPairs = Inf, cl=cl, k=7)
# print(paste0("time taken to run the DGCA pipeline with ", n0, " samples..."))
# ddcor_time = difftime(Sys.time(), ddcor_start, units = "secs")
# print(ddcor_time)
# 
# print("Calculating ROC for MI")
# observed = ddcor_res[ddcor_res["empPVals"] < 0.05,c("Gene1","Gene2")]
# observed = observed[order(observed[,"Gene1"],observed[,"Gene2"]),]
# observed_p = ddcor_res[ddcor_res["empPVals"] < 0.05,c("empPVals")]
# expected = as.data.frame(matrix(unlist(true_pairs), nrow=length(true_pairs),byrow = TRUE))
# colnames(expected) = colnames(observed)
# 
# pairs_categories = ddcor_res[,c("Gene1","Gene2","empPVals")]
# pairs_categories["obs_category"] = FALSE
# pairs_categories[pairs_categories["empPVals"] < 0.05,"obs_category"] = TRUE
# pairs_categories["exp_category"] = FALSE
# for (row in 1:nrow(pairs_categories)){
#   gene1 = pairs_categories[row,"Gene1"]
#   gene2 = pairs_categories[row,"Gene2"]
#   sorted = sort(c(gene1,gene2))
#   if (sum(expected[,"Gene1"]==sorted[1]&expected[,"Gene2"]==sorted[2])==1){
#     pairs_categories[row,"exp_category"] = TRUE
#   }
# }
# labels2 = unlist(pairs_categories[,"exp_category"])
# scores2 = as.numeric(unlist(pairs_categories[,"empPVals"]))
# print("ROC for MI")
# print(roc(labels2, scores2))
# roc_mi = roc(labels2, scores2)
# #========
# 
# #=========
# print("Using pearson...")
# ddcor_start = Sys.time()
# cl = startCluster(8)
# ddcor_res = ddcorAll(inputMat = t(X), design = design_mat, compare = labels, 
#                      corr_cutoff = 0.4, adjust = "perm", corrType = "pearson", nPerm = 15, 
#                      nPairs = Inf, cl=cl, k=5)
# print(paste0("time taken to run the DGCA pipeline with ", n0, " samples..."))
# ddcor_time = difftime(Sys.time(), ddcor_start, units = "secs")
# print(ddcor_time)
# 
# print("Calculating ROC for pearson")
# observed = ddcor_res[ddcor_res["pValDiff_adj"] < 0.05,c("Gene1","Gene2")]
# observed = observed[order(observed[,"Gene1"],observed[,"Gene2"]),]
# observed_p = ddcor_res[ddcor_res["empPVals"] < 0.05,c("pValDiff_adj")]
# expected = as.data.frame(matrix(unlist(true_pairs), nrow=length(true_pairs),byrow = TRUE))
# colnames(expected) = colnames(observed)
# 
# pairs_categories = ddcor_res[,c("Gene1","Gene2","pValDiff_adj")]
# pairs_categories["obs_category"] = FALSE
# pairs_categories[pairs_categories["pValDiff_adj"] < 0.05,"obs_category"] = TRUE
# pairs_categories["exp_category"] = FALSE
# for (row in 1:nrow(pairs_categories)){
#   gene1 = pairs_categories[row,"Gene1"]
#   gene2 = pairs_categories[row,"Gene2"]
#   sorted = sort(c(gene1,gene2))
#   if (sum(expected[,"Gene1"]==sorted[1]&expected[,"Gene2"]==sorted[2])==1){
#     pairs_categories[row,"exp_category"] = TRUE
#   }
# }
# labels3 = unlist(pairs_categories[,"exp_category"])
# scores3 = as.numeric(unlist(pairs_categories[,"pValDiff_adj"]))
# print("ROC for pearson")
# roc_pearson = roc(labels3, scores3)
# print(roc_pearson)
# #========
# 
# #=========
# print("Using spearman...")
# ddcor_start = Sys.time()
# cl = startCluster(8)
# ddcor_res = ddcorAll(inputMat = t(X), design = design_mat, compare = labels, 
#                      corr_cutoff = 0.4, adjust = "perm", corrType = "spearman", nPerm = 15, 
#                      nPairs = Inf, cl=cl, k=5)
# print(paste0("time taken to run the DGCA pipeline with ", n0, " samples..."))
# ddcor_time = difftime(Sys.time(), ddcor_start, units = "secs")
# print(ddcor_time)
# 
# print("Calculating ROC for spearman")
# observed = ddcor_res[ddcor_res["pValDiff_adj"] < 0.05,c("Gene1","Gene2")]
# observed = observed[order(observed[,"Gene1"],observed[,"Gene2"]),]
# observed_p = ddcor_res[ddcor_res["empPVals"] < 0.05,c("pValDiff_adj")]
# expected = as.data.frame(matrix(unlist(true_pairs), nrow=length(true_pairs),byrow = TRUE))
# colnames(expected) = colnames(observed)
# 
# pairs_categories = ddcor_res[,c("Gene1","Gene2","pValDiff_adj")]
# pairs_categories["obs_category"] = FALSE
# pairs_categories[pairs_categories["pValDiff_adj"] < 0.05,"obs_category"] = TRUE
# pairs_categories["exp_category"] = FALSE
# for (row in 1:nrow(pairs_categories)){
#   gene1 = pairs_categories[row,"Gene1"]
#   gene2 = pairs_categories[row,"Gene2"]
#   sorted = sort(c(gene1,gene2))
#   if (sum(expected[,"Gene1"]==sorted[1]&expected[,"Gene2"]==sorted[2])==1){
#     pairs_categories[row,"exp_category"] = TRUE
#   }
# }
# labels4 = unlist(pairs_categories[,"exp_category"])
# scores4 = as.numeric(unlist(pairs_categories[,"pValDiff_adj"]))
# print("ROC for spearman")
# roc_spearman = roc(labels4, scores4)
# print(roc_spearman)
# #========
# 
# plot.roc(labels3,scores3,print.auc = T, col = "dodgerblue",print.auc.y = 0.22,lty=1,lwd=3,print.auc.x=0.45, ci=T)
# plot.roc(labels4,scores4,add = T,print.auc = T,col = "darkgreen",print.auc.y=0.16,lty=2,lwd=3,print.auc.x=0.45,ci=T)
# plot.roc(labels1,scores1,add = T,print.auc = T,col = "firebrick",print.auc.y =0.10,lty=3,lwd=3,print.auc.x=0.45,ci=T)
# plot.roc(labels2,scores2,add = T,print.auc = T,col = "violetred",print.auc.y=0.04,lty=4, lwd=3,print.auc.x=0.45,ci=T)
# legend(0.35,0.55,legend=c("pearson","spearman","bicor","MI"), col = c("dodgerblue","darkgreen","firebrick","violetred"),lty=1:4,lwd = 1)
# 
# #====
# rocobj <- plot.roc(labels2,scores2,
#                      
#                      main="Confidence intervals", percent=TRUE,
#                      
#                      ci=TRUE, # compute AUC (of AUC by default)
#                      
#                      print.auc=TRUE) # print the AUC (will contain the CI)
# 
# ciobj <- ci.se(rocobj, # CI of sensitivity
#                
#                specificities=seq(0, 100, 5)) # over a select set of specificities
# 
# plot(ciobj, type="shape", col="#1c61b6AA") # plot as a blue shape
# 
# plot(ci(rocobj, of="thresholds", thresholds="best")) # add one threshold
