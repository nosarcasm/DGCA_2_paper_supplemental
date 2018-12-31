rm(list=ls())
set.seed(12345)
library(DGCA)
library(magrittr)
library(eclust)
library(factoextra)
library(pROC)

# simulation parameters
rho = 0.90; p = 600 ;SNR = 1 ; n = 400; n0 = n1 = 200 ; nActive = p*0.10 ; cluster_distance = "tom";
Ecluster_distance = "difftom"; rhoOther = 0.6; betaMean = 2;
alphaMean = 1; betaE = 3; distanceMethod = "euclidean"; clustMethod = "hclust";
cutMethod = "dynamic"; agglomerationMethod = "average"

#in this simulation its blocks 3 and 4 that are important
#leaveOut:  optional specification of modules that should be left out
#of the simulation, that is their genes will be simulated as unrelated
#("grey"). This can be useful when simulating several sets, in some which a module
#is present while in others it is absent.

#' @param n number of observations
#' @param p total number of predictors to simulate
#' @param exposed binary numeric vector of length \code{n} with 0 for unexposed
#'   and 1 for exposed
#' @param rho numeric value representing the expected correlation between green
#'   module and red module
#' @param ... arguments passed to the \code{\link[WGCNA]{simulateDatExpr}} function
#' @return \code{n x p} matrix of simulated data
#' 
s_modules <- function(n, p, rho, exposed, ...) {
    #Step 1: simulate the seed module eigengenes
    sMEturquoise <- stats::rnorm(n) #module 1
    
    #expected cor(sMEblue,sMEturquoise) = 0.60
    sMEblue <- sin((sMEturquoise)*pi/1.5+pi/3)*0.99+sqrt(1 - 0.99 ^ 2) * stats::rnorm(100) #module 2
    
    sMEblue<<-sMEblue
    sMEturquoise <<-sMEturquoise
    sMEyellow <- stats::rnorm(n) #module3
    
    sMEgreen <- stats::rnorm(n) #module 4
    
    #expected cor(e.continuous,seed.ME)=0.95
    temp0 <- rho[1] * sMEgreen + sqrt(1 - rho[1] ^ 2) * stats::rnorm(n)
    
    #expected cor(y.continuous,seed.ME) <- -0.95
    sMEred <- rho[1] * temp0 + sqrt(1 - rho[1] ^ 2) * stats::rnorm(n) #module 5
    
    datsME <- data.frame(sMEturquoise,sMEblue,sMEred,sMEgreen,sMEyellow) #module 6
    
    dat1 <- WGCNA::simulateDatExpr(eigengenes = datsME, nGenes = p, ...) #expression
    
    return(dat1)
}

d0 <- s_modules(n = n1, p = p, rho = rho, exposed = TRUE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.5,
                maxCor = 1,
                corPower = 0.5,
                propNegativeCor = 0.3,
                backgroundNoise = 0.2,
                signed = T)

d1 <- s_modules(n = n1, p = p, rho = rho, exposed = TRUE,
                modProportions = c(0.15,0.15,0.15,0.15,0.15,0.25),
                minCor = 0.4,
                maxCor = 1,
                corPower = 0.3,
                propNegativeCor = 0.3,
                backgroundNoise = 0.5,
                signed = TRUE ,
                leaveOut=c(2))

# module 2 is correlated with the turquoise module non-linearly (30% of data)
# module 4 is correlated to itself only linearly (15% of data)

#gene modules changed in 30% of data

truemodules <- d0$setLabels

X <- rbind(d0$datExpr, d1$datExpr) %>%
  magrittr::set_colnames(paste0("Gene", 1:p)) %>%
  magrittr::set_rownames(paste0("Subject",1:n))

true_pairs = combn(rownames(t(X)[(truemodules==2)|(truemodules==1),]),m=2,simplify=F,FUN=sort) #all 1x1, 1x2, and 2x2
true_pairs = setdiff(true_pairs,combn(rownames(t(X)[(truemodules==2),]),m=2,simplify=F,FUN=sort)) #remove the 2x2 so that we have half non-linear and half linear corrs
true_pairs = setdiff(true_pairs,combn(rownames(t(X)[(truemodules==1),]),m=2,simplify=F,FUN=sort)) #remove the 2x2 so that we have half non-linear and half linear corrs

group1 = rep("cond_a",nrow(d0$datExpr))
group2 = rep("cond_b",nrow(d0$datExpr))
conditions = c(group1,group2)
design_mat = model.matrix(~ conditions + 0)
labels = c("cond_a", "cond_b")
colnames(design_mat) = labels

#=========
print("Using bicor...")
ddcor_start = Sys.time()
cl = startCluster(8)
ddcor_res = ddcorAll(inputMat = t(X), design = design_mat, compare = labels, 
                     corr_cutoff = 0.4, adjust = "perm", corrType = "bicor", nPerm = 15, 
                     nPairs = Inf, cl=cl, k=5)
print(paste0("time taken to run the DGCA pipeline with ", n0, " samples..."))
ddcor_time = difftime(Sys.time(), ddcor_start, units = "secs")
print(ddcor_time)

print("Calculating ROC for bicor")
observed = ddcor_res[ddcor_res["pValDiff_adj"] < 0.05,c("Gene1","Gene2")]
observed = observed[order(observed[,"Gene1"],observed[,"Gene2"]),]
observed_p = ddcor_res[ddcor_res["empPVals"] < 0.05,c("pValDiff_adj")]
expected = as.data.frame(matrix(unlist(true_pairs), nrow=length(true_pairs),byrow = TRUE))
colnames(expected) = colnames(observed)

pairs_categories = ddcor_res[,c("Gene1","Gene2","pValDiff_adj")]
pairs_categories["obs_category"] = FALSE
pairs_categories[pairs_categories["pValDiff_adj"] < 0.05,"obs_category"] = TRUE
pairs_categories["exp_category"] = FALSE
for (row in 1:nrow(pairs_categories)){
  gene1 = pairs_categories[row,"Gene1"]
  gene2 = pairs_categories[row,"Gene2"]
  sorted = sort(c(gene1,gene2))
  if (sum(expected[,"Gene1"]==sorted[1]&expected[,"Gene2"]==sorted[2])==1){
    pairs_categories[row,"exp_category"] = TRUE
  }
}
labels1 = unlist(pairs_categories[,"exp_category"])
scores1 = as.numeric(unlist(pairs_categories[,"pValDiff_adj"]))
print("ROC for bicor")
print(roc(labels1, scores1))
roc_bicor = roc(labels1, scores1)
#========

#=========
print("Using MI...")

ddcor_start = Sys.time()
cl = startCluster(8)
ddcor_res = ddcorAll(inputMat = t(X), design = design_mat, compare = labels, 
                     corr_cutoff = 0.4, adjust = "perm", corrType = "mutualinformation", nPerm = 15, 
                     nPairs = Inf, cl=cl, k=7)
print(paste0("time taken to run the DGCA pipeline with ", n0, " samples..."))
ddcor_time = difftime(Sys.time(), ddcor_start, units = "secs")
print(ddcor_time)

print("Calculating ROC for MI")
observed = ddcor_res[ddcor_res["empPVals"] < 0.05,c("Gene1","Gene2")]
observed = observed[order(observed[,"Gene1"],observed[,"Gene2"]),]
observed_p = ddcor_res[ddcor_res["empPVals"] < 0.05,c("empPVals")]
expected = as.data.frame(matrix(unlist(true_pairs), nrow=length(true_pairs),byrow = TRUE))
colnames(expected) = colnames(observed)

pairs_categories = ddcor_res[,c("Gene1","Gene2","empPVals")]
pairs_categories["obs_category"] = FALSE
pairs_categories[pairs_categories["empPVals"] < 0.05,"obs_category"] = TRUE
pairs_categories["exp_category"] = FALSE
for (row in 1:nrow(pairs_categories)){
  gene1 = pairs_categories[row,"Gene1"]
  gene2 = pairs_categories[row,"Gene2"]
  sorted = sort(c(gene1,gene2))
  if (sum(expected[,"Gene1"]==sorted[1]&expected[,"Gene2"]==sorted[2])==1){
    pairs_categories[row,"exp_category"] = TRUE
  }
}
labels2 = unlist(pairs_categories[,"exp_category"])
scores2 = as.numeric(unlist(pairs_categories[,"empPVals"]))
print("ROC for MI")
print(roc(labels2, scores2))
roc_mi = roc(labels2, scores2)
#========

#=========
print("Using pearson...")
ddcor_start = Sys.time()
cl = startCluster(8)
ddcor_res = ddcorAll(inputMat = t(X), design = design_mat, compare = labels, 
                     corr_cutoff = 0.4, adjust = "perm", corrType = "pearson", nPerm = 15, 
                     nPairs = Inf, cl=cl, k=5)
print(paste0("time taken to run the DGCA pipeline with ", n0, " samples..."))
ddcor_time = difftime(Sys.time(), ddcor_start, units = "secs")
print(ddcor_time)

print("Calculating ROC for pearson")
observed = ddcor_res[ddcor_res["pValDiff_adj"] < 0.05,c("Gene1","Gene2")]
observed = observed[order(observed[,"Gene1"],observed[,"Gene2"]),]
observed_p = ddcor_res[ddcor_res["empPVals"] < 0.05,c("pValDiff_adj")]
expected = as.data.frame(matrix(unlist(true_pairs), nrow=length(true_pairs),byrow = TRUE))
colnames(expected) = colnames(observed)

pairs_categories = ddcor_res[,c("Gene1","Gene2","pValDiff_adj")]
pairs_categories["obs_category"] = FALSE
pairs_categories[pairs_categories["pValDiff_adj"] < 0.05,"obs_category"] = TRUE
pairs_categories["exp_category"] = FALSE
for (row in 1:nrow(pairs_categories)){
  gene1 = pairs_categories[row,"Gene1"]
  gene2 = pairs_categories[row,"Gene2"]
  sorted = sort(c(gene1,gene2))
  if (sum(expected[,"Gene1"]==sorted[1]&expected[,"Gene2"]==sorted[2])==1){
    pairs_categories[row,"exp_category"] = TRUE
  }
}
labels3 = unlist(pairs_categories[,"exp_category"])
scores3 = as.numeric(unlist(pairs_categories[,"pValDiff_adj"]))
print("ROC for pearson")
roc_pearson = roc(labels3, scores3)
print(roc_pearson)
#========

#=========
print("Using spearman...")
ddcor_start = Sys.time()
cl = startCluster(8)
ddcor_res = ddcorAll(inputMat = t(X), design = design_mat, compare = labels, 
                     corr_cutoff = 0.4, adjust = "perm", corrType = "spearman", nPerm = 15, 
                     nPairs = Inf, cl=cl, k=5)
print(paste0("time taken to run the DGCA pipeline with ", n0, " samples..."))
ddcor_time = difftime(Sys.time(), ddcor_start, units = "secs")
print(ddcor_time)

print("Calculating ROC for spearman")
observed = ddcor_res[ddcor_res["pValDiff_adj"] < 0.05,c("Gene1","Gene2")]
observed = observed[order(observed[,"Gene1"],observed[,"Gene2"]),]
observed_p = ddcor_res[ddcor_res["empPVals"] < 0.05,c("pValDiff_adj")]
expected = as.data.frame(matrix(unlist(true_pairs), nrow=length(true_pairs),byrow = TRUE))
colnames(expected) = colnames(observed)

pairs_categories = ddcor_res[,c("Gene1","Gene2","pValDiff_adj")]
pairs_categories["obs_category"] = FALSE
pairs_categories[pairs_categories["pValDiff_adj"] < 0.05,"obs_category"] = TRUE
pairs_categories["exp_category"] = FALSE
for (row in 1:nrow(pairs_categories)){
  gene1 = pairs_categories[row,"Gene1"]
  gene2 = pairs_categories[row,"Gene2"]
  sorted = sort(c(gene1,gene2))
  if (sum(expected[,"Gene1"]==sorted[1]&expected[,"Gene2"]==sorted[2])==1){
    pairs_categories[row,"exp_category"] = TRUE
  }
}
labels4 = unlist(pairs_categories[,"exp_category"])
scores4 = as.numeric(unlist(pairs_categories[,"pValDiff_adj"]))
print("ROC for spearman")
roc_spearman = roc(labels4, scores4)
print(roc_spearman)
#========

plot.roc(labels3,scores3,print.auc = T, col = "dodgerblue",print.auc.y = 0.22,lty=1,lwd=3,print.auc.x=0.45, ci=T)
plot.roc(labels4,scores4,add = T,print.auc = T,col = "darkgreen",print.auc.y=0.16,lty=2,lwd=3,print.auc.x=0.45,ci=T)
plot.roc(labels1,scores1,add = T,print.auc = T,col = "firebrick",print.auc.y =0.10,lty=3,lwd=3,print.auc.x=0.45,ci=T)
plot.roc(labels2,scores2,add = T,print.auc = T,col = "violetred",print.auc.y=0.04,lty=4, lwd=3,print.auc.x=0.45,ci=T)
legend(0.35,0.55,legend=c("pearson","spearman","bicor","MI"), col = c("dodgerblue","darkgreen","firebrick","violetred"),lty=1:4,lwd = 1)

#====
rocobj <- plot.roc(labels2,scores2,
                     
                     main="Confidence intervals", percent=TRUE,
                     
                     ci=TRUE, # compute AUC (of AUC by default)
                     
                     print.auc=TRUE) # print the AUC (will contain the CI)

ciobj <- ci.se(rocobj, # CI of sensitivity
               
               specificities=seq(0, 100, 5)) # over a select set of specificities

plot(ciobj, type="shape", col="#1c61b6AA") # plot as a blue shape

plot(ci(rocobj, of="thresholds", thresholds="best")) # add one threshold
