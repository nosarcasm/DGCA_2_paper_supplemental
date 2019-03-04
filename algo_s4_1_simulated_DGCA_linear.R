
library(genefilter)
library(MASS)
library(pROC)
library(docopt)

###############################
# steps to perform in the program
test = TRUE

# hardcode the parameters during a test
if(test){
	n_samples = 100
	seed = 12345

	ddcor_test = TRUE

	compare_ddcor_classes = TRUE
}

#reproducibility
set.seed(seed)

#color palette
#http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette = c("#000000", "darkgrey", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
	"#0072B2", "#CC79A7")

#################################
# functions

#Sets the superdiagonal and subdiagonal of a matrix to the number x.
subsup_diag <- function(matrix, x){
	diag(matrix[,-1]) = x
	diag(matrix[-1,]) = x
	return(matrix)
}

#changes the superdiagonal and subdiagonal of a particular submatrix of a larger matrix
add_sigma_to_matrix <- function(total_matrix, row_start, nrows, x){
	extract_mat = total_matrix[row_start:(row_start + nrows),
		row_start:(row_start + nrows)]
	mat = subsup_diag(extract_mat, x)
	total_matrix[row_start:(row_start + nrows),
		row_start:(row_start + nrows)] = mat
	return(total_matrix)
}

#updates the corrrelation matrix and adds the gene pairs that were modified to a class for subsequent extraction
update_corr_structure <- function(sigma, n_genes, variance, corr_factor, dcor = NULL,
	dcor_specific = NULL){
	total_dc_genes = sigma[["total_dc_genes"]]
	sigma[["sigma_tot"]] = add_sigma_to_matrix(sigma[["sigma_tot"]],
		row_start = total_dc_genes + 1,
		nrows = n_genes, x = variance * corr_factor)
	if(!is.null(dcor)){
		if(dcor == "goc"){
			sigma[["goc_gene_pairs"]] = c(sigma[["goc_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		}
		if(dcor == "loc"){
			sigma[["loc_gene_pairs"]] = c(sigma[["loc_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		}
	}
	if(!is.null(dcor_specific)){
		if(dcor_specific == "+/0"){
			sigma[["pos_none_gene_pairs"]] = c(sigma[["pos_none_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		}
		if(dcor_specific == "+/-"){
			sigma[["pos_neg_gene_pairs"]] = c(sigma[["pos_neg_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		}
		if(dcor_specific == "0/+"){
			sigma[["none_pos_gene_pairs"]] = c(sigma[["none_pos_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		}
		if(dcor_specific == "0/-"){
			sigma[["none_neg_gene_pairs"]] = c(sigma[["none_neg_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		}
		if(dcor_specific == "-/+"){
			sigma[["neg_pos_gene_pairs"]] = c(sigma[["neg_pos_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		}
		if(dcor_specific == "-/0"){
			sigma[["neg_none_gene_pairs"]] = c(sigma[["neg_none_gene_pairs"]],
				paste(letters_unique[(total_dc_genes + 1):(total_dc_genes + n_genes - 1)],
				letters_unique[(total_dc_genes + 2):(total_dc_genes + n_genes)], sep = " "))
		}
	}
	sigma[["total_dc_genes"]] = total_dc_genes + n_genes
	return(sigma)
}

#for gene name generation, i.e. "aa", "ab", "ac", etc
combine <- function (vecs, number){
  combn(vecs, number, paste, collapse = "")
}

letters_unique = combine(c(letters[1:26], letters[1:26]), 2)
letters_unique = sort(unique(letters_unique))

#################################
# parameters for simulation data

n_non_expr_genes = 300
n_housekeeping_genes = 100
n_activated_genes = 300
n_tot_genes = n_non_expr_genes + n_housekeeping_genes + n_activated_genes

low_mean = 2
high_mean = 50
low_var = 50
high_var = 100
positive_correlation = 0.5
negative_correlation = -0.5 #if this is set too low, then the matrix will be non-positive definite

######################################
#create overall correlation matrices and classes to store them
sigma_tot_a <- structure(list(
	sigma_tot = matrix(0, nrow = n_tot_genes, ncol = n_tot_genes),
	total_dc_genes = 0
	), class = "sigma")

sigma_tot_b <- structure(list(
	sigma_tot = matrix(0, nrow = n_tot_genes, ncol = n_tot_genes),
	total_dc_genes = 0,
	loc_gene_pairs = vector(),
	goc_gene_pairs = vector(),
	pos_none_gene_pairs = vector(),
	pos_neg_gene_pairs = vector(),
	none_pos_gene_pairs = vector(),
	none_neg_gene_pairs = vector(),
	neg_pos_gene_pairs = vector(),
	neg_none_gene_pairs = vector()
	), class = "sigma")

#A+, B+ = no differential correlation
n_apbp = 20
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_apbp,
	variance = high_var, corr_factor = positive_correlation, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_apbp,
	variance = high_var, corr_factor = positive_correlation, dcor = NULL)

#A+, B= = loss of correlation
n_apbe = 20
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_apbe,
	variance = high_var, corr_factor = positive_correlation, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_apbe,
	variance = high_var, corr_factor = 0, dcor = "loc", dcor_specific = "+/0")

#A+, B- = loss of correlation
n_apbm = 20
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_apbm,
	variance = high_var, corr_factor = positive_correlation, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_apbm,
	variance = high_var, corr_factor = negative_correlation, dcor = "loc",
	dcor_specific = "+/-")

#A=, B+ = gain of correlation
n_aebp = 20
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_aebp,
	variance = high_var, corr_factor = 0, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_aebp,
	variance = high_var, corr_factor = positive_correlation, dcor = "goc",
	dcor_specific = "0/+")

#A=, B- = loss of correlation
n_aebm = 20
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_aebm,
	variance = high_var, corr_factor = 0, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_aebm,
	variance = high_var, corr_factor = negative_correlation, dcor = "loc",
	dcor_specific = "0/-")

#A-, B+ = gain of correlation
n_ambp = 20
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_aebm,
	variance = high_var, corr_factor = negative_correlation, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_aebm,
	variance = high_var, corr_factor = positive_correlation, dcor = "goc",
	dcor_specific = "-/+")

#A-, B= = gain of correlation
n_ambe = 20
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_aebm,
	variance = high_var, corr_factor = negative_correlation, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_aebm,
	variance = high_var, corr_factor = 0, dcor = "goc",
	dcor_specific = "-/0")

#A-, B- = no differential correlation
n_ambm = 20
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_aebm,
	variance = high_var, corr_factor = negative_correlation, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_aebm,
	variance = high_var, corr_factor = negative_correlation, dcor = NULL)

#A=, B= = no differential correlation
n_aebe = n_activated_genes - sigma_tot_a[["total_dc_genes"]]
sigma_tot_a = update_corr_structure(sigma = sigma_tot_a, n_genes = n_aebm,
	variance = high_var, corr_factor = 0, dcor = NULL)
sigma_tot_b = update_corr_structure(sigma = sigma_tot_b, n_genes = n_aebm,
	variance = high_var, corr_factor = 0, dcor = NULL)

sigma_tot_a_mat = sigma_tot_a[["sigma_tot"]]
sigma_tot_b_mat = sigma_tot_b[["sigma_tot"]]

########################################
#set the means and variances of the genes

#set the variance
var_diag = c(rep(high_var, n_activated_genes),
	rep(low_var, n_housekeeping_genes),
	rep(high_var, n_non_expr_genes))

#set the mean for the high expression genes
mu = rnbinom(n = (n_activated_genes + n_housekeeping_genes), mu = high_mean, size = 0.5)

#set the mean for the low expression genes
mu = c(mu, rnbinom(n = n_non_expr_genes, mu = low_mean, size = 0.5))

#if mu is less than one, set it to one
mu[mu < 1] = 1

diag(sigma_tot_a_mat) = var_diag
diag(sigma_tot_b_mat) = var_diag

#generate the simulation matrix
sim_data_a = as.matrix(mvrnorm(n_samples, mu = mu, Sigma = sigma_tot_a_mat))
sim_data_b = as.matrix(mvrnorm(n_samples, mu = mu, Sigma = sigma_tot_b_mat))

true_cases = c(sigma_tot_b$loc_gene_pairs, sigma_tot_b$goc_gene_pairs)

####################################
#run the ddcor differential correlation pipeline on the simulated data

if(ddcor_test){

library(DGCA)

print("starting the ddcor pipeline")

sim_data_a = t(sim_data_a)
sim_data_b = t(sim_data_b)

sim_data_merge = cbind(sim_data_a, sim_data_b)
rownames(sim_data_merge) = letters_unique[1:nrow(sim_data_merge)]

conditions = c(rep("cond_a", ncol(sim_data_a)),
  rep("cond_b", ncol(sim_data_b)))
design_mat = model.matrix(~ conditions + 0)
labels = c("cond_a", "cond_b")
colnames(design_mat) = labels
npairs = (nrow(sim_data_merge)^2)/2 - nrow(sim_data_merge)

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
	ddcor_res$true = ddcor_res$combos %in% true_cases

	plot_groupwise_ROC_ddcor <- function(res_pval = ddcor_res$empPVals,
		pair_names = ddcor_res$combos, group_true, all_true = true_cases){

		group_true_vec = pair_names %in% group_true
		to_remove = pair_names %in% setdiff(all_true, group_true)
		tmp_df = data.frame(pvals = res_pval, true = group_true_vec)
		tmp_df = tmp_df[!to_remove, ]
		roc_list = roc(tmp_df$true, tmp_df$pvals, smooth = FALSE)
		return(roc_list)

	}

	print("starting ddcor ROC curve calculation")
	true_cases = c(sigma_tot_b$loc_gene_pairs, sigma_tot_b$goc_gene_pairs)
	full_roc_p = plot_groupwise_ROC_ddcor(res_pval = ddcor_res$pValDiff, group_true = true_cases)
	full_roc_emp = plot_groupwise_ROC_ddcor(res_pval = ddcor_res$empPVals, group_true = true_cases)
	full_roc_q_emp = plot_groupwise_ROC_ddcor(res_pval = ddcor_res$pValDiff_adj, group_true = true_cases)

	print("full AUC for ddcor with n_samples...")
	print(n_samples)
	print("p,emp p, emp q")
	print(full_roc_p)
	print(full_roc_emp)
	print(full_roc_q_emp)

	ddcor_line = c("ddcor", n_samples, seed, ddcor_time, as.numeric(full_roc_p$auc))
	write.table(ddcor_line, paste0(corrType,"_auc_results_p_mc.txt"), sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE,
		row.names = FALSE)

	ddcor_line = c("ddcor", n_samples, seed, ddcor_time, as.numeric(full_roc_emp$auc))
	write.table(ddcor_line, paste0(corrType,"_auc_results_perm_mc.txt"), sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE,
		row.names = FALSE)

	ddcor_line = c("ddcor", n_samples, seed, ddcor_time, as.numeric(full_roc_q_emp$auc))
	write.table(ddcor_line, paste0(corrType,"_auc_results_q_emp_mc.txt"), sep = "\t", quote = FALSE, append = TRUE, col.names = FALSE,
		row.names = FALSE)

	if(compare_ddcor_classes){
	  roc_loc_only = plot_groupwise_ROC_ddcor(res_pval = ddcor_res$pValDiff,group_true = sigma_tot_b$loc_gene_pairs)
	  roc_goc_only = plot_groupwise_ROC_ddcor(res_pval = ddcor_res$pValDiff,group_true = sigma_tot_b$goc_gene_pairs)
	  roc_p0_only = plot_groupwise_ROC_ddcor(res_pval = ddcor_res$pValDiff,group_true = sigma_tot_b$pos_none_gene_pairs)
	  roc_pn_only = plot_groupwise_ROC_ddcor(res_pval = ddcor_res$pValDiff,group_true = sigma_tot_b$pos_neg_gene_pairs)
	  roc_0p_only = plot_groupwise_ROC_ddcor(res_pval = ddcor_res$pValDiff,group_true = sigma_tot_b$none_pos_gene_pairs)
	  roc_0n_only = plot_groupwise_ROC_ddcor(res_pval = ddcor_res$pValDiff,group_true = sigma_tot_b$none_neg_gene_pairs)
	  roc_np_only = plot_groupwise_ROC_ddcor(res_pval = ddcor_res$pValDiff,group_true = sigma_tot_b$neg_pos_gene_pairs)
	  roc_n0_only = plot_groupwise_ROC_ddcor(res_pval = ddcor_res$pValDiff,group_true = sigma_tot_b$neg_none_gene_pairs)
	  
		#pdf(file = paste0("ddcor_classes_multcore_qemp_", n_samples, "_samples", seed, ".pdf"),
		#		width = 4, height = 4)
		plot(roc_goc_only, col = cbPalette[1])
		plot(roc_loc_only, add = TRUE, col = cbPalette[2])
		plot(roc_p0_only, add = TRUE, col = cbPalette[3])
		plot(roc_pn_only, add = TRUE, col = cbPalette[4])
		plot(roc_0p_only, add = TRUE, col = cbPalette[5])
		plot(roc_0n_only, add = TRUE, col = cbPalette[6])
		plot(roc_np_only, add = TRUE, col = cbPalette[7])
		plot(roc_n0_only, add = TRUE, col = cbPalette[8])
	  legend(0.2,0.5,legend=c("GOC","LOC","+/0","+/-","0/+","0/-","-/+","-/0"),lty=1,col = cbPalette[1:8],lwd=4)
		#dev.off()
	}
}

}



