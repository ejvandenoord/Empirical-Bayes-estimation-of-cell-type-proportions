rm(list=ls())
work_dir        = "G:\\My Drive\\SOP\\snRNAseq_brain\\paper\\panel\\github\\"

panel_data      = read.csv(paste0(work_dir,'panel_scaled.csv'),  stringsAsFactors=F,row.names=1,check.names=FALSE)
bulk_data       = read.csv(paste0(work_dir,'bulk_data.csv'),     stringsAsFactors=F,row.names=1,check.names=FALSE) # subject names in rows gene names in columns
covariate_data  = read.csv(paste0(work_dir,'covariate_data.csv'),stringsAsFactors=F,row.names=1,check.names=FALSE) # optional file: subject names in rows covariate names in columns

shrinkSD_factor = 1    # SD / shrink_factor give SD of the prior
n_iter          = 5000 # Number of iterations in Empirical Bayes estimation

# begin
source(paste0(work_dir,"functions.R"))
library(penalized)
library(rstanarm) 
library(openxlsx)
results_list = list()

# prepare the (bulk) data (covariate_data is optional: if passed to the function the covariates will be regressed out the bulk data)
temp       = prepare_data(bulk_data,panel_data) # covariates will NOT be regressed out the bulk data
#temp       = prepare_data(bulk_data,panel_data,covariate_data) # covariates will be regressed out the bulk data
bulk_data  = temp[[1]]
panel_data = temp[[2]]

temp       = constraint_estimation(bulk_data,panel_data) 
estimates  = temp[[1]]
qc         = temp[[2]]
results_list[[length(results_list)+1]] = estimates 
results_list[[length(results_list)+1]] = qc 

prior_mean = apply( estimates, 2, mean, na.rm=T)
prior_sd   = apply( estimates, 2, sd, na.rm=T)

temp       = eb_estimation(bulk_data,panel_data,prior_mean,prior_sd,shrinkSD_factor,n_iter)
estimates  = temp[[1]]
qc         = temp[[2]]
results_list[[length(results_list)+1]] = estimates 
results_list[[length(results_list)+1]] = qc 

names(results_list) = c('estimates_constraint','qc_constraint','estimates_eb','qc_eb')
write.xlsx(results_list,paste0(work_dir,'estimation_results.xlsx'),rowNames=T)




