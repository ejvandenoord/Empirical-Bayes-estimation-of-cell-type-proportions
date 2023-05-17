
prepare_data = function(bulk_data,panel_data,covariate_data=NULL) { 
  
  ind  = match(tolower(rownames(panel_data)), tolower(colnames(bulk_data)))
  
  genes_present = 
  print( paste0(sum(!is.na(ind))," of the ",nrow(panel_data)," panel genes are present in the bulk data") )

  bulk_data  = as.matrix( bulk_data[,ind[!is.na(ind)]] )
  panel_data = as.matrix( panel_data[!is.na(ind),]     )
  
  if (!is.null(covariate_data)) {
    print( "Regressing out covariates from the bulk data" )
    
    ind            = match(rownames(bulk_data), rownames(covariate_data))
    covariate_data = as.matrix(covariate_data[ind[!is.na(ind)],])
    bulk_data      = as.matrix(bulk_data[!is.na(ind),])

    dims                = dim( bulk_data )
    residuals           = matrix(NA,dims[1],dims[2])
    rownames(residuals) = rownames(bulk_data)
    colnames(residuals) = colnames(bulk_data)
    for(i in 1:dims[2]) { # i=1
      model  = lm(bulk_data[,i] ~ covariate_data )
      ind    = match(rownames(residuals),names(model$residuals))
      residuals[!is.na(ind),i] = model$residuals[ind[!is.na(ind)]]
    }
    
    bulk_data =  scale(residuals)
  } else {
    bulk_data = scale(bulk_data)
  }
  
  list(bulk_data,panel_data)
}


eb_estimation = function(bulk_data,panel,prior_mean,prior_sd,shrinkSD_factor,n_iter) { 
  
  log_dir = getwd()
  panel   = as.matrix(panel)
  
  n_samples       = nrow(bulk_data)
  n_sites         = ncol(bulk_data)
  n_cellTypes     = ncol(panel)
  
  props           = matrix(NA,n_samples,n_cellTypes)
  colnames(props) = colnames(panel)
  r2              = rep(0,n_samples)
  
  fmla  = bulk ~ -1 + panel
  scale = prior_sd/shrinkSD_factor
  scale[ scale==0] = 0.01
  
  print( "Start estimation: Empirical Bayes")
#  sink(paste0(log_dir,"/R.log"))  
  for(i in 1:n_samples) {  # i=1
#    if (i %% 10==0) {
#          closeAllConnections() 
#          print( paste0("Processing sample ",i," of ",n_samples) )
#          sink(paste0(work_dir,"/logs/R.log"))   
#        }  
    if (i %% 10==0) print( paste0("Processing sample ",i," of ",n_samples) )
    bulk = as.vector( t(bulk_data[i,]) )
    data = data.frame(bulk, panel )        
    
    mod_stan = stan_glm(fmla, data = data,
                        chains = 2, iter = n_iter, 
                        prior = normal(location = prior_mean, scale=scale),
                        refresh = 0)
    props[i,1:n_cellTypes] =  coef(mod_stan)    #  posterior median estimates      
    r2[i] = median( bayes_R2(  mod_stan  ) )
    
  }
#  closeAllConnections() 
  
  props[is.na(props)] = 0
  nprops              = data.frame(props)
  nprops[nprops < 0]  = 0
  estimates           = nprops/rowSums(nprops) # standardize
  rownames(estimates) = rownames(bulk_data)
  colnames(estimates) = colnames(props)
  
  colnames(nprops)    = paste0(colnames(props),"raw")
  nprops$nprops_max1  = apply(nprops[,colnames(nprops)], 1, max) 
  nprops$nprops_min0  = apply(nprops[,colnames(nprops)], 1, min) 
  nprops$nprops_max1  = round(nprops$nprops_max1 == 1) 
  nprops$nprops_min0  = round(nprops$nprops_min0 == 0) 
  
  results             = list()
  results$estimates   = estimates
  results$qc          = data.frame(nprops,r2,shrinkSD_factor) 
  rownames(results$qc) = rownames(bulk_data)
  results
  
  
  results
  
} 

constraint_estimation = function(bulk_data,panel) { 
  
  bulk_data  = as.matrix(bulk_data)
  panel      = as.matrix(panel)
  n_samples  = nrow(bulk_data)
  n_celltype = ncol(panel)
  
  cell_props           = matrix(NA,n_samples,n_celltype)
  colnames(cell_props) = colnames(panel)
  r2                   = rep(0,n_samples)
  
  print( "Start estimation: Use boundary contraint >= 0")
  bol_positive    = rep(T,n_celltype)
  for(i in 1:n_samples) { # i =1 
    if (i %% 10==0) print( paste0("Processing sample ",i," of ",n_samples) )
    mod = penalized(bulk_data[i,], ~ panel, ~-1,lambda1=0, lambda2=0, positive=bol_positive,trace=F)
    r2[i] = cor( fitted.values(mod), bulk_data[i,],use="pairwise.complete.obs")^2
    cell_props[i,1:n_celltype] = coef(mod)[1:n_celltype]           
  }  
  
  cell_props[is.na(cell_props)] = 0
  ncell_props              = data.frame(cell_props)
  ncell_props[ncell_props < 0]  = 0
  estimates           = ncell_props/rowSums(ncell_props) # standardize
  rownames(estimates) = rownames(bulk_data)
  colnames(estimates) = colnames(cell_props)
  
  colnames(ncell_props)    = paste0(colnames(cell_props),"raw")
  ncell_props$ncell_props_max1  = apply(ncell_props[,colnames(ncell_props)], 1, max) 
  ncell_props$ncell_props_min0  = apply(ncell_props[,colnames(ncell_props)], 1, min) 
  ncell_props$ncell_props_max1  = round(ncell_props$ncell_props_max1 == 1) 
  ncell_props$ncell_props_min0  = round(ncell_props$ncell_props_min0 == 0) 
  
  results             = list()
  results$estimates   = estimates
  results$qc          = data.frame(ncell_props,r2) 
  rownames(results$qc) =  rownames(bulk_data)
  results
  
} 

