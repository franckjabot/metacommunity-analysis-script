permutation.based.pathanalysis = function(model, data, nperm = 1000, verb = F)
{
  require(MASS) ; require(lavaan)
  
  npop = dim(data[[1]])[1]
  nvar = length(data)
  data_vec = NULL
  
  for(i in 1:nvar){
    temp = as.vector(as.dist(data[[i]]))
    data_vec = cbind(data_vec, temp)
  }
  colnames(data_vec) = names(data)
  
  fit = sem(model, data = data_vec, warn = F,std.ov=T)
  coeff = standardizedSolution(fit)$est.std[standardizedSolution(fit)$op == '~']
  
  coeff_perm = NULL
  pval_coef_perm = NULL
  
  data_perm = data_vec
  for (i in 1:nperm){
    check_conv = 0
    while (check_conv == 0){
      for (j in 1:nvar){
        rarray = sample(npop)
        data_perm[,j] = as.vector(as.dist(data[[j]][rarray,rarray]))
      }
      fit_perm = sem(model, data = data_perm, warn = F,std.ov=T)
      check_conv = round(sum(inspect(fit_perm, 'rsquare')), digit = 3)
    }
    coeff_perm = cbind(coeff_perm, standardizedSolution(fit_perm)$est.std[standardizedSolution(fit_perm)$op == '~'])
    if(verb == T){print(i)}
  }
  
  for (i in 1:length(coeff)){
    if (coeff[i]>=0){
      pval_coef_perm = c(pval_coef_perm, sum(coeff_perm[i,]>coeff[i])/(nperm + 1))
    } else {
      pval_coef_perm = c(pval_coef_perm, sum(coeff_perm[i,]<coeff[i])/(nperm + 1))
    }
  }
  coeffs_class = standardizedSolution(fit)[standardizedSolution(fit)$op == '~',]
  coeffs_class[,1] = paste(coeffs_class[,1], coeffs_class[,2], coeffs_class[,3])
  coeffs_class = coeffs_class[,-c(2,3,5:9)]
  coeffs = cbind(coeffs_class, pval_coef_perm)
  names(coeffs) = c('Path', 'Standardized estimate', 'Pvalue')
  return(coeffs)
}

 
