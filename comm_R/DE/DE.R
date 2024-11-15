






# df_DE=run_DE_limma(exp = exp,df_info = df_col,var_name = "group",comparisons = compare)

run_DE_limma=function(exp,df_info,var_name,comparisons,covariates=NULL,adjust_method="BH"){
  if (!var_name %in% colnames(df_info)) {stop(paste("Variable", var_name, "not found in df_info"))}
  
  comparison_levels=c(word(comparisons,1,sep="_vs_"),word(comparisons,2,sep="_vs_"))
  if(!all(comparison_levels %in% df_info[[var_name]])){stop("comparison level not in data")}
  
  #get df_info for analysis
  df_info_subset <- df_info[df_info[[var_name]] %in% comparison_levels, ]
  df_info_subset[[var_name]] <- factor(df_info_subset[[var_name]], levels = comparison_levels)
  
  #get exp for analysis 
  common_samples = intersect(colnames(exp), rownames(df_info_subset))
  if (length(common_samples) == 0) {stop("No matching samples between expression data and df_info")}
  exp_subset <- exp[, common_samples, drop = FALSE]
  df_info_subset <- df_info_subset[common_samples, , drop = FALSE]
  
  #get formula
  formula_terms <- c("0", var_name)
  if (!is.null(covariates)) {formula_terms <- c(formula_terms, covariates)}
  
  design_formula <- as.formula(paste("~", paste(formula_terms, collapse = " + ")))
  design = model.matrix(design_formula, data = df_info_subset)
  colnames(design) <- sort(levels(df_info_subset[[var_name]]))
  
  #run lmFit
  fit <- lmFit(exp_subset, design)
  
  #get contrast matrix
  contrast_formula <- paste0(comparison_levels[1], " - ", comparison_levels[2])
  contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  
  #run contrast
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  #get DE table
  results <- topTable(fit2, number = Inf, adjust.method = adjust_method) %>% 
    tibble::rownames_to_column(var="gene") %>% 
    mutate(contrast=comparisons)
  results
}