my_TSI <- function (data, os_names, score_types, se_names = NULL, metrics = NULL, 
          mean = NULL, var_ts = NULL, reliability = NULL, linked_obs_cor = NULL, 
          truncate = NULL, separated = rep(T, length(os_names)), ts_names = paste0("true_", 
                                                                                   os_names), mice_args) 
{
  p_to_impute = length(os_names)
  if (p_to_impute == 0) 
    stop("\nProvide the name of at least one observed score to impute as os_names.")
  args_to_test = setNames(list(list(se_names), list(metrics), 
                               list(score_types), list(mean), list(var_ts), list(separated)), 
                          c("se_names", "metrics", "score_types", "mean", "var_ts", 
                            "separated"))
  types_to_test = c("character", "character", "character", 
                    "numeric", "numeric", "logical")
  for (i in seq_len(length(args_to_test))) {
    if (!class(args_to_test[[i]][[1]]) %in% c("NULL", types_to_test[i])) 
      stop(paste0("\nExpected ", names(args_to_test)[i], 
                  " to be of type ", types_to_test[i], " or NULL; was ", 
                  class(args_to_test[[i]][[1]]), "."))
  }
  for (i in seq_len(length(args_to_test))) {
    if (!length(args_to_test[[i]][[1]]) %in% c(0, 1, p_to_impute)) 
      stop(paste0("\nThe number of elements in ", names(args_to_test)[i], 
                  " should be 1 or match the number of elements in os_names."))
  }
  if (length(metrics) == 1) 
    metrics = rep(metrics, p_to_impute)
  if (length(score_types) == 1) 
    score_types = rep(score_types, p_to_impute)
  if (length(mean) == 1) 
    mean = rep(mean, p_to_impute)
  if (length(var_ts) == 1) 
    var_ts = rep(var_ts, p_to_impute)
  if (length(separated) == 1) 
    separated = rep(separated, p_to_impute)
  if (!length(ts_names) == length(os_names)) 
    stop("Must provide as many true score names (ts_names) as observed score names (os_names),  or leave ts_names blank and the prefix 'true_' will be appended to os_names to name the resulting true scores")
  if (!all(score_types %in% c("CTT", "EAP", "ML", "crosswalk"))) 
    stop("Please specify score_type from the available types for true score imputation ('CTT', 'EAP', or 'ML')")
  if (!all(metrics %in% c(NULL, NA, "z", "T", "standard"))) 
    stop("Please specify metric from the available metrics for true score imputation ('z', 'T', 'standard')")
  # if (any(score_types %in% c("ML"))) 
  #   stop("Warning: Maximum likelihood (ML) scoring is not recommended due to poor simulation performance. Proceed with caution.")
  for (i in seq_len(length(score_types))) {
    is_null_or_na = function(x) if (!is.null(x)) 
      is.na(x)
    else is.null(x)
    is_null_metric = is_null_or_na(metrics[i])
    is_null_mean = is_null_or_na(mean[i])
    is_null_var_ts = is_null_or_na(var_ts[i])
    if ((is_null_metric & (is_null_mean | is_null_var_ts)) | 
        (!is_null_metric & (!is_null_mean | !is_null_var_ts))) 
      stop(paste0("Problem with variable ", i, ": Either assign a metric to each true score variable (e.g., 'T' for T scores) or assign BOTH a mean and var_ts for that variable."))
    if (score_types[i] %in% c("EAP", "ML")) {
      is_null_sename = is.null(se_names[i]) | is.na(se_names[i])
      if (is_null_sename) 
        stop(paste0("Problem with variable ", i, ": Each observed score (os_name) based on EAP or ML scoring must include a corresponding standard error (se_names)."))
    }
    else if (score_types[i] == "CTT") {
      is_null_reliability = is.null(reliability[i]) | is.na(reliability[i])
      if (is_null_reliability) 
        stop(paste0("Problem with variable ", i, ": Each observed score (os_name) based on CTT scoring must include a corresponding estimate of reliability."))
    }
    else if (score_types[i] == "crosswalk") {
      is_null_linked_obs_cor = is.null(linked_obs_cor[i]) | 
        is.na(linked_obs_cor[i])
      if (is_null_linked_obs_cor) 
        stop(paste0("Problem with variable ", i, ": Each observed score (os_name) based on crosswalk scoring must include a correlation between linked and observed scores."))
      is_null_truncate = is.null(truncate[i]) | is.na(truncate[i])
      if (is_null_truncate) 
        stop(paste0("Problem with variable ", i, ": Each observed score (os_name) based on crosswalk scoring must include a truncate between linked and observed scores."))
    }
  }
  required_variables = c(os_names, se_names)
  if (any(is.character(linked_obs_cor))) {
    required_variables = c(required_variables, linked_obs_cor[is.character(linked_obs_cor)])
  }
  missing_variables = required_variables[which(!required_variables %in% 
                                                 names(data))]
  if (length(missing_variables) > 0) 
    stop(paste0("The following os_names and/or se_names were not found in the data: ", 
                paste(missing_variables, collapse = ", "), "."))
  if (!is.data.frame(data)) 
    stop("Please provide data as data.frame")
  non_numeric_variables = names(data)[which(!sapply(data, class) %in% 
                                              c("integer", "numeric"))]
  if (length(non_numeric_variables) > 0) {
    if (any(required_variables %in% non_numeric_variables)) 
      stop(paste0("The following observed score and/or standard error variables are not numeric: ", 
                  paste(intersect(non_numeric_variables, required_variables), 
                        collapse = ", "), ". Please convert them to numeric prior to true score imputation."))
    warning(paste0("The following variables are not numeric and will be ignored during imputation: ", 
                   paste(non_numeric_variables, collapse = ", "), ". This implementation of true score imputation does not allow non-numeric variables in the imputation model."))
  }
  blocks = setNames(names(data), names(data))
  method = rep("pmm", ncol(data))
  predictor_matrix = matrix(1, ncol(data), ncol(data)) - diag(ncol(data))
  if (length(non_numeric_variables) > 0) {
    predictor_matrix[which(names(data) %in% non_numeric_variables), 
    ] = 0
    predictor_matrix[, which(names(data) %in% non_numeric_variables)] = 0
    method[which(names(data) %in% non_numeric_variables)] = ""
  }
  blots = list()
  for (i in seq_len(length(ts_names))) {
    blots[[ts_names[i]]] = list()
    blots[[ts_names[i]]]$os_name = os_names[i]
    if (score_types[i] %in% c("EAP", "ML")) {
      blots[[ts_names[i]]]$se_name = se_names[i]
    }
    else if (score_types[i] == "CTT") {
      blots[[ts_names[i]]]$reliability = reliability[i]
    }
    if (score_types[i] == "crosswalk") {
      blots[[ts_names[i]]]$linked_obs_cor = linked_obs_cor[i]
      blots[[ts_names[i]]]$truncate = truncate[i]
    }
    blots[[ts_names[i]]]$score_type = score_types[i]
    blots[[ts_names[i]]]$separated = separated[i]
    if (!is_null_or_na(metrics[i])) {
      if (metrics[i] == "z") {
        mean[i] = 0
        var_ts[i] = 1
      }
      else if (metrics[i] == "T") {
        mean[i] = 50
        var_ts[i] = 100
      }
      else if (metrics[i] == "standard") {
        mean[i] = 100
        var_ts[i] = 225
      }
    }
    blots[[ts_names[i]]]$mean = mean[i]
    blots[[ts_names[i]]]$var_ts = var_ts[i]
  }
  blots = lapply(blots, function(x) list(calibration = x))
  for (n in ts_names) {
    blocks[[n]] = n
    if (is.null(data[[n]])) {
      data[[n]] = NA
      method = c(method, "truescore")
      predictor_matrix = rbind(predictor_matrix, 1)
      predictor_matrix = cbind(predictor_matrix, 0)
    }
    else {
      which_n = which(names(data) == n)
      method[which_n] = "truescore"
      predictor_matrix[which_n, ] = 1
      predictor_matrix[, which_n] = 0
    }
  }
  #browser()
  if (length(non_numeric_variables) > 0) {
    predictor_matrix[, which(names(data) %in% non_numeric_variables)] = 0
  }
  do.call(mice, c(list(data = data, method = method, blocks = blocks, 
                       blots = blots, predictorMatrix = predictor_matrix, remove.constant = F, 
                       remove.collinear = F), mice_args))
}