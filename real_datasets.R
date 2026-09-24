suppressPackageStartupMessages({
  library(survival)
  library(SurvMetrics)
  library(riskRegression)
  library(caret)
  library(knitr)
})


dataset_name <- "CIBMTR"   # "CIBMTR", "flchain", "METABRIC", "TELCO", "SUPPORT", "GBSG","Breast_Cancer"


# METABRIC settings

# metabric_file <- "METABRIC.csv"
# metabric_time_col <- "Overall Survival (Months)"
# metabric_event_col <- "Overall Survival Status"
# metabric_id_cols <- c(
#   "Patient ID",
#   "Patient's Vital Status",
#   "Relapse Free Status",
#   "Relapse Free Survival (Months)",
#   "Cancer Type Detailed",
#   "Oncotree Code",
#   "Cancer Type",
#   "Sex"
# )
# metabric_categorical_cols <- NULL
# metabric_continuous_cols  <- NULL
# metabric_event_positive_labels <- c("DECEASED", "DIED", "DEAD")

# Use the pycox-extracted dataset instead of raw METABRIC CSV
metabric_file    <- "df_dataset.csv"
metabric_time_col  <- "duration"
metabric_event_col <- "event"

metabric_id_cols <- NULL

# x0-x3 continuous, x4-x7 binary categorical, x8 continuous (age)
metabric_continuous_cols  <- c("x0", "x1", "x2", "x3", "x8")
metabric_categorical_cols <- c("x4", "x5", "x6", "x7")

metabric_event_positive_labels <- c("1")



# SUPPORT settings

support_file <- "df_dataset_support.csv"
support_time_col <- "duration"
support_event_col <- "event"

support_id_cols <- NULL

# Based on the uploaded df_dataset_support.csv:
# x1-x6 behave like categorical/discrete variables,
# x0 and x7-x13 are treated as continuous.
support_continuous_cols  <- c("x0", "x7", "x8", "x9", "x10", "x11", "x12", "x13")
support_categorical_cols <- c("x1", "x2", "x3", "x4", "x5", "x6")

support_event_positive_labels <- c("1")

# TELCO settings

telco_file <- "customer_data (1).csv"
telco_time_col <- "tenure"
telco_event_col <- "Churn"

telco_id_cols <- c("customerID")

telco_drop_categorical_cols <- c(
  "MultipleLines",
  "OnlineSecurity",
  "OnlineBackup",
  "DeviceProtection",
  "TechSupport",
  "StreamingTV",
  "StreamingMovies"
)

telco_categorical_cols <- setdiff(
  c(
    "gender",
    "SeniorCitizen",
    "Partner",
    "Dependents",
    "PhoneService",
    "MultipleLines",
    "InternetService",
    "OnlineSecurity",
    "OnlineBackup",
    "DeviceProtection",
    "TechSupport",
    "StreamingTV",
    "StreamingMovies",
    "Contract",
    "PaperlessBilling",
    "PaymentMethod"
  ),
  telco_drop_categorical_cols
)

telco_continuous_cols <- c(
  "MonthlyCharges",
  "TotalCharges"
)

telco_event_positive_labels <- c("YES")


# GBSG settings
gbsg_file <- "df_dataset_gbsg.csv"
gbsg_time_col <- "duration"
gbsg_event_col <- "event"

gbsg_id_cols <- NULL

# Based on df_dataset_gbsg.csv:
# x0-x4 are continuous, x5-x6 are binary categorical
gbsg_continuous_cols  <- c( "x3", "x4","x5", "x6")
gbsg_categorical_cols <- c("x0", "x1", "x2")

gbsg_event_positive_labels <- c("1")





# Breast Cancer settings

breast_cancer_file <- "Breast_Cancer.csv"
breast_cancer_time_col <- "Survival Months"
breast_cancer_event_col <- "Status"

breast_cancer_id_cols <- NULL

breast_cancer_continuous_cols <- c(
  "Age",
  "Tumor Size",
  "Regional Node Examined",
  "Reginol Node Positive"
)

breast_cancer_categorical_cols <- c(
  "Race",
  "Marital Status",
  "T Stage ",
  "N Stage",
  "6th Stage",
  "differentiate",
  "Grade",
  "A Stage",
  "Estrogen Status",
  "Progesterone Status"
)

breast_cancer_event_positive_labels <- c("DEAD")




# Method toggles

include_km_greedy_separate <- FALSE
include_km_greedy_global   <- TRUE


use_joint_categorical_feature <- FALSE
joint_cat_name <- "cat_joint"

make_joint_categorical_feature <- function(df, cat_cols, new_name = "cat_joint") {
  if (length(cat_cols) == 0) return(df)
  
  joint_vals <- apply(
    df[, cat_cols, drop = FALSE],
    1,
    function(row) paste(row, collapse = " | ")
  )
  
  df[[new_name]] <- joint_vals
  df
}


# KM greedy global stopping settings
# km_greedy_min_rel_improve <- 0.001   # 0.01 = 1%, 0.001 = 0.1%, 0.0001 = 0.01%



# Repeated random split settings

R <- 20



# Output folder

plot_dir <- paste0("plots_all_methods_", dataset_name)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#  
# Helpers
#  
safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]+", "_", x)
}

save_pdf_plot <- function(filename, expr, width = 10, height = 7) {
  pdf(file = file.path(plot_dir, filename), width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  eval.parent(substitute(expr))
}

pretty_feature_name <- function(cc, original_name_lookup = NULL) {
  if (!is.null(original_name_lookup) && cc %in% names(original_name_lookup)) {
    return(original_name_lookup[[cc]])
  }
  cc
}

cox_aic <- function(cox_fit, n_for_bic_aicc) {
  ll <- as.numeric(cox_fit$loglik[2])
  beta <- coef(cox_fit)
  k <- sum(is.finite(beta) & !is.na(beta))
  
  aic  <- -2 * ll + 2 * k
  bic  <- -2 * ll + log(n_for_bic_aicc) * k
  denom <- max(1e-9, n_for_bic_aicc - k - 1)
  aicc <- aic + (2 * k * (k + 1)) / denom
  
  c(aic = aic, bic = bic, aicc = aicc)
}

cox_aic_value <- function(cox_fit, n_for_bic_aicc) {
  as.numeric(cox_aic(cox_fit, n_for_bic_aicc)["aic"])
}

cindex_like_lifelines <- function(time, event, risk_lp) {
  conc <- concordance(Surv(time, event) ~ I(-risk_lp))
  as.numeric(conc$concordance)
}

#  
# AUC / iAUC helpers
#  
iauc_weights_fS <- function(time, event, eval_times) {
  ok_eval <- is.finite(eval_times)
  eval_times <- eval_times[ok_eval]
  
  if (length(eval_times) < 2) return(rep(NA_real_, length(eval_times)))
  
  sf <- survfit(Surv(time, event) ~ 1)
  S_right <- as.numeric(summary(sf, times = eval_times, extend = TRUE)$surv)
  
  S_left <- c(1, head(S_right, -1))
  dF <- pmax(0, S_left - S_right)
  w_raw <- dF * S_right
  
  if (sum(w_raw, na.rm = TRUE) <= 0) {
    return(rep(NA_real_, length(eval_times)))
  }
  
  w_raw / sum(w_raw, na.rm = TRUE)
}

iauc_fS <- function(times, auc, time, event) {
  ok <- is.finite(times) & is.finite(auc)
  times <- times[ok]
  auc   <- auc[ok]
  
  if (length(times) < 2) return(NA_real_)
  
  o <- order(times)
  times <- times[o]
  auc   <- auc[o]
  
  w <- iauc_weights_fS(time = time, event = event, eval_times = times)
  if (all(!is.finite(w))) return(NA_real_)
  
  sum(w * auc, na.rm = TRUE)
}

#  
# BLL / IBLL helpers
#  
trapz_integral <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 2) return(NA_real_)
  
  o <- order(x)
  x <- x[o]
  y <- y[o]
  
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

km_censor_surv_at_times <- function(time, event, eval_times) {
  # G(t) = P(C > t), estimated by KM on censoring process
  # censoring event indicator = 1 - event
  sf_cens <- survfit(Surv(time, 1 - event) ~ 1)
  as.numeric(summary(sf_cens, times = eval_times, extend = TRUE)$surv)
}

bll_ipcw <- function(time, event, surv_mat, eval_times, eps = 1e-12) {
  # surv_mat: n x length(eval_times), entries S_hat(t | x_i)
  n <- length(time)
  m <- length(eval_times)
  
  if (!is.matrix(surv_mat)) {
    surv_mat <- as.matrix(surv_mat)
  }
  
  if (nrow(surv_mat) != n || ncol(surv_mat) != m) {
    stop("surv_mat must have dimensions length(time) x length(eval_times)")
  }
  
  # G(T_i) for event contributions
  G_Ti <- km_censor_surv_at_times(time, event, time)
  G_Ti <- pmax(G_Ti, eps)
  
  bll <- rep(NA_real_, m)
  
  for (k in seq_len(m)) {
    tt <- eval_times[k]
    
    G_t <- km_censor_surv_at_times(time, event, tt)
    G_t <- pmax(G_t, eps)
    
    S_t <- pmin(pmax(surv_mat[, k], eps), 1 - eps)
    
    event_term <- log(1 - S_t) * as.numeric(time <= tt & event == 1) / G_Ti
    cens_term  <- log(S_t)     * as.numeric(time >  tt)             / G_t
    
    bll[k] <- mean(event_term + cens_term, na.rm = TRUE)
  }
  
  data.frame(times = eval_times, BLL = bll)
}

ibll_ipcw <- function(times, bll) {
  ok <- is.finite(times) & is.finite(bll)
  times <- times[ok]
  bll   <- bll[ok]
  
  if (length(times) < 2) return(NA_real_)
  
  rng <- max(times) - min(times)
  if (!is.finite(rng) || rng <= 0) return(NA_real_)
  
  trapz_integral(times, bll) / rng
}


# Dataset loaders
# 

#  
# Read CIBMTR using dictionary
#  
read_cibmtr_from_dictionary <- function(
    cibmtr_file = "CIBMTR.csv",
    dict_file = "data_dictionary.csv",
    id_col = "ID",
    time_col = "efs_time",
    event_col = "efs"
) {
  dict <- read.csv(dict_file, stringsAsFactors = FALSE)
  df   <- read.csv(cibmtr_file, stringsAsFactors = FALSE)
  
  names(dict) <- trimws(names(dict))
  dict$variable <- trimws(dict$variable)
  dict$type     <- trimws(dict$type)
  
  dict <- dict[dict$variable %in% names(df), , drop = FALSE]
  
  cat_vars  <- dict$variable[dict$type == "Categorical"]
  cont_vars <- dict$variable[dict$type == "Numerical"]
  
  cat_vars  <- setdiff(cat_vars, c(id_col, time_col, event_col))
  cont_vars <- setdiff(cont_vars, c(id_col, time_col, event_col))
  
  for (v in cat_vars) df[[v]] <- as.character(df[[v]])
  for (v in cont_vars) df[[v]] <- as.numeric(df[[v]])
  
  df[[time_col]]  <- as.numeric(df[[time_col]])
  df[[event_col]] <- as.numeric(df[[event_col]])
  
  cat_new_names  <- paste0("cat_", seq_along(cat_vars))
  cont_new_names <- paste0("cont_", seq_along(cont_vars))
  
  names(cat_new_names)  <- cat_vars
  names(cont_new_names) <- cont_vars
  
  rename_map <- c(
    cat_new_names,
    cont_new_names,
    setNames("Time", time_col),
    setNames("Event", event_col)
  )
  
  common_to_rename <- intersect(names(rename_map), names(df))
  names(df)[match(common_to_rename, names(df))] <- rename_map[common_to_rename]
  
  if (id_col %in% names(df)) df[[id_col]] <- NULL
  
  final_cont <- unname(cont_new_names)
  final_cat  <- unname(cat_new_names)
  
  keep_cols <- c(final_cont, final_cat, "Time", "Event")
  keep_cols <- keep_cols[keep_cols %in% names(df)]
  df <- df[, keep_cols, drop = FALSE]
  
  reverse_cat_map  <- setNames(names(cat_new_names), unname(cat_new_names))
  reverse_cont_map <- setNames(names(cont_new_names), unname(cont_new_names))
  
  list(
    df = df,
    cont_cols = final_cont,
    cat_cols = final_cat,
    original_cont_names = cont_vars,
    original_cat_names = cat_vars,
    rename_map = rename_map,
    reverse_cat_map = reverse_cat_map,
    reverse_cont_map = reverse_cont_map
  )
}

#  
# Generic CSV dataset reader
#  
read_generic_survival_csv <- function(
    file,
    time_col,
    event_col,
    categorical_cols = NULL,
    continuous_cols = NULL,
    id_cols = NULL,
    positive_event_labels = c("1", "TRUE", "DECEASED", "DEAD", "EVENT", "YES")
) {
  df <- read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (!is.null(id_cols)) {
    drop_cols <- intersect(id_cols, names(df))
    if (length(drop_cols) > 0) {
      df <- df[, setdiff(names(df), drop_cols), drop = FALSE]
    }
  }
  
  if (!(time_col %in% names(df))) {
    stop("Time column not found in file: ", time_col)
  }
  if (!(event_col %in% names(df))) {
    stop("Event column not found in file: ", event_col)
  }
  
  feature_cols <- setdiff(names(df), c(time_col, event_col))
  
  if (is.null(categorical_cols) && is.null(continuous_cols)) {
    categorical_cols <- feature_cols[sapply(df[feature_cols], function(x) is.character(x) || is.factor(x))]
    continuous_cols  <- setdiff(feature_cols, categorical_cols)
  } else {
    if (is.null(categorical_cols)) categorical_cols <- character(0)
    if (is.null(continuous_cols))  continuous_cols  <- character(0)
  }
  
  categorical_cols <- intersect(categorical_cols, names(df))
  continuous_cols  <- intersect(continuous_cols, names(df))
  
  for (v in categorical_cols) df[[v]] <- as.character(df[[v]])
  for (v in continuous_cols)  df[[v]] <- as.numeric(df[[v]])
  
  df[[time_col]] <- as.numeric(df[[time_col]])
  
  if (is.character(df[[event_col]]) || is.factor(df[[event_col]])) {
    x <- toupper(trimws(as.character(df[[event_col]])))
    df[[event_col]] <- ifelse(x %in% toupper(positive_event_labels), 1, 0)
  } else {
    df[[event_col]] <- as.numeric(df[[event_col]])
  }
  
  cat_new_names  <- paste0("cat_", seq_along(categorical_cols))
  cont_new_names <- paste0("cont_", seq_along(continuous_cols))
  
  names(cat_new_names)  <- categorical_cols
  names(cont_new_names) <- continuous_cols
  
  rename_map <- c(
    cat_new_names,
    cont_new_names,
    setNames("Time", time_col),
    setNames("Event", event_col)
  )
  
  common_to_rename <- intersect(names(rename_map), names(df))
  names(df)[match(common_to_rename, names(df))] <- rename_map[common_to_rename]
  
  final_cont <- unname(cont_new_names)
  final_cat  <- unname(cat_new_names)
  
  keep_cols <- c(final_cont, final_cat, "Time", "Event")
  keep_cols <- keep_cols[keep_cols %in% names(df)]
  df <- df[, keep_cols, drop = FALSE]
  
  reverse_cat_map  <- setNames(names(cat_new_names), unname(cat_new_names))
  reverse_cont_map <- setNames(names(cont_new_names), unname(cont_new_names))
  
  list(
    df = df,
    cont_cols = final_cont,
    cat_cols = final_cat,
    original_cont_names = continuous_cols,
    original_cat_names = categorical_cols,
    rename_map = rename_map,
    reverse_cat_map = reverse_cat_map,
    reverse_cont_map = reverse_cont_map
  )
}

#  
# flchain loader
#  
read_flchain_dataset <- function() {
  data(flchain, package = "survival")
  df <- flchain
  
  id_like <- intersect(c("sample", "chapter"), names(df))
  if (length(id_like) > 0) {
    df <- df[, setdiff(names(df), id_like), drop = FALSE]
  }
  
  time_col  <- "futime"
  event_col <- "death"
  
  if (!(time_col %in% names(df))) stop("futime not found in flchain")
  if (!(event_col %in% names(df))) stop("death not found in flchain")
  
  feature_cols <- setdiff(names(df), c(time_col, event_col))
  categorical_cols <- feature_cols[sapply(df[feature_cols], function(x) is.factor(x) || is.character(x))]
  continuous_cols  <- setdiff(feature_cols, categorical_cols)
  
  for (v in categorical_cols) df[[v]] <- as.character(df[[v]])
  for (v in continuous_cols)  df[[v]] <- as.numeric(df[[v]])
  
  df[[time_col]]  <- as.numeric(df[[time_col]])
  df[[event_col]] <- as.numeric(df[[event_col]])
  
  cat_new_names  <- paste0("cat_", seq_along(categorical_cols))
  cont_new_names <- paste0("cont_", seq_along(continuous_cols))
  
  names(cat_new_names)  <- categorical_cols
  names(cont_new_names) <- continuous_cols
  
  rename_map <- c(
    cat_new_names,
    cont_new_names,
    setNames("Time", time_col),
    setNames("Event", event_col)
  )
  
  common_to_rename <- intersect(names(rename_map), names(df))
  names(df)[match(common_to_rename, names(df))] <- rename_map[common_to_rename]
  
  final_cont <- unname(cont_new_names)
  final_cat  <- unname(cat_new_names)
  
  keep_cols <- c(final_cont, final_cat, "Time", "Event")
  keep_cols <- keep_cols[keep_cols %in% names(df)]
  df <- df[, keep_cols, drop = FALSE]
  
  reverse_cat_map  <- setNames(names(cat_new_names), unname(cat_new_names))
  reverse_cont_map <- setNames(names(cont_new_names), unname(cont_new_names))
  
  list(
    df = df,
    cont_cols = final_cont,
    cat_cols = final_cat,
    original_cont_names = continuous_cols,
    original_cat_names = categorical_cols,
    rename_map = rename_map,
    reverse_cat_map = reverse_cat_map,
    reverse_cont_map = reverse_cont_map
  )
}

#  
# Master dataset switch
#  
read_dataset_by_name <- function(dataset_name) {
  ds <- toupper(trimws(dataset_name))
  
  if (ds == "CIBMTR") {
    return(
      read_cibmtr_from_dictionary(
        cibmtr_file = "CIBMTR.csv",
        dict_file = "data_dictionary.csv"
      )
    )
  }
  
  if (ds == "FLCHAIN") {
    return(read_flchain_dataset())
  }
  
  if (ds == "SUPPORT") {
    return(
      read_generic_survival_csv(
        file = support_file,
        time_col = support_time_col,
        event_col = support_event_col,
        categorical_cols = support_categorical_cols,
        continuous_cols = support_continuous_cols,
        id_cols = support_id_cols,
        positive_event_labels = support_event_positive_labels
      )
    )
  }
  
  if (ds == "METABRIC") {
    return(
      read_generic_survival_csv(
        file = metabric_file,
        time_col = metabric_time_col,
        event_col = metabric_event_col,
        categorical_cols = metabric_categorical_cols,
        continuous_cols = metabric_continuous_cols,
        id_cols = metabric_id_cols,
        positive_event_labels = metabric_event_positive_labels
      )
    )
  }
  
  if (ds == "GBSG") {
    return(
      read_generic_survival_csv(
        file = gbsg_file,
        time_col = gbsg_time_col,
        event_col = gbsg_event_col,
        categorical_cols = gbsg_categorical_cols,
        continuous_cols = gbsg_continuous_cols,
        id_cols = gbsg_id_cols,
        positive_event_labels = gbsg_event_positive_labels
      )
    )
  }
  
  if (ds == "TELCO") {
    return(
      read_generic_survival_csv(
        file = telco_file,
        time_col = telco_time_col,
        event_col = telco_event_col,
        categorical_cols = telco_categorical_cols,
        continuous_cols = telco_continuous_cols,
        id_cols = telco_id_cols,
        positive_event_labels = telco_event_positive_labels
      )
    )
  }
  
  if (ds %in% c("BREAST_CANCER", "BREAST CANCER", "BREASTCANCER")) {
    return(
      read_generic_survival_csv(
        file = breast_cancer_file,
        time_col = breast_cancer_time_col,
        event_col = breast_cancer_event_col,
        categorical_cols = breast_cancer_categorical_cols,
        continuous_cols = breast_cancer_continuous_cols,
        id_cols = breast_cancer_id_cols,
        positive_event_labels = breast_cancer_event_positive_labels
      )
    )
  }
  
  stop("Unknown dataset_name: ", dataset_name)
}

#  
# Diagnostic categorical plot
#  
plot_cat_distribution_one_pdf <- function(df, cat_cols, reverse_cat_map, file = "categorical_distributions.pdf") {
  pdf(file.path(plot_dir, file), width = 11, height = 7)
  on.exit(dev.off(), add = TRUE)
  
  for (cc in cat_cols) {
    tab <- table(df[[cc]], df$Event, useNA = "ifany")
    
    if (!("0" %in% colnames(tab))) tab <- cbind(tab, "0" = 0)
    if (!("1" %in% colnames(tab))) tab <- cbind(tab, "1" = 0)
    tab <- tab[, c("0", "1"), drop = FALSE]
    
    totals <- rowSums(tab)
    tab <- tab[order(totals, decreasing = TRUE), , drop = FALSE]
    
    original_name <- pretty_feature_name(cc, reverse_cat_map)
    
    par(mar = c(10, 5, 3, 2))
    barplot(
      t(tab),
      beside = FALSE,
      col = c("grey70", "tomato"),
      border = NA,
      las = 2,
      ylab = "Count",
      main = paste("Distribution of", original_name)
    )
    
    legend(
      "topright",
      legend = c("Censored", "Event"),
      fill = c("grey70", "tomato"),
      bty = "n"
    )
  }
}

#  
# One-hot encoding
#  
ohe <- function(cat_cols, train_df_noz, test_df_noz) {
  all_encoded_train <- train_df_noz
  all_encoded_test  <- test_df_noz
  all_censored_categories <- list()
  skipped_cols <- character(0)
  
  for (cc in cat_cols) {
    train_cats <- as.character(train_df_noz[[cc]])
    test_cats  <- as.character(test_df_noz[[cc]])
    
    all_levels <- unique(train_cats)
    
    has_event <- sapply(all_levels, function(lv) {
      any(train_df_noz$Event[train_cats == lv] == 1)
    })
    
    valid_levels    <- all_levels[has_event]
    censored_levels <- all_levels[!has_event]
    
    if (length(valid_levels) < 2) {
      skipped_cols <- c(skipped_cols, cc)
      all_encoded_train <- all_encoded_train[, setdiff(names(all_encoded_train), cc), drop = FALSE]
      all_encoded_test  <- all_encoded_test[, setdiff(names(all_encoded_test), cc), drop = FALSE]
      all_censored_categories[[cc]] <- censored_levels
      next
    }
    
    baseline    <- valid_levels[1]
    kept_levels <- setdiff(valid_levels, baseline)
    col_names   <- paste0(cc, "_", make.names(kept_levels))
    
    encode_cats <- function(cats) {
      mat <- matrix(
        0,
        nrow = length(cats),
        ncol = length(kept_levels),
        dimnames = list(NULL, col_names)
      )
      
      for (i in seq_along(cats)) {
        lv <- cats[i]
        if (lv %in% kept_levels) {
          mat[i, paste0(cc, "_", make.names(lv))] <- 1
        }
        # baseline / fully censored / unseen => all zero
      }
      
      mat
    }
    
    all_encoded_train <- cbind(
      all_encoded_train[, setdiff(names(all_encoded_train), cc), drop = FALSE],
      as.data.frame(encode_cats(train_cats))
    )
    
    all_encoded_test <- cbind(
      all_encoded_test[, setdiff(names(all_encoded_test), cc), drop = FALSE],
      as.data.frame(encode_cats(test_cats))
    )
    
    all_censored_categories[[cc]] <- censored_levels
  }
  
  list(
    train = all_encoded_train,
    test = all_encoded_test,
    cat_censored = all_censored_categories,
    skipped_cols = skipped_cols
  )
}

#  
# Random ordinal encoding
#  
ordinal_encode_random <- function(cat_cols, train_df_noz, test_df_noz) {
  train_out <- train_df_noz
  test_out  <- test_df_noz
  
  mappings <- list()
  fallback_codes <- list()
  
  for (cc in cat_cols) {
    train_cats <- as.character(train_df_noz[[cc]])
    test_cats  <- as.character(test_df_noz[[cc]])
    
    train_levels <- unique(train_cats)
    random_levels <- sample(train_levels, length(train_levels), replace = FALSE)
    code_map <- setNames(seq_along(random_levels), random_levels)
    
    fallback_code <- median(unname(code_map))
    
    train_out[[cc]] <- as.numeric(code_map[train_cats])
    
    test_codes <- as.numeric(code_map[test_cats])
    test_codes[is.na(test_codes)] <- fallback_code
    test_out[[cc]] <- test_codes
    
    mappings[[cc]] <- code_map
    fallback_codes[[cc]] <- fallback_code
  }
  
  list(
    train = train_out,
    test = test_out,
    mappings = mappings,
    fallback_codes = fallback_codes
  )
}

#  
# KM encoding helpers (cached)
#  
predict_surv_at_times <- function(sf, times) {
  out <- summary(sf, times = times, extend = TRUE)$surv
  as.numeric(out)
}

fit_marginal_km <- function(train_df) {
  survfit(Surv(Time, Event) ~ 1, data = train_df)
}

n_times_from_cardinality <- function(K, divisor) {
  max(1L, K %/% divisor)
}

make_fixed_times_from_cardinality <- function(train_df, cat_col, divisor) {
  K <- length(unique(as.character(train_df[[cat_col]])))
  n_times <- n_times_from_cardinality(K, divisor)
  probs <- seq(0.1, 0.9, length.out = n_times)
  unique(as.numeric(quantile(train_df$Time, probs = probs, na.rm = TRUE)))
}

make_fixed_times_list <- function(train_df, cat_cols, divisor) {
  out <- vector("list", length(cat_cols))
  names(out) <- cat_cols
  
  for (cc in cat_cols) {
    out[[cc]] <- make_fixed_times_from_cardinality(train_df, cc, divisor)
  }
  
  out
}

precompute_km_feature <- function(train_df, cat_col, candidate_times, fallback_km = NULL) {
  cats_train <- as.character(train_df[[cat_col]])
  levels_train <- unique(cats_train)
  
  if (is.null(fallback_km)) {
    fallback_km <- survfit(Surv(Time, Event) ~ 1, data = train_df)
  }
  
  fallback_vec <- predict_surv_at_times(fallback_km, candidate_times)
  
  km_lookup <- vector("list", length(levels_train))
  names(km_lookup) <- levels_train
  
  for (lv in levels_train) {
    sub <- train_df[cats_train == lv, , drop = FALSE]
    sf <- survfit(Surv(Time, Event) ~ 1, data = sub)
    km_lookup[[lv]] <- predict_surv_at_times(sf, candidate_times)
  }
  
  list(
    cat_col = cat_col,
    candidate_times = candidate_times,
    lookup = km_lookup,
    fallback = fallback_vec
  )
}

precompute_km_all_features <- function(train_df, cat_cols, candidate_times, fallback_km = NULL) {
  out <- vector("list", length(cat_cols))
  names(out) <- cat_cols
  
  if (is.null(fallback_km)) {
    fallback_km <- survfit(Surv(Time, Event) ~ 1, data = train_df)
  }
  
  for (cc in cat_cols) {
    out[[cc]] <- precompute_km_feature(
      train_df = train_df,
      cat_col = cc,
      candidate_times = candidate_times,
      fallback_km = fallback_km
    )
  }
  
  out
}

build_full_candidate_matrix_one <- function(df, precomp_feature) {
  cats <- as.character(df[[precomp_feature$cat_col]])
  uniq <- unique(cats)
  Tn <- length(precomp_feature$candidate_times)
  
  row_lookup <- lapply(uniq, function(lv) {
    vec <- precomp_feature$lookup[[lv]]
    if (is.null(vec)) vec <- precomp_feature$fallback
    vec
  })
  names(row_lookup) <- uniq
  
  x <- do.call(rbind, row_lookup[cats])
  x <- as.matrix(x)
  colnames(x) <- paste0(precomp_feature$cat_col, "_km_t", seq_len(Tn))
  x
}

build_full_candidate_matrices <- function(df, precomp_all, cat_cols) {
  out <- vector("list", length(cat_cols))
  names(out) <- cat_cols
  
  for (cc in cat_cols) {
    out[[cc]] <- build_full_candidate_matrix_one(df, precomp_all[[cc]])
  }
  
  out
}

times_to_indices <- function(times_vec, candidate_times, tol = 1e-10) {
  if (length(times_vec) == 0) return(integer(0))
  
  idx <- integer(length(times_vec))
  for (i in seq_along(times_vec)) {
    diffs <- abs(candidate_times - times_vec[i])
    j <- which.min(diffs)
    if (!is.finite(diffs[j]) || diffs[j] > tol) {
      stop("Could not match selected time to shared candidate grid.")
    }
    idx[i] <- j
  }
  unique(idx)
}

times_list_to_indices <- function(times_list, candidate_times, cat_cols, tol = 1e-10) {
  out <- vector("list", length(cat_cols))
  names(out) <- cat_cols
  
  for (cc in cat_cols) {
    out[[cc]] <- times_to_indices(times_list[[cc]], candidate_times, tol = tol)
  }
  
  out
}

assemble_selected_km_matrix <- function(full_mats, selected_idx, cat_cols) {
  parts <- list()
  
  for (cc in cat_cols) {
    idx <- selected_idx[[cc]]
    if (length(idx) > 0) {
      parts[[length(parts) + 1]] <- full_mats[[cc]][, idx, drop = FALSE]
    }
  }
  
  if (length(parts) == 0) {
    first_cc <- cat_cols[1]
    matrix(nrow = nrow(full_mats[[first_cc]]), ncol = 0)
  } else {
    do.call(cbind, parts)
  }
}

km_encode_multiple_cats <- function(df, full_mats, selected_idx, cat_cols) {
  assemble_selected_km_matrix(full_mats = full_mats, selected_idx = selected_idx, cat_cols = cat_cols)
}

km_encode_greedy_separate <- function(df, full_mats, selected_idx, cat_cols) {
  assemble_selected_km_matrix(full_mats = full_mats, selected_idx = selected_idx, cat_cols = cat_cols)
}

km_encode_greedy <- function(df, full_mats, selected_idx, cat_cols) {
  assemble_selected_km_matrix(full_mats = full_mats, selected_idx = selected_idx, cat_cols = cat_cols)
}

#  
# Evaluation
#  
evaluate_cox_model <- function(cox_fit, df_train, df_test, time_interest) {
  b2 <- Score(
    object = list(cox = cox_fit),
    formula = Surv(Time, Event) ~ 1,
    data = df_test,
    times = time_interest,
    metrics = c("Brier", "auc"),
    cens.model = "km",
    cens.data = df_train,
    conf.int = FALSE,
    null.model = FALSE,
    summary = "ibs"
  )
  
  brier_values <- as.data.frame(b2$Brier$score)
  auc_values   <- as.data.frame(b2$AUC$score)
  ibs <- tail(b2$Brier$score$IBS, 1)
  
  lp <- as.numeric(predict(cox_fit, newdata = df_test, type = "lp"))
  cindex <- cindex_like_lifelines(df_test$Time, df_test$Event, lp)
  
  n_events <- max(2L, as.integer(sum(df_train$Event)))
  ic <- cox_aic(cox_fit, n_events)
  
  iauc <- iauc_fS(
    times = auc_values$times,
    auc = auc_values$AUC,
    time = df_test$Time,
    event = df_test$Event
  )
  
  # Survival predictions S_hat(t | x) = 1 - risk(t | x)
  risk_mat <- predictRisk(
    object = cox_fit,
    newdata = df_test,
    times = time_interest
  )
  
  surv_mat <- 1 - as.matrix(risk_mat)
  surv_mat <- pmin(pmax(surv_mat, 1e-12), 1 - 1e-12)
  
  bll_values <- bll_ipcw(
    time = df_test$Time,
    event = df_test$Event,
    surv_mat = surv_mat,
    eval_times = time_interest
  )
  
  ibll <- ibll_ipcw(
    times = bll_values$times,
    bll = bll_values$BLL
  )
  
  list(
    brier_values = brier_values,
    auc_values = auc_values,
    bll_values = bll_values,
    ibs = as.numeric(ibs),
    cindex = cindex,
    iauc = iauc,
    ibll = ibll,
    aic = as.numeric(ic["aic"]),
    bic = as.numeric(ic["bic"]),
    aicc = as.numeric(ic["aicc"])
  )
}

safe_eval_cox <- function(fit_obj, df_train, df_test, time_interest) {
  out <- list(
    fit_success = FALSE,
    score_success = FALSE,
    has_na_coef = NA,
    brier_values = NULL,
    auc_values = NULL,
    bll_values = NULL,
    ibs = NA_real_,
    cindex = NA_real_,
    iauc = NA_real_,
    ibll = NA_real_,
    aic = NA_real_,
    bic = NA_real_,
    aicc = NA_real_
  )
  
  if (inherits(fit_obj, "try-error") || is.null(fit_obj)) return(out)
  
  out$fit_success <- TRUE
  
  cf <- try(coef(fit_obj), silent = TRUE)
  if (!inherits(cf, "try-error")) out$has_na_coef <- anyNA(cf)
  
  ev <- try(evaluate_cox_model(fit_obj, df_train, df_test, time_interest), silent = TRUE)
  
  if (!inherits(ev, "try-error")) {
    out$score_success <- TRUE
    out$brier_values <- ev$brier_values
    out$auc_values   <- ev$auc_values
    out$bll_values   <- ev$bll_values
    out$ibs    <- ev$ibs
    out$cindex <- ev$cindex
    out$iauc   <- ev$iauc
    out$ibll   <- ev$ibll
    out$aic    <- ev$aic
    out$bic    <- ev$bic
    out$aicc   <- ev$aicc
  }
  
  out
}

# 
# km_greedy_select_times <- function(train_df_noz, cont_cols, cat_cols, candidate_times,
#                                    candidate_pool_idx, full_mats_train,
#                                    min_rel_improve = 0.0001, max_total = Inf) {
#   selected_idx <- vector("list", length(cat_cols))
#   names(selected_idx) <- cat_cols
#   for (cc in cat_cols) selected_idx[[cc]] <- integer(0)
#   
#   selection_log <- data.frame(
#     step = integer(0),
#     cat = character(0),
#     time = numeric(0),
#     loss = numeric(0),
#     improve = numeric(0),
#     rel_improve = numeric(0)
#   )
#   
#   max_times_per_cat <- setNames(
#     sapply(cat_cols, function(cc) length(unique(as.character(train_df_noz[[cc]]))) %/% 2),
#     cat_cols
#   )
#   
#   base_df <- train_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE]
#   
#   repeat {
#     best_cc <- NULL
#     best_j <- NA_integer_
#     best_loss <- Inf
#     
#     candidate_grid <- list()
#     for (cc in cat_cols) {
#       if (length(selected_idx[[cc]]) >= max_times_per_cat[[cc]]) next
#       remaining_idx <- candidate_pool_idx[!(candidate_pool_idx %in% selected_idx[[cc]])]
#       if (length(remaining_idx) == 0) next
#       candidate_grid[[cc]] <- remaining_idx
#     }
#     
#     total_trials <- sum(lengths(candidate_grid))
#     if (total_trials == 0) break
#     
#     X_current <- assemble_selected_km_matrix(full_mats_train, selected_idx, cat_cols)
#     
#     for (cc in names(candidate_grid)) {
#       for (j in candidate_grid[[cc]]) {
#         x_new <- full_mats_train[[cc]][, j, drop = FALSE]
#         X_trial <- if (ncol(X_current) == 0) x_new else cbind(X_current, x_new)
#         
#         df_trial <- cbind(base_df, as.data.frame(X_trial))
#         
#         fit_trial <- try(
#           coxph(
#             Surv(Time, Event) ~ .,
#             data = df_trial,
#             ties = "efron",
#             singular.ok = TRUE
#           ),
#           silent = TRUE
#         )
#         if (inherits(fit_trial, "try-error")) next
#         
#         loss_trial <- -as.numeric(fit_trial$loglik[2])
#         
#         if (is.finite(loss_trial) && loss_trial < best_loss) {
#           best_loss <- loss_trial
#           best_cc <- cc
#           best_j <- j
#         }
#       }
#     }
#     
#     if (is.null(best_cc) || !is.finite(best_loss)) break
#     
#     # first selected term is always accepted
#     if (all(lengths(selected_idx) == 0)) {
#       selected_idx[[best_cc]] <- c(selected_idx[[best_cc]], best_j)
#       selection_log <- rbind(
#         selection_log,
#         data.frame(
#           step = sum(lengths(selected_idx)),
#           cat = best_cc,
#           time = candidate_times[best_j],
#           loss = best_loss,
#           improve = NA_real_,
#           rel_improve = NA_real_
#         )
#       )
#       current_best_loss <- best_loss
#       
#       if (sum(lengths(selected_idx)) >= max_total) break
#       next
#     }
#     
#     abs_improve <- current_best_loss - best_loss
#     rel_improve <- abs_improve / max(abs(current_best_loss), 1e-12)
#     
#     if (rel_improve >= min_rel_improve) {
#       selected_idx[[best_cc]] <- c(selected_idx[[best_cc]], best_j)
#       selection_log <- rbind(
#         selection_log,
#         data.frame(
#           step = sum(lengths(selected_idx)),
#           cat = best_cc,
#           time = candidate_times[best_j],
#           loss = best_loss,
#           improve = abs_improve,
#           rel_improve = rel_improve
#         )
#       )
#       current_best_loss <- best_loss
#     } else {
#       break
#     }
#     
#   }
#   
#   list(
#     selected_idx = selected_idx,
#     selected_times = lapply(selected_idx, function(idx) candidate_times[idx]),
#     selection_log = selection_log
#   )
# }



km_greedy_select_times <- function(train_df_noz, cont_cols, cat_cols, candidate_times,
                                   candidate_pool_idx, full_mats_train) {
  selected_idx <- setNames(lapply(cat_cols, function(cc) integer(0)), cat_cols)
  
  selection_log <- data.frame(
    step    = integer(0),
    cat     = character(0),
    time    = numeric(0),
    aic     = numeric(0),
    improve = numeric(0)
  )
  
  max_times_per_cat <- setNames(
    sapply(cat_cols, function(cc) length(unique(as.character(train_df_noz[[cc]])))),
    cat_cols
  )
  
  base_df <- train_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE]
  
  base_fit <- try(
    coxph(
      Surv(Time, Event) ~ .,
      data = base_df,
      ties = "efron",
      singular.ok = FALSE
    ),
    silent = TRUE
  )
  
  if (inherits(base_fit, "try-error")) {
    return(list(
      selected_idx   = selected_idx,
      selected_times = lapply(selected_idx, function(idx) candidate_times[idx]),
      selection_log  = selection_log
    ))
  }
  
  n_events <- max(2L, as.integer(sum(train_df_noz$Event)))
  current_best_aic <- cox_aic_value(base_fit, n_events)
  
  repeat {
    candidate_grid <- lapply(setNames(cat_cols, cat_cols), function(cc) {
      if (length(selected_idx[[cc]]) >= max_times_per_cat[[cc]]) return(integer(0))
      candidate_pool_idx[!(candidate_pool_idx %in% selected_idx[[cc]])]
    })
    
    if (sum(lengths(candidate_grid)) == 0) break
    
    X_current <- assemble_selected_km_matrix(full_mats_train, selected_idx, cat_cols)
    
    best_cc  <- NULL
    best_j   <- NA_integer_
    best_aic <- Inf
    
    for (cc in cat_cols) {
      for (j in candidate_grid[[cc]]) {
        x_new <- full_mats_train[[cc]][, j, drop = FALSE]
        X_trial <- if (ncol(X_current) == 0) x_new else cbind(X_current, x_new)
        
        df_trial <- cbind(base_df, as.data.frame(X_trial))
        
        fit_trial <- try(
          coxph(
            Surv(Time, Event) ~ .,
            data = df_trial,
            ties = "efron",
            singular.ok = FALSE
          ),
          silent = TRUE
        )
        
        if (inherits(fit_trial, "try-error")) next
        
        aic_trial <- cox_aic_value(fit_trial, n_events)
        
        if (is.finite(aic_trial) && aic_trial < best_aic) {
          best_aic <- aic_trial
          best_cc  <- cc
          best_j   <- j
        }
      }
    }
    
    if (is.null(best_cc) || !is.finite(best_aic)) break
    
    aic_improve <- current_best_aic - best_aic
    
    if (is.finite(aic_improve) && aic_improve > 0) {
      selected_idx[[best_cc]] <- c(selected_idx[[best_cc]], best_j)
      
      selection_log <- rbind(
        selection_log,
        data.frame(
          step    = sum(lengths(selected_idx)),
          cat     = best_cc,
          time    = candidate_times[best_j],
          aic     = best_aic,
          improve = aic_improve
        )
      )
      
      current_best_aic <- best_aic
    } else {
      break
    }
  }
  
  list(
    selected_idx   = selected_idx,
    selected_times = lapply(selected_idx, function(idx) candidate_times[idx]),
    selection_log  = selection_log
  )
}



#  
# Plot helpers for selection logs
#  
plot_selection_log_global <- function(selection_log, reverse_cat_map, rep_id, prefix = dataset_name) {
  save_pdf_plot(
    sprintf("%s_km_greedy_global_aic_improvement_iteration_%02d.pdf", prefix, rep_id),
    {
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par), add = TRUE)
      
      par(mar = c(5, 4, 3, 2))
      
      if (nrow(selection_log) > 0) {
        cum_improve_plot <- cumsum(selection_log$improve)
        
        plot(
          selection_log$step,
          cum_improve_plot,
          type = "b",
          xlab = "Global step",
          ylab = "Cumulative AIC decrease",
          xaxt = "n",
          main = "KM greedy global AIC improvement"
        )
        axis(1, at = selection_log$step)
        
        labels <- paste0(
          sapply(selection_log$cat, pretty_feature_name, original_name_lookup = reverse_cat_map),
          "\n",
          "t=", round(selection_log$time, 3)
        )
        
        par(xpd = TRUE)
        text(selection_log$step, cum_improve_plot, labels = labels, pos = 3, cex = 0.75)
      } else {
        plot.new()
        text(0.5, 0.5, "No KM greedy global selections made")
      }
    },
    width = 10,
    height = 7
  )
}

# ============================================================
# Read and preprocess data
# ============================================================
dataset_obj <- read_dataset_by_name(dataset_name)

df <- dataset_obj$df
cont_cols <- dataset_obj$cont_cols
cat_cols  <- dataset_obj$cat_cols
reverse_cat_map <- dataset_obj$reverse_cat_map
reverse_cont_map <- dataset_obj$reverse_cont_map


if (length(cont_cols) == 0) {
  stop("No continuous columns detected. At least one continuous predictor is needed.")
}
if (length(cat_cols) == 0) {
  warning("No categorical columns detected. KM and categorical encoding methods may fail or be empty.")
}

# Drop missing continuous rows
if (length(cont_cols) > 0) {
  df <- df[stats::complete.cases(df[, cont_cols, drop = FALSE]), , drop = FALSE]
}

# Fill missing categorical with explicit "NA"
for (v in cat_cols) {
  x <- as.character(df[[v]])
  x[is.na(x) | trimws(x) == ""] <- "NA"
  df[[v]] <- x
}

if (use_joint_categorical_feature && length(cat_cols) > 0) {
  
  original_cat_cols <- cat_cols
  
  df <- make_joint_categorical_feature(
    df,
    original_cat_cols,
    new_name = joint_cat_name
  )
  
  df <- df[, setdiff(names(df), original_cat_cols), drop = FALSE]
  
  cat_cols <- joint_cat_name
  
  reverse_cat_map <- setNames(
    paste0(
      "Joint categorical feature: ",
      paste(
        sapply(original_cat_cols, pretty_feature_name,
               original_name_lookup = reverse_cat_map),
        collapse = " + "
      )
    ),
    joint_cat_name
  )
}

# Remove rows with missing Time / Event
df <- df[!is.na(df$Time) & !is.na(df$Event), , drop = FALSE]
df$Event <- as.numeric(df$Event)
df <- df[df$Time >= 0, , drop = FALSE]
df$Event[df$Event != 1] <- 0



#### SELECT TRAIN SIZE
#train_frac <- 500/nrow(df)
#train_frac <- 0.7
train_frac <- 0.1

# Remove exact linear-combination columns among continuous features
if (length(cont_cols) > 1) {
  X_cont <- as.matrix(df[, cont_cols, drop = FALSE])
  lc_cont <- findLinearCombos(X_cont)
  
  if (!is.null(lc_cont$linearCombos)) {
    cat("\nLinear combination groups among continuous columns:\n")
    print(lapply(lc_cont$linearCombos, function(idx) cont_cols[idx]))
  }
  
  if (!is.null(lc_cont$remove)) {
    cols_to_remove <- cont_cols[lc_cont$remove]
    cat("\nRemoving linearly dependent continuous columns:\n")
    print(setNames(cols_to_remove, sapply(cols_to_remove, pretty_feature_name, original_name_lookup = reverse_cont_map)))
    df <- df[, setdiff(names(df), cols_to_remove), drop = FALSE]
    cont_cols <- setdiff(cont_cols, cols_to_remove)
  } else {
    cat("\nNo exact linear combinations found among continuous columns.\n")
  }
}

cat("\n========================================\n")
cat("Dataset:", dataset_name, "\n")
cat("Dataset dimensions:", nrow(df), "x", ncol(df), "\n")
cat("Number of events:", sum(df$Event, na.rm = TRUE), "\n")
cat("Event fraction:", mean(df$Event, na.rm = TRUE), "\n")

cat("\nContinuous columns:\n")
print(setNames(cont_cols, sapply(cont_cols, pretty_feature_name, original_name_lookup = reverse_cont_map)))

cat("\nCategorical columns:\n")
print(setNames(cat_cols, sapply(cat_cols, pretty_feature_name, original_name_lookup = reverse_cat_map)))

cat("\nAny missing in Time? ", any(is.na(df$Time)), "\n", sep = "")
cat("Any missing in Event? ", any(is.na(df$Event)), "\n", sep = "")

if (length(cat_cols) > 0) {
  plot_cat_distribution_one_pdf(
    df,
    cat_cols,
    reverse_cat_map,
    file = sprintf("%s_categorical_distributions.pdf", dataset_name)
  )
}

# ============================================================
# Repeated random splits
# ============================================================
results <- vector("list", R)

for (rep in seq_len(R)) {
  cat("\n========================================\n")
  cat("Repetition:", rep, "\n")
  cat("========================================\n")
  
  n <- nrow(df)
  train_idx <- sample(seq_len(n), size = floor(train_frac * n), replace = FALSE)
  
  train_df <- df[train_idx, , drop = FALSE]
  test_df  <- df[-train_idx, , drop = FALSE]
  
  train_df_noz <- train_df
  test_df_noz  <- test_df
  
  #  
  # 1) Continuous only / baseline
  #  
  df_train_cont_only <- train_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE]
  df_test_cont_only  <- test_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE]
  
  cox_cont_only <- coxph(
    Surv(Time, Event) ~ .,
    data = df_train_cont_only,
    ties = "efron",
    singular.ok = FALSE,
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  
  #  
  # 2) Random ordinal encoding
  #  
  if (length(cat_cols) > 0) {
    ordinal_result <- ordinal_encode_random(cat_cols, train_df_noz, test_df_noz)
    df_train_ordinal <- ordinal_result$train
    df_test_ordinal  <- ordinal_result$test
  } else {
    df_train_ordinal <- train_df_noz
    df_test_ordinal  <- test_df_noz
  }
  
  cox_ordinal <- coxph(
    Surv(Time, Event) ~ .,
    data = df_train_ordinal,
    ties = "efron",
    singular.ok = FALSE,
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
  
  #  
  # 3) One-hot encoding
  #  
  if (length(cat_cols) > 0) {
    ohe_result <- ohe(cat_cols, train_df_noz, test_df_noz)
    df_train_ohe <- ohe_result$train
    df_test_ohe  <- ohe_result$test
    
    cat("Fully censored categories in train:\n")
    print(unlist(ohe_result$cat_censored))
    
    cat("Skipped OHE columns:\n")
    print(setNames(ohe_result$skipped_cols, sapply(ohe_result$skipped_cols, pretty_feature_name, original_name_lookup = reverse_cat_map)))
  } else {
    ohe_result <- list(cat_censored = list(), skipped_cols = character(0))
    df_train_ohe <- train_df_noz
    df_test_ohe  <- test_df_noz
  }
  
  cox_ohe <- tryCatch(
    coxph(
      Surv(Time, Event) ~ .,
      data = df_train_ohe,
      ties = "efron",
      singular.ok = FALSE,
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    error = function(e) NULL
  )
  
  cat("Categories with NA coefficients in OHE:\n")
  if (!is.null(cox_ohe)) {
    print(names(coef(cox_ohe))[is.na(coef(cox_ohe))])
  } else {
    print("OHE Cox fit failed")
  }
  
  #  
  # KM methods only if categorical features exist
  #  
  if (length(cat_cols) > 0) {
    fallback_km <- fit_marginal_km(train_df_noz)
    
    fixed_times_k10 <- make_fixed_times_list(
      train_df = train_df_noz,
      cat_cols = cat_cols,
      divisor = 10
    )
    
    fixed_times_k20_raw <- make_fixed_times_list(
      train_df = train_df_noz,
      cat_cols = cat_cols,
      divisor = 20
    )
    
    fixed_times_k20 <- fixed_times_k20_raw
    for (cc in cat_cols) {
      if (length(fixed_times_k20_raw[[cc]]) == length(fixed_times_k10[[cc]])) {
        fixed_times_k20[[cc]] <- numeric(0)
      }
    }
    do_km_k20 <- any(sapply(fixed_times_k20, length) > 0)
    
    
    
    event_times <- sort(unique(train_df_noz$Time[train_df_noz$Event == 1]))
    candidate_times_search <- unique(as.numeric(quantile(
      event_times, probs = seq(0.1, 0.9, length.out = 70)
    )))
    
    # candidate_times_search <- unique(as.numeric(quantile(
    #   train_df_noz$Time,
    #   probs = seq(0.1, 0.9, length.out = 30),
    #   na.rm = TRUE
    # )))
    
    all_times_union <- sort(unique(c(
      unlist(fixed_times_k10, use.names = FALSE),
      unlist(fixed_times_k20_raw, use.names = FALSE),
      candidate_times_search
    )))
    
    
    km_precomp <- precompute_km_all_features(
      train_df = train_df_noz,
      cat_cols = cat_cols,
      candidate_times = all_times_union,
      fallback_km = fallback_km
    )
    
    full_mats_train <- build_full_candidate_matrices(
      df = train_df_noz,
      precomp_all = km_precomp,
      cat_cols = cat_cols
    )
    
    full_mats_test <- build_full_candidate_matrices(
      df = test_df_noz,
      precomp_all = km_precomp,
      cat_cols = cat_cols
    )
    
    fixed_idx_k10 <- times_list_to_indices(
      times_list = fixed_times_k10,
      candidate_times = all_times_union,
      cat_cols = cat_cols
    )
    
    fixed_idx_k20 <- times_list_to_indices(
      times_list = fixed_times_k20,
      candidate_times = all_times_union,
      cat_cols = cat_cols
    )
    
    candidate_pool_idx <- times_to_indices(
      times_vec = candidate_times_search,
      candidate_times = all_times_union
    )
    
    
    #  
    # 4) KM fixed K/10
    # #  
    # X_train_km_k10 <- km_encode_multiple_cats(
    #   df = train_df_noz,
    #   full_mats = full_mats_train,
    #   selected_idx = fixed_idx_k10,
    #   cat_cols = cat_cols
    # )
    # 
    # X_test_km_k10 <- km_encode_multiple_cats(
    #   df = test_df_noz,
    #   full_mats = full_mats_test,
    #   selected_idx = fixed_idx_k10,
    #   cat_cols = cat_cols
    # )
    # 
    # df_train_km_k10 <- cbind(
    #   train_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE],
    #   as.data.frame(X_train_km_k10)
    # )
    # 
    # df_test_km_k10 <- cbind(
    #   test_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE],
    #   as.data.frame(X_test_km_k10)
    # )
    # 
    # cox_km_k10 <- coxph(
    #   Surv(Time, Event) ~ .,
    #   data = df_train_km_k10,
    #   ties = "efron",
    #   singular.ok = TRUE,
    #   x = TRUE,
    #   y = TRUE,
    #   model = TRUE
    # )
    
    #  
    # 5) KM fixed K/20
    #  
    # if (do_km_k20) {
    #   X_train_km_k20 <- km_encode_multiple_cats(
    #     df = train_df_noz,
    #     full_mats = full_mats_train,
    #     selected_idx = fixed_idx_k20,
    #     cat_cols = cat_cols
    #   )
    #   
    #   X_test_km_k20 <- km_encode_multiple_cats(
    #     df = test_df_noz,
    #     full_mats = full_mats_test,
    #     selected_idx = fixed_idx_k20,
    #     cat_cols = cat_cols
    #   )
    #   
    #   df_train_km_k20 <- cbind(
    #     train_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE],
    #     as.data.frame(X_train_km_k20)
    #   )
    #   
    #   df_test_km_k20 <- cbind(
    #     test_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE],
    #     as.data.frame(X_test_km_k20)
    #   )
    #   
    #   cox_km_k20 <- coxph(
    #     Surv(Time, Event) ~ .,
    #     data = df_train_km_k20,
    #     ties = "efron",
    #     singular.ok = TRUE,
    #     x = TRUE,
    #     y = TRUE,
    #     model = TRUE
    #   )
    # }
    
    #  
    # 6) KM greedy separate-by-feature (AIC-based)
    #  
    # if (include_km_greedy_separate) {
    #   greedy_separate_out <- greedy_km_separate_all_cats(
    #     train_df_noz = train_df_noz,
    #     cont_cols = cont_cols,
    #     cat_cols = cat_cols,
    #     candidate_times = all_times_union,
    #     candidate_pool_idx = candidate_pool_idx,
    #     full_mats_train = full_mats_train,
    #     min_aic_drop = km_greedy_min_aic_drop
    #   )
    #   
    #   selected_idx_greedy_separate <- greedy_separate_out$selected_idx
    #   selected_times_greedy_separate <- greedy_separate_out$selected_times
    #   selection_log_greedy_separate <- greedy_separate_out$selection_log
    #   
    #   X_train_km_greedy_separate <- km_encode_greedy_separate(
    #     df = train_df_noz,
    #     full_mats = full_mats_train,
    #     selected_idx = selected_idx_greedy_separate,
    #     cat_cols = cat_cols
    #   )
    #   
    #   X_test_km_greedy_separate <- km_encode_greedy_separate(
    #     df = test_df_noz,
    #     full_mats = full_mats_test,
    #     selected_idx = selected_idx_greedy_separate,
    #     cat_cols = cat_cols
    #   )
    #   
    #   df_train_km_greedy_separate <- cbind(
    #     train_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE],
    #     as.data.frame(X_train_km_greedy_separate)
    #   )
    #   
    #   df_test_km_greedy_separate <- cbind(
    #     test_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE],
    #     as.data.frame(X_test_km_greedy_separate)
    #   )
    #   
    #   cox_km_greedy_separate <- coxph(
    #     Surv(Time, Event) ~ .,
    #     data = df_train_km_greedy_separate,
    #     ties = "efron",
    #     singular.ok = TRUE,
    #     x = TRUE,
    #     y = TRUE,
    #     model = TRUE
    #   )
    # } else {
    #   selected_idx_greedy_separate <- setNames(vector("list", length(cat_cols)), cat_cols)
    #   selected_times_greedy_separate <- setNames(vector("list", length(cat_cols)), cat_cols)
    #   for (cc in cat_cols) {
    #     selected_idx_greedy_separate[[cc]] <- integer(0)
    #     selected_times_greedy_separate[[cc]] <- numeric(0)
    #   }
    #   selection_log_greedy_separate <- data.frame(
    #     cat = character(0), step = integer(0), time = numeric(0),
    #     aic = numeric(0), improve = numeric(0)
    #   )
    # }
    # 
    #  
    # 7) KM greedy  
    #  
    if (include_km_greedy_global) {
      # greedy_out <- km_greedy_select_times(
      #   train_df_noz = train_df_noz,
      #   cont_cols = cont_cols,
      #   cat_cols = cat_cols,
      #   candidate_times = all_times_union,
      #   candidate_pool_idx = candidate_pool_idx,
      #   full_mats_train = full_mats_train,
      #   min_rel_improve = km_greedy_min_rel_improve
      # )
      # 
      greedy_out <- km_greedy_select_times(
        train_df_noz       = train_df_noz,
        cont_cols          = cont_cols,
        cat_cols           = cat_cols,
        candidate_times    = all_times_union,
        candidate_pool_idx = candidate_pool_idx,
        full_mats_train    = full_mats_train
      )
      
      selected_idx_greedy <- greedy_out$selected_idx
      selected_times_greedy <- greedy_out$selected_times
      selection_log_greedy <- greedy_out$selection_log
      
      X_train_km_greedy <- km_encode_greedy(
        df = train_df_noz,
        full_mats = full_mats_train,
        selected_idx = selected_idx_greedy,
        cat_cols = cat_cols
      )
      
      X_test_km_greedy <- km_encode_greedy(
        df = test_df_noz,
        full_mats = full_mats_test,
        selected_idx = selected_idx_greedy,
        cat_cols = cat_cols
      )
      
      df_train_km_greedy <- cbind(
        train_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE],
        as.data.frame(X_train_km_greedy)
      )
      
      df_test_km_greedy <- cbind(
        test_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE],
        as.data.frame(X_test_km_greedy)
      )
      
      cox_km_greedy <- coxph(
        Surv(Time, Event) ~ .,
        data = df_train_km_greedy,
        ties = "efron",
        singular.ok = FALSE,
        x = TRUE,
        y = TRUE,
        model = TRUE
      )
    } else {
      selected_idx_greedy <- setNames(vector("list", length(cat_cols)), cat_cols)
      selected_times_greedy <- setNames(vector("list", length(cat_cols)), cat_cols)
      for (cc in cat_cols) {
        selected_idx_greedy[[cc]] <- integer(0)
        selected_times_greedy[[cc]] <- numeric(0)
      }
      selection_log_greedy <- data.frame(
        step    = integer(0),
        cat     = character(0),
        time    = numeric(0),
        aic     = numeric(0),
        improve = numeric(0)
      )
    }
  } else {
    # do_km_k20 <- FALSE
    # selected_idx_greedy_separate <- list()
    # selected_times_greedy_separate <- list()
    # selection_log_greedy_separate <- data.frame(
    #   cat = character(0), step = integer(0), time = numeric(0),
    #   aic = numeric(0), improve = numeric(0)
    # )
    selected_idx_greedy <- list()
    selected_times_greedy <- list()
    selection_log_greedy <- data.frame(
      step    = integer(0),
      cat     = character(0),
      time    = numeric(0),
      aic     = numeric(0),
      improve = numeric(0)
    )
  }
  
  #  
  # Evaluation time grid
  #  
  time_interest <- unique(seq(
    quantile(df_test_ohe$Time, 0.05, na.rm = TRUE),
    quantile(df_test_ohe$Time, 0.95, na.rm = TRUE),
    length.out = 100
  ))
  
  #  
  # Evaluate all methods
  #  
  eval_cont_only <- evaluate_cox_model(cox_cont_only, df_train_cont_only, df_test_cont_only, time_interest)
  eval_ordinal   <- evaluate_cox_model(cox_ordinal, df_train_ordinal, df_test_ordinal, time_interest)
  eval_ohe <- safe_eval_cox(cox_ohe, df_train_ohe, df_test_ohe, time_interest)
  
  if (length(cat_cols) > 0) {
    # eval_km_k10 <- evaluate_cox_model(cox_km_k10, df_train_km_k10, df_test_km_k10, time_interest)
    # 
    # if (do_km_k20) {
    #   eval_km_k20 <- evaluate_cox_model(cox_km_k20, df_train_km_k20, df_test_km_k20, time_interest)
    # }
    
    if (include_km_greedy_separate) {
      eval_km_greedy_separate <- evaluate_cox_model(
        cox_km_greedy_separate,
        df_train_km_greedy_separate,
        df_test_km_greedy_separate,
        time_interest
      )
    }
    
    if (include_km_greedy_global) {
      eval_km_greedy <- evaluate_cox_model(
        cox_km_greedy,
        df_train_km_greedy,
        df_test_km_greedy,
        time_interest
      )
    }
    
    # n_selected_greedy_separate <- sapply(cat_cols, function(cc) length(selected_times_greedy_separate[[cc]]))
    # counts_text_greedy_separate <- paste0(
    #   sapply(names(n_selected_greedy_separate), pretty_feature_name, original_name_lookup = reverse_cat_map),
    #   "=",
    #   n_selected_greedy_separate,
    #   collapse = ", "
    # )
    
    n_selected_greedy <- sapply(cat_cols, function(cc) length(selected_times_greedy[[cc]]))
    counts_text_greedy <- paste0(
      sapply(names(n_selected_greedy), pretty_feature_name, original_name_lookup = reverse_cat_map),
      "=",
      n_selected_greedy,
      collapse = ", "
    )
  } else {
    eval_km_k10 <- list(cindex = NA_real_, ibs = NA_real_, iauc = NA_real_, aic = NA_real_, bic = NA_real_, aicc = NA_real_,
                        brier_values = eval_cont_only$brier_values, auc_values = eval_cont_only$auc_values)
    counts_text_greedy_separate <- ""
    counts_text_greedy <- ""
  }
  
  censor_prop_train <- 1 - mean(train_df_noz$Event)
  
  #  
  # Per-iteration Brier plot
  #  
  save_pdf_plot(
    sprintf("%s_all_methods_brier_iteration_%02d.pdf", dataset_name, rep),
    {
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par), add = TRUE)
      
      par(mar = c(12, 5, 2, 2))
      par(mgp = c(2, 0.8, 0))
      
      brier_all <- c(
        eval_cont_only$brier_values$Brier,
        eval_ordinal$brier_values$Brier,
        eval_ohe$brier_values$Brier
      )
      # if (length(cat_cols) > 0) brier_all <- c(brier_all, eval_km_k10$brier_values$Brier)
      # if (length(cat_cols) > 0 && do_km_k20) brier_all <- c(brier_all, eval_km_k20$brier_values$Brier)
      if (length(cat_cols) > 0 && include_km_greedy_separate) brier_all <- c(brier_all, eval_km_greedy_separate$brier_values$Brier)
      if (length(cat_cols) > 0 && include_km_greedy_global) brier_all <- c(brier_all, eval_km_greedy$brier_values$Brier)
      
      plot(
        eval_cont_only$brier_values$times, eval_cont_only$brier_values$Brier,
        type = "l", lwd = 2, lty = 1, col = 8,
        xlab = "Time",
        ylab = paste0("Brier score (train censoring=", round(censor_prop_train, 3), ")"),
        ylim = range(brier_all, na.rm = TRUE)
      )
      lines(eval_ordinal$brier_values$times, eval_ordinal$brier_values$Brier, lwd = 2, lty = 2, col = 7)
      lines(eval_ohe$brier_values$times, eval_ohe$brier_values$Brier, lwd = 2, lty = 3, col = 1)
      
      legend_labels <- c(
        paste0("Continuous only (IBS=", round(eval_cont_only$ibs, 3), ")"),
        paste0("ORD (IBS=", round(eval_ordinal$ibs, 3), ")"),
        paste0("OH (IBS=", round(eval_ohe$ibs, 3), ")")
      )
      legend_lty <- c(1, 2, 3)
      legend_col <- c(8, 7, 1)
      
      # if (length(cat_cols) > 0) {
      #   lines(eval_km_k10$brier_values$times, eval_km_k10$brier_values$Brier, lwd = 2, lty = 4, col = 2)
      #   legend_labels <- c(legend_labels, paste0("KM fixed (K/10) (IBS=", round(eval_km_k10$ibs, 3), ")"))
      #   legend_lty <- c(legend_lty, 4)
      #   legend_col <- c(legend_col, 2)
      # }
      # 
      # if (length(cat_cols) > 0 && do_km_k20) {
      #   lines(eval_km_k20$brier_values$times, eval_km_k20$brier_values$Brier, lwd = 2, lty = 5, col = 4)
      #   legend_labels <- c(legend_labels, paste0("KM fixed (K/20) (IBS=", round(eval_km_k20$ibs, 3), ")"))
      #   legend_lty <- c(legend_lty, 5)
      #   legend_col <- c(legend_col, 4)
      # }
      
      if (length(cat_cols) > 0 && include_km_greedy_separate) {
        lines(eval_km_greedy_separate$brier_values$times, eval_km_greedy_separate$brier_values$Brier, lwd = 2, lty = 6, col = 3)
        legend_labels <- c(
          legend_labels,
          paste0("KM greedy separate (", counts_text_greedy_separate, ", IBS=", round(eval_km_greedy_separate$ibs, 3), ")")
        )
        legend_lty <- c(legend_lty, 6)
        legend_col <- c(legend_col, 3)
      }
      
      if (length(cat_cols) > 0 && include_km_greedy_global) {
        lines(eval_km_greedy$brier_values$times, eval_km_greedy$brier_values$Brier, lwd = 2, lty = 7, col = 6)
        legend_labels <- c(
          legend_labels,
          paste0("KM greedy global (", counts_text_greedy, ", IBS=", round(eval_km_greedy$ibs, 3), ")")
        )
        legend_lty <- c(legend_lty, 7)
        legend_col <- c(legend_col, 6)
      }
      
      legend(
        "bottom",
        inset = c(0, -0.32),
        legend = legend_labels,
        lwd = 2,
        lty = legend_lty,
        col = legend_col,
        bty = "n",
        xpd = TRUE,
        ncol = 2,
        cex = 0.85
      )
    },
    width = 11,
    height = 9
  )
  
  #  
  # Per-iteration AUC plot
  #  
  save_pdf_plot(
    sprintf("%s_all_methods_auc_iteration_%02d.pdf", dataset_name, rep),
    {
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par), add = TRUE)
      
      par(mar = c(12, 5, 2, 2))
      par(mgp = c(2, 0.8, 0))
      
      auc_all <- c(
        eval_cont_only$auc_values$AUC,
        eval_ordinal$auc_values$AUC,
        eval_ohe$auc_values$AUC
      )
      # if (length(cat_cols) > 0) auc_all <- c(auc_all, eval_km_k10$auc_values$AUC)
      # if (length(cat_cols) > 0 && do_km_k20) auc_all <- c(auc_all, eval_km_k20$auc_values$AUC)
      if (length(cat_cols) > 0 && include_km_greedy_separate) auc_all <- c(auc_all, eval_km_greedy_separate$auc_values$AUC)
      if (length(cat_cols) > 0 && include_km_greedy_global) auc_all <- c(auc_all, eval_km_greedy$auc_values$AUC)
      
      plot(
        eval_cont_only$auc_values$times, eval_cont_only$auc_values$AUC,
        type = "l", lwd = 2, lty = 1, col = 8,
        xlab = "Time",
        ylab = paste0("AUC (train censoring=", round(censor_prop_train, 3), ")"),
        ylim = range(auc_all, na.rm = TRUE)
      )
      lines(eval_ordinal$auc_values$times, eval_ordinal$auc_values$AUC, lwd = 2, lty = 2, col = 7)
      lines(eval_ohe$auc_values$times, eval_ohe$auc_values$AUC, lwd = 2, lty = 3, col = 1)
      
      legend_labels <- c(
        paste0("Continuous only (iAUC=", round(eval_cont_only$iauc, 3), ")"),
        paste0("ORD (iAUC=", round(eval_ordinal$iauc, 3), ")"),
        paste0("OH (iAUC=", round(eval_ohe$iauc, 3), ")")
      )
      legend_lty <- c(1, 2, 3)
      legend_col <- c(8, 7, 1)
      
      # if (length(cat_cols) > 0) {
      #   lines(eval_km_k10$auc_values$times, eval_km_k10$auc_values$AUC, lwd = 2, lty = 4, col = 2)
      #   legend_labels <- c(legend_labels, paste0("KM fixed Kdiv10 (iAUC=", round(eval_km_k10$iauc, 3), ")"))
      #   legend_lty <- c(legend_lty, 4)
      #   legend_col <- c(legend_col, 2)
      # }
      
      # if (length(cat_cols) > 0 && do_km_k20) {
      #   lines(eval_km_k20$auc_values$times, eval_km_k20$auc_values$AUC, lwd = 2, lty = 5, col = 4)
      #   legend_labels <- c(legend_labels, paste0("KM fixed Kdiv20 (iAUC=", round(eval_km_k20$iauc, 3), ")"))
      #   legend_lty <- c(legend_lty, 5)
      #   legend_col <- c(legend_col, 4)
      # }
      
      if (length(cat_cols) > 0 && include_km_greedy_separate) {
        lines(eval_km_greedy_separate$auc_values$times, eval_km_greedy_separate$auc_values$AUC, lwd = 2, lty = 6, col = 3)
        legend_labels <- c(
          legend_labels,
          paste0("KM greedy separate (", counts_text_greedy_separate, ", iAUC=", round(eval_km_greedy_separate$iauc, 3), ")")
        )
        legend_lty <- c(legend_lty, 6)
        legend_col <- c(legend_col, 3)
      }
      
      if (length(cat_cols) > 0 && include_km_greedy_global) {
        lines(eval_km_greedy$auc_values$times, eval_km_greedy$auc_values$AUC, lwd = 2, lty = 7, col = 6)
        legend_labels <- c(
          legend_labels,
          paste0("KM greedy global (", counts_text_greedy, ", iAUC=", round(eval_km_greedy$iauc, 3), ")")
        )
        legend_lty <- c(legend_lty, 7)
        legend_col <- c(legend_col, 6)
      }
      
      legend(
        "bottom",
        inset = c(0, -0.32),
        legend = legend_labels,
        lwd = 2,
        lty = legend_lty,
        col = legend_col,
        bty = "n",
        xpd = TRUE,
        ncol = 2,
        cex = 0.85
      )
    },
    width = 11,
    height = 9
  )
  
  if (length(cat_cols) > 0 && include_km_greedy_separate) {
    plot_selection_log_separate(selection_log_greedy_separate, cat_cols, reverse_cat_map, rep_id = rep, prefix = dataset_name)
  }
  if (length(cat_cols) > 0 && include_km_greedy_global) {
    plot_selection_log_global(selection_log_greedy, reverse_cat_map, rep_id = rep, prefix = dataset_name)
  }
  
  
  save_pdf_plot(
    sprintf("%s_all_methods_bll_iteration_%02d.pdf", dataset_name, rep),
    {
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par), add = TRUE)
      
      par(mar = c(12, 5, 2, 2))
      par(mgp = c(2, 0.8, 0))
      
      bll_all <- c(
        eval_cont_only$bll_values$BLL,
        eval_ordinal$bll_values$BLL,
        eval_ohe$bll_values$BLL
      )
      # if (length(cat_cols) > 0) bll_all <- c(bll_all, eval_km_k10$bll_values$BLL)
      # if (length(cat_cols) > 0 && do_km_k20) bll_all <- c(bll_all, eval_km_k20$bll_values$BLL)
      if (length(cat_cols) > 0 && include_km_greedy_separate) bll_all <- c(bll_all, eval_km_greedy_separate$bll_values$BLL)
      if (length(cat_cols) > 0 && include_km_greedy_global) bll_all <- c(bll_all, eval_km_greedy$bll_values$BLL)
      
      plot(
        eval_cont_only$bll_values$times, eval_cont_only$bll_values$BLL,
        type = "l", lwd = 2, lty = 1, col = 8,
        xlab = "Time",
        ylab = paste0("Binomial log-likelihood (train censoring=", round(censor_prop_train, 3), ")"),
        ylim = range(bll_all, na.rm = TRUE)
      )
      lines(eval_ordinal$bll_values$times, eval_ordinal$bll_values$BLL, lwd = 2, lty = 2, col = 7)
      lines(eval_ohe$bll_values$times, eval_ohe$bll_values$BLL, lwd = 2, lty = 3, col = 1)
      
      legend_labels <- c(
        paste0("Continuous only (IBLL=", round(eval_cont_only$ibll, 3), ")"),
        paste0("ORD (IBLL=", round(eval_ordinal$ibll, 3), ")"),
        paste0("OH (IBLL=", round(eval_ohe$ibll, 3), ")")
      )
      legend_lty <- c(1, 2, 3)
      legend_col <- c(8, 7, 1)
      
      # if (length(cat_cols) > 0) {
      #   lines(eval_km_k10$bll_values$times, eval_km_k10$bll_values$BLL, lwd = 2, lty = 4, col = 2)
      #   legend_labels <- c(legend_labels, paste0("KM fixed (K/10) (IBLL=", round(eval_km_k10$ibll, 3), ")"))
      #   legend_lty <- c(legend_lty, 4)
      #   legend_col <- c(legend_col, 2)
      # }
      # 
      # if (length(cat_cols) > 0 && do_km_k20) {
      #   lines(eval_km_k20$bll_values$times, eval_km_k20$bll_values$BLL, lwd = 2, lty = 5, col = 4)
      #   legend_labels <- c(legend_labels, paste0("KM fixed (K/20) (IBLL=", round(eval_km_k20$ibll, 3), ")"))
      #   legend_lty <- c(legend_lty, 5)
      #   legend_col <- c(legend_col, 4)
      # }
      
      if (length(cat_cols) > 0 && include_km_greedy_separate) {
        lines(eval_km_greedy_separate$bll_values$times, eval_km_greedy_separate$bll_values$BLL, lwd = 2, lty = 6, col = 3)
        legend_labels <- c(
          legend_labels,
          paste0("KM greedy separate (", counts_text_greedy_separate, ", IBLL=", round(eval_km_greedy_separate$ibll, 3), ")")
        )
        legend_lty <- c(legend_lty, 6)
        legend_col <- c(legend_col, 3)
      }
      
      if (length(cat_cols) > 0 && include_km_greedy_global) {
        lines(eval_km_greedy$bll_values$times, eval_km_greedy$bll_values$BLL, lwd = 2, lty = 7, col = 6)
        legend_labels <- c(
          legend_labels,
          paste0("KM greedy global (", counts_text_greedy, ", IBLL=", round(eval_km_greedy$ibll, 3), ")")
        )
        legend_lty <- c(legend_lty, 7)
        legend_col <- c(legend_col, 6)
      }
      
      legend(
        "bottom",
        inset = c(0, -0.32),
        legend = legend_labels,
        lwd = 2,
        lty = legend_lty,
        col = legend_col,
        bty = "n",
        xpd = TRUE,
        ncol = 2,
        cex = 0.85
      )
    },
    width = 11,
    height = 9
  )
  

  
  results[[rep]] <- list(
    rep = rep,
    
    # selected_times_greedy_separate = selected_times_greedy_separate,
    selected_times_greedy = selected_times_greedy,
    # selection_log_greedy_separate = selection_log_greedy_separate,
    selection_log_greedy = selection_log_greedy,
    
    cindex_cont_only = eval_cont_only$cindex,
    ibs_cont_only = eval_cont_only$ibs,
    iauc_cont_only = eval_cont_only$iauc,
    aic_cont_only = eval_cont_only$aic,
    bic_cont_only = eval_cont_only$bic,
    aicc_cont_only = eval_cont_only$aicc,
    
    cindex_ordinal = eval_ordinal$cindex,
    ibs_ordinal = eval_ordinal$ibs,
    iauc_ordinal = eval_ordinal$iauc,
    aic_ordinal = eval_ordinal$aic,
    bic_ordinal = eval_ordinal$bic,
    aicc_ordinal = eval_ordinal$aicc,
    
    cindex_ohe = eval_ohe$cindex,
    ibs_ohe = eval_ohe$ibs,
    iauc_ohe = eval_ohe$iauc,
    aic_ohe = eval_ohe$aic,
    bic_ohe = eval_ohe$bic,
    aicc_ohe = eval_ohe$aicc,
    
    ibll_cont_only = eval_cont_only$ibll,
    ibll_ordinal   = eval_ordinal$ibll,
    ibll_ohe       = eval_ohe$ibll,
    # ibll_km_k10 = if (length(cat_cols) > 0) eval_km_k10$ibll else NA_real_,
    # ibll_km_k20 = if (length(cat_cols) > 0 && do_km_k20) eval_km_k20$ibll else NA_real_,
    ibll_km_greedy_separate = if (length(cat_cols) > 0 && include_km_greedy_separate) eval_km_greedy_separate$ibll else NA_real_,
    ibll_km_greedy = if (length(cat_cols) > 0 && include_km_greedy_global) eval_km_greedy$ibll else NA_real_,
    
    # cindex_km_k10 = if (length(cat_cols) > 0) eval_km_k10$cindex else NA_real_,
    # ibs_km_k10 = if (length(cat_cols) > 0) eval_km_k10$ibs else NA_real_,
    # iauc_km_k10 = if (length(cat_cols) > 0) eval_km_k10$iauc else NA_real_,
    # aic_km_k10 = if (length(cat_cols) > 0) eval_km_k10$aic else NA_real_,
    # bic_km_k10 = if (length(cat_cols) > 0) eval_km_k10$bic else NA_real_,
    # aicc_km_k10 = if (length(cat_cols) > 0) eval_km_k10$aicc else NA_real_,
    
    # cindex_km_k20 = if (length(cat_cols) > 0 && do_km_k20) eval_km_k20$cindex else NA_real_,
    # ibs_km_k20 = if (length(cat_cols) > 0 && do_km_k20) eval_km_k20$ibs else NA_real_,
    # iauc_km_k20 = if (length(cat_cols) > 0 && do_km_k20) eval_km_k20$iauc else NA_real_,
    # aic_km_k20 = if (length(cat_cols) > 0 && do_km_k20) eval_km_k20$aic else NA_real_,
    # bic_km_k20 = if (length(cat_cols) > 0 && do_km_k20) eval_km_k20$bic else NA_real_,
    # aicc_km_k20 = if (length(cat_cols) > 0 && do_km_k20) eval_km_k20$aicc else NA_real_,
    # 
    cindex_km_greedy_separate = if (length(cat_cols) > 0 && include_km_greedy_separate) eval_km_greedy_separate$cindex else NA_real_,
    ibs_km_greedy_separate = if (length(cat_cols) > 0 && include_km_greedy_separate) eval_km_greedy_separate$ibs else NA_real_,
    iauc_km_greedy_separate = if (length(cat_cols) > 0 && include_km_greedy_separate) eval_km_greedy_separate$iauc else NA_real_,
    aic_km_greedy_separate = if (length(cat_cols) > 0 && include_km_greedy_separate) eval_km_greedy_separate$aic else NA_real_,
    bic_km_greedy_separate = if (length(cat_cols) > 0 && include_km_greedy_separate) eval_km_greedy_separate$bic else NA_real_,
    aicc_km_greedy_separate = if (length(cat_cols) > 0 && include_km_greedy_separate) eval_km_greedy_separate$aicc else NA_real_,
    
    cindex_km_greedy = if (length(cat_cols) > 0 && include_km_greedy_global) eval_km_greedy$cindex else NA_real_,
    ibs_km_greedy = if (length(cat_cols) > 0 && include_km_greedy_global) eval_km_greedy$ibs else NA_real_,
    iauc_km_greedy = if (length(cat_cols) > 0 && include_km_greedy_global) eval_km_greedy$iauc else NA_real_,
    aic_km_greedy = if (length(cat_cols) > 0 && include_km_greedy_global) eval_km_greedy$aic else NA_real_,
    bic_km_greedy = if (length(cat_cols) > 0 && include_km_greedy_global) eval_km_greedy$bic else NA_real_,
    aicc_km_greedy = if (length(cat_cols) > 0 && include_km_greedy_global) eval_km_greedy$aicc else NA_real_
  )
}

# ============================================================
# Results summary table
# ============================================================
df_res <- data.frame(
  rep = sapply(results, `[[`, "rep"),
  
  cindex_cont_only = sapply(results, `[[`, "cindex_cont_only"),
  ibs_cont_only = sapply(results, `[[`, "ibs_cont_only"),
  iauc_cont_only = sapply(results, `[[`, "iauc_cont_only"),
  aic_cont_only = sapply(results, `[[`, "aic_cont_only"),
  bic_cont_only = sapply(results, `[[`, "bic_cont_only"),
  aicc_cont_only = sapply(results, `[[`, "aicc_cont_only"),
  
  cindex_ordinal = sapply(results, `[[`, "cindex_ordinal"),
  ibs_ordinal = sapply(results, `[[`, "ibs_ordinal"),
  iauc_ordinal = sapply(results, `[[`, "iauc_ordinal"),
  aic_ordinal = sapply(results, `[[`, "aic_ordinal"),
  bic_ordinal = sapply(results, `[[`, "bic_ordinal"),
  aicc_ordinal = sapply(results, `[[`, "aicc_ordinal"),
  
  cindex_ohe = sapply(results, `[[`, "cindex_ohe"),
  ibs_ohe = sapply(results, `[[`, "ibs_ohe"),
  iauc_ohe = sapply(results, `[[`, "iauc_ohe"),
  aic_ohe = sapply(results, `[[`, "aic_ohe"),
  bic_ohe = sapply(results, `[[`, "bic_ohe"),
  aicc_ohe = sapply(results, `[[`, "aicc_ohe"),
  
  # cindex_km_k10 = sapply(results, `[[`, "cindex_km_k10"),
  # ibs_km_k10 = sapply(results, `[[`, "ibs_km_k10"),
  # iauc_km_k10 = sapply(results, `[[`, "iauc_km_k10"),
  # aic_km_k10 = sapply(results, `[[`, "aic_km_k10"),
  # bic_km_k10 = sapply(results, `[[`, "bic_km_k10"),
  # aicc_km_k10 = sapply(results, `[[`, "aicc_km_k10"),
  # 
  # cindex_km_k20 = sapply(results, `[[`, "cindex_km_k20"),
  # ibs_km_k20 = sapply(results, `[[`, "ibs_km_k20"),
  # iauc_km_k20 = sapply(results, `[[`, "iauc_km_k20"),
  # aic_km_k20 = sapply(results, `[[`, "aic_km_k20"),
  # bic_km_k20 = sapply(results, `[[`, "bic_km_k20"),
  # aicc_km_k20 = sapply(results, `[[`, "aicc_km_k20"),
  
  cindex_km_greedy_separate = sapply(results, `[[`, "cindex_km_greedy_separate"),
  ibs_km_greedy_separate = sapply(results, `[[`, "ibs_km_greedy_separate"),
  iauc_km_greedy_separate = sapply(results, `[[`, "iauc_km_greedy_separate"),
  aic_km_greedy_separate = sapply(results, `[[`, "aic_km_greedy_separate"),
  bic_km_greedy_separate = sapply(results, `[[`, "bic_km_greedy_separate"),
  aicc_km_greedy_separate = sapply(results, `[[`, "aicc_km_greedy_separate"),
  
  cindex_km_greedy = sapply(results, `[[`, "cindex_km_greedy"),
  ibs_km_greedy = sapply(results, `[[`, "ibs_km_greedy"),
  iauc_km_greedy = sapply(results, `[[`, "iauc_km_greedy"),
  aic_km_greedy = sapply(results, `[[`, "aic_km_greedy"),
  bic_km_greedy = sapply(results, `[[`, "bic_km_greedy"),
  aicc_km_greedy = sapply(results, `[[`, "aicc_km_greedy"),
  
  ibll_cont_only = sapply(results, `[[`, "ibll_cont_only"),
  ibll_ordinal   = sapply(results, `[[`, "ibll_ordinal"),
  ibll_ohe       = sapply(results, `[[`, "ibll_ohe"),
  # ibll_km_k10 = sapply(results, `[[`, "ibll_km_k10"),
  # ibll_km_k20 = sapply(results, `[[`, "ibll_km_k20"),
  ibll_km_greedy_separate = sapply(results, `[[`, "ibll_km_greedy_separate"),
  ibll_km_greedy = sapply(results, `[[`, "ibll_km_greedy")
)

save_metric_boxplot <- function(filename, metric_cols, metric_names, ylab_text) {
  keep <- sapply(metric_cols, function(col) {
    col %in% names(df_res) && any(!is.na(df_res[[col]]))
  })
  
  metric_cols  <- metric_cols[keep]
  metric_names <- metric_names[keep]
  
  save_pdf_plot(filename, {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    par(mar = c(12, 6, 3, 2), mgp = c(4.2, 1, 0))
    
    boxplot(
      df_res[, metric_cols, drop = FALSE],
      names = metric_names,
      ylab = ylab_text,
      las = 2,
      cex.axis = 1.4,
      cex.lab = 1.4
    )
  }, width = 11, height = 7)
}

metric_cols_common <- c(
  "cont_only",
  "ordinal",
  "ohe",
  # "km_k10",
  # "km_k20",
  "km_greedy_separate",
  "km_greedy"
)

metric_names_common <- c(
  "CONT ONLY",
  "ORD",
  "OH",
  # "KM fixed(K/10)",
  # "KM fixed(K/20)",
  "KM greedy separate",
  "KM"
)

save_metric_boxplot(
  sprintf("%s_boxplot_cindex_all_methods.pdf", dataset_name),
  paste0("cindex_", metric_cols_common),
  metric_names_common,
  "C-index"
)

save_metric_boxplot(
  sprintf("%s_boxplot_ibll_all_methods.pdf", dataset_name),
  paste0("ibll_", metric_cols_common),
  metric_names_common,
  "IBLL"
)

save_metric_boxplot(
  sprintf("%s_boxplot_ibs_all_methods.pdf", dataset_name),
  paste0("ibs_", metric_cols_common),
  metric_names_common,
  "IBS"
)

save_metric_boxplot(
  sprintf("%s_boxplot_iauc_all_methods.pdf", dataset_name),
  paste0("iauc_", metric_cols_common),
  metric_names_common,
  "iAUC"
)

save_metric_boxplot(
  sprintf("%s_boxplot_aic_all_methods.pdf", dataset_name),
  paste0("aic_", metric_cols_common),
  metric_names_common,
  "AIC"
)

save_metric_boxplot(
  sprintf("%s_boxplot_bic_all_methods.pdf", dataset_name),
  paste0("bic_", metric_cols_common),
  metric_names_common,
  "BIC"
)

save_metric_boxplot(
  sprintf("%s_boxplot_aicc_all_methods.pdf", dataset_name),
  paste0("aicc_", metric_cols_common),
  metric_names_common,
  "AICc"
)

collapse_times <- function(x) {
  if (length(x) == 0) return("")
  paste(round(x, 3), collapse = ", ")
}

if (length(cat_cols) > 0) {
  selected_times_wide <- do.call(
    rbind,
    lapply(seq_along(results), function(r) {
      row <- data.frame(rep = r)
      
      for (cc in cat_cols) {
        pretty_name <- safe_filename(pretty_feature_name(cc, reverse_cat_map))
        # row[[paste0("greedy_separate_", pretty_name)]] <- collapse_times(results[[r]]$selected_times_greedy_separate[[cc]])
        row[[paste0("greedy_global_", pretty_name)]]   <- collapse_times(results[[r]]$selected_times_greedy[[cc]])
      }
      
      row
    })
  )
  
  save_pdf_plot(
    sprintf("%s_counts_greedy_global_distribution.pdf", dataset_name),
    {
      counts_greedy <- do.call(
        rbind,
        lapply(results, function(res) {
          sapply(res$selected_times_greedy, length)
        })
      )
      
      counts_greedy <- as.data.frame(counts_greedy)
      
      par(mfrow = c(length(cat_cols), 1), mar = c(2, 4, 2, 1), oma = c(4, 0, 0, 0))
      
      for (cc in cat_cols) {
        tab <- table(counts_greedy[[cc]])
        barplot(
          tab,
          xlab = "",
          ylab = pretty_feature_name(cc, reverse_cat_map),
          main = "Global greedy selected time counts"
        )
      }
      
      mtext("Number of selected time points", side = 1, outer = TRUE, line = 2.5)
    },
    width = 8,
    height = max(6, 2.5 * length(cat_cols))
  )
  
  if (include_km_greedy_separate) {
    save_pdf_plot(
      sprintf("%s_counts_greedy_separate_distribution.pdf", dataset_name),
      {
        counts_sep <- do.call(
          rbind,
          lapply(results, function(res) {
            sapply(res$selected_times_greedy_separate, length)
          })
        )
        
        counts_sep <- as.data.frame(counts_sep)
        
        par(mfrow = c(length(cat_cols), 1), mar = c(2, 4, 2, 1), oma = c(4, 0, 0, 0))
        
        for (cc in cat_cols) {
          tab <- table(counts_sep[[cc]])
          barplot(
            tab,
            xlab = "",
            ylab = pretty_feature_name(cc, reverse_cat_map),
            main = "Separate greedy selected time counts"
          )
        }
        
        mtext("Number of selected time points", side = 1, outer = TRUE, line = 2.5)
      },
      width = 8,
      height = max(6, 2.5 * length(cat_cols))
    )
  }
} else {
  selected_times_wide <- data.frame(rep = seq_along(results))
}

write.csv(
  df_res,
  file = file.path(plot_dir, sprintf("results_summary_all_methods_%s.csv", dataset_name)),
  row.names = FALSE
)

write.csv(
  selected_times_wide,
  file = file.path(plot_dir, sprintf("selected_times_all_methods_wide_%s.csv", dataset_name)),
  row.names = FALSE
)

cat("\nAll PDFs and CSVs saved in folder:", normalizePath(plot_dir), "\n")

# ============================================================
# Summary table
# ============================================================
make_summary_table <- function(metrics, methods, method_labels) {
  rows <- lapply(seq_along(metrics), function(m) {
    metric <- metrics[m]
    fmt <- if (metric == "bic") "%.2f (%.2f, %.2f)" else "%.3f (%.3f, %.3f)"
    cols <- paste0(metric, "_", methods)
    sapply(cols, function(col) {
      x <- df_res[[col]]
      med <- median(x, na.rm = TRUE)
      q25 <- quantile(x, 0.25, na.rm = TRUE)
      q75 <- quantile(x, 0.75, na.rm = TRUE)
      sprintf(fmt, med, q25, q75)
    })
  })
  
  tab <- do.call(rbind, rows)
  rownames(tab) <- toupper(metrics)
  colnames(tab) <- method_labels
  as.data.frame(t(tab))
}

summary_tab <- make_summary_table(
  metrics       = c("cindex","iauc", "ibs", "ibll", "bic","aic"),
  methods       = metric_cols_common,
  method_labels = metric_names_common
)

summary_tab$Method <- rownames(summary_tab)
summary_tab <- summary_tab[, c("Method", "CINDEX","IAUC", "IBS", "IBLL", "BIC","AIC")]
colnames(summary_tab) <- c("Method","C-index", "iAUC", "IBS", "IBLL", "BIC","AIC")

write.csv2(
  summary_tab,
  file = file.path(plot_dir, sprintf("%s_summary.csv", dataset_name)),
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# Saving details about columns
# ============================================================
cont_names <- sapply(
  cont_cols,
  pretty_feature_name,
  original_name_lookup = reverse_cont_map
)

cat_names <- if (length(cat_cols) > 0) {
  sapply(cat_cols, function(cc) {
    nm <- pretty_feature_name(cc, reverse_cat_map)
    k  <- length(setdiff(unique(df[[cc]]), "NA"))
    paste0(nm, " (", k, ")")
  })
} else {
  character(0)
}

half <- ceiling(length(cat_names) / 2)

cat_1 <- if (length(cat_names) > 0) cat_names[1:half] else character(0)
cat_2 <- if (length(cat_names) > 0) cat_names[(half + 1):length(cat_names)] else character(0)

if (length(cat_2) < half) {
  cat_2 <- c(cat_2, rep("", half - length(cat_2)))
}

max_len <- max(length(cont_names), length(cat_1), 1)

cont_names <- c(cont_names, rep("", max_len - length(cont_names)))
cat_1 <- c(cat_1, rep("", max_len - length(cat_1)))
cat_2 <- c(cat_2, rep("", max_len - length(cat_2)))

feature_table <- data.frame(
  "Continuous covariates" = cont_names,
  "Categorical covariates" = cat_1,
  " " = cat_2,
  check.names = FALSE,
  row.names = NULL
)

latex_code <- kable(
  feature_table,
  format = "latex",
  booktabs = FALSE,
  align = "|l|ll",
  caption = paste("Summary of continuous and categorical covariates in", dataset_name),
  label = paste0("tab:features_", safe_filename(dataset_name)),
  row.names = FALSE
)

writeLines(latex_code, file.path(plot_dir, sprintf("feature_table_%s.tex", dataset_name)))

# ============================================================
# Gain of KM over OHE
# Positive gain = KM is better
# ============================================================

df_res$gain_cindex_km_vs_ohe <-
  df_res$cindex_km_greedy - df_res$cindex_ohe

df_res$gain_ibs_km_vs_ohe <-
  df_res$ibs_ohe - df_res$ibs_km_greedy

df_res$gain_ibll_km_vs_ohe <-
  df_res$ibll_km_greedy - df_res$ibll_ohe

save_gain_boxplot <- function(filename, gain_values, ylab_text) {
  
  save_pdf_plot(filename, {
    
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE)
    
    par(mar = c(5, 5, 3, 2))
    
    boxplot(
      gain_values,
      names = "KM vs OHE",
      ylab = ylab_text,
      cex.axis = 1.4,
      cex.lab = 1.4
    )
    
    abline(h = 0, lty = 2)
    
  }, width = 7, height = 7)
}

save_gain_boxplot(
  sprintf("%s_boxplot_gain_cindex_km_over_ohe.pdf", dataset_name),
  df_res$cindex_km_greedy - df_res$cindex_ohe,
  "C-index gain"
)

save_gain_boxplot(
  sprintf("%s_boxplot_gain_ibs_km_over_ohe.pdf", dataset_name),
  df_res$ibs_ohe - df_res$ibs_km_greedy,
  "IBS gain"
)

save_gain_boxplot(
  sprintf("%s_boxplot_gain_ibll_km_over_ohe.pdf", dataset_name),
  df_res$ibll_km_greedy - df_res$ibll_ohe,
  "IBLL gain"
)





#### after
for (cc in names(selected_idx_greedy)) {
  
  idx <- selected_idx_greedy[[cc]]
  
  if (length(idx) == 0) next
  
  cat("\n\n====================================\n")
  cat(
    "Feature:",
    pretty_feature_name(cc, reverse_cat_map),
    "\n"
  )
  cat("====================================\n")
  
  for (j in idx) {
    
    cat(
      "\nSelected time:",
      round(km_precomp[[cc]]$candidate_times[j], 3),
      "\n"
    )
    
    vals <- sapply(
      km_precomp[[cc]]$lookup,
      function(x) x[j]
    )
    
    print(
      data.frame(
        category = names(vals),
        KM_value = round(vals, 4),
        row.names = NULL
      )
    )
  }
}