suppressPackageStartupMessages({
  library(survival)
  # library(SurvMetrics)
  library(riskRegression)
  library(compareC)
  library(viridisLite)
})


#  
# Output folder
#  
plot_dir <- "plots_km_greedy_vs_ohe_sim"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

save_pdf_plot <- function(filename, expr, width = 10, height = 7) {
  pdf(file = file.path(plot_dir, filename), width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  eval.parent(substitute(expr))
}

#  
# Simulation settings
#  

R               <- 50
N_TRAIN_GRID    <- c(1000)
K_GRID          <- c(20)
CENSORING_GRID  <- c(0.4)

N_TEST          <- 5000
TIME_GRID_N     <- 100

# DGP parameters
LAM             <- 0.08
BETA_SD         <- 0.4
HEAD_PROP       <- 1
HEAD_FRAC       <- 1
ADD_CONT        <- TRUE   # include 2 continuous covariates

# Representative regime for Brier curve
REPREG_N        <- 500
REPREG_K        <- 25
REPREG_CENSOR   <- 0.5

#  
# Helpers
#  
cox_aic <- function(cox_fit, n_for_bic_aicc) {
  ll    <- as.numeric(cox_fit$loglik[2])
  beta  <- coef(cox_fit)
  k     <- sum(is.finite(beta) & !is.na(beta))
  aic   <- -2 * ll + 2 * k
  bic   <- -2 * ll + log(n_for_bic_aicc) * k
  denom <- max(1e-9, n_for_bic_aicc - k - 1)
  aicc  <- aic + (2 * k * (k + 1)) / denom
  c(aic = aic, bic = bic, aicc = aicc, k = k)
}

cox_aic_value <- function(cox_fit, n_for_bic_aicc) {
  as.numeric(cox_aic(cox_fit, n_for_bic_aicc)["aic"])
}

cindex_like_lifelines <- function(time, event, risk_lp) {
  as.numeric(concordance(Surv(time, event) ~ I(-risk_lp))$concordance)
}

iauc_weights_fS <- function(time, event, eval_times) {
  eval_times <- eval_times[is.finite(eval_times)]
  if (length(eval_times) < 2) return(rep(NA_real_, length(eval_times)))
  sf      <- survfit(Surv(time, event) ~ 1)
  S_right <- as.numeric(summary(sf, times = eval_times, extend = TRUE)$surv)
  S_left  <- c(1, head(S_right, -1))
  dF      <- pmax(0, S_left - S_right)
  w_raw   <- dF * S_right
  if (sum(w_raw, na.rm = TRUE) <= 0) return(rep(NA_real_, length(eval_times)))
  w_raw / sum(w_raw, na.rm = TRUE)
}

iauc_fS <- function(times, auc, time, event) {
  ok    <- is.finite(times) & is.finite(auc)
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
  
  auc_values <- as.data.frame(b2$AUC$score)
  ibs        <- as.numeric(tail(b2$Brier$score$IBS, 1))
  lp         <- as.numeric(predict(cox_fit, newdata = df_test, type = "lp"))
  cindex     <- cindex_like_lifelines(df_test$Time, df_test$Event, lp)
  iauc       <- iauc_fS(auc_values$times, auc_values$AUC, df_test$Time, df_test$Event)
  ic         <- cox_aic(cox_fit, max(2L, as.integer(sum(df_train$Event))))
  
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
    ibs = ibs,
    cindex = cindex,
    iauc = iauc,
    ibll = ibll,
    aic = as.numeric(ic["aic"]),
    bic = as.numeric(ic["bic"]),
    aicc = as.numeric(ic["aicc"]),
    k = as.numeric(ic["k"]),
    score_obj = b2,
    bll_values = bll_values
  )
}

safe_eval_cox <- function(fit_obj, df_train, df_test, time_interest) {
  out <- list(
    fit_success = FALSE,
    score_success = FALSE,
    has_na_coef = NA,
    ibs = NA_real_,
    cindex = NA_real_,
    iauc = NA_real_,
    ibll = NA_real_,
    aic = NA_real_,
    bic = NA_real_,
    aicc = NA_real_,
    k = NA_real_,
    score_obj = NULL,
    bll_values = NULL
  )
  
  if (inherits(fit_obj, "try-error") || is.null(fit_obj)) return(out)
  if (is.null(df_train) || is.null(df_test)) return(out)
  
  out$fit_success <- TRUE
  
  cf <- try(coef(fit_obj), silent = TRUE)
  if (!inherits(cf, "try-error")) out$has_na_coef <- anyNA(cf)
  
  ev <- try(evaluate_cox_model(fit_obj, df_train, df_test, time_interest), silent = TRUE)
  if (!inherits(ev, "try-error")) {
    out$score_success <- TRUE
    out$ibs        <- ev$ibs
    out$cindex     <- ev$cindex
    out$iauc       <- ev$iauc
    out$ibll       <- ev$ibll
    out$aic        <- ev$aic
    out$bic        <- ev$bic
    out$aicc       <- ev$aicc
    out$k          <- ev$k
    out$score_obj  <- ev$score_obj
    out$bll_values <- ev$bll_values
  }
  
  out
}

#  
# KM jump diagnostics
#  
km_jumps_per_category <- function(df, cat_col) {
  cats <- unique(as.character(df[[cat_col]]))
  
  out <- lapply(cats, function(lv) {
    sub <- df[as.character(df[[cat_col]]) == lv, , drop = FALSE]
    event_times <- unique(sub$Time[sub$Event == 1])
    
    data.frame(
      category = lv,
      n = nrow(sub),
      events = sum(sub$Event),
      jumps = length(event_times)
    )
  })
  
  do.call(rbind, out)
}

#  
# DGP: one sparse, high-cardinality, imbalanced categorical feature
#  
category_labels_z1 <- function(K) {
  paste0("C_{1,", seq_len(K), "}")
}

make_probs_long_tail <- function(K, head_prop = 0.7, head_frac = 0.2) {
  n_head <- max(2, floor(K * head_frac))
  n_tail <- K - n_head
  
  p_head <- rep(head_prop / n_head, n_head)
  
  if (n_tail > 0) {
    tail_raw <- 1 / (seq_len(n_tail)^1.2)
    tail_raw <- tail_raw / sum(tail_raw)
    p_tail <- (1 - head_prop) * tail_raw
    p <- c(p_head, p_tail)
  } else {
    p <- p_head
  }
  
  p / sum(p)
}

generate_dataset_onecat <- function(n,
                                    K,
                                    lam = 0.08,
                                    censor_prop = 0.4,
                                    beta_sd = 0.6,
                                    beta_g = NULL,
                                    head_prop = 0.5,
                                    head_frac = 0.4,
                                    add_cont = TRUE) {
  probs <- make_probs_long_tail(K = K, head_prop = head_prop, head_frac = head_frac)
  g <- sample(seq_len(K), size = n, replace = TRUE, prob = probs)
  
  if (is.null(beta_g)) {
    beta_g <- rnorm(K, mean = 0, sd = beta_sd)
    beta_g <- beta_g - mean(beta_g)
  }
  
  cont_1 <- rnorm(n)
  cont_2 <- rnorm(n)
  
  lp <- beta_g[g]
  if (add_cont) {
    lp <- lp + 0.5 * cont_1 - 0.4 * cont_2
  }
  
  rate <- lam * exp(lp)
  T_event <- rexp(n, rate = rate)
  
  if (censor_prop <= 0) {
    T_cens <- rep(Inf, n)
  } else {
    lambda_c <- (censor_prop / (1 - censor_prop)) * lam
    T_cens <- rexp(n, rate = lambda_c)
  }
  
  data.frame(
    cont_1 = cont_1,
    cont_2 = cont_2,
    cat_1  = category_labels_z1(K)[g],
    Time   = pmin(T_event, T_cens),
    Event  = as.integer(T_event <= T_cens)
  )
}

#  
# One-hot encoding
#  
ohe <- function(cat_cols, train_df_noz, test_df_noz) {
  all_encoded_train       <- train_df_noz
  all_encoded_test        <- test_df_noz
  all_censored_categories <- list()
  
  for (cc in cat_cols) {
    train_cats <- as.character(train_df_noz[[cc]])
    test_cats  <- as.character(test_df_noz[[cc]])
    all_levels <- unique(train_cats)
    
    has_event       <- sapply(all_levels, function(lv) any(train_df_noz$Event[train_cats == lv] == 1))
    valid_levels    <- all_levels[has_event]
    censored_levels <- all_levels[!has_event]
    
    K_valid <- length(valid_levels)
    if (K_valid < 2) stop(sprintf("Column '%s' has fewer than 2 non-fully-censored categories.", cc))
    
    baseline    <- valid_levels[1]
    kept_levels <- setdiff(valid_levels, baseline)
    K_minus_1   <- length(kept_levels)
    col_names   <- make.names(paste0(cc, "_", kept_levels))
    
    encode_cats <- function(cats) {
      mat <- matrix(0, nrow = length(cats), ncol = K_minus_1,
                    dimnames = list(NULL, col_names))
      for (i in seq_along(cats)) {
        lv <- cats[i]
        if (lv %in% kept_levels) {
          mat[i, make.names(paste0(cc, "_", lv))] <- 1
        }
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
  
  list(train = all_encoded_train, test = all_encoded_test, cat_censored = all_censored_categories)
}

#  
# KM greedy helpers
#  
predict_surv_at_times <- function(sf, times) {
  as.numeric(summary(sf, times = times, extend = TRUE)$surv)
}

precompute_km_all_features <- function(train_df, cat_cols, candidate_times) {
  fallback_km  <- survfit(Surv(Time, Event) ~ 1, data = train_df)
  fallback_vec <- predict_surv_at_times(fallback_km, candidate_times)
  
  lapply(setNames(cat_cols, cat_cols), function(cc) {
    cats_train   <- as.character(train_df[[cc]])
    levels_train <- unique(cats_train)
    km_lookup    <- lapply(setNames(levels_train, levels_train), function(lv) {
      sub <- train_df[cats_train == lv, , drop = FALSE]
      if (nrow(sub) == 0 || sum(sub$Event) == 0) {
        fallback_vec
      } else {
        predict_surv_at_times(survfit(Surv(Time, Event) ~ 1, data = sub), candidate_times)
      }
    })
    list(cat_col = cc, candidate_times = candidate_times, lookup = km_lookup, fallback = fallback_vec)
  })
}

build_full_candidate_matrices <- function(df, precomp_all, cat_cols) {
  lapply(setNames(cat_cols, cat_cols), function(cc) {
    pc   <- precomp_all[[cc]]
    cats <- as.character(df[[cc]])
    uniq <- unique(cats)
    row_lookup <- lapply(setNames(uniq, uniq), function(lv) {
      v <- pc$lookup[[lv]]
      if (is.null(v)) v <- pc$fallback
      v
    })
    x <- do.call(rbind, row_lookup[cats])
    colnames(x) <- paste0(cc, "_km_t", seq_len(ncol(x)))
    as.matrix(x)
  })
}

assemble_selected_km_matrix <- function(full_mats, selected_idx, cat_cols) {
  parts <- Filter(Negate(is.null), lapply(cat_cols, function(cc) {
    idx <- selected_idx[[cc]]
    if (length(idx) > 0) full_mats[[cc]][, idx, drop = FALSE] else NULL
  }))
  if (length(parts) == 0) {
    matrix(nrow = nrow(full_mats[[cat_cols[1]]]), ncol = 0)
  } else {
    do.call(cbind, parts)
  }
}

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
    coxph(Surv(Time, Event) ~ ., data = base_df, ties = "efron", singular.ok = FALSE),
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
        x_new   <- full_mats_train[[cc]][, j, drop = FALSE]
        X_trial <- if (ncol(X_current) == 0) x_new else cbind(X_current, x_new)
        
        df_trial <- cbind(base_df, as.data.frame(X_trial))
        
        fit_trial <- try(
          coxph(Surv(Time, Event) ~ ., data = df_trial, ties = "efron", singular.ok = FALSE),
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
# Category-distribution plotting helper
#  
plot_category_distribution <- function(df_sub, filename) {
  save_pdf_plot(filename, {
    par(mar = c(6, 5, 4, 2))
    
    K_here <- unique(df_sub$K)
    n_head <- max(2, floor(HEAD_FRAC * K_here))
    cols <- ifelse(df_sub$category_index <= n_head, "grey40", "grey75")
    
    barplot(
      df_sub$observed_prop,
      names.arg = df_sub$category_index,
      col = cols,
      border = NA,
      xlab = expression("Category index " * i * " in " * C[1*i]),
      ylab = "Frequency (proportion)",
      las = 2,
      cex.names = 0.65
    )
    
    legend(
      "topright",
      legend = c("Head categories", "Tail categories"),
      fill = c("grey40", "grey75"),
      border = NA,
      bty = "n"
    )
  }, width = 11, height = 6)
}


### Intepretability plots
plot_pca_km_encoding_by_beta <- function(X_train_grd, beta_g, cat_vec, K,
                                         n_train, cens, r, plot_dir) {
  if (ncol(X_train_grd) < 2) {
    message("Skipping PCA: fewer than 2 KM greedy features selected.")
    return(invisible(NULL))
  }
  
  all_labs <- category_labels_z1(K)
  
  # one row per category: mean encoded KM features
  cat_means <- aggregate(X_train_grd, by = list(cat = cat_vec), FUN = mean)
  rownames(cat_means) <- cat_means$cat
  cat_means$cat <- NULL
  cat_means <- as.matrix(cat_means)
  
  pca_res <- prcomp(cat_means, center = TRUE, scale. = TRUE)
  scores  <- pca_res$x[, 1:2, drop = FALSE]
  
  
  var_exp <- summary(pca_res)$importance[2, 1:2] * 100
  
  cat_index <- match(rownames(cat_means), all_labs)
  obs_beta  <- beta_g[cat_index]
  
  cor_pc1 <- cor(scores[,1], obs_beta)
  
  fname <- sprintf("pca_km_encoding_true_beta_K%d_n%d_c%s_r%d.pdf",
                   K, n_train, gsub("\\.", "", as.character(cens)), r)
  
  pdf(file.path(plot_dir, fname), width = 7, height = 6)
  
  pal <- viridisLite::viridis(100)
  
  z <- (obs_beta - min(obs_beta, na.rm = TRUE)) /
    (max(obs_beta, na.rm = TRUE) - min(obs_beta, na.rm = TRUE))
  
  cols <- pal[pmax(1, pmin(100, round(1 + 99 * z)))]
  
  
  plot(scores[, 1], scores[, 2],
       col = cols,
       pch = 16,
       cex = 0.9,
       xlab = sprintf("PC1 (%.1f%% var)", var_exp[1]),
       ylab = sprintf("PC2 (%.1f%% var)", var_exp[2]),
       main = "")
  
  # Draw a horizontal line at y = 0
  abline(h = 0, col = "grey", lty = 2)
  
  # Draw a vertical line at x = 0
  abline(v = 0, col = "grey", lty = 2)
  
  text(scores[,1], scores[,2],
       labels = cat_index,
       cex = 0.6,
       pos = 3)
  
  legend("topright",
         legend = c(
           expression("low " * beta),
           expression("high " * beta)
         ),
         col = c(viridisLite::viridis(100)[1],
                 viridisLite::viridis(100)[100]),
         pch = 16,
         bty = "n")
  
  dev.off()
  
  fname2 <- sprintf(
    "pc1_vs_beta_K%d_n%d_c%s_r%d.pdf",
    K, n_train,
    gsub("\\.", "", as.character(cens)),
    r
  )
  
  pdf(file.path(plot_dir, fname2), width = 6, height = 5)
  
  plot(scores[,1],
       obs_beta,
       pch = 16,
       # col = cols,
       xlab = sprintf("PC1 (%.1f%% variance)", var_exp[1]),
       ylab = expression(True~beta[k]),
       main = sprintf("Correlation = %.3f", cor_pc1))
  
  abline(
    lm(obs_beta ~ scores[,1]),
    col = "red",
    lwd = 2
  )
  
  dev.off()
  invisible(pca_res)
  
  cor_pc1
}


plot_ohe_coef_vs_beta <- function(cox_ohe, beta_g, K,
                                  n_train, cens, r, plot_dir) {
  
  coef_hat <- rep(0, K)
  
  cf <- coef(cox_ohe)
  
  coef_names <- names(cf)
  
  for (k in 2:K) {
    
    nm <- make.names(
      paste0("cat_1_", category_labels_z1(K)[k])
    )
    
    if (nm %in% coef_names) {
      coef_hat[k] <- cf[nm]
    }
  }
  
  coef_hat <- coef_hat - mean(coef_hat)
  beta_g_c <- beta_g - mean(beta_g)
  
  cor_beta <- cor(coef_hat, beta_g_c)
  
  fname <- sprintf(
    "ohe_coef_vs_true_beta_K%d_n%d_c%s_r%d.pdf",
    K,
    n_train,
    gsub("\\.", "", as.character(cens)),
    r
  )
  
  pdf(file.path(plot_dir, fname), width = 6, height = 5)
  
  plot(
    coef_hat,
    beta_g_c,
    pch = 16,
    xlab = expression(hat(beta)[OHE]),
    ylab = expression(True~beta[k]),
    main = sprintf("Correlation = %.3f", cor_beta)
  )
  
  abline(
    lm(beta_g_c ~ coef_hat),
    col = "red",
    lwd = 2
  )
  
  dev.off()
  
  invisible(cor_beta)
  
  cor_beta
}

plot_true_beta <- function(beta_g, K, n_train, cens, plot_dir) {
  
  fname <- sprintf(
    "true_beta_structure_K%d_n%d_c%s.pdf",
    K, n_train,
    gsub("\\.", "", as.character(cens))
  )
  
  pdf(file.path(plot_dir, fname), width = 7, height = 4)
  
  plot(seq_len(K),
       beta_g,
       type = "b",
       pch = 16,
       xlab = "Category index",
       ylab = expression(beta[k]))
  
  dev.off()
}

#  
# Main simulation
#  
all_results <- list()
idx <- 0L
rep_curve_store <- list()
best_rep_gain <- -Inf
compareC_store <- list()
compare_idx <- 0L
cat_dist_store <- list()

selected_time_store <- list()
selected_idx_store <- 0L

for (n_train in N_TRAIN_GRID) {
  for (K in K_GRID) {
    i_K <- seq_len(K)
    
    beta_g_K <-
      0.7 * sin(2 * pi * i_K / K) +
      0.35 * cos(4 * pi * i_K / K) +
      rnorm(K, mean = 0, sd = 0.15)
    
    
    beta_g_K <- beta_g_K - mean(beta_g_K)
    beta_g_K <- beta_g_K / sd(beta_g_K) * BETA_SD
    
    for (cens in CENSORING_GRID) {
      cat("Running: n =", n_train, "| K =", K, "| censor =", cens, "\n")
      
      store_correlation <- numeric(0)
      store_correlation_ohe <- numeric(0)
      
      test_df_noz <- generate_dataset_onecat(
        n = N_TEST,
        K = K,
        beta_g = beta_g_K,
        lam = LAM,
        censor_prop = cens,
        beta_sd = BETA_SD,
        head_prop = HEAD_PROP,
        head_frac = HEAD_FRAC,
        add_cont = ADD_CONT
      )
      
      for (r in seq_len(R)) {
        train_df_noz <- generate_dataset_onecat(
          n = n_train,
          K = K,
          lam = LAM,
          censor_prop = cens,
          beta_sd = BETA_SD,
          beta_g = beta_g_K,
          head_prop = HEAD_PROP,
          head_frac = HEAD_FRAC,
          add_cont = ADD_CONT
        )
        
        if (r == 1) {
          all_labs <- category_labels_z1(K)
          tab <- table(factor(train_df_noz$cat_1, levels = all_labs))
          true_probs <- make_probs_long_tail(
            K = K,
            head_prop = HEAD_PROP,
            head_frac = HEAD_FRAC
          )
          
          plot_true_beta(
            beta_g_K,
            K,
            n_train,
            cens,
            plot_dir
          )
          
          cat_dist_store[[paste0("K", K, "_c", cens)]] <- data.frame(
            n_train = n_train,
            K = K,
            censor_prop = cens,
            category_index = seq_len(K),
            category = all_labs,
            count = as.numeric(tab),
            observed_prop = as.numeric(tab) / n_train,
            true_prob = true_probs,
            expected_count = n_train * true_probs
          )
        }
        
        cont_cols <- grep("^cont_", names(train_df_noz), value = TRUE)
        cat_cols  <- "cat_1"
        
        time_interest <- seq(
          quantile(test_df_noz$Time, 0.1),
          quantile(test_df_noz$Time, 0.9),
          length.out = TIME_GRID_N
        )
        
        # OHE
        ohe_result   <- tryCatch(ohe(cat_cols, train_df_noz, test_df_noz), error = function(e) NULL)
        df_train_ohe <- if (!is.null(ohe_result)) ohe_result$train else NULL
        df_test_ohe  <- if (!is.null(ohe_result)) ohe_result$test  else NULL
        
        cox_ohe <- if (!is.null(df_train_ohe)) {
          tryCatch(
            coxph(Surv(Time, Event) ~ ., data = df_train_ohe,
                  ties = "efron", singular.ok = FALSE,
                  x = TRUE, y = TRUE, model = TRUE),
            error = function(e) NULL
          )
        } else NULL
        
        eval_ohe <- safe_eval_cox(cox_ohe, df_train_ohe, df_test_ohe, time_interest)
        
        #  KM greedy
        event_times <- sort(unique(train_df_noz$Time[train_df_noz$Event == 1]))
        candidate_times <- sample(unique(as.numeric(quantile(
          event_times, probs = seq(0.1, 0.9, length.out = 70), na.rm = TRUE
        ))))
        candidate_times <- candidate_times[is.finite(candidate_times)]
        
        if (length(candidate_times) == 0) {
          candidate_times <- unique(as.numeric(quantile(
            train_df_noz$Time, probs = seq(0.1, 0.9, length.out = 70), na.rm = TRUE
          )))
          candidate_times <- candidate_times[is.finite(candidate_times)]
        }
        
        km_precomp      <- precompute_km_all_features(train_df_noz, cat_cols, candidate_times)
        full_mats_train <- build_full_candidate_matrices(train_df_noz, km_precomp, cat_cols)
        full_mats_test  <- build_full_candidate_matrices(test_df_noz,  km_precomp, cat_cols)
        
        greedy_out <- km_greedy_select_times(
          train_df_noz       = train_df_noz,
          cont_cols          = cont_cols,
          cat_cols           = cat_cols,
          candidate_times    = candidate_times,
          candidate_pool_idx = seq_along(candidate_times),
          full_mats_train    = full_mats_train
        )
        
        selected_times_here <- unlist(greedy_out$selected_times)
        
        if (length(selected_times_here) > 0) {
          selected_idx_store <- selected_idx_store + 1L
          
          selected_time_store[[selected_idx_store]] <-
            cbind(
              data.frame(
                rep = r,
                n_train = n_train,
                K = K,
                censor_prop = cens
              ),
              greedy_out$selection_log
            )
          
        }
        
        
        
        X_train_grd <- assemble_selected_km_matrix(full_mats_train, greedy_out$selected_idx, cat_cols)
        X_test_grd  <- assemble_selected_km_matrix(full_mats_test,  greedy_out$selected_idx, cat_cols)
        
        if (ncol(X_train_grd) >= 2) {
          cor_pc1 = plot_pca_km_encoding_by_beta(
            X_train_grd = X_train_grd,
            beta_g      = beta_g_K,
            cat_vec     = as.character(train_df_noz$cat_1),
            K           = K,
            n_train     = n_train,
            cens        = cens,
            r           = r,
            plot_dir    = plot_dir
          )
          
          store_correlation <- c(store_correlation, cor_pc1)
        }
        
        if (!is.null(cox_ohe)) {
          
          cor_ohe <- plot_ohe_coef_vs_beta(
            cox_ohe   = cox_ohe,
            beta_g    = beta_g_K,
            K         = K,
            n_train   = n_train,
            cens      = cens,
            r         = r,
            plot_dir  = plot_dir
          )
          
          store_correlation_ohe <- c(store_correlation_ohe, cor_ohe)
          
        }
        
        df_train_grd <- cbind(
          train_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE],
          as.data.frame(X_train_grd)
        )
        
        df_test_grd <- cbind(
          test_df_noz[, c(cont_cols, "Time", "Event"), drop = FALSE],
          as.data.frame(X_test_grd)
        )
        
        cox_grd <- tryCatch(
          coxph(Surv(Time, Event) ~ ., data = df_train_grd,
                ties = "efron", singular.ok = FALSE,
                x = TRUE, y = TRUE, model = TRUE),
          error = function(e) NULL
        )
        
        eval_grd <- safe_eval_cox(cox_grd, df_train_grd, df_test_grd, time_interest)
        
        grp_sizes <- table(train_df_noz$cat_1)
        km_jump_df <- km_jumps_per_category(train_df_noz, "cat_1")
        
        mean_km_jumps   <- mean(km_jump_df$jumps, na.rm = TRUE)
        median_km_jumps <- median(km_jump_df$jumps, na.rm = TRUE)
        min_km_jumps    <- min(km_jump_df$jumps, na.rm = TRUE)
        max_km_jumps    <- max(km_jump_df$jumps, na.rm = TRUE)
        total_km_jumps  <- sum(km_jump_df$jumps, na.rm = TRUE)
        
        mean_events_per_cat   <- mean(km_jump_df$events, na.rm = TRUE)
        median_events_per_cat <- median(km_jump_df$events, na.rm = TRUE)
        
        lp_ohe <- try(
          as.numeric(predict(cox_ohe, newdata = df_test_ohe, type = "lp")),
          silent = TRUE
        )
        
        lp_km <- try(
          as.numeric(predict(cox_grd, newdata = df_test_grd, type = "lp")),
          silent = TRUE
        )
        
        if (!inherits(lp_ohe, "try-error") &&
            !inherits(lp_km, "try-error")) {
          compare_idx <- compare_idx + 1L
          
          compareC_store[[compare_idx]] <- data.frame(
            rep = r,
            n_train = n_train,
            K = K,
            censor_prop = cens,
            Time = df_test_ohe$Time,
            Event = df_test_ohe$Event,
            lp_ohe = lp_ohe,
            lp_km = lp_km
          )
        }
        
        idx <- idx + 1L
        all_results[[idx]] <- data.frame(
          rep = r,
          n_train = n_train,
          K = K,
          censor_prop = cens,
          events_train = sum(train_df_noz$Event),
          actual_event_frac = mean(train_df_noz$Event),
          median_group_size = median(as.numeric(grp_sizes)),
          min_group_size = min(as.numeric(grp_sizes)),
          n_over_k = n_train / K,
          
          mean_km_jumps = mean_km_jumps,
          median_km_jumps = median_km_jumps,
          min_km_jumps = min_km_jumps,
          max_km_jumps = max_km_jumps,
          total_km_jumps = total_km_jumps,
          mean_events_per_cat = mean_events_per_cat,
          median_events_per_cat = median_events_per_cat,
          
          ohe_fit_success   = eval_ohe$fit_success,
          ohe_score_success = eval_ohe$score_success,
          ohe_has_na_coef   = eval_ohe$has_na_coef,
          
          km_fit_success    = eval_grd$fit_success,
          km_score_success  = eval_grd$score_success,
          km_has_na_coef    = eval_grd$has_na_coef,
          
          ibs_ohe = eval_ohe$ibs,
          ibll_ohe = eval_ohe$ibll,
          cindex_ohe = eval_ohe$cindex,
          iauc_ohe = eval_ohe$iauc,
          aic_ohe = eval_ohe$aic,
          bic_ohe = eval_ohe$bic,
          aicc_ohe = eval_ohe$aicc,
          k_ohe = eval_ohe$k,
          
          ibs_km_greedy = eval_grd$ibs,
          ibll_km_greedy = eval_grd$ibll,
          cindex_km_greedy = eval_grd$cindex,
          iauc_km_greedy = eval_grd$iauc,
          aic_km_greedy = eval_grd$aic,
          bic_km_greedy = eval_grd$bic,
          aicc_km_greedy = eval_grd$aicc,
          k_km_greedy = eval_grd$k,
          
          n_selected_total = sum(lengths(greedy_out$selected_times))
        )
        
        if (n_train == REPREG_N && K == REPREG_K && cens == REPREG_CENSOR) {
          ibs_gain_here <- eval_ohe$ibs - eval_grd$ibs
          if (is.finite(ibs_gain_here) && ibs_gain_here > best_rep_gain &&
              !is.null(eval_ohe$score_obj) && !is.null(eval_grd$score_obj)) {
            best_rep_gain <- ibs_gain_here
            rep_curve_store <- list(
              score_ohe = eval_ohe$score_obj,
              score_km_greedy = eval_grd$score_obj,
              bll_ohe = eval_ohe$bll_values,
              bll_km_greedy = eval_grd$bll_values
            )
          }
        }
      }
    }
  }
}

df_res <- do.call(rbind, all_results)

df_selected_times <- do.call(rbind, selected_time_store)

write.csv(
  df_selected_times,
  file.path(plot_dir, "selected_times.csv"),
  row.names = FALSE
)

df_res$ohe_fit_failure   <- 1 - as.numeric(df_res$ohe_fit_success)
df_res$ohe_score_failure <- 1 - as.numeric(df_res$ohe_score_success)
df_res$km_fit_failure    <- 1 - as.numeric(df_res$km_fit_success)
df_res$km_score_failure  <- 1 - as.numeric(df_res$km_score_success)

df_res$ibs_gain    <- df_res$ibs_ohe - df_res$ibs_km_greedy
# BLL / IBLL are log-likelihood style quantities, so higher / less negative is better.
df_res$ibll_gain   <- df_res$ibll_km_greedy - df_res$ibll_ohe
df_res$cindex_gain <- df_res$cindex_km_greedy - df_res$cindex_ohe
df_res$iauc_gain   <- df_res$iauc_km_greedy - df_res$iauc_ohe
df_res$aic_gain    <- df_res$aic_ohe - df_res$aic_km_greedy
df_res$bic_gain    <- df_res$bic_ohe - df_res$bic_km_greedy
df_res$aicc_gain   <- df_res$aicc_ohe - df_res$aicc_km_greedy
df_res$k_reduction <- df_res$k_ohe - df_res$k_km_greedy


df_res$selected_per_mean_jump <- df_res$n_selected_total / df_res$mean_km_jumps
df_res$selected_per_total_jump <- df_res$n_selected_total / df_res$total_km_jumps

df_res$selected_per_mean_jump[!is.finite(df_res$selected_per_mean_jump)] <- NA_real_
df_res$selected_per_total_jump[!is.finite(df_res$selected_per_total_jump)] <- NA_real_

saveRDS(df_res, file.path(plot_dir, "df_res_main.rds"))
write.csv(df_res, file.path(plot_dir, "df_res_main.csv"), row.names = FALSE)

# Save one-replicate category distribution data for each (K, censoring) pair.
df_cat_dist <- do.call(rbind, cat_dist_store)
write.csv(
  df_cat_dist,
  file.path(plot_dir, "category_distribution_one_rep_per_K_censor.csv"),
  row.names = FALSE
)

failure_summary <- aggregate(
  cbind(
    ohe_fit_success,
    ohe_score_success,
    ohe_has_na_coef,
    km_fit_success,
    km_score_success,
    km_has_na_coef
  ) ~ n_train + K + censor_prop,
  data = df_res,
  FUN = function(x) mean(x, na.rm = TRUE)
)

write.csv(failure_summary, file.path(plot_dir, "failure_summary.csv"), row.names = FALSE)

#  
# Plots
#  

# Category imbalance plots: one training sample replicate for each (K, censoring) pair.
for (nm in names(cat_dist_store)) {
  plot_category_distribution(
    cat_dist_store[[nm]],
    paste0("category_distribution_", nm, ".pdf")
  )
}

#  
# Mean / median gain vs K by censoring
#  

metric_summary <- function(df, gain_col) {
  do.call(rbind, lapply(
    split(df, list(df$K, df$censor_prop), drop = TRUE),
    function(sub) {
      x <- sub[[gain_col]]
      x <- x[is.finite(x)]
      if (length(x) < 2) return(NULL)
      
      se <- sd(x) / sqrt(length(x))
      ci <- mean(x) + c(-1, 1) * qt(0.975, df = length(x) - 1) * se
      
      data.frame(
        K = unique(sub$K),
        censor_prop = unique(sub$censor_prop),
        n = length(x),
        mean = mean(x),
        median = median(x),
        ci_low = ci[1],
        ci_high = ci[2]
      )
    }
  ))
}

plot_gain_vs_K_by_censoring <- function(summary_df, value_col, ylab_text, title_text, filename) {
  save_pdf_plot(filename, {
    cps <- sort(unique(summary_df$censor_prop))
    cols <- seq_along(cps)
    
    ylim <- range(summary_df[[value_col]], na.rm = TRUE)
    ylim <- ylim + c(-1, 1) * diff(ylim) * 0.08
    
    plot(NA,
         xlim = range(summary_df$K),
         ylim = ylim,
         xlab = "Number of categories K",
         ylab = ylab_text,
         main = title_text)
    
    abline(h = 0, lty = 2, col = "grey50")
    
    for (i in seq_along(cps)) {
      cp <- cps[i]
      sub <- summary_df[summary_df$censor_prop == cp, ]
      sub <- sub[order(sub$K), ]
      
      lines(sub$K, sub[[value_col]], type = "b",
            pch = 16, lwd = 2, col = cols[i])
    }
    
    legend("topleft",
           legend = paste("censor =", cps),
           col = cols,
           pch = 16,
           lwd = 2,
           bty = "n")
  })
}

plot_gain_ci_ribbons <- function(summary_df, ylab_text, title_text, filename) {
  save_pdf_plot(filename, {
    cps <- sort(unique(summary_df$censor_prop))
    cols <- seq_along(cps)
    
    ylim <- range(c(summary_df$ci_low, summary_df$ci_high), na.rm = TRUE)
    ylim <- ylim + c(-1, 1) * diff(ylim) * 0.08
    
    plot(NA,
         xlim = range(summary_df$K),
         ylim = ylim,
         xlab = "K",
         ylab = ylab_text,
         main = title_text)
    
    abline(h = 0, lty = 2, col = "grey50")
    
    for (i in seq_along(cps)) {
      cp <- cps[i]
      sub <- summary_df[summary_df$censor_prop == cp, ]
      sub <- sub[order(sub$K), ]
      
      polygon(
        c(sub$K, rev(sub$K)),
        c(sub$ci_high, rev(sub$ci_low)),
        col = adjustcolor(cols[i], alpha.f = 0.15),
        border = NA
      )
      
      lines(sub$K, sub$mean, type = "b",
            pch = 16, lwd = 2, col = cols[i])
    }
    
    legend("topleft",
           legend = paste("censor =", cps),
           col = cols,
           pch = 16,
           lwd = 2,
           bty = "n")
  })
}

save_pdf_plot("failure_rate_by_censoring.pdf", {
  agg_fail <- aggregate(
    cbind(ohe_score_failure, km_score_failure) ~ censor_prop,
    data = df_res,
    FUN = mean, na.rm = TRUE
  )
  ylim <- range(c(agg_fail$ohe_score_failure, agg_fail$km_score_failure), na.rm = TRUE)
  plot(agg_fail$censor_prop, agg_fail$ohe_score_failure,
       type = "b", pch = 16, lwd = 2, col = "#E41A1C",
       ylim = ylim,
       xlab = "Training censoring proportion",
       ylab = "Scoring failure rate")
  lines(agg_fail$censor_prop, agg_fail$km_score_failure,
        type = "b", pch = 17, lwd = 2, col = "#377EB8", lty = 2)
  legend("topleft", legend = c("OH", "KM"),
         col = c("#E41A1C", "#377EB8"),
         pch = c(16, 17), lty = c(1, 2), bty = "n")
})

save_pdf_plot("ibs_vs_n_over_k.pdf", {
  plot(df_res$n_over_k, df_res$ibs_gain,
       pch = 16, col = adjustcolor("steelblue", 0.35),
       xlab = "n / K",
       ylab = "IBS gain (OH - KM)")
  abline(h = 0, lty = 2, col = "grey50")
  lines(lowess(df_res$n_over_k, df_res$ibs_gain, f = 0.5), lwd = 2)
})

save_pdf_plot("ibll_vs_n_over_k.pdf", {
  plot(df_res$n_over_k, df_res$ibll_gain,
       pch = 16, col = adjustcolor("steelblue", 0.35),
       xlab = "n / K",
       ylab = "IBLL gain (KM - OH)")
  abline(h = 0, lty = 2, col = "grey50")
  lines(lowess(df_res$n_over_k, df_res$ibll_gain, f = 0.5), lwd = 2)
})

save_pdf_plot("cindex_vs_n_over_k.pdf", {
  plot(df_res$n_over_k, df_res$cindex_gain,
       pch = 16, col = adjustcolor("darkgreen", 0.35),
       xlab = "n / K",
       ylab = "C-index gain (KM - OH)")
  abline(h = 0, lty = 2, col = "grey50")
  lines(lowess(df_res$n_over_k, df_res$cindex_gain, f = 0.5), lwd = 2)
})

agg_aic <- aggregate(cbind(aic_ohe, aic_km_greedy) ~ K, data = df_res, FUN = mean, na.rm = TRUE)
save_pdf_plot("aic_vs_K.pdf", {
  ylim <- range(c(agg_aic$aic_ohe, agg_aic$aic_km_greedy), na.rm = TRUE)
  plot(agg_aic$K, agg_aic$aic_ohe, type = "b", pch = 16, lwd = 2,
       col = "#E41A1C", ylim = ylim,
       xlab = "Number of categories K", ylab = "Mean AIC")
  lines(agg_aic$K, agg_aic$aic_km_greedy, type = "b", pch = 17, lwd = 2,
        col = "#377EB8", lty = 2)
  legend("topleft", legend = c("OH", "KM"),
         col = c("#E41A1C", "#377EB8"), pch = c(16, 17), lty = c(1, 2), bty = "n")
})

agg_k <- aggregate(cbind(k_ohe, k_km_greedy) ~ K, data = df_res, FUN = mean, na.rm = TRUE)
save_pdf_plot("nparams_vs_K.pdf", {
  ylim <- range(c(agg_k$k_ohe, agg_k$k_km_greedy), na.rm = TRUE)
  plot(agg_k$K, agg_k$k_ohe, type = "b", pch = 16, lwd = 2,
       col = "#E41A1C", ylim = ylim,
       xlab = "Number of categories K", ylab = "Mean number of fitted coefficients")
  lines(agg_k$K, agg_k$k_km_greedy, type = "b", pch = 17, lwd = 2,
        col = "#377EB8", lty = 2)
  legend("topleft", legend = c("OH", "KM"),
         col = c("#E41A1C", "#377EB8"), pch = c(16, 17), lty = c(1, 2), bty = "n")
})

save_pdf_plot("boxplot_ibs_gain_by_censoring.pdf", {
  boxplot(split(df_res$ibs_gain, df_res$censor_prop),
          xlab = "Training censoring proportion",
          ylab = "IBS gain (OH - KM)",
          main = "")
  abline(h = 0, lty = 2, col = "grey50")
})

save_pdf_plot("boxplot_ibll_gain_by_censoring.pdf", {
  boxplot(split(df_res$ibll_gain, df_res$censor_prop),
          xlab = "Training censoring proportion",
          ylab = "IBLL gain (KM - OH)",
          main = "")
  abline(h = 0, lty = 2, col = "grey50")
})

save_pdf_plot("boxplot_aic_gain_by_censoring.pdf", {
  boxplot(split(df_res$aic_gain, df_res$censor_prop),
          xlab = "Training censoring proportion",
          ylab = "AIC gain (OH - KM)",
          main = "")
  abline(h = 0, lty = 2, col = "grey50")
})

save_pdf_plot("boxplot_n_selected_greedy_total.pdf", {
  boxplot(split(df_res$n_selected_total, df_res$censor_prop),
          xlab = "Training censoring proportion",
          ylab = "Total selected KM time points",
          main = "")
})

if (length(rep_curve_store) > 0) {
  brier_ohe <- as.data.frame(rep_curve_store$score_ohe$Brier$score)
  brier_km  <- as.data.frame(rep_curve_store$score_km_greedy$Brier$score)
  
  brier_ohe <- brier_ohe[is.finite(brier_ohe$times), , drop = FALSE]
  brier_km  <- brier_km[is.finite(brier_km$times), , drop = FALSE]
  
  save_pdf_plot("representative_brier_curves.pdf", {
    ylim <- range(c(brier_ohe$Brier, brier_km$Brier), na.rm = TRUE)
    plot(brier_ohe$times, brier_ohe$Brier, type = "l", lwd = 2,
         col = "#E41A1C", ylim = ylim,
         xlab = "Time", ylab = "Brier score",
         main = paste0("Representative regime: n=", REPREG_N,
                       ", K=", REPREG_K, ", censor=", REPREG_CENSOR))
    lines(brier_km$times, brier_km$Brier, lwd = 2, col = "#377EB8", lty = 2)
    legend("topright", legend = c("OH", "KM"),
           col = c("#E41A1C", "#377EB8"),
           lwd = 2, lty = c(1, 2), bty = "n")
  })
  
  if (!is.null(rep_curve_store$bll_ohe) && !is.null(rep_curve_store$bll_km_greedy)) {
    bll_ohe <- rep_curve_store$bll_ohe
    bll_km  <- rep_curve_store$bll_km_greedy
    
    save_pdf_plot("representative_bll_curves.pdf", {
      ylim <- range(c(bll_ohe$BLL, bll_km$BLL), na.rm = TRUE)
      plot(bll_ohe$times, bll_ohe$BLL, type = "l", lwd = 2,
           col = "#E41A1C", ylim = ylim,
           xlab = "Time", ylab = "Binomial log-likelihood",
           main = paste0("Representative regime: n=", REPREG_N,
                         ", K=", REPREG_K, ", censor=", REPREG_CENSOR))
      lines(bll_km$times, bll_km$BLL, lwd = 2, col = "#377EB8", lty = 2)
      legend("topright", legend = c("OH", "KM"),
             col = c("#E41A1C", "#377EB8"),
             lwd = 2, lty = c(1, 2), bty = "n")
    })
  }
}

save_pdf_plot("boxplot_cindex_gain_by_censoring.pdf", {
  boxplot(split(df_res$cindex_gain, df_res$censor_prop),
          xlab = "Training censoring proportion",
          ylab = "C-index gain (KM - OH)",
          main = "")
  abline(h = 0, lty = 2, col = "grey50")
})

save_pdf_plot("boxplot_iauc_gain_by_censoring.pdf", {
  boxplot(split(df_res$iauc_gain, df_res$censor_prop),
          xlab = "Training censoring proportion",
          ylab = "iAUC gain (KM - OH)",
          main = "")
  abline(h = 0, lty = 2, col = "grey50")
})

cat("\nAll PDFs saved in:", normalizePath(plot_dir), "\n")
print(head(df_res))

#  
# Robust significance plots
#  
make_sig_summary <- function(df, gain_col) {
  out <- do.call(rbind, lapply(
    split(df, list(df$K, df$censor_prop), drop = TRUE),
    function(sub) {
      x <- sub[[gain_col]]
      x <- x[is.finite(x)]
      
      if (length(x) < 5) return(NULL)
      
      t_res <- t.test(x, mu = 0)
      w_res <- wilcox.test(x, mu = 0, exact = FALSE)
      
      data.frame(
        K = unique(sub$K),
        censor_prop = unique(sub$censor_prop),
        n = length(x),
        mean_gain = mean(x),
        median_gain = median(x),
        ci_low = t_res$conf.int[1],
        ci_high = t_res$conf.int[2],
        p_t = t_res$p.value,
        p_w = w_res$p.value
      )
    }
  ))
  
  out$p_adj_w <- p.adjust(out$p_w, method = "BH")
  out$stars_w <- cut(
    out$p_adj_w,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", ""),
    right = TRUE
  )
  out$stars_w <- as.character(out$stars_w)
  out
}

sig_ibs <- make_sig_summary(df_res, "ibs_gain")
sig_ibll <- make_sig_summary(df_res, "ibll_gain")
sig_cindex <- make_sig_summary(df_res, "cindex_gain")

write.csv(sig_ibs, file.path(plot_dir, "sig_ibs.csv"), row.names = FALSE)
write.csv(sig_ibll, file.path(plot_dir, "sig_ibll.csv"), row.names = FALSE)
write.csv(sig_cindex, file.path(plot_dir, "sig_cindex.csv"), row.names = FALSE)

plot_gain_with_sig <- function(sig_df, ylab_text, filename) {
  save_pdf_plot(filename, {
    old_par <- par(no.readonly = TRUE); on.exit(par(old_par), add = TRUE)
    
    cps <- sort(unique(sig_df$censor_prop))
    par(mfrow = c(length(cps), 1),
        mar = c(4, 5, 3, 2),
        oma = c(2, 0, 2, 0))
    
    for (cp in cps) {
      sub <- sig_df[sig_df$censor_prop == cp, , drop = FALSE]
      sub <- sub[order(sub$K), , drop = FALSE]
      sub <- sub[is.finite(sub$mean_gain), , drop = FALSE]
      
      if (nrow(sub) == 0) {
        plot.new()
        title(main = "")
        next
      }
      
      y_for_range <- c(sub$mean_gain, sub$ci_low, sub$ci_high)
      y_for_range <- y_for_range[is.finite(y_for_range)]
      
      if (length(y_for_range) == 0) {
        ylim <- c(-1e-4, 1e-4)
      } else if (diff(range(y_for_range)) == 0) {
        center <- mean(y_for_range)
        ylim <- c(center - 1e-4, center + 1e-4)
      } else {
        pad <- 0.15 * diff(range(y_for_range))
        ylim <- c(min(y_for_range) - pad, max(y_for_range) + 2 * pad)
      }
      
      plot(sub$K, sub$mean_gain,
           type = "b", pch = 16, lwd = 2,
           ylim = ylim,
           xlab = "K",
           ylab = ylab_text,
           main = "")
      
      abline(h = 0, lty = 2, col = "grey50")
      
      ok_ci <- is.finite(sub$ci_low) & is.finite(sub$ci_high)
      if (any(ok_ci)) {
        arrows(sub$K[ok_ci], sub$ci_low[ok_ci],
               sub$K[ok_ci], sub$ci_high[ok_ci],
               angle = 90, code = 3, length = 0.05)
      }
      
      rng <- diff(range(ylim))
      if (!is.finite(rng) || rng == 0) rng <- 1e-4
      star_y <- ifelse(is.finite(sub$ci_high), sub$ci_high, sub$mean_gain) + 0.05 * rng
      
      text(sub$K, star_y,
           labels = sub$stars_w,
           cex = 1.3)
    }
  }, width = 8, height = 3.5 * length(sort(unique(sig_df$censor_prop))))
}

plot_gain_with_sig(
  sig_df = sig_ibs,
  ylab_text = "IBS gain (OH - KM)",
  filename = "ibs_gain_mean_ci_significance.pdf"
)

plot_gain_with_sig(
  sig_df = sig_ibll,
  ylab_text = "IBLL gain (KM - OH)",
  filename = "ibll_gain_mean_ci_significance.pdf"
)

plot_gain_with_sig(
  sig_df = sig_cindex,
  ylab_text = "C-index gain (KM - OH)",
  filename = "cindex_gain_mean_ci_significance.pdf"
)

#  
# KM jump diagnostics plots
#  

save_pdf_plot("km_jumps_vs_selected_times.pdf", {
  ok <- is.finite(df_res$mean_km_jumps) & is.finite(df_res$n_selected_total)
  
  plot(
    df_res$mean_km_jumps[ok],
    df_res$n_selected_total[ok],
    pch = 16,
    col = adjustcolor("purple", 0.35),
    xlab = "Mean KM jumps per category",
    ylab = "Total selected KM greedy time points",
    main = ""
  )
  
  if (sum(ok) >= 3) {
    lines(lowess(df_res$mean_km_jumps[ok], df_res$n_selected_total[ok]), lwd = 2)
  }
})

save_pdf_plot("km_jumps_vs_selected_by_censoring.pdf", {
  ok_all <- is.finite(df_res$mean_km_jumps) & is.finite(df_res$n_selected_total)
  df_plot <- df_res[ok_all, , drop = FALSE]
  
  cps <- sort(unique(df_plot$censor_prop))
  
  plot(
    NA,
    xlim = range(df_plot$mean_km_jumps, na.rm = TRUE),
    ylim = range(df_plot$n_selected_total, na.rm = TRUE),
    xlab = "Mean KM jumps per category",
    ylab = "Total selected KM greedy time points",
    main = ""
  )
  
  for (i in seq_along(cps)) {
    cp <- cps[i]
    sub <- df_plot[df_plot$censor_prop == cp, , drop = FALSE]
    
    points(
      sub$mean_km_jumps,
      sub$n_selected_total,
      pch = 15 + i,
      col = i
    )
    
    if (nrow(sub) >= 3) {
      lines(
        lowess(sub$mean_km_jumps, sub$n_selected_total),
        col = i,
        lwd = 2
      )
    }
  }
  
  legend(
    "topleft",
    legend = paste("Censoring =", cps),
    col = seq_along(cps),
    pch = 15 + seq_along(cps),
    lwd = 2,
    bty = "n"
  )
})

save_pdf_plot("selected_per_mean_jump_vs_K.pdf", {
  ok <- is.finite(df_res$selected_per_mean_jump)
  df_plot <- df_res[ok, , drop = FALSE]
  
  boxplot(
    selected_per_mean_jump ~ K,
    data = df_plot,
    xlab = "Number of categories K",
    ylab = "Selected times / mean KM jumps",
    main = ""
  )
  
  abline(h = 1, lty = 2, col = "grey50")
})

save_pdf_plot("selected_per_mean_jump_by_censoring.pdf", {
  ok <- is.finite(df_res$selected_per_mean_jump)
  df_plot <- df_res[ok, , drop = FALSE]
  
  boxplot(
    selected_per_mean_jump ~ censor_prop,
    data = df_plot,
    xlab = "Training censoring proportion",
    ylab = "Selected times / mean KM jumps",
    main = ""
  )
  
  abline(h = 1, lty = 2, col = "grey50")
})

save_pdf_plot("mean_km_jumps_by_K_and_censoring.pdf", {
  df_plot <- df_res[is.finite(df_res$mean_km_jumps), , drop = FALSE]
  
  Ks <- sort(unique(df_plot$K))
  cps <- sort(unique(df_plot$censor_prop))
  offset <- seq(-0.27, 0.27, length.out = length(cps))
  
  plot(
    NA,
    xlim = c(0.5, length(Ks) + 0.5),
    ylim = range(df_plot$mean_km_jumps, na.rm = TRUE),
    xaxt = "n",
    xlab = "Number of categories K",
    ylab = "Mean KM jumps per category",
    main = ""
  )
  
  axis(1, at = seq_along(Ks), labels = Ks)
  
  for (i in seq_along(cps)) {
    cp <- cps[i]
    
    boxplot(
      mean_km_jumps ~ factor(K, levels = Ks),
      data = df_plot[df_plot$censor_prop == cp, , drop = FALSE],
      at = seq_along(Ks) + offset[i],
      add = TRUE,
      boxwex = 0.12,
      xaxt = "n",
      yaxt = "n",
      outline = FALSE,
      border = i
    )
  }
  
  legend(
    "topright",
    legend = paste("Censoring =", cps),
    lwd = 2,
    col = seq_along(cps),
    bty = "n"
  )
})

#  
# Create summaries
#  

ibs_sum    <- metric_summary(df_res, "ibs_gain")
ibll_sum   <- metric_summary(df_res, "ibll_gain")
cindex_sum <- metric_summary(df_res, "cindex_gain")
aic_sum    <- metric_summary(df_res, "aic_gain")
iauc_sum   <- metric_summary(df_res, "iauc_gain")

#  
# Mean gain plots
#  

plot_gain_vs_K_by_censoring(
  ibs_sum,
  value_col = "mean",
  ylab_text = "Mean IBS gain (OH - KM)",
  title_text = "",
  filename = "ibs_gain_vs_K_by_censoring.pdf"
)

plot_gain_vs_K_by_censoring(
  ibll_sum,
  value_col = "mean",
  ylab_text = "Mean IBLL gain (KM - OH)",
  title_text = "",
  filename = "ibll_gain_vs_K_by_censoring.pdf"
)

plot_gain_vs_K_by_censoring(
  iauc_sum,
  value_col = "mean",
  ylab_text = "Mean iAUC gain (KM - OH)",
  title_text = "",
  filename = "iauc_gain_vs_K_by_censoring.pdf"
)

plot_gain_vs_K_by_censoring(
  cindex_sum,
  value_col = "mean",
  ylab_text = "Mean C-index gain (KM - OH)",
  title_text = "",
  filename = "cindex_gain_vs_K_by_censoring.pdf"
)

plot_gain_vs_K_by_censoring(
  aic_sum,
  value_col = "mean",
  ylab_text = "Mean AIC gain (OH - KM)",
  title_text = "",
  filename = "aic_gain_vs_K_by_censoring.pdf"
)

#  
# Median gain plots
#  

plot_gain_vs_K_by_censoring(
  ibs_sum,
  value_col = "median",
  ylab_text = "Median IBS gain (OH - KM)",
  title_text = "",
  filename = "ibs_gain_median_vs_K_by_censoring.pdf"
)

plot_gain_vs_K_by_censoring(
  ibll_sum,
  value_col = "median",
  ylab_text = "Median IBLL gain (KM - OH)",
  title_text = "",
  filename = "ibll_gain_median_vs_K_by_censoring.pdf"
)

plot_gain_vs_K_by_censoring(
  cindex_sum,
  value_col = "median",
  ylab_text = "Median C-index gain (KM - OH)",
  title_text = "",
  filename = "cindex_gain_median_vs_K_by_censoring.pdf"
)

#  
# 95% CI ribbon plots
#  

plot_gain_ci_ribbons(
  ibs_sum,
  ylab_text = "IBS gain (OH - KM)",
  title_text = "",
  filename = "ibs_gain_with_ci_ribbons.pdf"
)

plot_gain_ci_ribbons(
  ibll_sum,
  ylab_text = "IBLL gain (KM - OH)",
  title_text = "",
  filename = "ibll_gain_with_ci_ribbons.pdf"
)

plot_gain_ci_ribbons(
  cindex_sum,
  ylab_text = "C-index gain (KM - OH)",
  title_text = "",
  filename = "cindex_gain_with_ci_ribbons.pdf"
)

#  
# Mean AIC vs K by censoring, OHE and KM greedy
#  

aic_model_sum <- aggregate(
  cbind(aic_ohe, aic_km_greedy) ~ K + censor_prop,
  data = df_res,
  FUN = mean,
  na.rm = TRUE
)

save_pdf_plot("aic_vs_K_by_censoring.pdf", {
  cps <- sort(unique(aic_model_sum$censor_prop))
  cols <- seq_along(cps)
  
  ylim <- range(c(aic_model_sum$aic_ohe, aic_model_sum$aic_km_greedy), na.rm = TRUE)
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.08
  
  plot(NA,
       xlim = range(aic_model_sum$K),
       ylim = ylim,
       xlab = "Number of categories K",
       ylab = "Mean AIC",
       main = "")
  
  for (i in seq_along(cps)) {
    cp <- cps[i]
    sub <- aic_model_sum[aic_model_sum$censor_prop == cp, ]
    sub <- sub[order(sub$K), ]
    
    lines(sub$K, sub$aic_ohe,
          type = "b", pch = 16, lwd = 2, col = cols[i])
    
    lines(sub$K, sub$aic_km_greedy,
          type = "b", pch = 17, lwd = 2, lty = 2, col = cols[i])
  }
  
  legend("topleft",
         legend = paste("censor =", cps),
         col = cols,
         pch = 16,
         lwd = 2,
         bty = "n")
  
  legend("topright",
         legend = c("OH", "KM"),
         pch = c(16, 17),
         lty = c(1, 2),
         lwd = 2,
         col = "black",
         bty = "n")
})

cat("\nCategory distribution PDFs saved in:", normalizePath(plot_dir), "\n")
cat("Category distribution CSV:", file.path(normalizePath(plot_dir), "category_distribution_one_rep_per_K_censor.csv"), "\n")

agg_k <- aggregate(
  cbind(k_ohe, k_km_greedy) ~ K + censor_prop,
  data = df_res,
  FUN = mean,
  na.rm = TRUE
)

save_pdf_plot("nparams_vs_K_by_censoring.pdf", {
  cps  <- sort(unique(agg_k$censor_prop))
  cols <- seq_along(cps)
  
  ylim <- range(c(agg_k$k_ohe, agg_k$k_km_greedy), na.rm = TRUE)
  ylim <- ylim + c(-1, 1) * diff(ylim) * 0.08
  
  plot(NA,
       xlim = range(agg_k$K),
       ylim = ylim,
       xlab = "Number of categories K",
       ylab = "Mean number of fitted coefficients")
  
  for (i in seq_along(cps)) {
    cp  <- cps[i]
    sub <- agg_k[agg_k$censor_prop == cp, ]
    sub <- sub[order(sub$K), ]
    
    lines(sub$K, sub$k_ohe,
          type = "b", pch = 16, lwd = 2,
          col = cols[i], lty = 1)
    
    lines(sub$K, sub$k_km_greedy,
          type = "b", pch = 17, lwd = 2,
          col = cols[i], lty = 2)
  }
  
  legend("topleft",
         legend = paste("censor =", cps),
         col    = cols,
         lwd    = 2,
         pch    = 16,
         bty    = "n")
  
  legend("bottomright",
         legend = c("OH", "KM"),
         col    = "black",
         pch    = c(16, 17),
         lty    = c(1, 2),
         lwd    = 2,
         bty    = "n")
})

for (model in c("ohe", "km_greedy")) {
  k_col   <- paste0("k_", model)
  fname   <- paste0("nparams_vs_K_by_censoring_", model, ".pdf")
  
  save_pdf_plot(fname, {
    cps  <- sort(unique(agg_k$censor_prop))
    cols <- seq_along(cps)
    
    ylim <- range(agg_k[[k_col]], na.rm = TRUE)
    ylim <- ylim + c(-1, 1) * diff(ylim) * 0.08
    
    plot(NA,
         xlim = range(agg_k$K),
         ylim = ylim,
         xlab = "Number of categories K",
         ylab = "Mean number of fitted coefficients",
         main = "")
    
    for (i in seq_along(cps)) {
      sub <- agg_k[agg_k$censor_prop == cps[i], ]
      sub <- sub[order(sub$K), ]
      lines(sub$K, sub[[k_col]],
            type = "b", pch = 16, lwd = 2, col = cols[i])
    }
    
    legend("topleft",
           legend = paste("censor =", cps),
           col    = cols,
           pch    = 16,
           lwd    = 2,
           bty    = "n")
  })
}

save_pdf_plot("first_selected_time_histogram.pdf", {
  hist(
    df_selected_times$time[df_selected_times$step == 1],
    breaks = 10,
    col = "grey70",
    border = "white",
    xlab = "First selected KM time",
    main = "",
    xlim = range(candidate_times)
  )
})

save_pdf_plot("selected_times_aic_improvement.pdf", {
  
  gl <- df_selected_times
  
  if (is.null(gl) || nrow(gl) == 0) {
    plot.new()
    text(0.5, 0.5, "No greedy selections made", cex = 1.5)
    
  } else {
    
    reps <- sort(unique(gl$rep))
    
    par(
      mfrow = c(length(reps), 1),
      mar = c(5, 5, 3, 2),
      oma = c(0, 0, 2, 0)
    )
    
    for (rr in reps) {
      
      sub <- gl[gl$rep == rr, , drop = FALSE]
      
      plot(
        NA,
        xlim = range(candidate_times),
        ylim = c(0, max(c(gl$improve, 1), na.rm = TRUE) * 1.25),
        xlab = "Candidate time",
        ylab = expression(Delta * "AIC  (positive = improvement)")
      )
      
      if (nrow(sub) > 0) {
        
        for (s in seq_len(nrow(sub))) {
          
          x_s <- sub$time[s]
          y_s <- sub$improve[s]
          
          segments(
            x_s, 0, x_s, y_s,
            lwd = 2,
            col = adjustcolor("steelblue", 0.8)
          )
          
          points(
            x_s, y_s,
            pch = 21,
            bg = "steelblue",
            col = "white",
            cex = 1.6
          )
          
          text(
            x_s, y_s,
            labels = sub$step[s],
            pos = 3,
            cex = 0.8,
            col = "steelblue4"
          )
        }
        
        legend(
          "topright",
          legend = paste0(
            "Step ", sub$step,
            ":  t=", round(sub$time, 2)
          ),
          pch = 21,
          pt.bg = "steelblue",
          col = "steelblue4",
          bty = "n",
          cex = 0.85
        )
        
      } else {
        text(
          mean(candidate_times),
          0.5,
          "No selections for this replicate",
          cex = 1.2,
          col = "grey50"
        )
      }
    }
  }
  
}, width = 11, height = 4 * length(unique(df_selected_times$rep)))


save_pdf_plot("correlation_recovery_boxplot.pdf", {
  
  boxplot(
    list(
      KM = abs(store_correlation),
      OHE = abs(store_correlation_ohe)
    ),
    ylab = "|corr|"
  )
  
}, width = 6, height = 5)

cat("\nAll done. Outputs saved in:", normalizePath(plot_dir), "\n")
