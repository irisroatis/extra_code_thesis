suppressPackageStartupMessages({
  library(survival)
  library(riskRegression)
})



N_CAT_FEATURES <- 1   # set to 2 to add a second categorical feature

# DGP parameters
N_TRAIN      <- 1000
N_TEST       <- 10000
LAM          <- 0.02
CENSOR_PROP  <- 0.3
BETA_CONT    <- c(0.5, -0.4)   # continuous covariate effects
TAU_VEC      <- c(25)                # change-point(s); one per cat feature
ALPHA_VEC    <- c(1.0)               # latent effect magnitude(s)

if (N_CAT_FEATURES == 2) {
  TAU_VEC   <- c(25, 40)
  ALPHA_VEC <- c(1.0, 1.5)
}

C_PER_CAT <- 50   # number of categories per categorical feature

# Greedy search grid
N_CANDIDATE_TIMES <- 60
MIN_AIC_IMPROVE   <- 0    # use AIC stopping (any improvement)

# Output
PLOT_DIR <- "plots_interpretability"
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

set.seed(4)

# =============================================================================
# HELPERS
# =============================================================================

save_pdf <- function(filename, expr, width = 10, height = 7) {
  pdf(file.path(PLOT_DIR, filename), width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  eval.parent(substitute(expr))
}

viridis_ramp <- colorRampPalette(
  c("#440154", "#31688e", "#35b779", "#fde725")
)

scale_to_color <- function(x, n_colors = 256, pal = viridis_ramp) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(pal(n_colors)[ceiling(n_colors / 2)], length(x)))
  idx <- pmax(1, pmin(n_colors, round(1 + (n_colors - 1) * (x - rng[1]) / diff(rng))))
  pal(n_colors)[idx]
}

add_colorbar <- function(zlim, pal = viridis_ramp, label = "z",
                         n = 200, x_frac_start = 1.03, x_frac_end = 1.08) {
  usr <- par("usr")
  xw  <- diff(usr[1:2])
  yw  <- diff(usr[3:4])
  x1  <- usr[2] + x_frac_start * xw
  x2  <- usr[2] + x_frac_end   * xw
  ys  <- seq(usr[3], usr[4], length.out = n)
  cols <- pal(n - 1)
  rect(x1, head(ys, -1), x2, tail(ys, -1),
       col = cols, border = cols, xpd = TRUE)
  text(x2, usr[3], round(zlim[1], 2), pos = 4, xpd = TRUE, cex = 0.75)
  text(x2, usr[4], round(zlim[2], 2), pos = 4, xpd = TRUE, cex = 0.75)
  text((x1 + x2) / 2, usr[4] + 0.06 * yw,
       label, xpd = TRUE, cex = 0.8, font = 2)
}

cox_aic_value <- function(fit, n_events) {
  ll <- as.numeric(fit$loglik[2])
  k  <- sum(is.finite(coef(fit)) & !is.na(coef(fit)))
  -2 * ll + 2 * k
}

predict_surv_at <- function(sf, times) {
  as.numeric(summary(sf, times = times, extend = TRUE)$surv)
}

# =============================================================================
# DATA GENERATING PROCESS  (non-proportional hazards)
# =============================================================================

simulate_data <- function(n, C_vec, lam, censor_prop, tau_vec, beta_cont, alpha_vec,
                          binning = "quantile") {
  stopifnot(length(C_vec) == length(tau_vec),
            length(tau_vec) == length(alpha_vec))
  
  n_cont <- length(beta_cont)
  n_cat  <- length(tau_vec)
  
  # --- continuous covariates ---
  cont_df <- as.data.frame(
    setNames(lapply(seq_len(n_cont), function(i) rnorm(n)),
             paste0("cont_", seq_len(n_cont)))
  )
  lp_cont <- as.numeric(as.matrix(cont_df) %*% beta_cont)
  
  # --- categorical covariates (binned latent) ---
  latent_mat <- matrix(0, n, n_cat)
  cat_df     <- data.frame(matrix(NA_character_, n, n_cat))
  names(cat_df) <- paste0("cat_", seq_len(n_cat))
  
  for (i in seq_len(n_cat)) {
    z <- rnorm(n)
    if (binning == "quantile") {
      edges <- quantile(z, probs = seq(0, 1, length.out = C_vec[i] + 1))
    } else {
      edges <- seq(min(z), max(z), length.out = C_vec[i] + 1)
    }
    edges[1] <- -Inf; edges[length(edges)] <- Inf
    g <- findInterval(z, edges, rightmost.closed = FALSE, all.inside = TRUE)
    cat_df[[paste0("cat_", i)]] <- paste0("C", i, "_", g)
    latent_mat[, i] <- z
  }
  colnames(latent_mat) <- paste0("latent_", seq_len(n_cat))
  
  # --- piecewise-constant hazard (sign flip at each tau) ---
  tau_sorted <- sort(unique(tau_vec))
  starts <- c(0, tau_sorted)
  ends   <- c(tau_sorted, Inf)
  
  rate_mat <- matrix(0, n, length(starts))
  for (k in seq_along(starts)) {
    signs    <- ifelse(ends[k] <= tau_vec, -1, 1)
    lp_cat   <- latent_mat %*% (alpha_vec * signs)
    rate_mat[, k] <- lam * exp(lp_cont + lp_cat)
  }
  
  # inversion sampling
  u       <- runif(n)
  target  <- -log(u)
  T_event <- numeric(n)
  
  for (i in seq_len(n)) {
    cumhaz <- 0
    for (k in seq_along(starts)) {
      r <- rate_mat[i, k]
      if (is.infinite(ends[k])) {
        T_event[i] <- starts[k] + (target[i] - cumhaz) / r
        break
      }
      h <- r * (ends[k] - starts[k])
      if (target[i] <= cumhaz + h) {
        T_event[i] <- starts[k] + (target[i] - cumhaz) / r
        break
      }
      cumhaz <- cumhaz + h
    }
  }
  
  lambda_c <- (censor_prop / (1 - censor_prop)) * lam
  T_cens   <- rexp(n, lambda_c)
  Time     <- pmin(T_event, T_cens)
  Event    <- as.integer(T_event <= T_cens)
  
  cbind(cont_df, cat_df, as.data.frame(latent_mat),
        data.frame(Time = Time, Event = Event))
}

# =============================================================================
# ONE-HOT ENCODING  (strict: fully-censored categories → all zeros)
# =============================================================================

ohe <- function(cat_cols, train_df, test_df) {
  tr <- train_df; te <- test_df
  
  for (cc in cat_cols) {
    tr_cats <- as.character(train_df[[cc]])
    te_cats <- as.character(test_df[[cc]])
    lvls    <- unique(tr_cats)
    has_ev  <- sapply(lvls, function(lv) any(train_df$Event[tr_cats == lv] == 1))
    valid   <- lvls[has_ev]
    if (length(valid) < 2) stop(paste("Too few valid levels in", cc))
    
    base   <- valid[1]
    kept   <- setdiff(valid, base)
    cnames <- make.names(paste0(cc, "_", kept))
    
    encode <- function(cats) {
      m <- matrix(0, length(cats), length(kept),
                  dimnames = list(NULL, cnames))
      for (i in seq_along(cats))
        if (cats[i] %in% kept)
          m[i, make.names(paste0(cc, "_", cats[i]))] <- 1
      m
    }
    
    tr <- cbind(tr[, setdiff(names(tr), cc), drop = FALSE],
                as.data.frame(encode(tr_cats)))
    te <- cbind(te[, setdiff(names(te), cc), drop = FALSE],
                as.data.frame(encode(te_cats)))
  }
  list(train = tr, test = te)
}

# =============================================================================
# KM PRECOMPUTATION
# =============================================================================

precompute_km <- function(train_df, cat_cols, candidate_times) {
  fallback <- survfit(Surv(Time, Event) ~ 1, data = train_df)
  fb_vec   <- predict_surv_at(fallback, candidate_times)
  
  lapply(setNames(cat_cols, cat_cols), function(cc) {
    lvls   <- unique(as.character(train_df[[cc]]))
    lookup <- lapply(setNames(lvls, lvls), function(lv) {
      sub <- train_df[as.character(train_df[[cc]]) == lv, , drop = FALSE]
      if (nrow(sub) == 0 || sum(sub$Event) == 0) return(fb_vec)
      predict_surv_at(survfit(Surv(Time, Event) ~ 1, data = sub), candidate_times)
    })
    list(candidate_times = candidate_times, lookup = lookup, fallback = fb_vec)
  })
}

build_km_matrices <- function(df, precomp, cat_cols) {
  lapply(setNames(cat_cols, cat_cols), function(cc) {
    pc   <- precomp[[cc]]
    cats <- as.character(df[[cc]])
    uniq <- unique(cats)
    rl   <- lapply(setNames(uniq, uniq), function(lv) {
      v <- pc$lookup[[lv]]; if (is.null(v)) v <- pc$fallback; v
    })
    m <- do.call(rbind, rl[cats])
    colnames(m) <- paste0(cc, "_km_t", seq_len(ncol(m)))
    as.matrix(m)
  })
}

assemble_km_matrix <- function(full_mats, selected_idx, cat_cols) {
  parts <- Filter(Negate(is.null), lapply(cat_cols, function(cc) {
    idx <- selected_idx[[cc]]
    if (length(idx) > 0) full_mats[[cc]][, idx, drop = FALSE] else NULL
  }))
  if (length(parts) == 0)
    return(matrix(nrow = nrow(full_mats[[cat_cols[1]]]), ncol = 0))
  do.call(cbind, parts)
}

# =============================================================================
# GREEDY AIC SELECTION
# =============================================================================

km_greedy_aic <- function(train_df, cont_cols, cat_cols,
                          candidate_times, full_mats) {
  selected  <- setNames(lapply(cat_cols, function(cc) integer(0)), cat_cols)
  base_df   <- train_df[, c(cont_cols, "Time", "Event"), drop = FALSE]
  n_events  <- max(2L, sum(train_df$Event))
  
  base_fit  <- coxph(Surv(Time, Event) ~ ., data = base_df,
                     ties = "efron", singular.ok = FALSE)
  best_aic  <- cox_aic_value(base_fit, n_events)
  
  log_df <- data.frame(step = integer(), cat = character(),
                       time_idx = integer(), time_val = numeric(),
                       aic = numeric(), delta_aic = numeric())
  
  max_per_cat <- setNames(
    sapply(cat_cols, function(cc) length(unique(as.character(train_df[[cc]])))),
    cat_cols
  )
  
  repeat {
    pool <- lapply(setNames(cat_cols, cat_cols), function(cc) {
      if (length(selected[[cc]]) >= max_per_cat[[cc]]) return(integer(0))
      setdiff(seq_along(candidate_times), selected[[cc]])
    })
    if (sum(lengths(pool)) == 0) break
    
    X_cur    <- assemble_km_matrix(full_mats, selected, cat_cols)
    best_cc  <- NULL; best_j <- NA_integer_; best_new_aic <- Inf
    
    for (cc in cat_cols) {
      for (j in pool[[cc]]) {
        xn  <- full_mats[[cc]][, j, drop = FALSE]
        Xtr <- if (ncol(X_cur) == 0) xn else cbind(X_cur, xn)
        fit <- try(coxph(Surv(Time, Event) ~ .,
                         data = cbind(base_df, as.data.frame(Xtr)),
                         ties = "efron", singular.ok = FALSE), silent = TRUE)
        if (inherits(fit, "try-error")) next
        a <- cox_aic_value(fit, n_events)
        if (is.finite(a) && a < best_new_aic) {
          best_new_aic <- a; best_cc <- cc; best_j <- j
        }
      }
    }
    
    if (is.null(best_cc) || !is.finite(best_new_aic)) break
    delta <- best_aic - best_new_aic
    if (delta <= MIN_AIC_IMPROVE) break
    
    selected[[best_cc]] <- c(selected[[best_cc]], best_j)
    log_df <- rbind(log_df, data.frame(
      step      = sum(lengths(selected)),
      cat       = best_cc,
      time_idx  = best_j,
      time_val  = candidate_times[best_j],
      aic       = best_new_aic,
      delta_aic = delta
    ))
    best_aic <- best_new_aic
  }
  
  list(selected = selected, log = log_df, final_aic = best_aic)
}

# =============================================================================
# SIMULATE DATA
# =============================================================================

cat_cols_all  <- paste0("cat_",     seq_len(N_CAT_FEATURES))
lat_cols_all  <- paste0("latent_",  seq_len(N_CAT_FEATURES))
cont_cols_all <- paste0("cont_",    seq_along(BETA_CONT))
C_VEC         <- rep(C_PER_CAT, N_CAT_FEATURES)

df_full <- simulate_data(
  n           = N_TRAIN + N_TEST,
  C_vec       = C_VEC,
  lam         = LAM,
  censor_prop = CENSOR_PROP,
  tau_vec     = TAU_VEC,
  beta_cont   = BETA_CONT,
  alpha_vec   = ALPHA_VEC
)

train_idx <- seq_len(N_TRAIN)
test_idx  <- (N_TRAIN + 1):(N_TRAIN + N_TEST)

train_all <- df_full[train_idx, , drop = FALSE]
test_all  <- df_full[test_idx,  , drop = FALSE]

# strip latent columns from the modelling frames
noz_cols   <- setdiff(names(train_all), lat_cols_all)
train_noz  <- train_all[, noz_cols, drop = FALSE]
test_noz   <- test_all[,  noz_cols, drop = FALSE]

cat("Events in train:", sum(train_noz$Event),
    "| Event fraction:", round(mean(train_noz$Event), 3), "\n")

# =============================================================================
# FIT MODELS
# =============================================================================

# ---- latent (oracle) ----
lat_train <- train_all[, c(cont_cols_all, lat_cols_all, "Time", "Event")]
lat_test  <- test_all[,  c(cont_cols_all, lat_cols_all, "Time", "Event")]
cox_lat   <- coxph(Surv(Time, Event) ~ ., data = lat_train,
                   ties = "efron", x = TRUE, y = TRUE, model = TRUE)

# ---- OHE ----
ohe_res   <- ohe(cat_cols_all, train_noz, test_noz)
cox_ohe   <- coxph(Surv(Time, Event) ~ ., data = ohe_res$train,
                   ties = "efron", singular.ok = FALSE,
                   x = TRUE, y = TRUE, model = TRUE)

# ---- KM greedy ----
event_times     <- sort(unique(train_noz$Time[train_noz$Event == 1]))
candidate_times <- unique(as.numeric(quantile(
  event_times, probs = seq(0.05, 0.95, length.out = N_CANDIDATE_TIMES),
  na.rm = TRUE
)))

km_precomp      <- precompute_km(train_noz, cat_cols_all, candidate_times)
full_mats_train <- build_km_matrices(train_noz, km_precomp, cat_cols_all)
full_mats_test  <- build_km_matrices(test_noz,  km_precomp, cat_cols_all)

greedy_out <- km_greedy_aic(
  train_df        = train_noz,
  cont_cols       = cont_cols_all,
  cat_cols        = cat_cols_all,
  candidate_times = candidate_times,
  full_mats       = full_mats_train
)

X_tr_grd <- assemble_km_matrix(full_mats_train, greedy_out$selected, cat_cols_all)
X_te_grd <- assemble_km_matrix(full_mats_test,  greedy_out$selected, cat_cols_all)

grd_train <- cbind(train_noz[, c(cont_cols_all, "Time", "Event")],
                   as.data.frame(X_tr_grd))
grd_test  <- cbind(test_noz[,  c(cont_cols_all, "Time", "Event")],
                   as.data.frame(X_te_grd))

cox_grd <- coxph(Surv(Time, Event) ~ ., data = grd_train,
                 ties = "efron", singular.ok = FALSE,
                 x = TRUE, y = TRUE, model = TRUE)

# quick summary
cat("\nGreedy selection log:\n")
print(greedy_out$log)


save_pdf("interp_1_greedy_selection_path.pdf", {
  gl <- greedy_out$log
  
  if (nrow(gl) == 0) {
    plot.new()
    text(0.5, 0.5, "No greedy selections made", cex = 1.5)
  } else {
    # one panel per cat feature
    n_panels <- length(cat_cols_all)
    par(mfrow = c(n_panels, 1), mar = c(5, 5, 3, 2), oma = c(0, 0, 2, 0))
    
    for (i in seq_along(cat_cols_all)) {
      cc  <- cat_cols_all[i]
      tau <- TAU_VEC[i]
      sub <- gl[gl$cat == cc, , drop = FALSE]
      
      # background: AIC improvement landscape at each candidate time
      # (computed once, for the final selected set + each candidate)
      plot(NA,
           xlim = range(candidate_times),
           ylim = c(0, max(c(gl$delta_aic, 1), na.rm = TRUE) * 1.25),
           xlab = "Candidate time",
           ylab = expression(Delta * "AIC  (positive = improvement)"))
      
      abline(v = tau, lty = 2, lwd = 2, col = "firebrick")
      mtext(bquote(tau == .(tau)), side = 3, at = tau,
            col = "firebrick", cex = 0.85, line = 0.2)
      
      if (nrow(sub) > 0) {
        # draw selected times as vertical segments from 0
        for (s in seq_len(nrow(sub))) {
          x_s <- sub$time_val[s]
          y_s <- sub$delta_aic[s]
          segments(x_s, 0, x_s, y_s,
                   lwd = 2,
                   col = adjustcolor("steelblue", 0.8))
          points(x_s, y_s, pch = 21,
                 bg = "steelblue", col = "white", cex = 1.6)
          text(x_s, y_s, labels = sub$step[s],
               pos = 3, cex = 0.8, col = "steelblue4")
        }
        
        # annotate whether each selected time is pre- or post-tau
        sides <- ifelse(sub$time_val < tau, "pre", "post")
        legend("topright",
               legend = paste0("Step ", sub$step, ":  t=",
                               round(sub$time_val, 2)),
               pch    = 21,
               pt.bg  = "steelblue",
               col    = "steelblue4",
               bty    = "n", cex = 0.85)
      } else {
        text(mean(candidate_times), 0.5, "No selections for this feature",
             cex = 1.2, col = "grey50")
      }
    }
    
  }
}, width = 11, height = 4 * length(cat_cols_all))



save_pdf("interp_2_km_curves_selected_times.pdf", {
  par(mfrow = c(1, length(cat_cols_all)),
      mar   = c(5, 5, 4, 7))
  
  for (i in seq_along(cat_cols_all)) {
    cc       <- cat_cols_all[i]
    lat_col  <- lat_cols_all[i]
    tau      <- TAU_VEC[i]
    pc       <- km_precomp[[cc]]
    sel_idx  <- greedy_out$selected[[cc]]
    sel_times <- candidate_times[sel_idx]
    
    # mean latent per category (from full training set)
    cats   <- names(pc$lookup)
    z_mean <- sapply(cats, function(lv)
      mean(train_all[[lat_col]][as.character(train_all[[cc]]) == lv], na.rm = TRUE))
    
    pt_cols <- scale_to_color(z_mean)
    
    # time axis for plotting
    t_axis <- seq(min(candidate_times), max(candidate_times), length.out = 300)
    
    ylim <- c(0, 1)
    plot(NA, xlim = range(t_axis), ylim = ylim,
         xlab = "Time", ylab = "S(t | category)",
         main = paste0(cc, "\nKM curves coloured by mean latent z"))
    
    abline(v = tau, lty = 2, lwd = 2, col = "firebrick")
    mtext(bquote(tau == .(tau)), side = 3, at = tau,
          col = "firebrick", cex = 0.8, line = 0.1)
    
    for (ci in seq_along(cats)) {
      lv  <- cats[ci]
      vec <- pc$lookup[[lv]]
      # interpolate onto t_axis
      s_t <- approx(pc$candidate_times, vec, xout = t_axis,
                    method = "constant", rule = 2)$y
      lines(t_axis, s_t, col = adjustcolor(pt_cols[ci], 0.45), lwd = 0.9)
    }
    
    # overlay selected times
    if (length(sel_times) > 0) {
      abline(v = sel_times, lty = 3, lwd = 1.5, col = "steelblue")
      text(sel_times, 0.02,
           labels = paste0("t", seq_along(sel_times)),
           col = "steelblue4", cex = 0.75, srt = 90, adj = 0)
    }
    
    # colorbar for latent z
    add_colorbar(
      zlim  = range(z_mean),
      label = bquote(bar(z)[.(sub("cat_", "", cc))])
    )
  }
}, width = 7 * length(cat_cols_all), height = 7)


save_pdf("interp_3_km_embedding_2d.pdf", {
  n_panels <- sum(sapply(cat_cols_all, function(cc)
    length(greedy_out$selected[[cc]]) >= 2))
  
  if (n_panels == 0) {
    plot.new()
    text(0.5, 0.5, "Fewer than 2 dimensions selected for all features", cex = 1.2)
  } else {
    par(mfrow = c(1, max(n_panels, 1)), mar = c(5, 5, 4, 7))
    
    for (i in seq_along(cat_cols_all)) {
      cc      <- cat_cols_all[i]
      lat_col <- lat_cols_all[i]
      sel_idx <- greedy_out$selected[[cc]]
      
      if (length(sel_idx) < 2) {
        plot.new()
        text(0.5, 0.5, paste("Only", length(sel_idx), "dim selected for", cc),
             cex = 1.1, col = "grey50")
        next
      }
      
      idx12  <- sel_idx[1:2]
      t1     <- candidate_times[idx12[1]]
      t2     <- candidate_times[idx12[2]]
      pc     <- km_precomp[[cc]]
      cats   <- names(pc$lookup)
      
      coords <- t(sapply(cats, function(lv) {
        v <- pc$lookup[[lv]]; if (is.null(v)) v <- pc$fallback
        c(v[idx12[1]], v[idx12[2]])
      }))
      
      z_mean  <- sapply(cats, function(lv)
        mean(train_all[[lat_col]][as.character(train_all[[cc]]) == lv], na.rm = TRUE))
      pt_cols <- scale_to_color(z_mean)
      
      ylim <- c(0, 1)
      
      plot(
        coords[, 1], coords[, 2],
        pch = 19,
        cex = 1.2,
        col = pt_cols,
        # ylim = ylim,
        xlab = bquote(
          hat(S)(
            t[1] == .(round(t1, 1))
            ~ "|" ~
              x[1] == C[1 * "," * j]
          )
        ),
        ylab = bquote(
          hat(S)(
            t[2] == .(round(t2, 1))
            ~ "|" ~
              x[1] == C[1 * "," * j]
          )
        )
      )
      
      add_colorbar(
        zlim  = range(z_mean),
        label = bquote(bar(z)[.(sub("cat_", "", cc))])
      )
    }
  }
}, width = 7 * max(n_panels, 1), height = 7)


km_lp_per_category <- function(cox_grd, cat_cols, cat_col_focus, train_noz,
                               full_mats_train, selected, candidate_times,
                               cont_cols) {
  cats <- unique(as.character(train_noz[[cat_col_focus]]))
  
  # median continuous covariates (to isolate the categorical contribution)
  med_cont <- lapply(setNames(cont_cols, cont_cols), function(cc)
    median(train_noz[[cc]], na.rm = TRUE))
  
  # one representative row per category
  lp_vals <- sapply(cats, function(lv) {
    # build the KM-encoded row for this category
    km_row <- lapply(cat_cols, function(cc) {
      idx <- selected[[cc]]
      if (length(idx) == 0) return(NULL)
      pc  <- km_precomp[[cc]]
      # use this cat's KM vec if cc == cat_col_focus, else use marginal
      vec <- if (cc == cat_col_focus) {
        v <- pc$lookup[[lv]]; if (is.null(v)) v <- pc$fallback; v
      } else {
        pc$fallback   # marginal for other features
      }
      setNames(as.list(vec[idx]),
               paste0(cc, "_km_t", idx))
    })
    km_df <- as.data.frame(do.call(c, Filter(Negate(is.null), km_row)))
    
    row <- cbind(as.data.frame(med_cont), km_df,
                 data.frame(Time = 1, Event = 0))
    as.numeric(predict(cox_grd, newdata = row, type = "lp"))
  })
  
  setNames(lp_vals, cats)
}

save_pdf("interp_4_coefficient_comparison.pdf", {
  n_rows <- length(cat_cols_all)
  # 3 columns: pre-tau, post-tau, KM-lp vs latent
  par(mfrow = c(n_rows, 3), mar = c(5, 5, 4, 2), oma = c(0, 0, 3, 0))
  
  for (i in seq_along(cat_cols_all)) {
    cc      <- cat_cols_all[i]
    lat_col <- lat_cols_all[i]
    tau     <- TAU_VEC[i]
    alpha   <- ALPHA_VEC[i]
    
    cats <- unique(as.character(train_noz[[cc]]))
    
    # mean latent per category
    z_mean <- sapply(cats, function(lv)
      mean(train_all[[lat_col]][as.character(train_all[[cc]]) == lv], na.rm = TRUE))
    
    # true DGP effects
    true_pre  <- -alpha * z_mean   # sign = -1 before tau
    true_post <-  alpha * z_mean   # sign = +1 after  tau
    
    # OHE coefficients (0 for baseline / censored)
    prefix  <- paste0(cc, "_")
    ohe_map <- setNames(rep(0, length(cats)), cats)
    for (nm in names(coef(cox_ohe))) {
      if (startsWith(nm, prefix)) {
        lv_raw <- sub(prefix, "", nm, fixed = TRUE)
        # reverse make.names: find best match
        match_lv <- cats[make.names(cats) == lv_raw]
        if (length(match_lv) == 1) ohe_map[match_lv] <- coef(cox_ohe)[nm]
      }
    }
    
    # KM greedy per-category LP
    km_lp <- km_lp_per_category(
      cox_grd         = cox_grd,
      cat_cols        = cat_cols_all,
      cat_col_focus   = cc,
      train_noz       = train_noz,
      full_mats_train = full_mats_train,
      selected        = greedy_out$selected,
      candidate_times = candidate_times,
      cont_cols       = cont_cols_all
    )
    km_lp_vec <- km_lp[cats]
    
    pt_cols <- scale_to_color(z_mean)
    
    # ---- panel A: OHE beta vs true pre-tau effect ----
    plot(true_pre, ohe_map[cats],
         pch = 19, col = pt_cols, cex = 1.1,
         xlab = bquote("True pre-\u03c4 effect  (-\u03b1 " * bar(z) * ")"),
         ylab = "OHE fitted \u03b2",
         main = paste0(cc, "\nOHE vs true (pre-\u03c4)"))
    abline(lm(ohe_map[cats] ~ true_pre), lty = 2, col = "grey50")
    abline(0, 1, lty = 3, col = "tomato")
    
    # ---- panel B: OHE beta vs true post-tau effect ----
    plot(true_post, ohe_map[cats],
         pch = 19, col = pt_cols, cex = 1.1,
         xlab = bquote("True post-\u03c4 effect  (+\u03b1 " * bar(z) * ")"),
         ylab = "OHE fitted \u03b2",
         main = paste0(cc, "\nOHE vs true (post-\u03c4)"))
    abline(lm(ohe_map[cats] ~ true_post), lty = 2, col = "grey50")
    abline(0, 1, lty = 3, col = "tomato")
    
    # ---- panel C: KM greedy LP vs latent z (time-averaged) ----
    plot(z_mean, km_lp_vec,
         pch = 19, col = pt_cols, cex = 1.1,
         xlab = bquote("Mean latent " * bar(z)),
         ylab = "KM greedy LP (marginal)",
         main = paste0(cc, "\nKM greedy LP vs latent z"))
    abline(lm(km_lp_vec ~ z_mean), lty = 2, col = "grey50")
    # note: KM encoding is non-linear in z because of the sign flip;
    # we expect a *non-monotone* or attenuated relationship — that is the point
  }
  
  mtext("Coefficient comparison: OHE (pre/post-\u03c4) vs KM greedy LP",
        outer = TRUE, cex = 1.1)
}, width = 14, height = 5 * length(cat_cols_all))



martingale_by_category <- function(cox_fit, df_model, cat_col, train_full, lat_col) {
  mr   <- residuals(cox_fit, type = "martingale")
  cats <- as.character(df_model[[cat_col]])
  
  uniq  <- unique(cats)
  z_mean <- sapply(uniq, function(lv)
    mean(train_full[[lat_col]][as.character(train_full[[cat_col]]) == lv], na.rm = TRUE))
  
  mr_mean <- sapply(uniq, function(lv) mean(mr[cats == lv], na.rm = TRUE))
  mr_sd   <- sapply(uniq, function(lv) sd(mr[cats == lv],   na.rm = TRUE))
  n_cat   <- sapply(uniq, function(lv) sum(cats == lv))
  
  data.frame(
    cat     = uniq,
    z_mean  = z_mean,
    mr_mean = mr_mean,
    mr_sd   = mr_sd,
    n       = n_cat
  )
}

save_pdf("interp_5_martingale_residuals.pdf", {
  n_rows <- length(cat_cols_all)
  par(mfrow = c(n_rows, 2), mar = c(5, 5, 4, 2), oma = c(0, 0, 3, 0))
  
  for (i in seq_along(cat_cols_all)) {
    cc      <- cat_cols_all[i]
    lat_col <- lat_cols_all[i]
    
    # OHE model uses ohe_res$train which has the cat column removed
    # we need to attach the original cat labels back for grouping
    ohe_tr_with_cat <- cbind(
      ohe_res$train,
      setNames(data.frame(as.character(train_noz[[cc]])), cc)
    )
    
    mr_ohe <- martingale_by_category(
      cox_fit    = cox_ohe,
      df_model   = ohe_tr_with_cat,
      cat_col    = cc,
      train_full = train_all,
      lat_col    = lat_col
    )
    
    # KM greedy model: grd_train also has cat removed, re-attach
    grd_tr_with_cat <- cbind(
      grd_train,
      setNames(data.frame(as.character(train_noz[[cc]])), cc)
    )
    
    mr_grd <- martingale_by_category(
      cox_fit    = cox_grd,
      df_model   = grd_tr_with_cat,
      cat_col    = cc,
      train_full = train_all,
      lat_col    = lat_col
    )
    
    ylim <- range(c(mr_ohe$mr_mean, mr_grd$mr_mean,
                    mr_ohe$mr_mean + mr_ohe$mr_sd,
                    mr_grd$mr_mean + mr_grd$mr_sd,
                    mr_ohe$mr_mean - mr_ohe$mr_sd,
                    mr_grd$mr_mean - mr_grd$mr_sd), na.rm = TRUE)
    
    pt_cols_ohe <- scale_to_color(mr_ohe$z_mean)
    pt_cols_grd <- scale_to_color(mr_grd$z_mean)
    
    for (mod_name in c("OHE", "KM greedy")) {
      df_mr   <- if (mod_name == "OHE") mr_ohe else mr_grd
      pt_cols <- if (mod_name == "OHE") pt_cols_ohe else pt_cols_grd
      
      plot(df_mr$z_mean, df_mr$mr_mean,
           pch = 19, cex = 1.1, col = pt_cols,
           ylim = ylim,
           xlab = bquote("Mean latent " * bar(z)[.(sub("cat_", "", cc))]),
           ylab = "Mean martingale residual",
           main = paste0(cc, " — ", mod_name))
      
      # error bars: ±1 SD / sqrt(n)
      se <- df_mr$mr_sd / sqrt(pmax(df_mr$n, 1))
      arrows(df_mr$z_mean, df_mr$mr_mean - se,
             df_mr$z_mean, df_mr$mr_mean + se,
             angle = 90, code = 3, length = 0.04,
             col = adjustcolor("grey30", 0.5))
      
      abline(h = 0, lty = 2, col = "firebrick", lwd = 1.5)
      lines(lowess(df_mr$z_mean, df_mr$mr_mean, f = 0.6),
            lwd = 2, col = "steelblue")
    }
  }
  
  mtext(
    "Martingale residuals by category  |  well-specified \u2192 residuals centered at 0\n(error bars = \u00b11 SE, blue = lowess smooth)",
    outer = TRUE, cex = 0.95
  )
}, width = 12, height = 5 * length(cat_cols_all))



save_pdf("interp_6_survival_reversal_at_selected_times.pdf", {
  n_panels_total <- sum(sapply(cat_cols_all, function(cc)
    length(greedy_out$selected[[cc]]) >= 2))
  
  if (n_panels_total == 0) {
    plot.new()
    text(0.5, 0.5, "Need at least 2 selected times per feature", cex = 1.2)
  } else {
    par(mfrow = c(length(cat_cols_all), 2),
        mar   = c(5, 5, 4, 7))
    
    for (i in seq_along(cat_cols_all)) {
      cc      <- cat_cols_all[i]
      lat_col <- lat_cols_all[i]
      tau     <- TAU_VEC[i]
      sel_idx <- greedy_out$selected[[cc]]
      pc      <- km_precomp[[cc]]
      cats    <- names(pc$lookup)
      
      z_mean  <- sapply(cats, function(lv)
        mean(train_all[[lat_col]][as.character(train_all[[cc]]) == lv], na.rm = TRUE))
      pt_cols <- scale_to_color(z_mean)
      
      if (length(sel_idx) < 2) {
        plot.new(); plot.new()
        next
      }
      
      for (which_dim in 1:2) {
        j   <- sel_idx[which_dim]
        t_j <- candidate_times[j]
        s_j <- sapply(cats, function(lv) {
          v <- pc$lookup[[lv]]; if (is.null(v)) v <- pc$fallback; v[j]
        })
        
        side_label <- if (t_j < tau) "pre-\u03c4" else "post-\u03c4"
        
        plot(z_mean, s_j,
             pch = 19, cex = 1.2, col = pt_cols,
             xlab = bquote("Mean latent " * bar(z)),
             ylab = bquote(hat(S)(.(round(t_j, 1)) * " | cat")),
             main = paste0(cc, "  |  dim ", which_dim,
                           "  (t = ", round(t_j, 1), ",  ", side_label, ")"))
        
        abline(lm(s_j ~ z_mean), lty = 2, col = "grey50", lwd = 1.5)
        abline(v = 0, lty = 3, col = "grey70")
        
        add_colorbar(
          zlim  = range(z_mean),
          label = bquote(bar(z)[.(sub("cat_", "", cc))])
        )
      }
    }
  }
}, width = 14, height = 7 * length(cat_cols_all))


save_pdf("interp_0_true_survival_switch.pdf", {
  
  i     <- 1
  tau   <- TAU_VEC[i]
  alpha <- ALPHA_VEC[i]
  
  # latent values to illustrate
  z_vals <- c(-1, -0.5, 0, 0.5, 1)
  
  # colours
  cols <- viridis_ramp(length(z_vals))
  
  # plotting grid
  t_grid <- seq(0, max(candidate_times), length.out = 500)
  
  # true survival curves
  surv_mat <- sapply(z_vals, function(z) {
    
    h_pre  <- LAM * exp(-alpha * z)
    h_post <- LAM * exp( alpha * z)
    
    ifelse(
      t_grid <= tau,
      exp(-h_pre * t_grid),
      exp(-h_pre * tau - h_post * (t_grid - tau))
    )
  })
  
  matplot(
    t_grid,
    surv_mat,
    type = "l",
    lty = 1,
    lwd = 2.5,
    col = cols,
    ylim = c(0, 1),
    xlab = expression(t),
    ylab = expression(S(t ~ "|" ~ u[1])),
    main = " "
  )
  
  # change point
  abline(v = tau, lty = 2, lwd = 2, col = "firebrick")
  
  text(
    tau,
    0.98,
    labels = expression(t == tau),
    pos = 4,
    col = "firebrick",
    cex = 0.9
  )
  
  legend(
    "topright",
    legend = as.expression(lapply(z_vals, function(v) bquote(u[1] == .(v)))),
    col = cols,
    lwd = 2.5,
    lty = 1,
    bty = "n",
    title = expression(u[1])
  )
  
})

cat("\nAll interpretability plots saved in:", normalizePath(PLOT_DIR), "\n")
cat("Files produced:\n")
cat("  interp_1_greedy_selection_path.pdf        — AIC steps with tau marked\n")
cat("  interp_2_km_curves_selected_times.pdf     — S(t|cat) curves + selected times\n")
cat("  interp_3_km_embedding_2d.pdf              — 2D embedding coloured by latent z\n")
cat("  interp_4_coefficient_comparison.pdf       — OHE beta vs KM LP vs true effects\n")
cat("  interp_5_martingale_residuals.pdf         — residuals by category\n")
cat("  interp_6_survival_reversal.pdf            — S reversal at selected times\n")