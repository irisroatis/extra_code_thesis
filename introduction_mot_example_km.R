# motivating example code for KM



# Required package: survival
# install.packages("survival")

suppressPackageStartupMessages(library(survival))

set.seed(5)

#  parameter settings

N              <- 200
LAMBDA0        <- 0.025
BETA_X         <- 0.50
BETA_X1_C12   <- log(2.00)  # C_{1,2} vs C_{1,1}
BETA_X1_C13   <- log(1.20)  # C_{1,3} vs C_{1,1}
BETA_X2_C22   <- log(0.1)  # C_{2,2} vs C_{2,1}
BETA_X2_C23   <- log(1.50)  # C_{2,3} vs C_{2,1}
CENSOR_RATE    <- 0.012
CANDIDATE_TIMES <- c(20, 35)
MIN_AIC_IMPROVE <- 0


# name of each category

CAT1_LABELS <- expression(
  C[1*","*1],
  C[1*","*2],
  C[1*","*3]
)

CAT2_LABELS <- expression(
  C[2*","*1],
  C[2*","*2],
  C[2*","*3]
)

CATEGORY_LABELS <- list(
  cat1 = CAT1_LABELS,
  cat2 = CAT2_LABELS
)

OUT_DIR <- "km_greedy_two_three_level_ph_output"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

CATEGORY_COLOURS <- list(
  cat1 = c("C_{1,1}" = "#1F77B4", "C_{1,2}" = "#D62728", "C_{1,3}" = "#FF7F0E"),
  cat2 = c("C_{2,1}" = "#2CA02C", "C_{2,2}" = "#9467BD", "C_{2,3}" = "#8C564B")
)


predict_surv_at <- function(sf, times) {
  as.numeric(summary(sf, times = times, extend = TRUE)$surv)
}

# define aic
cox_aic <- function(fit) {
  ll <- as.numeric(fit$loglik[2])
  k  <- sum(is.finite(coef(fit)) & !is.na(coef(fit)))
  -2 * ll + 2 * k
}

safe_coxph <- function(df) {
  try(
    coxph(
      Surv(Time, Event) ~ .,
      data = df,
      ties = "efron",
      singular.ok = FALSE,
      x = TRUE,
      y = TRUE,
      model = TRUE
    ),
    silent = TRUE
  )
}

save_plot <- function(filename, expr, width = 10, height = 7) {
  png(file.path(OUT_DIR, filename), width = width, height = height,
      units = "in", res = 180)
  on.exit(dev.off(), add = TRUE)
  eval.parent(substitute(expr))
}


# PH data generating
# h(t|x,x1,x2) = lambda0 * exp(beta_x*x + beta_12*I(x1=C_{1,2}) +
#                                      beta_13*I(x1=C_{1,3}) + beta_22*I(x2=C_{2,2}) +
#                                      beta_23*I(x2=C_{2,3}))

simulate_ph_data <- function(n) {
  x    <- rnorm(n)
  cat1 <- factor(sample(c("C_{1,1}", "C_{1,2}", "C_{1,3}"), n, replace = TRUE),
                 levels = c("C_{1,1}", "C_{1,2}", "C_{1,3}"))
  cat2 <- factor(sample(c("C_{2,1}", "C_{2,2}", "C_{2,3}"), n, replace = TRUE),
                 levels = c("C_{2,1}", "C_{2,2}", "C_{2,3}"))

  lp <- BETA_X * x +
    BETA_X1_C12 * as.integer(cat1 == "C_{1,2}") +
    BETA_X1_C13 * as.integer(cat1 == "C_{1,3}") +
    BETA_X2_C22 * as.integer(cat2 == "C_{2,2}") +
    BETA_X2_C23 * as.integer(cat2 == "C_{2,3}")

  event_time <- -log(runif(n)) / (LAMBDA0 * exp(lp))
  censor_time <- rexp(n, rate = CENSOR_RATE)

  data.frame(
    x = x,
    cat1 = cat1,
    cat2 = cat2,
    Time = pmin(event_time, censor_time),
    Event = as.integer(event_time <= censor_time),
    true_lp = lp,
    event_time = event_time,
    censor_time = censor_time
  )
}

dat <- simulate_ph_data(N)
model_dat <- dat[, c("x", "cat1", "cat2", "Time", "Event")]
cat_cols <- c("cat1", "cat2")
cont_cols <- "x"

cat("N:", nrow(dat), "\n")
cat("Events:", sum(dat$Event), "\n")
cat("Event fraction:", round(mean(dat$Event), 3), "\n")


#  Precompute category-level KM values
precompute_km <- function(train_df, cat_cols, candidate_times) {
  marginal_fit <- survfit(Surv(Time, Event) ~ 1, data = train_df)
  fallback <- predict_surv_at(marginal_fit, candidate_times)

  out <- vector("list", length(cat_cols))
  names(out) <- cat_cols

  for (cc in cat_cols) {
    levels_cc <- levels(train_df[[cc]])
    lookup <- vector("list", length(levels_cc))
    names(lookup) <- levels_cc

    for (lv in levels_cc) {
      sub <- train_df[train_df[[cc]] == lv, , drop = FALSE]
      if (nrow(sub) == 0 || sum(sub$Event) == 0) {
        lookup[[lv]] <- fallback
      } else {
        lookup[[lv]] <- predict_surv_at(
          survfit(Surv(Time, Event) ~ 1, data = sub),
          candidate_times
        )
      }
    }

    out[[cc]] <- list(
      candidate_times = candidate_times,
      lookup = lookup,
      fallback = fallback
    )
  }

  out
}

build_km_matrices <- function(df, precomp, cat_cols) {
  out <- vector("list", length(cat_cols))
  names(out) <- cat_cols

  for (cc in cat_cols) {
    pc <- precomp[[cc]]
    labels <- as.character(df[[cc]])

    mat <- t(vapply(labels, function(lv) {
      value <- pc$lookup[[lv]]
      if (is.null(value)) value <- pc$fallback
      value
    }, numeric(length(pc$candidate_times))))

    colnames(mat) <- paste0(cc, "_KM_t", pc$candidate_times)
    out[[cc]] <- mat
  }

  out
}

assemble_km_matrix <- function(full_mats, selected, cat_cols) {
  pieces <- list()
  for (cc in cat_cols) {
    idx <- selected[[cc]]
    if (length(idx) > 0) {
      pieces[[length(pieces) + 1]] <- full_mats[[cc]][, idx, drop = FALSE]
    }
  }

  if (length(pieces) == 0) {
    return(matrix(numeric(0), nrow = nrow(full_mats[[cat_cols[1]]]), ncol = 0))
  }
  do.call(cbind, pieces)
}

km_precomp <- precompute_km(model_dat, cat_cols, CANDIDATE_TIMES)
full_km <- build_km_matrices(model_dat, km_precomp, cat_cols)




# Category-level encoding table
encoding_table <- do.call(rbind, lapply(cat_cols, function(cc) {
  do.call(rbind, lapply(names(km_precomp[[cc]]$lookup), function(lv) {
    vals <- km_precomp[[cc]]$lookup[[lv]]
    data.frame(
      feature = cc,
      level = lv,
      S_t1 = vals[1],
      S_t2 = vals[2],
      stringsAsFactors = FALSE
    )
  }))
}))

write.csv(encoding_table,
          file.path(OUT_DIR, "category_km_encodings.csv"),
          row.names = FALSE)


#  forward selection (using AIC)

km_greedy_aic <- function(train_df, cont_cols, cat_cols,
                          candidate_times, full_mats,
                          min_aic_improve = 0) {
  selected <- setNames(lapply(cat_cols, function(x) integer(0)), cat_cols)
  base_df <- train_df[, c(cont_cols, "Time", "Event"), drop = FALSE]

  base_fit <- safe_coxph(base_df)
  if (inherits(base_fit, "try-error")) stop("Continuous-only Cox model failed.")
  current_aic <- cox_aic(base_fit)

  selection_log <- data.frame()
  evaluation_log <- data.frame()
  step <- 1L

  cat("STEP 0: CONTINUOUS-ONLY MODEL\n")
  cat("Model: Surv(Time, Event) ~ ", paste(cont_cols, collapse = " + "), "\n", sep = "")
  cat(sprintf("AIC = %.4f\n", current_aic))

  repeat {
    current_X <- assemble_km_matrix(full_mats, selected, cat_cols)
    selected_names <- if (ncol(current_X) == 0) character(0) else colnames(current_X)
    current_model <- paste(c(cont_cols, selected_names), collapse = " + ")
    step_results <- data.frame()

    for (cc in cat_cols) {
      remaining <- setdiff(seq_along(candidate_times), selected[[cc]])

      for (j in remaining) {
        candidate_col <- full_mats[[cc]][, j, drop = FALSE]
        X_trial <- if (ncol(current_X) == 0) candidate_col else cbind(current_X, candidate_col)

        trial_df <- cbind(base_df, as.data.frame(X_trial, check.names = FALSE))
        fit <- safe_coxph(trial_df)

        if (inherits(fit, "try-error")) {
          trial_aic <- Inf
          delta <- NA_real_
          status <- "singular/redundant"
        } else {
          trial_aic <- cox_aic(fit)
          delta <- current_aic - trial_aic
          status <- "estimable"
        }

        step_results <- rbind(step_results, data.frame(
          step = step,
          current_model = current_model,
          feature = cc,
          time_idx = j,
          time = candidate_times[j],
          candidate = paste0(cc, " @ t=", candidate_times[j]),
          candidate_column = colnames(candidate_col),
          aic_before = current_aic,
          trial_aic = trial_aic,
          delta_aic = delta,
          status = status,
          stringsAsFactors = FALSE
        ))
      }
    }

    if (nrow(step_results) == 0) break

    step_results$chosen <- FALSE
    step_results$accepted <- FALSE
    finite_rows <- which(is.finite(step_results$trial_aic))

    cat(sprintf("STEP %d: TEST EVERY REMAINING CANDIDATE\n", step))
    cat("Current model: ", current_model, "\n", sep = "")
    cat(sprintf("Current AIC: %.4f\n\n", current_aic))

    if (length(finite_rows) == 0) {
      evaluation_log <- rbind(evaluation_log, step_results)
      print(step_results[, c("candidate", "aic_before", "trial_aic", "delta_aic", "status")], row.names = FALSE)
      cat("Decision: all remaining candidates are singular/redundant. Stop.\n")
      break
    }

    winner_row <- finite_rows[which.min(step_results$trial_aic[finite_rows])]
    step_results$chosen[winner_row] <- TRUE
    best <- step_results[winner_row, ]
    accepted <- is.finite(best$delta_aic) && best$delta_aic > min_aic_improve
    step_results$accepted[winner_row] <- accepted
    evaluation_log <- rbind(evaluation_log, step_results)

    display_table <- step_results[, c("candidate", "aic_before", "trial_aic", "delta_aic", "status", "chosen", "accepted")]
    display_table$aic_before <- round(display_table$aic_before, 4)
    display_table$trial_aic <- ifelse(is.finite(display_table$trial_aic), round(display_table$trial_aic, 4), Inf)
    display_table$delta_aic <- round(display_table$delta_aic, 4)
    print(display_table, row.names = FALSE)

    if (!accepted) {
      cat(sprintf(
        "\nDecision: %s has the smallest trial AIC (%.4f), but its improvement %.4f is not greater than %.4f. Stop.\n",
        best$candidate, best$trial_aic, best$delta_aic, min_aic_improve
      ))
      break
    }

    cat(sprintf(
      "\nDecision: SELECT %s because it gives the smallest AIC: %.4f.\n",
      best$candidate, best$trial_aic
    ))
    cat(sprintf(
      "AIC reduction: %.4f - %.4f = %.4f.\n",
      current_aic, best$trial_aic, best$delta_aic
    ))

    selected[[best$feature]] <- c(selected[[best$feature]], best$time_idx)
    selection_log <- rbind(selection_log, data.frame(
      step = step,
      feature = best$feature,
      time_idx = best$time_idx,
      time = best$time,
      selected_column = best$candidate,
      aic_before = current_aic,
      aic_after = best$trial_aic,
      delta_aic = best$delta_aic,
      stringsAsFactors = FALSE
    ))

    current_aic <- best$trial_aic
    step <- step + 1L
  }

  cat("FINAL GREEDY MODEL\n")
  final_names <- unlist(lapply(cat_cols, function(cc) {
    idx <- selected[[cc]]
    if (length(idx) == 0) character(0) else colnames(full_mats[[cc]])[idx]
  }))
  cat("Model: ", paste(c(cont_cols, final_names), collapse = " + "), "\n", sep = "")
  cat(sprintf("Final AIC = %.4f\n", current_aic))

  list(
    selected = selected,
    selection_log = selection_log,
    evaluation_log = evaluation_log,
    base_aic = cox_aic(base_fit),
    final_aic = current_aic
  )
}

greedy <- km_greedy_aic(
  train_df = model_dat,
  cont_cols = cont_cols,
  cat_cols = cat_cols,
  candidate_times = CANDIDATE_TIMES,
  full_mats = full_km,
  min_aic_improve = MIN_AIC_IMPROVE
)

write.csv(greedy$selection_log,
          file.path(OUT_DIR, "greedy_selection_log.csv"),
          row.names = FALSE)
write.csv(greedy$evaluation_log,
          file.path(OUT_DIR, "greedy_all_candidate_evaluations.csv"),
          row.names = FALSE)

#  table showing every AIC comparison
step_by_step_aic <- greedy$evaluation_log[, c(
  "step", "current_model", "candidate", "aic_before", "trial_aic",
  "delta_aic", "status", "chosen", "accepted"
)]
write.csv(step_by_step_aic,
          file.path(OUT_DIR, "greedy_step_by_step_aic.csv"),
          row.names = FALSE)
capture.output(
  {
    cat(sprintf("Step 0: continuous-only AIC = %.4f\n\n", greedy$base_aic))
    for (s in sort(unique(step_by_step_aic$step))) {
      z <- step_by_step_aic[step_by_step_aic$step == s, ]
      cat(sprintf("Step %d; current model: %s\n", s, z$current_model[1]))
      print(z[, c("candidate", "aic_before", "trial_aic", "delta_aic", "status", "chosen", "accepted")], row.names = FALSE)
      cat("\n")
    }
    cat(sprintf("Final AIC = %.4f\n", greedy$final_aic))
  },
  file = file.path(OUT_DIR, "greedy_step_by_step_aic.txt")
)

cat("\nGreedy selection log:\n")
print(greedy$selection_log)

# Fit final and comparator models

X_selected <- assemble_km_matrix(full_km, greedy$selected, cat_cols)
km_model_df <- cbind(
  model_dat[, c("x", "Time", "Event")],
  as.data.frame(X_selected, check.names = FALSE)
)

fit_continuous <- coxph(Surv(Time, Event) ~ x, data = model_dat,
                        ties = "efron", x = TRUE, y = TRUE, model = TRUE)
fit_ohe <- coxph(Surv(Time, Event) ~ x + cat1 + cat2, data = model_dat,
                 ties = "efron", x = TRUE, y = TRUE, model = TRUE)
fit_km <- coxph(Surv(Time, Event) ~ ., data = km_model_df,
                ties = "efron", singular.ok = FALSE,
                x = TRUE, y = TRUE, model = TRUE)

model_comparison <- data.frame(
  model = c("Continuous only", "One-hot Cox", "KM greedy Cox"),
  n_parameters = c(length(coef(fit_continuous)), length(coef(fit_ohe)), length(coef(fit_km))),
  AIC = c(cox_aic(fit_continuous), cox_aic(fit_ohe), cox_aic(fit_km)),
  concordance = c(summary(fit_continuous)$concordance[1],
                  summary(fit_ohe)$concordance[1],
                  summary(fit_km)$concordance[1])
)
write.csv(model_comparison,
          file.path(OUT_DIR, "model_comparison.csv"),
          row.names = FALSE)

capture.output(
  list(
    data_summary = c(N = N, events = sum(dat$Event), event_fraction = mean(dat$Event)),
    true_coefficients = c(
      beta_x = BETA_X,
      beta_x1_C12 = BETA_X1_C12,
      beta_x1_C13 = BETA_X1_C13,
      beta_x2_C22 = BETA_X2_C22,
      beta_x2_C23 = BETA_X2_C23
    ),
    encoding_table = encoding_table,
    greedy_selection = greedy$selection_log,
    final_km_model = summary(fit_km),
    one_hot_model = summary(fit_ohe),
    model_comparison = model_comparison
  ),
  file = file.path(OUT_DIR, "experiment_summary.txt")
)


# Plot 1: true PH survival curves

save_plot("01_true_ph_survival_curves.png", {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
  t_grid <- seq(0, 60, length.out = 300)

  # Hold x=0 and the other categorical feature at its reference level
  s_11 <- exp(-LAMBDA0 * t_grid)
  s_12 <- exp(-LAMBDA0 * exp(BETA_X1_C12) * t_grid)
  s_13 <- exp(-LAMBDA0 * exp(BETA_X1_C13) * t_grid)
  plot(t_grid, s_11, type = "l", lwd = 3, lty = 1,
       col = CATEGORY_COLOURS$cat1["C_{1,1}"], ylim = c(0, 1),
       xlab = "Time", ylab = "True survival probability",
       main = expression(x[1]*": proportional hazards"))
  lines(t_grid, s_12, lwd = 3, lty = 2, col = CATEGORY_COLOURS$cat1["C_{1,2}"])
  lines(t_grid, s_13, lwd = 3, lty = 3, col = CATEGORY_COLOURS$cat1["C_{1,3}"])
  abline(v = CANDIDATE_TIMES, lty = 3)
  legend(
    "topright",
    legend = c(
      expression(C[1*","*1] * " (reference)"),
      CAT1_LABELS[[2]],
      CAT1_LABELS[[3]]
    ),
    lty = 1:3,
    lwd = 3,
    col = unname(CATEGORY_COLOURS$cat1),
    text.col = "black",
    bty = "n"
  )

  s_21 <- exp(-LAMBDA0 * t_grid)
  s_22 <- exp(-LAMBDA0 * exp(BETA_X2_C22) * t_grid)
  s_23 <- exp(-LAMBDA0 * exp(BETA_X2_C23) * t_grid)
  plot(t_grid, s_21, type = "l", lwd = 3, lty = 1,
       col = CATEGORY_COLOURS$cat2["C_{2,1}"], ylim = c(0, 1),
       xlab = "Time", ylab = "True survival probability",
       main = expression(x[2]*": proportional hazards"))
  lines(t_grid, s_22, lwd = 3, lty = 2, col = CATEGORY_COLOURS$cat2["C_{2,2}"])
  lines(t_grid, s_23, lwd = 3, lty = 3, col = CATEGORY_COLOURS$cat2["C_{2,3}"])
  abline(v = CANDIDATE_TIMES, lty = 3)
  legend(
    "topright",
    legend = c(
      expression(C[2*","*1] * " (reference)"),
      CAT2_LABELS[[2]],
      CAT2_LABELS[[3]]
    ),
    lty = 1:3,
    lwd = 3,
    col = unname(CATEGORY_COLOURS$cat2),
    text.col = "black",
    bty = "n"
  )
}, width = 12, height = 5.5)



# Plot 2a and 2b: empirical KM curves


plot_empirical_km <- function(cc, filename) {
  
  save_plot(filename, {
    
    par(
      mar = c(5, 5, 4, 2),
      col.axis = "black",
      col.lab = "black",
      col.main = "black"
    )
    
    feature_idx <- match(cc, cat_cols)
    
    form <- as.formula(
      paste0("Surv(Time, Event) ~ ", cc)
    )
    
    sf <- survfit(
      form,
      data = model_dat
    )
    
    levels_cc <- levels(model_dat[[cc]])
    
    cols <- unname(
      CATEGORY_COLOURS[[cc]][levels_cc]
    )
    
    plot(
      sf,
      lwd = 2,
      lty = seq_along(levels_cc),
      col = cols,
      mark.time = FALSE,
      conf.int = FALSE,
      xlim = c(0, 60),
      ylim = c(0, 1),
      xaxt = "n",
      xlab = "Time",
      ylab = "KM curves",
      #  main = bquote(
      #    x[.(feature_idx)] * ": KM curves"
      #  ),
      main = " ",
      col.axis = "black",
      col.lab = "black",
      col.main = "black"
    )
    
    # Draw the ordinary numerical x-axis ticks, excluding
    # the candidate times 10 and 30.
    regular_ticks <- pretty(c(0, 60))
    regular_ticks <- setdiff(
      regular_ticks,
      CANDIDATE_TIMES
    )
    
    axis(
      side = 1,
      at = regular_ticks,
      labels = regular_ticks,
      col = "black",
      col.axis = "black"
    )
    
    # Display the candidate times as t_1 and t_2.
    axis(
      side = 1,
      at = CANDIDATE_TIMES,
      labels = expression(t[1], t[2]),
      col = "black",
      col.axis = "black"
    )
    
    abline(
      v = CANDIDATE_TIMES,
      lty = 3,
      col = "black"
    )
    
    for (i in seq_along(levels_cc)) {
      
      vals <- km_precomp[[cc]]$lookup[[levels_cc[i]]]
      
      points(
        CANDIDATE_TIMES,
        vals,
        pch = 19,
        cex = 1.2,
        col = cols[i]
      )
      
      for (j in seq_along(CANDIDATE_TIMES)) {
        
        label <- bquote(
          hat(S)(
            t[.(j)] ~ "|" ~
              x[.(feature_idx)] ==
              C[.(feature_idx) * "," * .(i)]
          )
        )
        
        label_position <- if (i == 1) {
          3
        } else if (i == 2) {
          1
        } else {
          1
        }
        
        text(
          CANDIDATE_TIMES[j],
          vals[j],
          labels = label,
          pos = label_position,
          offset = 1,
          cex = 0.9,
          col = "black",
          xpd = NA
        )
      }
    }
    
    legend(
      "topright",
      legend = CATEGORY_LABELS[[cc]],
      lty = seq_along(levels_cc),
      col = cols,
      lwd = 2,
      text.col = "black",
      bty = "n"
    )
    
  }, width = 8, height = 6)
}


# Separate KM plot for x_1
plot_empirical_km(
  cc = "cat1",
  filename = "02a_empirical_km_curves_x1.png"
)


# Separate KM plot for x_2
plot_empirical_km(
  cc = "cat2",
  filename = "02b_empirical_km_curves_x2.png"
)



# Plot 3: two-time KM category embedding

save_plot("03_two_time_km_embedding.png", {
  par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

  for (cc in cat_cols) {
    feature_idx <- match(cc, cat_cols)
    sub <- encoding_table[encoding_table$feature == cc, ]
    cols <- unname(CATEGORY_COLOURS[[cc]][sub$level])
    category_labels <- CATEGORY_LABELS[[cc]]

    plot(
      sub$S_t1,
      sub$S_t2,
      pch = 19,
      cex = 2,
      col = cols,
      xlim = range(c(sub$S_t1, 0, 1)),
      ylim = range(c(sub$S_t2, 0, 1)),
      xlab = bquote(hat(S)(t[1] ~ "|" ~ x[.(feature_idx)])),
      ylab = bquote(hat(S)(t[2] ~ "|" ~ x[.(feature_idx)])),
      main = bquote(x[.(feature_idx)] * ": two candidate dimensions")
    )

    text(
      sub$S_t1,
      sub$S_t2,
      labels = category_labels,
      pos = 3,
      cex = 1.1,
      col = "black"
    )

    if (nrow(sub) > 1) {
      pairs_idx <- combn(seq_len(nrow(sub)), 2)
      apply(pairs_idx, 2, function(idx) {
        segments(
          sub$S_t1[idx[1]],
          sub$S_t2[idx[1]],
          sub$S_t1[idx[2]],
          sub$S_t2[idx[2]],
          lty = 3
        )
      })
    }

    mtext(
      "Three levels can require two KM coordinates for a full-rank encoding",
      side = 3,
      line = 0.3,
      cex = 0.75
    )
  }
}, width = 12, height = 5.5)


# Plot 4: full AIC landscape by greedy step

save_plot("04_greedy_aic_selection.png", {
  eval_log <- greedy$evaluation_log
  n_steps <- length(unique(eval_log$step))
  par(mfrow = c(n_steps, 1), mar = c(7, 5, 4, 2))

  candidate_cols <- c(
    "cat1 @ t=10" = "#1F77B4",
    "cat1 @ t=30" = "#6BAED6",
    "cat2 @ t=10" = "#2CA02C",
    "cat2 @ t=30" = "#98DF8A"
  )

  for (s in sort(unique(eval_log$step))) {
    sub <- eval_log[eval_log$step == s, ]
    finite_aic <- sub$trial_aic[is.finite(sub$trial_aic)]
    plot_vals <- sub$trial_aic
    # Put redundant candidates just above the visible plotting range.
    if (length(finite_aic) > 0) {
      plot_vals[!is.finite(plot_vals)] <- max(c(finite_aic, sub$aic_before), na.rm = TRUE) + 2
    }

    lower <- min(c(finite_aic, sub$aic_before), na.rm = TRUE) - 1
    upper <- max(c(plot_vals, sub$aic_before), na.rm = TRUE) + 2
    bar_cols <- unname(candidate_cols[sub$candidate])

    bp <- barplot(
      plot_vals,
      names.arg = sub$candidate,
      col = bar_cols,
      border = NA,
      las = 2,
      ylim = c(lower, upper),
      ylab = "AIC after adding candidate",
      main = paste0("Greedy step ", s, ": current AIC = ", sprintf("%.2f", sub$aic_before[1]))
    )
    abline(h = sub$aic_before[1], lty = 2, lwd = 2)

    for (i in seq_len(nrow(sub))) {
      if (sub$status[i] != "estimable") {
        lab <- "redundant"
      } else {
        lab <- sprintf("AIC %.2f\nDelta %.2f", sub$trial_aic[i], sub$delta_aic[i])
      }
      text(bp[i], plot_vals[i], labels = lab, pos = 3, cex = 0.78)
      if (isTRUE(sub$accepted[i])) {
        points(bp[i], plot_vals[i], pch = 8, cex = 1.8)
      }
    }
    legend("topright", c("Current-model AIC", "Accepted candidate"),
           lty = c(2, NA), pch = c(NA, 8), lwd = c(2, NA), bty = "n")
  }
}, width = 11, height = max(6, 5 * length(unique(greedy$evaluation_log$step))))


# Plot 5: selected design matrix
save_plot("05_selected_km_design.png", {
  mat <- matrix(0L, nrow = length(cat_cols), ncol = length(CANDIDATE_TIMES),
                dimnames = list(cat_cols, paste0("t=", CANDIDATE_TIMES)))

  for (i in seq_len(nrow(greedy$selection_log))) {
    r <- match(greedy$selection_log$feature[i], cat_cols)
    c <- match(greedy$selection_log$time[i], CANDIDATE_TIMES)
    mat[r, c] <- 1L
  }

  image(x = seq_len(ncol(mat)), y = seq_len(nrow(mat)),
        z = t(mat[nrow(mat):1, , drop = FALSE]),
        axes = FALSE, xlab = "Candidate time", ylab = "Categorical feature",
        main = "Final KM dimensions selected by greedy AIC")
  axis(1, at = seq_len(ncol(mat)), labels = colnames(mat))
  axis(
    2,
    at = seq_len(nrow(mat)),
    labels = expression(x[2], x[1]),
    las = 1
  )
  box()

  mat_plot <- mat[nrow(mat):1, , drop = FALSE]
  for (r in seq_len(nrow(mat_plot))) {
    for (c in seq_len(ncol(mat_plot))) {
      text(c, r, if (mat_plot[r, c] == 1) "SELECTED" else "REDUNDANT", cex = 0.9)
    }
  }
}, width = 8, height = 5.5)


# Plot 6: model comparison

save_plot("06_model_comparison.png", {
  par(mfrow = c(1, 2), mar = c(7, 5, 4, 2))

  barplot(model_comparison$AIC,
          names.arg = model_comparison$model,
          las = 2, ylab = "AIC", main = "Model fit")

  barplot(model_comparison$concordance,
          names.arg = model_comparison$model,
          las = 2, ylim = c(0.5, 1),
          ylab = "Harrell concordance", main = "Discrimination")
  abline(h = 0.5, lty = 2)
}, width = 12, height = 5.5)




cat("\nFinished. Outputs written to:\n", normalizePath(OUT_DIR), "\n")
cat("\nKey interpretation:\n")
cat("- The DGP is exactly proportional hazards.\n")
cat("- Each categorical feature has three levels and therefore three KM values at each time.\n")
cat("- A single KM column supplies one numeric contrast among the three levels.\n")
cat("- A second KM time can add an independent contrast and need not be redundant.\n")
cat("- Greedy AIC may retain up to two KM times per categorical feature in this experiment.\n")
