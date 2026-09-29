#======== PROJECT COM2LIFE ========
## Modelisation GLMM on alpha diversity estimates in "estimates_alpha_diversity.R" 
## First load the libraries. Do a correlation matrix to chosse environmeental variables 
##### Here, we choose 5 variables 

# Load metadata object
metadata_alpha <- readRDS(here::here("data",
                                     "metadata_alpha.rds"))

# ---------------------------------------------------------------- 0. CONFIG
set.seed(42)
cfg <- list(
  responses    = c("taxo_q0", "taxo_q1", "phylo_q0", "phylo_q1"),
  predictors   = c("Temperature_median", "Chla_median", "Salinity_median",
                   "pH_median", "Oxygen_median"),
  factor_term  = "species",
  random       = "(1|lake)",
  candidate_df = 1:4,      # spline df tried per predictor (1 = linear)
  delta_aic    = 2,        # keep the simplest basis within 2 AIC of the best
  alpha        = 0.05,
  out_dir      = here::here("output", "glmm")
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

metadata_alpha <- readRDS(here::here("data", "metadata_alpha.rds")) |>
  mutate(species = factor(origin), lake = factor(lake))

# z-score predictors: much better convergence for glmmTMB and comparable effect sizes
metadata_alpha <- metadata_alpha |>
  mutate(across(all_of(cfg$predictors), ~ as.numeric(scale(.x))))

# ------------------------------------------- 1. CORRELATION & COLLINEARITY
cand_vars <- c("Temperature_median", "Oxygen_median", "pH_median", "Chla_median",
               "Salinity_median", "TPC", "TPN", "TPC.TPN", "NH4", "PO4", "NO3_N02")

cor_pmat <- function(mat) {
  mat <- as.matrix(mat); n <- ncol(mat)
  p <- matrix(0, n, n, dimnames = list(colnames(mat), colnames(mat)))
  for (i in 1:(n - 1)) for (j in (i + 1):n)
    p[i, j] <- p[j, i] <- suppressWarnings(
      cor.test(mat[, i], mat[, j], method = "spearman")$p.value)
  p
}
cor_mat <- cor(metadata_alpha[, cand_vars], method = "spearman", use = "pairwise.complete.obs")
pdf(file.path(cfg$out_dir, "correlation_matrix.pdf"), width = 8, height = 8)
corrplot::corrplot(cor_mat, type = "upper", order = "hclust", method = "number",
                   p.mat = cor_pmat(metadata_alpha[, cand_vars]),
                   sig.level = cfg$alpha, insig = "blank")
dev.off()

vif_tab <- performance::check_collinearity(
  lm(reformulate(c(cfg$predictors, cfg$factor_term), "taxo_q0"), data = metadata_alpha))
print(vif_tab)   # VIF < 5 (conservative) or < 10 is acceptable

# ---------------------------------------------------------- 2. HELPERS
term_str <- function(x, df) if (df == 1) x else sprintf("splines::ns(%s, df = %d)", x, df)

make_formula <- function(resp, dfs) {
  rhs <- c(unname(mapply(term_str, names(dfs), dfs)), cfg$factor_term, cfg$random)
  reformulate(rhs, response = resp)
}
fit_glmm <- function(f, data)
  glmmTMB(f, data = data, family = gaussian,
          control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))
ok_fit   <- function(m) !is.null(m) && isTRUE(m$sdr$pdHess)

lin_dfs <- setNames(rep(1L, length(cfg$predictors)), cfg$predictors)

# ------------------------------------------- 3. SPLINE BASIS SELECTION
# Coordinate-wise search: for each predictor try each df (others fixed),
# keep the simplest basis within delta_aic of the best AIC.
select_basis <- function(resp, data) {
  dfs <- lin_dfs; tab <- list()
  for (p in cfg$predictors) {
    aic <- sapply(cfg$candidate_df, function(k) {
      d <- dfs; d[p] <- k
      m <- tryCatch(suppressWarnings(fit_glmm(make_formula(resp, d), data)),
                    error = function(e) NULL)
      if (ok_fit(m)) AIC(m) else NA_real_
    })
    keep   <- which(aic <= min(aic, na.rm = TRUE) + cfg$delta_aic)[1]
    dfs[p] <- cfg$candidate_df[keep]
    tab[[p]] <- data.frame(response = resp, predictor = p, df = cfg$candidate_df,
                           AIC = aic, dAIC = aic - min(aic, na.rm = TRUE),
                           selected = cfg$candidate_df == dfs[p])
  }
  list(dfs = dfs, table = do.call(rbind, tab))
}

# ---------------------------------------------------- 4. PERFORMANCE
safe <- function(expr) tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) NA_real_)
num1 <- function(x, nm) if (is.list(x) && !is.null(x[[nm]])) as.numeric(x[[nm]][1]) else NA_real_

perf_table <- function(m, label, resp) {
  vc_lake <- if (!is.null(VarCorr(m)$cond$lake)) as.numeric(VarCorr(m)$cond$lake[1, 1]) else NA_real_
  r2  <- safe(performance::r2_nakagawa(m))
  icc <- safe(performance::icc(m))
  data.frame(response = resp, model = label, AIC = AIC(m), BIC = BIC(m),
             logLik = as.numeric(logLik(m)),
             R2_marginal    = num1(r2, "R2_marginal"),
             R2_conditional = num1(r2, "R2_conditional"),
             ICC   = num1(icc, "ICC_adjusted"),
             RMSE  = safe(performance::performance_rmse(m)), sigma = sigma(m),
             var_lake = vc_lake,
             converged = ok_fit(m),                              # positive-definite Hessian
             singular  = !is.na(vc_lake) && vc_lake < 1e-6)      # lake variance ~ 0
}

# Leave-one-lake-out CV: how well does the model predict a NEW lake?
cv_lolo <- function(f, data, label, resp) {
  out <- lapply(levels(data$lake), function(l) {
    tr <- droplevels(data[data$lake != l, ]); te <- data[data$lake == l, ]
    te <- te[te$species %in% levels(tr$species), ]
    if (!nrow(te)) return(NULL)
    m <- tryCatch(suppressWarnings(fit_glmm(f, tr)), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    data.frame(obs = te[[resp]],
               pred = predict(m, newdata = te, re.form = NA, allow.new.levels = TRUE))
  }) |> bind_rows()
  data.frame(response = resp, model = label,
             CV_RMSE = sqrt(mean((out$obs - out$pred)^2)),
             CV_MAE  = mean(abs(out$obs - out$pred)),
             CV_R2   = 1 - sum((out$obs - out$pred)^2) / sum((out$obs - mean(out$obs))^2))
}

# ------------------------------------------- 5. RESIDUAL DIAGNOSTICS
run_diagnostics <- function(m, data, resp) {
  res <- DHARMa::simulateResiduals(m, n = 1000, plot = FALSE)
  tests <- data.frame(
    response = resp,
    test = c("KS uniformity", "Dispersion", "Outliers", "RE normality (Shapiro)"),
    p_value = c(DHARMa::testUniformity(res, plot = FALSE)$p.value,
                DHARMa::testDispersion(res, plot = FALSE)$p.value,
                DHARMa::testOutliers(res, plot = FALSE)$p.value,
                shapiro.test(ranef(m)$cond$lake[, 1])$p.value))
  tests$pass <- tests$p_value > cfg$alpha
  
  pdf(file.path(cfg$out_dir, paste0("diagnostics_", resp, ".pdf")), width = 10, height = 8)
  plot(res)
  for (p in cfg$predictors) DHARMa::plotResiduals(res, data[[p]], xlab = p)
  DHARMa::plotResiduals(res, data$species, xlab = "species")
  qqnorm(ranef(m)$cond$lake[, 1], main = "Random intercepts (lake)"); qqline(ranef(m)$cond$lake[, 1])
  dev.off()
  tests
}

# ------------------------------------------- 6. SIGNIFICANCE (LRT)
term_tests <- function(m, data, resp) {
  labs <- attr(terms(lme4::nobars(formula(m))), "term.labels")
  out <- lapply(labs, function(lab) {
    m0 <- tryCatch(update(m, formula = as.formula(paste(". ~ . -", lab)), data = data),
                   error = function(e) NULL)
    if (is.null(m0)) return(NULL)
    a <- anova(m0, m)
    data.frame(response = resp, term = lab,
               predictor = sub("^splines::ns\\(([A-Za-z0-9_.]+),.*$", "\\1", lab),
               df = a[2, "Chi Df"], Chisq = a[2, "Chisq"], p_value = a[2, "Pr(>Chisq)"],
               dAIC_if_dropped = AIC(m0) - AIC(m))
  })
  bind_rows(out) |> mutate(p_BH = p.adjust(p_value, "BH"))
}

# ---------------------------------------------------------- 8. MAIN LOOP
results <- lapply(cfg$responses, function(resp) {
  message("==== ", resp, " ====")
  dat <- droplevels(metadata_alpha[complete.cases(
    metadata_alpha[, c(resp, cfg$predictors, "species", "lake")]), ])
  
  # 8.1 spline basis
  sel <- select_basis(resp, dat)
  f_list <- list(null     = reformulate(c(cfg$factor_term, cfg$random), response = resp),
                 linear   = make_formula(resp, lin_dfs),
                 selected = make_formula(resp, sel$dfs))
  if (all(sel$dfs == 1)) f_list$selected <- NULL       # selected == linear
  mods  <- lapply(f_list, fit_glmm, data = dat)
  
  # best model = simplest converged model within delta_aic of the lowest AIC
  # (order in `mods` is null < linear < selected, i.e. simplest first)
  aic_all <- sapply(mods, function(m) if (ok_fit(m)) AIC(m) else NA_real_)
  best    <- names(mods)[which(aic_all <= min(aic_all, na.rm = TRUE) + cfg$delta_aic)[1]]
  final   <- mods[[best]]
  message("Best model for ", resp, ": ", best)
  
  # 8.2 model comparison + performance
  lrt  <- do.call(anova, unname(mods)); rownames(lrt) <- names(mods)
  perf <- bind_rows(Map(perf_table, mods, names(mods), resp))
  cv   <- bind_rows(Map(cv_lolo, f_list, names(f_list), MoreArgs = list(data = dat, resp = resp)))
  
  # 8.3 diagnostics, significance
  diag  <- run_diagnostics(final, dat, resp)
  terms <- term_tests(final, dat, resp)
  
  # 8.4 effect plots (replaces the 6 hand-written plot_model calls)
  in_model <- intersect(c(cfg$predictors, cfg$factor_term),
                        all.vars(lme4::nobars(formula(final))))
  eff <- lapply(in_model,
                function(t) plot(ggeffects::ggpredict(final, terms = t)) + labs(title = NULL))
  ggsave(file.path(cfg$out_dir, paste0("effects_", resp, ".png")),
         patchwork::wrap_plots(eff), width = 12, height = 8, dpi = 200)
  
  list(response = resp, model = final, best = best, basis = sel$table, dfs = sel$dfs,
       lrt = lrt, perf = perf, cv = cv, diag = diag, terms = terms)
})
names(results) <- cfg$responses

# ------------------------------------------- 9. COMBINE & EXPORT
collect <- function(name) bind_rows(lapply(results, `[[`, name))
tabs <- list(
  basis_selection = collect("basis"),
  performance     = collect("perf"),
  cross_validation = collect("cv"),
  diagnostics     = collect("diag"),
  term_tests      = collect("terms") |>
    mutate(p_BH_all = p.adjust(p_value, "BH"),                     # across all 4 models
           sig_BH   = p_BH_all < cfg$alpha,
           sig_bonf = p_value < cfg$alpha / length(cfg$responses)) # stricter threshold
)
for (n in names(tabs)) write.csv(tabs[[n]], file.path(cfg$out_dir, paste0(n, ".csv")), row.names = FALSE)
lapply(tabs, print)

best_models <- lapply(results, `[[`, "model")
best_names  <- sapply(results, `[[`, "best")
write.csv(data.frame(response = names(best_names), best_model = unname(best_names)),
          file.path(cfg$out_dir, "best_model_per_response.csv"), row.names = FALSE)

# one table with the best model of each response
sjPlot::tab_model(best_models, show.aic = TRUE,
                  dv.labels = paste0(cfg$responses, " (", best_names, ")"),
                  file = file.path(cfg$out_dir, "glmm_results.html"))
