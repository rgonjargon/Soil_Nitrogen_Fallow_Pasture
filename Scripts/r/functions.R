# =============================================================================
# Fallow vs pasture analysis — shared functions (sourced by 1_analysis.qmd)
# =============================================================================
# Exports: clean_formula, fit_models, best_model_equations,
#          check_model_performance, best_model_summaries, plot_model_fits

# -----------------------------------------------------------------------------
# clean_formula: format formula for display (RHS only, × for *, ² for poly)
# -----------------------------------------------------------------------------
clean_formula <- function(s) {
  s <- sub("^[^~]*~\\s*", "", s)
  s <- gsub(" \\* ", " × ", s)
  s <- gsub("poly\\(([^,]+), 2, raw = FALSE\\)", "\\1\u00b2", s)
  s <- gsub("poly\\(([^,]+), 2\\)", "\\1\u00b2", s)
  s <- gsub("I\\(sqrt\\(([^)]+)\\)\\)", "sqrt(\\1)", s)
  s
}

# -----------------------------------------------------------------------------
# fit_models: fit Linear, Quadratic, Square Root; return fits + AICc table
# -----------------------------------------------------------------------------
fit_models <- function(data, response_vars, time, treatment = NULL, depth = NULL,
                       random_effects = NULL) {

  fit_fn <- function(f, d = data) {
    if (!is.null(random_effects)) lmer(f, data = d, REML = FALSE) else lm(f, data = d)
  }

  safe_fit <- function(expr, label, resp) {
    tryCatch(expr, error = function(e) {
      warning(label, " model failed for '", resp, "': ", conditionMessage(e))
      NULL
    })
  }

  make_formula <- function(resp, term) {
    rhs <- if (!is.null(random_effects)) paste(term, "+", random_effects) else term
    as.formula(paste(resp, "~", rhs))
  }

  aicc_fn <- function(m) {
    if (requireNamespace("MuMIn", quietly = TRUE)) return(MuMIn::AICc(m))
    k <- length(coef(m)) + 1
    n <- nobs(m)
    AIC(m) + (2 * k * (k + 1)) / (n - k - 1)
  }

  r2_fn <- function(m) {
    if (requireNamespace("MuMIn", quietly = TRUE))
      return(round(MuMIn::r.squaredGLMM(m)[1, "R2m"], 2))
    round(summary(m)$r.squared, 2)
  }

  fits      <- list()
  aicc_rows <- list()

  for (resp in response_vars) {

    # Decide whether to include depth (need ≥2 levels in non-NA rows)
    use_depth <- FALSE
    if (!is.null(depth) && depth %in% names(data)) {
      valid_rows   <- !is.na(data[[resp]])
      n_depth_lvls <- length(unique(na.omit(data[[depth]][valid_rows])))
      if (n_depth_lvls < 2) {
        warning("'", resp, "': depth has < 2 levels — dropped; fitting time * treatment only.")
      } else {
        use_depth <- TRUE
      }
    }

    # Model terms by case
    if (!is.null(treatment) && use_depth) {
      lin_term  <- paste(time, "*", treatment, "*", depth)
      quad_term <- paste0("poly(", time, ", 2, raw = FALSE) * ", treatment, " * ", depth)
      sqrt_term <- paste0("I(sqrt(", time, ")) * ", treatment, " * ", depth)
    } else if (!is.null(treatment)) {
      lin_term  <- paste(time, "*", treatment)
      quad_term <- paste0("poly(", time, ", 2, raw = FALSE) * ", treatment)
      sqrt_term <- paste0("I(sqrt(", time, ")) * ", treatment)
    } else if (use_depth) {
      lin_term  <- paste(time, "*", depth)
      quad_term <- paste0("poly(", time, ", 2, raw = FALSE) * ", depth)
      sqrt_term <- paste0("I(sqrt(", time, ")) * ", depth)
    } else {
      lin_term  <- time
      quad_term <- paste0("poly(", time, ", 2, raw = FALSE)")
      sqrt_term <- paste0("I(sqrt(", time, "))")
    }

    f1 <- make_formula(resp, lin_term)
    f2 <- make_formula(resp, quad_term)
    f3 <- make_formula(resp, sqrt_term)

    # Complete-case data for this response (avoids NA/singularity when depth dropped)
    vars_needed <- c(resp, time, if (!is.null(treatment)) treatment, if (use_depth && !is.null(depth)) depth)
    vars_needed <- vars_needed[vars_needed %in% names(data)]
    data_fit    <- data[complete.cases(data[, vars_needed, drop = FALSE]), , drop = FALSE]
    if (nrow(data_fit) < 4L) {
      warning("'", resp, "': too few complete cases (n = ", nrow(data_fit), ") — skipping.")
      aicc_rows[[resp]] <- data.frame(
        Response = c(resp, rep("", 2)), Model = c("Linear", "Quadratic", "Square Root"),
        K = NA, n = NA, R2 = NA, AICc = NA, Delta_i = NA, omega_i = NA, Evid_ratio = NA,
        Formula = clean_formula(c(deparse(f1), deparse(f2), deparse(f3))),
        stringsAsFactors = FALSE
      )
      next
    }

    m1 <- safe_fit(fit_fn(f1, data_fit), "Linear", resp)
    m2 <- safe_fit(fit_fn(f2, data_fit), "Quadratic", resp)
    m3 <- safe_fit(fit_fn(f3, data_fit), "Square Root", resp)

    fits[[resp]] <- list(linear = m1, quadratic = m2, `square root` = m3)
    models       <- list(linear = m1, quadratic = m2, `square root` = m3)
    model_names  <- c("Linear", "Quadratic", "Square Root")
    formulas     <- clean_formula(c(deparse(f1), deparse(f2), deparse(f3)))

    if (all(sapply(models, is.null))) {
      warning("All models failed for '", resp, "' — returning NA rows.")
      aicc_rows[[resp]] <- data.frame(
        Response = c(resp, rep("", 2)), Model = model_names,
        K = NA, n = NA, R2 = NA, AICc = NA, Delta_i = NA, omega_i = NA, Evid_ratio = NA,
        Formula = formulas, stringsAsFactors = FALSE
      )
      next
    }

    aicc_vals <- sapply(models, function(m) if (is.null(m)) NA else aicc_fn(m))
    n_obs     <- nobs(Filter(Negate(is.null), models)[[1]])

    k_vals <- sapply(models, function(m) {
      if (is.null(m)) return(NA)
      if (inherits(m, "lmerMod")) length(fixef(m)) + 1 else length(coef(m)) + 1
    })
    r2_vals <- sapply(models, function(m) {
      if (is.null(m)) return(NA)
      tryCatch(if (inherits(m, "lmerMod")) r2_fn(m) else round(summary(m)$r.squared, 2), error = function(e) NA)
    })

    # AICc weights and evidence ratios
    ord   <- order(aicc_vals, na.last = TRUE)
    delta <- aicc_vals[ord] - min(aicc_vals, na.rm = TRUE)
    w     <- exp(-0.5 * delta) / sum(exp(-0.5 * delta), na.rm = TRUE)
    evid  <- w[1] / w

    n_models <- length(models)
    aicc_rows[[resp]] <- data.frame(
      Response   = c(resp, rep("", n_models - 1)),
      Model      = model_names[ord],
      K          = k_vals[ord],
      n          = n_obs,
      R2         = r2_vals[ord],
      AICc       = round(aicc_vals[ord], 1),
      Delta_i    = round(delta, 1),
      omega_i    = round(w, 2),
      Evid_ratio = round(evid, 2),
      Formula    = formulas[ord],
      stringsAsFactors = FALSE
    )
  }

  aicc_table <- do.call(rbind, aicc_rows)
  rownames(aicc_table) <- NULL
  colnames(aicc_table) <- c("Response/Hypothesis", "Model", "K", "n",
                             "R2", "AICc", "\u0394i", "\u03c9i", "Evid. ratio", "Model Formula")

  list(fits = fits, aicc_table = aicc_table)
}


# -----------------------------------------------------------------------------
# best_model_equations: coefficient table for best model only (one col per coef)
# -----------------------------------------------------------------------------
best_model_equations <- function(model_fits, digits = 3) {
  coef_list <- lapply(names(model_fits$fits), function(resp) {
    models   <- model_fits$fits[[resp]]
    aicc_sub <- model_fits$aicc_table[model_fits$aicc_table[[1]] == resp, ]
    best_nm  <- tolower(aicc_sub$Model[1])
    m        <- models[[best_nm]]
    if (is.null(m)) return(NULL)

    coefs <- if (inherits(m, "lmerMod")) fixef(m) else coef(m)

    row <- data.frame(
      Response   = resp,
      Best_Model = aicc_sub$Model[1],
      Formula    = aicc_sub[["Model Formula"]][1],
      R2         = aicc_sub$R2[1],

      stringsAsFactors = FALSE
    )
    # Add each coefficient as its own column
    for (nm in names(coefs)) {
      row[[nm]] <- round(coefs[[nm]], digits)
    }
    row
  })

  coef_list <- Filter(Negate(is.null), coef_list)

  # Bind rows — missing coefficients (different model structures) become NA
  all_cols <- unique(unlist(lapply(coef_list, names)))
  out <- do.call(rbind, lapply(coef_list, function(r) {
    missing <- setdiff(all_cols, names(r))
    for (col in missing) r[[col]] <- NA
    r[, all_cols]
  }))

  rownames(out) <- NULL
  out
}


# -----------------------------------------------------------------------------
# check_model_performance: diagnostic plots for best model only, per response
# -----------------------------------------------------------------------------
check_model_performance <- function(model_fits) {
  for (resp in names(model_fits$fits)) {
    models   <- model_fits$fits[[resp]]
    aicc_sub <- model_fits$aicc_table[model_fits$aicc_table[[1]] == resp, ]
    best_nm  <- tolower(aicc_sub$Model[1])
    m        <- models[[best_nm]]
    if (is.null(m)) next

    cat("\n", strrep("=", 60), "\n")
    cat(" Response:", resp, "| Best model:", aicc_sub$Model[1], "\n")
    cat(strrep("=", 60), "\n")

    tryCatch({
      # panel = TRUE gives one combined figure (posterior predictive, linearity, homogeneity, collinearity, normality)
      cm <- performance::check_model(m, panel = TRUE)
      p  <- plot(cm)
      if (inherits(p, "list")) {
        gridExtra::grid.arrange(grobs = p, ncol = 1)
      } else {
        print(p)
      }
      invisible(NULL)  # prevent knitr from auto-printing the same figure again
    }, error = function(e) {
      if (grepl("contrasts can be applied only to factors", conditionMessage(e)))
        message("  Skipped — factor has < 2 levels")
      else stop(e)
    })
  }
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# best_model_summaries: print summary() of best model only, per response
# -----------------------------------------------------------------------------
best_model_summaries <- function(model_fits) {
  for (resp in names(model_fits$fits)) {
    models   <- model_fits$fits[[resp]]
    aicc_sub <- model_fits$aicc_table[model_fits$aicc_table[[1]] == resp, ]
    best_nm  <- tolower(aicc_sub$Model[1])
    m        <- models[[best_nm]]
    if (is.null(m)) next

    cat("\n")
    cat(strrep("=", 70), "\n")
    cat("  Response variable:", resp, "\n")
    cat("  Best model:       ", aicc_sub$Model[1], "\n")
    cat(strrep("=", 70), "\n\n")

    print(summary(m))
    cat("\n")
  }
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# plot_model_fits: one panel per response; best model fit + raw data; Fallow=br, Pasture=gr
# -----------------------------------------------------------------------------
plot_model_fits <- function(model_fits, time, treatment = NULL, depth = NULL,
                            raw_data, response_labels = NULL,
                            depth_order = NULL, ncol = 2, force_model = NULL,
                            save_dir = NULL, save_prefix = "", save_dpi = 600,
                            save_width = 5, save_height = 5) {
  library(ggplot2)
  library(RColorBrewer)
  plots <- list()

  for (resp in names(model_fits$fits)) {

    models     <- model_fits$fits[[resp]]
    aicc_sub   <- model_fits$aicc_table[model_fits$aicc_table[[1]] == resp, ]
    best_model <- aicc_sub$Model[1]
    best_nm    <- if (!is.null(force_model)) tolower(force_model) else tolower(best_model)
    m_best     <- models[[best_nm]]
    if (is.null(m_best)) next

    time_seq     <- seq(min(raw_data[[time]], na.rm = TRUE),
                        max(raw_data[[time]], na.rm = TRUE),
                        length.out = 200)
    treat_levels <- if (!is.null(treatment)) levels(raw_data[[treatment]]) else "all"

    has_depth <- !is.null(depth) && depth %in% names(raw_data)

    # Depth ordering — user-supplied or sorted
    depth_levels <- if (has_depth) {
      present <- as.character(unique(na.omit(raw_data[[depth]][!is.na(raw_data[[resp]])])))
      if (!is.null(depth_order)) intersect(depth_order, present) else sort(present)
    } else "all"

    multi_depth <- has_depth && length(depth_levels) > 1

    # Linetype: all solid (depth shown by colour only)
    depth_linetypes <- if (multi_depth) setNames(rep("solid", length(depth_levels)), depth_levels) else NULL

    # Colour: Fallow = brown, Pasture = green; darker = shallower depth
    n_trt       <- length(treat_levels)
    trt_colors  <- c("Fallow" = "#8B4513", "Pasture" = "#228B22")  # saddle brown, forest green
    base_hues   <- setNames(character(length(treat_levels)), treat_levels)
    for (trt in treat_levels) {
      base_hues[trt] <- if (trt %in% names(trt_colors)) trt_colors[trt]
                       else brewer.pal(max(3, n_trt), "Set1")[match(trt, treat_levels)]
    }

    if (multi_depth) {
      combo_levels <- as.vector(outer(treat_levels, depth_levels, paste, sep = " | "))
      pal <- setNames(vector("character", length(combo_levels)), combo_levels)
      for (trt in treat_levels) {
        # Ramp dark (shallow) → light (deep): wider spread so depths are more distinct
        dark_end <- colorspace::darken(base_hues[trt], 0.15)
        light_end <- colorspace::lighten(base_hues[trt], 0.7)
        ramp   <- colorRamp(c(dark_end, light_end))
        shades <- apply(ramp(seq(0, 1, length.out = length(depth_levels))), 1,
                        function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
        for (di in seq_along(depth_levels))
          pal[paste(trt, depth_levels[di], sep = " | ")] <- shades[di]
      }
    } else {
      pal <- base_hues
    }

    # Predictions on a fine time grid for each treatment × depth
    combos <- expand.grid(trt = treat_levels, dep = depth_levels, stringsAsFactors = FALSE)

    pred_list <- lapply(seq_len(nrow(combos)), function(i) {
      trt <- combos$trt[i]; dep <- combos$dep[i]
      nd  <- data.frame(time_seq); names(nd) <- time
      if (!is.null(treatment)) nd[[treatment]] <- factor(trt, levels = treat_levels)
      if (has_depth && dep != "all") nd[[depth]] <- factor(dep, levels = depth_levels)

      pred <- tryCatch({
        p <- predict(m_best, newdata = nd, re.form = NA, se.fit = TRUE)
        data.frame(fit = p$fit, se = p$se.fit)
      }, error = function(e) {
        fit <- tryCatch(predict(m_best, newdata = nd), error = function(e2) rep(NA, nrow(nd)))
        data.frame(fit = fit, se = NA)
      })

      nd[[resp]]      <- pred$fit
      nd$ymin         <- pred$fit - 1.96 * pred$se
      nd$ymax         <- pred$fit + 1.96 * pred$se
      nd$trt_group    <- if (!is.null(treatment)) trt else "Predicted"
      nd$depth_group  <- if (has_depth && dep != "all") dep else NA
      nd$colour_group <- if (multi_depth) paste(trt, dep, sep = " | ") else nd$trt_group
      nd
    })

    pred_df <- do.call(rbind, pred_list)
    pred_df$colour_group <- factor(pred_df$colour_group, levels = names(pal))
    if (multi_depth) pred_df$depth_group <- factor(pred_df$depth_group, levels = depth_levels)

    # Raw data grouping
    plot_data <- raw_data[!is.na(raw_data[[resp]]), ]
    plot_data$colour_group <- if (multi_depth) {
      factor(paste(plot_data[[treatment]], as.character(plot_data[[depth]]), sep = " | "), levels = names(pal))
    } else if (!is.null(treatment)) {
      factor(plot_data[[treatment]], levels = treat_levels)
    } else "Predicted"

    if (multi_depth) plot_data$depth_group <- factor(plot_data[[depth]], levels = depth_levels)

    y_label   <- if (!is.null(response_labels) && resp %in% names(response_labels)) response_labels[[resp]] else gsub("_", " ", resp)
    y_label   <- gsub(",", "", y_label)  # remove commas from variable names
    col_label <- if (multi_depth) paste(gsub("_"," ",treatment), "\u00d7", gsub("_"," ",depth)) else if (!is.null(treatment)) gsub("_"," ",treatment) else NULL
    if (!is.null(col_label)) col_label <- tools::toTitleCase(col_label)

    # Month axis: breaks every 12 months
    x_range      <- range(raw_data[[time]], na.rm = TRUE)
    month_breaks <- seq(floor(x_range[1] / 12) * 12, ceiling(x_range[2] / 12) * 12, by = 12)

    # Formula form for subtitle and equation LHS: y, sqrt(y), or y²
    model_form <- switch(best_nm,
      "linear"       = "y",
      "square root"  = "sqrt(y)",
      "quadratic"    = "y\u00b2",
      "y")
    eq_lhs    <- switch(best_nm,
      "linear"       = "y = ",
      "square root"  = "sqrt(y) = ",
      "quadratic"    = "y\u00b2 = ",
      "y = ")
    model_name <- tools::toTitleCase(best_nm)

    # Subtitle: equations in order Fallow (depths low→high), then Pasture (depths low→high)
    trt_order <- c(intersect(c("Fallow", "Pasture"), treat_levels), setdiff(treat_levels, c("Fallow", "Pasture")))
    eq_order  <- if (multi_depth) as.vector(t(outer(trt_order, depth_levels, paste, sep = " | "))) else treat_levels
    eq_lines  <- character(0)
    for (grp in eq_order) {
      if (!grp %in% names(pal)) next
      pg <- pred_df[pred_df$colour_group == grp, ]
      if (nrow(pg) >= 2) {
        pg <- pg[order(pg[[time]]), ]
        x1 <- pg[[time]][1]; x2 <- pg[[time]][nrow(pg)]
        y1 <- pg[[resp]][1]; y2 <- pg[[resp]][nrow(pg)]
        slope     <- (y2 - y1) / (x2 - x1)
        intercept <- y1 - slope * x1
        slope_str     <- format(round(slope, 3), trim = TRUE)
        intercept_str <- format(round(intercept, 3), trim = TRUE)
        pm <- if (slope >= 0) " + " else " \u2212 "
        slope_str <- if (slope >= 0) slope_str else format(round(-slope, 3), trim = TRUE)
        eq_lines <- c(eq_lines, paste0(grp, ": ", eq_lhs, intercept_str, pm, slope_str, "\u00b7", sub("months", "Months", time, ignore.case = TRUE)))
      }
    }
    # Subtitle: Best model name (Linear / Square Root / Quadratic) and formula form, then R², then formulas
    r2_val    <- aicc_sub$R2[1]
    r2_line   <- paste0("R\u00b2 = ", r2_val)
    best_line <- paste0(if (!is.null(force_model)) "Forced model: " else "Best model: ", model_name)
    subtitle_text <- paste(c(best_line, r2_line, eq_lines), collapse = "\n")

    p <- ggplot() +
      # Ribbon
      { if (!all(is.na(pred_df$ymin)))
          geom_ribbon(data = pred_df,
                      aes(x = .data[[time]], ymin = ymin, ymax = ymax,
                          fill  = colour_group,
                          group = colour_group),
                      alpha = 0.12, show.legend = FALSE)
      } +
      # Raw points
      geom_point(data = plot_data,
                 aes(x      = .data[[time]],
                     y      = .data[[resp]],
                     colour = colour_group,
                     shape  = if (multi_depth) depth_group else NULL),
                 alpha = 0.45, size = 1.8) +
      # Predicted lines — no alpha (fully opaque); alpha only on points
      { if (multi_depth)
          geom_line(data = pred_df,
                    aes(x        = .data[[time]],
                        y        = .data[[resp]],
                        colour   = colour_group,
                        linetype = depth_group),
                    linewidth = 1, alpha = 1)
        else
          geom_line(data = pred_df,
                    aes(x      = .data[[time]],
                        y      = .data[[resp]],
                        colour = colour_group),
                    linewidth = 1, alpha = 1)
      } +
      # Legend: multi_depth = 6 combos with key linetypes; no depth = Fallow/Pasture stacked
      { if (multi_depth) {
          combo_ltys <- setNames(
            depth_linetypes[sub(".* \\| ", "", names(pal))],
            names(pal)
          )
          scale_colour_manual(
            values = pal, name = col_label, drop = FALSE,
            guide  = guide_legend(
              keywidth     = unit(1, "cm"),
              override.aes = list(linetype = unname(combo_ltys), linewidth = 1)
            )
          )
        } else {
          # No depth: order Fallow then Pasture, stack legend vertically (one on top of the other)
          trt_order_leg <- c(intersect(c("Fallow", "Pasture"), names(pal)), setdiff(names(pal), c("Fallow", "Pasture")))
          pal_ordered   <- pal[trt_order_leg]
          scale_colour_manual(
            values = pal_ordered, name = col_label, drop = FALSE,
            guide  = guide_legend(nrow = 2)
          )
        }
      } +
      scale_fill_manual(values = pal, drop = FALSE, guide = "none") +
      scale_x_continuous(breaks = month_breaks, expand = expansion(mult = c(0.02, 0.06))) +
      scale_y_continuous(expand = expansion(mult = 0.12)) +
      { if (multi_depth) scale_linetype_manual(values = depth_linetypes, guide = "none") } +
      { if (multi_depth) scale_shape_discrete(guide = "none") } +
      labs(
        title    = "",
        subtitle = subtitle_text,
        x = sub("months", "Months", gsub("_", " ", time), ignore.case = TRUE),
        y = y_label
      ) +
      theme_classic(base_size = 11) +
      theme(plot.title            = element_text(face = "bold", size = 11),
            plot.subtitle         = element_text(size = 8, colour = "grey40", hjust = 0, lineheight = 1,
                                                margin = margin(b = 10)),
            legend.position       = "bottom",
            legend.box.spacing    = unit(10, "pt"),
            legend.title.position = "top",
            legend.key.size       = unit(0.5, "cm"),
            legend.text           = element_text(size = 7),
            legend.title          = element_text(size = 8))
    # Aspect: square if one panel; 0.5 if multi-panel so panel heights match in grid
    n_plots <- length(model_fits$fits)
    p <- p + theme(aspect.ratio = if (n_plots == 1L) 1 else 0.5)
    if (!multi_depth) p <- p + theme(legend.justification = "left")
    plots[[resp]] <- p
  }

  # Equalise panel height: pad subtitle to same number of lines
  n_subtitle_lines <- vapply(plots, function(p) {
    s <- p$labels$subtitle
    if (is.null(s) || !is.character(s)) return(0L)
    length(strsplit(s, "\n", fixed = TRUE)[[1L]])
  }, integer(1L))
  max_subtitle_lines <- max(n_subtitle_lines)
  for (i in seq_along(plots)) {
    if (n_subtitle_lines[i] < max_subtitle_lines) {
      s <- plots[[i]]$labels$subtitle
      plots[[i]] <- plots[[i]] + labs(subtitle = paste0(s, strrep("\n", max_subtitle_lines - n_subtitle_lines[i])))
    }
  }

  # Save each response plot to disk when save_dir is set (same dimensions as Quarto panels, 600 dpi)
  if (!is.null(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    for (resp in names(plots)) {
      fname <- paste0(gsub("[^a-zA-Z0-9._-]", "_", resp), ".png")
      if (nzchar(save_prefix)) fname <- paste0(save_prefix, "_", fname)
      ggplot2::ggsave(
        filename = file.path(save_dir, fname),
        plot     = plots[[resp]],
        width    = save_width,
        height   = save_height,
        dpi      = save_dpi
      )
    }
  }

  # Each plot keeps its own legend (no shared legend)
  if (length(plots) == 0L) return(invisible(plots))
  nrow_ <- ceiling(length(plots) / ncol)
  gridExtra::grid.arrange(grobs = plots, ncol = ncol, nrow = nrow_, widths = rep(1, ncol))

  invisible(plots)
}
