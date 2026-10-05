source(file.path(script_dir, "upstream_ml_utils.R"))

fit_penalized_nested <- function(x, y, groups, alpha_grid, model_name) {
  x <- as.matrix(x)
  s <- apply(x, 2, stats::sd)
  keep <- names(s)[is.finite(s) & s > 1e-8]
  stopifnot(length(keep) > 0L)
  xk <- x[, keep, drop = FALSE]
  foldid <- make_inner_foldid(y, groups, "GSE103940")
  stopifnot(all(vapply(split(foldid, groups), function(z) length(unique(z)) == 1L,
                       logical(1))))
  weights <- unname(class_weights(y)[ifelse(y == 1L, "positive", "negative")])
  fits <- vector("list", length(alpha_grid))
  scores <- numeric(length(alpha_grid))
  selected_n <- integer(length(alpha_grid))
  for (i in seq_along(alpha_grid)) {
    fits[[i]] <- suppressWarnings(glmnet::cv.glmnet(
      xk, as.integer(y), family = "binomial", alpha = alpha_grid[i],
      foldid = foldid, type.measure = "deviance", standardize = TRUE,
      weights = weights, nlambda = 40L, lambda.min.ratio = 0.01
    ))
    lambda_index <- which.min(abs(fits[[i]]$lambda - fits[[i]]$lambda.1se))
    scores[i] <- fits[[i]]$cvm[lambda_index]
    cf <- as.matrix(stats::coef(fits[[i]], s = "lambda.1se"))
    selected_n[i] <- sum(cf[-1, 1] != 0)
  }
  best <- which(scores == min(scores, na.rm = TRUE))[1]
  fit <- fits[[best]]
  cf <- as.matrix(stats::coef(fit, s = "lambda.1se"))
  selected <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
  list(
    model = fit, keep = keep, selected = selected,
    tuning = data.frame(model = model_name, alpha = alpha_grid[best],
                        lambda = fit$lambda.1se, lambda_min = fit$lambda.min,
                        inner_deviance = scores[best], selected_n = length(selected)),
    alpha_path = data.frame(alpha = alpha_grid, inner_deviance = scores,
                            selected_n = selected_n)
  )
}

predict_penalized_nested <- function(object, newx) {
  as.numeric(stats::predict(object$model,
                            newx = as.matrix(newx)[, object$keep, drop = FALSE],
                            s = "lambda.1se", type = "response"))
}

fit_svm_rfe_strict_nested <- function(x, y, groups, seed) {
  x <- as.matrix(x)
  foldid <- make_inner_foldid(y, groups, "GSE103940")
  costs <- c(0.01, 0.1, 1, 10)
  p <- ncol(x)
  sizes <- sort(unique(pmin(p, c(2L, 5L, 10L, 20L, p))))
  rows <- list(); counter <- 0L
  for (cost in costs) {
    for (size in sizes) {
      fold_scores <- numeric(max(foldid))
      for (fold in seq_len(max(foldid))) {
        tr <- foldid != fold; va <- foldid == fold
        prep <- fit_preprocessor(x[tr, , drop = FALSE])
        xtr <- apply_preprocessor(x[tr, , drop = FALSE], prep)
        xva <- apply_preprocessor(x[va, , drop = FALSE], prep)
        ranking <- svm_rfe_ranking(xtr, y[tr], cost)
        features <- head(ranking, min(size, length(ranking)))
        set.seed(seed + counter * 100L + fold)
        model <- svm_fit_core(xtr[, features, drop = FALSE], y[tr], cost,
                              probability = FALSE)
        pred <- stats::predict(model, xva[, features, drop = FALSE])
        fold_scores[fold] <- balanced_accuracy_labels(y[va],
                                                       as.integer(pred == "positive"))
      }
      counter <- counter + 1L
      rows[[counter]] <- data.frame(cost = cost, subset_size = size,
                                    inner_balanced_accuracy = mean(fold_scores))
    }
  }
  tuning <- do.call(rbind, rows)
  tuning <- tuning[order(-tuning$inner_balanced_accuracy, tuning$subset_size,
                         tuning$cost), ]
  best <- tuning[1, ]
  prep <- fit_preprocessor(x)
  xs <- apply_preprocessor(x, prep)
  ranking <- svm_rfe_ranking(xs, y, best$cost)
  selected <- head(ranking, best$subset_size)
  set.seed(seed + 10000L)
  model <- svm_fit_core(xs[, selected, drop = FALSE], y, best$cost,
                        probability = TRUE)
  list(model = model, prep = prep, selected = selected,
       tuning = best, ranking = ranking)
}

final_methods <- c(LASSO = "LASSO_probability",
                   ElasticNet = "ElasticNet_probability",
                   RF = "RF_probability", SVM = "SVM_probability",
                   Ensemble = "Ensemble_probability")

run_final_loso <- function(x, y, subjects, seed = 103940L,
                           keep_details = TRUE) {
  x <- as.matrix(x)
  y <- as.integer(y)
  subjects <- as.character(subjects)
  stopifnot(nrow(x) == length(y), length(y) == length(subjects),
            all(y %in% 0:1), !anyNA(x),
            all(tapply(y, subjects, function(z) setequal(z, 0:1))))
  held_subjects <- sort(unique(subjects))
  prediction_rows <- tuning_rows <- selection_rows <- audit_rows <- vector("list", length(held_subjects))
  for (i in seq_along(held_subjects)) {
    held <- held_subjects[i]
    test <- subjects == held
    train <- !test
    stopifnot(sum(test) == 2L, setequal(y[test], 0:1),
              !held %in% subjects[train])
    fold_hash <- sum(utf8ToInt(held) * seq_along(utf8ToInt(held)))
    fold_seed <- seed + as.integer(fold_hash * 100L)
    lasso <- fit_penalized_nested(x[train, , drop = FALSE], y[train],
                                  subjects[train], 1, "LASSO")
    elastic <- fit_penalized_nested(x[train, , drop = FALSE], y[train],
                                    subjects[train], c(0.1, 0.3, 0.5, 0.7, 0.9),
                                    "Elastic Net")
    rf <- fit_rf(x[train, , drop = FALSE], y[train], fold_seed + 20000L)
    svm <- fit_svm_rfe_strict_nested(x[train, , drop = FALSE], y[train],
                                     subjects[train], fold_seed + 40000L)
    prediction_rows[[i]] <- data.frame(
      row_id = rownames(x)[test], subject_id = subjects[test], observed = y[test],
      LASSO_probability = predict_penalized_nested(lasso, x[test, , drop = FALSE]),
      ElasticNet_probability = predict_penalized_nested(elastic, x[test, , drop = FALSE]),
      RF_probability = predict_rf(rf, x[test, , drop = FALSE]),
      SVM_probability = predict_svm(svm, x[test, , drop = FALSE]),
      outer_fold = held, stringsAsFactors = FALSE
    )
    if (keep_details) {
      tuning_rows[[i]] <- rbind(
        data.frame(outer_fold = held, method = "LASSO", parameter = "alpha/lambda.1se",
                   value = paste(lasso$tuning$alpha, signif(lasso$tuning$lambda, 7), sep = "/"),
                   secondary = paste0("selected_n=", length(lasso$selected))),
        data.frame(outer_fold = held, method = "Elastic Net", parameter = "alpha/lambda.1se",
                   value = paste(elastic$tuning$alpha, signif(elastic$tuning$lambda, 7), sep = "/"),
                   secondary = paste0("selected_n=", length(elastic$selected))),
        data.frame(outer_fold = held, method = "RF", parameter = "mtry",
                   value = rf$tuning$mtry,
                   secondary = paste0("OOB_balanced_accuracy=", signif(rf$tuning$OOB_balanced_accuracy, 6))),
        data.frame(outer_fold = held, method = "SVM-RFE", parameter = "cost",
                   value = svm$tuning$cost,
                   secondary = paste0("subset_size=", svm$tuning$subset_size,
                                      ";inner_balanced_accuracy=",
                                      signif(svm$tuning$inner_balanced_accuracy, 6)))
      )
      make_selection <- function(method, genes) {
        data.frame(outer_fold = rep(held, length(genes)),
                   method = rep(method, length(genes)), gene = genes)
      }
      selection_rows[[i]] <- rbind(make_selection("LASSO", lasso$selected),
                                   make_selection("Elastic Net", elastic$selected),
                                   make_selection("SVM-RFE", svm$selected))
      inner_fold <- make_inner_foldid(y[train], subjects[train], "GSE103940")
      audit_rows[[i]] <- data.frame(
        outer_fold = held, test_subjects = held,
        training_subjects = paste(sort(unique(subjects[train])), collapse = ";"),
        train_test_subject_overlap = length(intersect(subjects[train], subjects[test])),
        test_samples = sum(test), training_samples = sum(train),
        test_pair_intact = setequal(y[test], 0:1),
        inner_subject_split = any(vapply(split(inner_fold, subjects[train]),
                                         function(z) length(unique(z)) > 1L, logical(1))),
        preprocessing_scope = "outer-training subjects only",
        feature_selection_scope = "outer-training subjects only",
        tuning_scope = "grouped inner CV within outer-training subjects",
        stringsAsFactors = FALSE
      )
    }
  }
  predictions <- do.call(rbind, prediction_rows)
  predictions$Ensemble_probability <- rowMeans(
    predictions[, c("LASSO_probability", "ElasticNet_probability",
                    "RF_probability", "SVM_probability")])
  predictions <- predictions[match(rownames(x), predictions$row_id), ]
  list(predictions = predictions,
       tuning = if (keep_details) do.call(rbind, tuning_rows) else NULL,
       selections = if (keep_details) do.call(rbind, selection_rows) else NULL,
       leakage_audit = if (keep_details) do.call(rbind, audit_rows) else NULL)
}

final_performance <- function(predictions) {
  do.call(rbind, lapply(names(final_methods), function(method) {
    met <- classification_metrics(predictions$observed,
                                  predictions[[final_methods[[method]]]])
    data.frame(dataset = "GSE103940", method = method, met, row.names = NULL)
  }))
}

paired_scores <- function(predictions) {
  rows <- list()
  k <- 0L
  for (subject in sort(unique(predictions$subject_id))) {
    z <- predictions[predictions$subject_id == subject, ]
    stopifnot(nrow(z) == 2L, setequal(z$observed, 0:1))
    for (method in names(final_methods)) {
      k <- k + 1L
      column <- final_methods[[method]]
      high <- z[[column]][z$observed == 1L]
      plain <- z[[column]][z$observed == 0L]
      rows[[k]] <- data.frame(subject_id = subject, method = method,
                              plain_probability = plain,
                              high_altitude_probability = high,
                              score_difference_high_minus_plain = high - plain,
                              correctly_ranked = high > plain)
    }
  }
  do.call(rbind, rows)
}

subject_bootstrap <- function(predictions, n_boot = 5000L, seed = 20260917L) {
  subjects <- sort(unique(predictions$subject_id))
  set.seed(seed)
  rows <- vector("list", n_boot * length(final_methods))
  k <- 0L
  for (b in seq_len(n_boot)) {
    sampled <- sample(subjects, length(subjects), replace = TRUE)
    idx <- unlist(lapply(sampled, function(s) which(predictions$subject_id == s)),
                  use.names = FALSE)
    for (method in names(final_methods)) {
      k <- k + 1L
      met <- classification_metrics(predictions$observed[idx],
                                    predictions[[final_methods[[method]]]][idx])
      rows[[k]] <- data.frame(iteration = b, method = method, met[, 1:4])
    }
  }
  do.call(rbind, rows)
}
