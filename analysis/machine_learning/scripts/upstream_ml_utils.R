safe_auc <- function(y, prob) {
  y <- as.integer(y)
  if (length(unique(y)) < 2L || any(!is.finite(prob))) return(NA_real_)
  pos <- prob[y == 1L]
  neg <- prob[y == 0L]
  mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
}

classification_metrics <- function(y, prob, threshold = 0.5) {
  y <- as.integer(y)
  pred <- as.integer(prob >= threshold)
  tp <- sum(pred == 1L & y == 1L)
  tn <- sum(pred == 0L & y == 0L)
  fp <- sum(pred == 1L & y == 0L)
  fn <- sum(pred == 0L & y == 1L)
  sensitivity <- if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_
  specificity <- if ((tn + fp) > 0L) tn / (tn + fp) else NA_real_
  data.frame(
    ROC_AUC = safe_auc(y, prob),
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    sensitivity = sensitivity,
    specificity = specificity,
    TP = tp, TN = tn, FP = fp, FN = fn
  )
}

fit_preprocessor <- function(x) {
  x <- as.matrix(x)
  s <- apply(x, 2, sd)
  keep <- is.finite(s) & s > 1e-8
  stopifnot(any(keep))
  x <- x[, keep, drop = FALSE]
  center <- colMeans(x)
  scale <- apply(x, 2, sd)
  list(keep = names(keep)[keep], center = center, scale = scale)
}

apply_preprocessor <- function(x, prep) {
  x <- as.matrix(x)[, prep$keep, drop = FALSE]
  x <- sweep(x, 2, prep$center, "-")
  sweep(x, 2, prep$scale, "/")
}

make_inner_foldid <- function(y, groups, task) {
  y <- as.integer(y)
  if (task == "GSE103940") {
    ug <- sort(unique(groups))
    k <- min(5L, length(ug))
    group_fold <- setNames(rep(seq_len(k), length.out = length(ug)), ug)
    return(unname(group_fold[as.character(groups)]))
  }
  k <- min(3L, min(table(y)))
  stopifnot(k >= 2L)
  foldid <- integer(length(y))
  for (class_value in 0:1) {
    idx <- which(y == class_value)
    idx <- idx[order(as.character(groups[idx]))]
    foldid[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  foldid
}

class_weights <- function(y) {
  tab <- table(factor(y, levels = 0:1))
  w <- sum(tab) / (2 * tab)
  setNames(as.numeric(w), c("negative", "positive"))
}

factor_outcome <- function(y) {
  factor(ifelse(as.integer(y) == 1L, "positive", "negative"),
         levels = c("negative", "positive"))
}

fit_lasso <- function(x, y, groups, task) {
  prep <- fit_preprocessor(x)
  xs <- apply_preprocessor(x, prep)
  foldid <- make_inner_foldid(y, groups, task)
  fit <- suppressWarnings(glmnet::cv.glmnet(
    xs, as.integer(y), family = "binomial", alpha = 1,
    foldid = foldid, type.measure = "deviance", standardize = FALSE,
    weights = unname(class_weights(y)[ifelse(y == 1L, "positive", "negative")]),
    nlambda = 40, lambda.min.ratio = 0.01
  ))
  cf <- as.matrix(stats::coef(fit, s = "lambda.1se"))
  selected <- setdiff(rownames(cf)[cf[, 1] != 0], "(Intercept)")
  list(model = fit, prep = prep, selected = selected,
       tuning = data.frame(lambda = fit$lambda.1se,
                           lambda_min = fit$lambda.min,
                           selected_n = length(selected)))
}

predict_lasso <- function(object, newx) {
  xs <- apply_preprocessor(newx, object$prep)
  as.numeric(stats::predict(object$model, newx = xs, s = "lambda.1se",
                            type = "response"))
}

rf_importance_vector <- function(model) {
  imp <- randomForest::importance(model, type = 1, scale = FALSE)
  if (is.matrix(imp)) imp <- imp[, 1]
  imp <- as.numeric(imp)
  names(imp) <- rownames(randomForest::importance(model, type = 1, scale = FALSE))
  imp
}

fit_rf <- function(x, y, seed, ntree = 500L, fixed_mtry = NULL) {
  prep <- fit_preprocessor(x)
  xs <- apply_preprocessor(x, prep)
  p <- ncol(xs)
  candidates <- if (is.null(fixed_mtry)) {
    sort(unique(pmax(1L, pmin(p, c(1L, floor(sqrt(p)), floor(p / 3), floor(p / 2))))))
  } else {
    as.integer(max(1L, min(p, fixed_mtry)))
  }
  fits <- vector("list", length(candidates))
  scores <- numeric(length(candidates))
  for (i in seq_along(candidates)) {
    set.seed(seed + i)
    fits[[i]] <- randomForest::randomForest(
      x = xs, y = factor_outcome(y), ntree = ntree, mtry = candidates[i],
      importance = TRUE, classwt = class_weights(y),
      sampsize = rep(min(table(factor_outcome(y))), 2L), keep.forest = TRUE
    )
    class_error <- fits[[i]]$confusion[, "class.error"]
    scores[i] <- 1 - mean(class_error, na.rm = TRUE)
  }
  best <- which(scores == max(scores, na.rm = TRUE))[1]
  model <- fits[[best]]
  list(model = model, prep = prep, selected = character(0),
       importance = rf_importance_vector(model),
       tuning = data.frame(mtry = candidates[best], OOB_balanced_accuracy = scores[best],
                           ntree = ntree))
}

predict_rf <- function(object, newx) {
  xs <- apply_preprocessor(newx, object$prep)
  as.numeric(stats::predict(object$model, newdata = xs, type = "prob")[, "positive"])
}

svm_fit_core <- function(x, y, cost, probability = FALSE) {
  e1071::svm(x = as.matrix(x), y = factor_outcome(y), kernel = "linear",
             cost = cost, scale = FALSE, probability = probability,
             class.weights = class_weights(y))
}

svm_weight <- function(model) {
  w <- drop(t(model$coefs) %*% model$SV)
  names(w) <- colnames(model$SV)
  w
}

svm_rfe_ranking <- function(x, y, cost) {
  remaining <- colnames(x)
  eliminated <- character(0)
  while (length(remaining) > 2L) {
    model <- svm_fit_core(x[, remaining, drop = FALSE], y, cost, probability = FALSE)
    importance <- abs(svm_weight(model))
    remove_n <- max(1L, floor(length(remaining) * 0.20))
    remove_n <- min(remove_n, length(remaining) - 2L)
    remove <- names(sort(importance, decreasing = FALSE))[seq_len(remove_n)]
    eliminated <- c(eliminated, remove)
    remaining <- setdiff(remaining, remove)
  }
  model <- svm_fit_core(x[, remaining, drop = FALSE], y, cost, probability = FALSE)
  final_order <- names(sort(abs(svm_weight(model)), decreasing = FALSE))
  rev(c(eliminated, final_order))
}

balanced_accuracy_labels <- function(y, pred) {
  sens <- mean(pred[y == 1L] == 1L)
  spec <- mean(pred[y == 0L] == 0L)
  mean(c(sens, spec))
}

fit_svm_rfe <- function(x, y, groups, task, seed) {
  prep <- fit_preprocessor(x)
  xs <- apply_preprocessor(x, prep)
  p <- ncol(xs)
  costs <- c(0.01, 0.1, 1, 10)
  sizes <- sort(unique(pmin(p, c(2L, 5L, 10L, 20L, p))))
  foldid <- make_inner_foldid(y, groups, task)
  rows <- list()
  rankings <- list()
  counter <- 0L
  for (cost in costs) {
    ranking <- svm_rfe_ranking(xs, y, cost)
    rankings[[as.character(cost)]] <- ranking
    for (size in sizes) {
      features <- head(ranking, size)
      fold_scores <- numeric(max(foldid))
      for (fold in seq_len(max(foldid))) {
        tr <- foldid != fold
        va <- foldid == fold
        inner_prep <- fit_preprocessor(x[tr, features, drop = FALSE])
        xtr <- apply_preprocessor(x[tr, features, drop = FALSE], inner_prep)
        xva <- apply_preprocessor(x[va, features, drop = FALSE], inner_prep)
        set.seed(seed + counter + fold)
        model <- svm_fit_core(xtr, y[tr], cost, probability = FALSE)
        pred <- stats::predict(model, xva)
        pred <- as.integer(pred == "positive")
        fold_scores[fold] <- balanced_accuracy_labels(y[va], pred)
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
  selected <- head(rankings[[as.character(best$cost)]], best$subset_size)
  set.seed(seed + 10000L)
  model <- svm_fit_core(xs[, selected, drop = FALSE], y, best$cost, probability = TRUE)
  list(model = model, prep = prep, selected = selected,
       tuning = best, ranking = rankings[[as.character(best$cost)]])
}

predict_svm <- function(object, newx) {
  xs <- apply_preprocessor(newx, object$prep)
  pred <- stats::predict(object$model, xs[, object$selected, drop = FALSE],
                         probability = TRUE)
  probabilities <- attr(pred, "probabilities")
  if (is.null(probabilities) || !"positive" %in% colnames(probabilities)) {
    stop("SVM probability output is unavailable")
  }
  as.numeric(probabilities[, "positive"])
}

run_loso <- function(x, y, subjects, task, seed, keep_details = TRUE) {
  x <- as.matrix(x)
  y <- as.integer(y)
  stopifnot(nrow(x) == length(y), length(y) == length(subjects),
            all(y %in% 0:1), !anyNA(x))
  if (task == "GSE103940") {
    fold_subjects <- lapply(sort(unique(as.character(subjects))), identity)
    fold_labels <- vapply(fold_subjects, as.character, character(1))
  } else {
    negative <- sort(unique(as.character(subjects)[y == 0L]))
    positive <- sort(unique(as.character(subjects)[y == 1L]))
    pair_grid <- expand.grid(negative = negative, positive = positive,
                             stringsAsFactors = FALSE)
    fold_subjects <- lapply(seq_len(nrow(pair_grid)), function(i) {
      c(pair_grid$negative[i], pair_grid$positive[i])
    })
    fold_labels <- paste(pair_grid$negative, pair_grid$positive, sep = "__")
  }
  prediction_rows <- list()
  tuning_rows <- list()
  selection_rows <- list()
  for (i in seq_along(fold_subjects)) {
    held_out <- fold_subjects[[i]]
    fold_label <- paste(sort(held_out), collapse = "__")
    test <- as.character(subjects) %in% held_out
    train <- !test
    stopifnot(!any(as.character(subjects[train]) %in% held_out),
              length(unique(y[train])) == 2L)
    fold_hash <- sum(utf8ToInt(fold_label) * seq_along(utf8ToInt(fold_label)))
    fold_seed <- seed + as.integer(fold_hash * 100L)

    lasso <- fit_lasso(x[train, , drop = FALSE], y[train], subjects[train], task)
    rf <- fit_rf(x[train, , drop = FALSE], y[train], fold_seed + 20000L)
    svm <- fit_svm_rfe(x[train, , drop = FALSE], y[train], subjects[train], task,
                       fold_seed + 40000L)

    prediction_rows[[i]] <- data.frame(
      row_id = rownames(x)[test], subject_id = as.character(subjects[test]),
      observed = y[test],
      LASSO_probability = predict_lasso(lasso, x[test, , drop = FALSE]),
      RF_probability = predict_rf(rf, x[test, , drop = FALSE]),
      SVM_probability = predict_svm(svm, x[test, , drop = FALSE]),
      outer_fold = fold_label, stringsAsFactors = FALSE
    )
    if (keep_details) {
      tuning_rows[[i]] <- rbind(
        data.frame(outer_fold = fold_label, method = "LASSO",
                   parameter = "lambda.1se", value = lasso$tuning$lambda,
                   secondary = paste0("selected_n=", lasso$tuning$selected_n)),
        data.frame(outer_fold = fold_label, method = "RF", parameter = "mtry",
                   value = rf$tuning$mtry,
                   secondary = paste0("OOB_balanced_accuracy=",
                                      signif(rf$tuning$OOB_balanced_accuracy, 6))),
        data.frame(outer_fold = fold_label, method = "SVM-RFE",
                   parameter = "cost", value = svm$tuning$cost,
                   secondary = paste0("subset_size=", svm$tuning$subset_size,
                                      ";inner_balanced_accuracy=",
                                      signif(svm$tuning$inner_balanced_accuracy, 6)))
      )
      selection_rows[[i]] <- rbind(
        data.frame(outer_fold = rep(fold_label, length(lasso$selected)),
                   method = rep("LASSO", length(lasso$selected)),
                   gene = lasso$selected),
        data.frame(outer_fold = rep(fold_label, length(svm$selected)),
                   method = rep("SVM-RFE", length(svm$selected)),
                   gene = svm$selected)
      )
    }
  }
  fold_predictions <- do.call(rbind, prediction_rows)
  fold_predictions$Ensemble_probability <- rowMeans(
    fold_predictions[, c("LASSO_probability", "RF_probability", "SVM_probability")]
  )
  if (anyDuplicated(fold_predictions$row_id)) {
    probabilities <- c("LASSO_probability", "RF_probability", "SVM_probability",
                       "Ensemble_probability")
    predictions <- stats::aggregate(
      fold_predictions[, probabilities],
      fold_predictions[, c("row_id", "subject_id", "observed")], mean
    )
    predictions$outer_fold <- "mean_of_balanced_leave_pair_out_predictions"
    predictions <- predictions[match(rownames(x), predictions$row_id), ]
  } else {
    predictions <- fold_predictions
  }
  list(
    predictions = predictions,
    fold_predictions = fold_predictions,
    tuning = if (keep_details) do.call(rbind, tuning_rows) else NULL,
    selections = if (keep_details) do.call(rbind, selection_rows) else NULL
  )
}

performance_from_predictions <- function(predictions, dataset) {
  methods <- c(LASSO = "LASSO_probability", RF = "RF_probability",
               SVM = "SVM_probability", Ensemble = "Ensemble_probability")
  do.call(rbind, lapply(names(methods), function(method) {
    result <- classification_metrics(predictions$observed,
                                     predictions[[methods[[method]]]])
    data.frame(dataset = dataset, method = method, result, row.names = NULL)
  }))
}

bootstrap_oof_metrics <- function(predictions, dataset, task, n_boot = 2000L,
                                  seed = 9001L) {
  methods <- c(LASSO = "LASSO_probability", RF = "RF_probability",
               SVM = "SVM_probability", Ensemble = "Ensemble_probability")
  subjects <- unique(predictions$subject_id)
  output <- list()
  counter <- 0L
  set.seed(seed)
  for (b in seq_len(n_boot)) {
    if (task == "GSE103940") {
      sampled <- sample(subjects, length(subjects), replace = TRUE)
    } else {
      subject_label <- setNames(predictions$observed, predictions$subject_id)
      negative <- names(subject_label)[subject_label == 0L]
      positive <- names(subject_label)[subject_label == 1L]
      sampled <- c(sample(negative, length(negative), replace = TRUE),
                   sample(positive, length(positive), replace = TRUE))
    }
    idx <- unlist(lapply(sampled, function(s) which(predictions$subject_id == s)),
                  use.names = FALSE)
    for (method in names(methods)) {
      counter <- counter + 1L
      met <- classification_metrics(predictions$observed[idx],
                                    predictions[[methods[[method]]]][idx])
      output[[counter]] <- data.frame(iteration = b, dataset = dataset,
                                      method = method, met[, 1:4])
    }
  }
  do.call(rbind, output)
}

paired_permute_labels <- function(y, subjects, flip) {
  out <- as.integer(y)
  ug <- sort(unique(as.character(subjects)))
  stopifnot(length(flip) == length(ug))
  for (i in seq_along(ug)) {
    idx <- which(as.character(subjects) == ug[i])
    stopifnot(length(idx) == 2L, setequal(out[idx], 0:1))
    if (flip[i]) out[idx] <- 1L - out[idx]
  }
  out
}

run_permutation_set <- function(x, y_matrix, subjects, task, seed,
                                observed_y, observed_performance,
                                is_observed = rep(FALSE, ncol(y_matrix)),
                                progress_every = 25L) {
  cores <- as.integer(Sys.getenv("ML_PERMUTATION_CORES", "1"))
  cores <- max(1L, cores)
  run_one <- function(j) {
    fit <- run_loso(x, y_matrix[, j], subjects, task,
                    seed = seed, keep_details = FALSE)
    perf <- performance_from_predictions(fit$predictions, task)
    perf$permutation_id <- colnames(y_matrix)[j]
    perf$is_observed_assignment <- is_observed[j]
    perf
  }
  message(task, ": running ", ncol(y_matrix), " assignments with ", cores,
          " process(es); RNG is fixed by held-out subject/pair, not assignment index")
  if (cores == 1L) {
    output <- lapply(seq_len(ncol(y_matrix)), run_one)
  } else {
    output <- parallel::mclapply(seq_len(ncol(y_matrix)), run_one,
                                 mc.cores = cores, mc.preschedule = FALSE,
                                 mc.set.seed = FALSE)
  }
  do.call(rbind, output)
}

stability_analysis <- function(x, y, subjects, task, genes, repeats = 200L,
                               rf_null_permutations = 20L, seed = 7001L) {
  x <- as.matrix(x)
  subject_label <- tapply(y, subjects, function(z) unique(z))
  if (task == "GSE103940") {
    stopifnot(all(vapply(subject_label, length, integer(1)) == 2L))
  }
  selected_rows <- list()
  rf_importance_rows <- list()
  counter <- 0L
  rf_counter <- 0L
  set.seed(seed)
  for (b in seq_len(repeats)) {
    if (b %% 20L == 0L || b == 1L) message(task, " stability ", b, "/", repeats)
    if (task == "GSE103940") {
      ug <- sort(unique(as.character(subjects)))
      sampled_subjects <- sample(ug, ceiling(0.80 * length(ug)), replace = FALSE)
    } else {
      labels <- setNames(as.integer(y), as.character(subjects))
      negative <- names(labels)[labels == 0L]
      positive <- names(labels)[labels == 1L]
      sampled_subjects <- c(sample(negative, 4L, replace = FALSE),
                            sample(positive, 4L, replace = FALSE))
    }
    idx <- as.character(subjects) %in% sampled_subjects
    xb <- x[idx, , drop = FALSE]
    yb <- as.integer(y[idx])
    gb <- subjects[idx]
    rep_seed <- seed + b * 100000L

    lasso <- fit_lasso(xb, yb, gb, task)
    svm <- fit_svm_rfe(xb, yb, gb, task, rep_seed + 30000L)
    rf <- fit_rf(xb, yb, rep_seed + 50000L)

    null_importance <- matrix(NA_real_, nrow = length(rf$importance),
                              ncol = rf_null_permutations,
                              dimnames = list(names(rf$importance), NULL))
    for (q in seq_len(rf_null_permutations)) {
      if (task == "GSE103940") {
        flips <- sample(c(FALSE, TRUE), length(unique(gb)), replace = TRUE)
        y_null <- paired_permute_labels(yb, gb, flips)
      } else {
        y_null <- sample(yb, replace = FALSE)
      }
      rf_null <- fit_rf(xb, y_null, rep_seed + 60000L + q * 100L,
                        fixed_mtry = rf$tuning$mtry)
      null_importance[, q] <- rf_null$importance[rownames(null_importance)]
    }
    null_q95 <- apply(null_importance, 1, stats::quantile, probs = 0.95,
                      na.rm = TRUE, names = FALSE, type = 8)
    rf_selected <- names(rf$importance)[rf$importance > 0 &
                                         rf$importance > null_q95]

    for (method in c("LASSO", "RF", "SVM")) {
      selected <- switch(method, LASSO = lasso$selected,
                         RF = rf_selected, SVM = svm$selected)
      counter <- counter + 1L
      selected_rows[[counter]] <- data.frame(
        dataset = task, resample = b, method = method, gene = genes,
        selected = genes %in% selected, stringsAsFactors = FALSE
      )
    }
    rf_counter <- rf_counter + 1L
    rf_importance_rows[[rf_counter]] <- data.frame(
      dataset = task, resample = b, gene = names(rf$importance),
      observed_importance = rf$importance,
      null_q95 = null_q95[names(rf$importance)],
      selected = names(rf$importance) %in% rf_selected,
      mtry = rf$tuning$mtry, stringsAsFactors = FALSE
    )
  }
  list(selections = do.call(rbind, selected_rows),
       rf_importance = do.call(rbind, rf_importance_rows))
}

summarize_stability <- function(stability, all_targets, dataset,
                                eligible_genes, threshold = 0.60) {
  x <- stability$selections
  frequencies <- stats::aggregate(selected ~ method + gene, x, mean)
  wide <- reshape(frequencies, idvar = "gene", timevar = "method",
                  direction = "wide")
  colnames(wide) <- sub("selected\\.", "", colnames(wide))
  result <- data.frame(dataset = dataset, gene = all_targets,
                       expression_eligible = all_targets %in% eligible_genes,
                       stringsAsFactors = FALSE)
  result <- merge(result, wide, by = "gene", all.x = TRUE, sort = FALSE)
  result <- result[match(all_targets, result$gene), ]
  names(result)[names(result) == "LASSO"] <- "LASSO_frequency"
  names(result)[names(result) == "RF"] <- "RF_frequency"
  names(result)[names(result) == "SVM"] <- "SVM_frequency"
  for (column in c("LASSO_frequency", "RF_frequency", "SVM_frequency")) {
    result[[column]][!result$expression_eligible] <- NA_real_
  }
  result$methods_supported <- rowSums(
    cbind(result$LASSO_frequency >= threshold,
          result$RF_frequency >= threshold,
          result$SVM_frequency >= threshold), na.rm = TRUE
  )
  result$stable_any_method <- result$methods_supported >= 1L
  result
}
