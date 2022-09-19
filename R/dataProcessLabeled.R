#' Summarize heavy and light features separately
#' @export
dataProcessLabeled = function(raw, impute = FALSE) {
  input = MSstats::MSstatsPrepareForDataProcess(raw, 2, NULL)
  input = MSstats::MSstatsNormalize(input, "NONE", NULL, NULL)
  input = MSstats::MSstatsMergeFractions(input)
  input = MSstats::MSstatsHandleMissing(input, "TMP", impute,
                                        "0", 0.9999)
  input = MSstats::MSstatsSelectFeatures(input, "all", 3,
                                         3)
  processed = MSstats::getProcessed(input)

  input = MSstatsPrepareForLabeledSummarization(input, "labeledTMP", TRUE, "0",
                                                FALSE)

  input_split = split(input, input$PROTEIN)

  summarized = MSstatsSummarizeLabeledTMP(input_split, impute, "0", FALSE)

  output = MSstatsSummarizationLabeledOutput(input, summarized,
                                             processed, "LabeledTMP", T, "0")
  output
}


MSstatsSummarizationLabeledOutput = function(input, summarized, processed,
                                             method, impute, censored_symbol) {
  LABEL = TotalGroupMeasurements = GROUP = Protein = RUN = NULL

  input = MSstats:::.finalizeInput(input, summarized, "TMP", impute, censored_symbol)
  summarized = lapply(summarized, function(x) x[[1]])
  summarized = data.table::rbindlist(summarized)
  if (inherits(summarized, "try-error")) {
    msg = paste("*** error : can't summarize per subplot with ",
                method, ".")
    getOption("MSstatsLog")("ERROR", msg)
    getOption("MSstatsMsg")("ERROR", msg)
    rqall = NULL
    rqmodelqc = NULL
    workpred = NULL
  } else {
    input[, TotalGroupMeasurements := uniqueN(.SD),
          by = c("PROTEIN", "GROUP", "LABEL"),
          .SDcols = c("FEATURE", "originalRUN")]
    cols = intersect(c("PROTEIN", "originalRUN", "RUN", "GROUP",
                       "GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "LABEL",
                       "TotalGroupMeasurements",
                       "NumMeasuredFeature", "MissingPercentage",
                       "more50missing", "NumImputedFeature"),
                     colnames(input))
    merge_col = ifelse(is.element("RUN", colnames(summarized)),
                       "RUN", "SUBJECT_ORIGINAL")
    input_merge = input[, cols, with = FALSE]
    input_merge = input_merge[, colnames(input_merge) != "GROUP", with = FALSE]
    rqall = merge(summarized, input_merge, by.x = c(merge_col, "Protein", "LABEL"),
                  by.y = c(merge_col, "PROTEIN", "LABEL"))
    data.table::setnames(rqall, c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL"),
                         c("GROUP", "SUBJECT"), skip_absent = TRUE)

    rqall$GROUP = factor(as.character(rqall$GROUP))
    rqall$Protein = factor(rqall$Protein)
    rqmodelqc = summarized$ModelQC
  }

  if (is.element("RUN", colnames(rqall)) & !is.null(rqall)) {
    rqall = rqall[order(Protein, as.numeric(as.character(RUN))), ]
    rownames(rqall) = NULL
  }
  output_cols = intersect(c("PROTEIN", "PEPTIDE", "TRANSITION", "FEATURE",
                            "LABEL", "GROUP", "RUN", "SUBJECT", "FRACTION",
                            "originalRUN", "censored", "INTENSITY", "ABUNDANCE",
                            "newABUNDANCE", "predicted", "feature_quality",
                            "is_outlier", "remove"), colnames(input))
  input = input[, output_cols, with = FALSE]

  if (is.element("remove", colnames(processed))) {
    processed = processed[(remove),
                          intersect(output_cols,
                                    colnames(processed)), with = FALSE]
    input = rbind(input, processed, fill = TRUE)
  }
  list(FeatureLevelData = as.data.frame(input),
       ProteinLevelData = as.data.frame(rqall),
       SummaryMethod = method)
}


.finalizeInput = function(input, summarized, method, impute, censored_symbol) {
  if (method == "TMP") {
    input = .finalizeTMP(input, censored_symbol, impute, summarized)
  } else if (method == "linear") {
    input = .finalizeLinear(input, censored_symbol)
  } else {
    input = .finalizeLabeledTMP(input, censored_symbol)
  }
  input
}

.finalizeLabeledTMP = function(input, censored_symbol, impute, summarized) {
  NonMissingStats = NumMeasuredFeature = MissingPercentage = LABEL = NULL
  total_features = more50missing = nonmissing_orig = censored = NULL
  INTENSITY = newABUNDANCE = NumImputedFeature = NULL

  survival_predictions = lapply(summarized, function(x) x[[2]])
  predicted_survival = data.table::rbindlist(survival_predictions)
  if (impute) {
    cols = intersect(colnames(input), c("newABUNDANCE",
                                        "cen", "RUN",
                                        "FEATURE", "ref"))
    input = merge(input[, colnames(input) != "newABUNDANCE", with = FALSE],
                  predicted_survival,
                  by = setdiff(cols, "newABUNDANCE"),
                  all.x = TRUE)
  }
  input[, NonMissingStats := .getNonMissingFilterStats(.SD, censored_symbol)]
  input[, NumMeasuredFeature := sum(NonMissingStats),
        by = c("PROTEIN", "RUN")]
  input[, MissingPercentage := 1 - (NumMeasuredFeature / total_features)]
  input[, more50missing := MissingPercentage >= 0.5]
  if (!is.null(censored_symbol)) {
    if (is.element("censored", colnames(input))) {
      input[, nonmissing_orig := LABEL == "L" & !censored]
    } else {
      input[, nonmissing_orig := LABEL == "L" & !is.na(INTENSITY)]
    }
    input[, nonmissing_orig := ifelse(is.na(newABUNDANCE), TRUE, nonmissing_orig)]
    if (impute) {
      input[, NumImputedFeature := sum(LABEL == "L" & !nonmissing_orig),
            by = c("PROTEIN", "RUN")]
    } else {
      input[, NumImputedFeature := 0]
    }
  }
  input
}


MSstatsSummarizeLabeledTMP = function(proteins_list, impute, censored_symbol,
                                      remove50missing) {
  num_proteins = length(proteins_list)
  summarized_results = vector("list", num_proteins)

  pb = utils::txtProgressBar(min = 0, max = num_proteins, style = 3)
  for (protein_id in seq_len(num_proteins)) {
    single_protein = proteins_list[[protein_id]]
    summarized_results[[protein_id]] = MSstatsSummarizeSingleLabeledTMP(
      single_protein, impute, censored_symbol, remove50missing)
    setTxtProgressBar(pb, protein_id)
  }
  close(pb)

  summarized_results
}


MSstatsSummarizeSingleLabeledTMP = function(single_protein, impute, censored_symbol,
                                            remove50missing) {
  newABUNDANCE = n_obs = n_obs_run = RUN = FEATURE = LABEL = NULL
  predicted = censored = NULL
  cols = intersect(colnames(single_protein), c("newABUNDANCE", "cen", "RUN",
                                               "FEATURE", "ref"))
  single_protein = single_protein[(n_obs > 1 & !is.na(n_obs)) &
                                    (n_obs_run > 0 & !is.na(n_obs_run))]
  if (nrow(single_protein) == 0) {
    return(list(NULL, NULL))
  }
  single_protein[, RUN := factor(RUN)]
  single_protein[, FEATURE := factor(FEATURE)]
  if (impute & any(single_protein[["censored"]])) {
    survival_fit = MSstats:::.fitSurvival(single_protein[LABEL == "L", cols,
                                                         with = FALSE])
    single_protein[, predicted := predict(survival_fit,
                                          newdata = .SD)]
    single_protein[, predicted := ifelse(censored & (LABEL == "L"), predicted, NA)]
    single_protein[, newABUNDANCE := ifelse(censored & LABEL == "L",
                                            predicted, newABUNDANCE)]
    survival = single_protein[, c(cols, "predicted"), with = FALSE]
  } else {
    survival = single_protein[, cols, with = FALSE]
    survival[, predicted := NA]
  }

  single_protein = MSstats:::.isSummarizable(single_protein, remove50missing)
  if (is.null(single_protein)) {
    return(list(NULL, NULL))
  } else {
    single_protein = single_protein[!is.na(newABUNDANCE), ]
    result = lapply(split(single_protein, single_protein$LABEL),
                    function(x) {
                      if (nrow(x) > 0) {
                        label = unique(x$LABEL)
                        result = MSstats:::.fitTukey(x)
                        result[, LABEL := label]
                        result[, Protein := unique(x$PROTEIN)]
                        result
                      } else {
                        NULL
                      }
                    })
    result = rbindlist(result)
  }
  list(result, survival)
}


.runTukeyLabeled = function(single_protein, censored_symbol, remove50missing) {
  Protein = RUN = newABUNDANCE = NULL

  if (nlevels(single_protein$FEATURE) > 1) {
    tmp_result = .fitTukey(single_protein)
  } else {
    tmp_result = single_protein[, list(RUN, LABEL, LogIntensities = newABUNDANCE)]
  }
  tmp_result[, Protein := unique(single_protein$PROTEIN)]
  tmp_result
}


MSstatsPrepareForLabeledSummarization = function(input, method, impute, censored_symbol,
                                                 remove_uninformative_feature_outlier) {
  ABUNDANCE = feature_quality = is_outlier = PROTEIN = NULL

  label = data.table::uniqueN(input$LABEL) == 2
  if (label) {
    input$ref = factor(ifelse(input$LABEL == "L",
                              input$RUN, 0))
  }

  if (is.element("remove", colnames(input))) {
    input = input[!(remove)]
  }

  if (remove_uninformative_feature_outlier &
      is.element("feature_quality", colnames(input))) {
    input[, ABUNDANCE := ifelse(feature_quality == "Uninformative",
                                NA, ABUNDANCE)]
    input[, ABUNDANCE := ifelse(is_outlier, NA, ABUNDANCE)]
    msg = "** Filtered out uninformative features and outliers."
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
  }

  input = .prepareSummary(input, method, impute, censored_symbol)
  input[, PROTEIN := factor(PROTEIN)]
  input
}


.prepareSummary = function(input, method, impute, censored_symbol) {
  if (method == "TMP") {
    input = .prepareTMP(input, impute, censored_symbol)
  } else if (method == "linear") {
    input = .prepareLinear(input, FALSE, censored_symbol)
  } else if (method == "labeledTMP") {
    input = .prepareLabeledTMP(input, impute, censored_symbol)
  }
  input
}


.prepareLabeledTMP = function(input, impute, censored_symbol) {
  censored = feature_quality = newABUNDANCE = cen = nonmissing = n_obs = NULL
  n_obs_run = total_features = FEATURE = prop_features = NULL
  remove50missing = ABUNDANCE = NULL

  if (impute & !is.null(censored_symbol)) {
    if (is.element("feature_quality", colnames(input))) {
      input[, censored := ifelse(feature_quality == "Informative",
                                 censored, FALSE)]
    }
    if (censored_symbol == "0") {
      input[, newABUNDANCE := ifelse(censored, 0, ABUNDANCE)]
    } else if (censored_symbol == "NA") {
      input[, newABUNDANCE := ifelse(censored, NA, ABUNDANCE)]
    }
    input[, cen := ifelse(censored, 0, 1)]
  } else {
    input[, newABUNDANCE := ABUNDANCE]
  }

  input[, nonmissing := .getNonMissingLabeledFilter(input, impute, censored_symbol)]
  input[, n_obs := sum(nonmissing), by = c("PROTEIN", "FEATURE")]
  input[, nonmissing := ifelse(n_obs <= 1, FALSE, nonmissing)]
  input[, n_obs_run := sum(nonmissing), by = c("PROTEIN", "RUN")]

  input[, total_features := uniqueN(FEATURE), by = "PROTEIN"]
  input[, prop_features := sum(nonmissing) / total_features,
        by = c("PROTEIN", "RUN")]

  if (is.element("cen", colnames(input))) {
    if (any(input[["cen"]] == 0)) {
      MSstats:::.setCensoredByThreshold(input, censored_symbol, remove50missing)
    }
  }

  input
}


.getNonMissingLabeledFilter = function(input, impute, censored_symbol) {
  if (impute) {
    if (!is.null(censored_symbol)) {
      if (censored_symbol == "0") {
        nonmissing_filter = !is.na(input$newABUNDANCE) & input$newABUNDANCE != 0
      } else if (censored_symbol == "NA") {
        nonmissing_filter = !is.na(input$newABUNDANCE)
      }
    }
  } else {
    nonmissing_filter = !is.na(input$newABUNDANCE) & input$newABUNDANCE != 0
  }
  nonmissing_filter
}

#' Plot protein-level summaries separately for heavy and light features
#' @import ggplot2
#' @export
plotLabeledProfiles = function(summarized, proteins) {
  protein_level = summarized[["ProteinLevelData"]][summarized[["ProteinLevelData"]][["Protein"]] %in% proteins, ]
  feature_level = summarized[["FeatureLevelData"]][summarized[["FeatureLevelData"]][["PROTEIN"]] %in% proteins, ]

  protein_level = as.data.table(protein_level)
  feature_level = as.data.table(feature_level)
  setnames(protein_level,
           c("RUN", "Protein", "LABEL", "LogIntensities", "GROUP"),
           c("Run", "ProteinName", "Label", "LogIntensity", "Condition"))
  setnames(feature_level,
           c("RUN", "PROTEIN", "FEATURE", "LABEL", "newABUNDANCE", "GROUP"),
           c("Run", "ProteinName", "PSM", "Label", "LogIntensity", "Condition"))
  feature_level = feature_level[LogIntensity > 0]

  ggplot(feature_level,
         aes(x = Run, y = LogIntensity, group = paste(PSM, Label))) +
    geom_point(color = "grey") +
    geom_line(color = "grey") +
    geom_point(aes(x = Run, y = LogIntensity),
               data = protein_level, size = 1.2, inherit.aes = FALSE) +
    geom_line(aes(x = Run, y = LogIntensity, group = paste(ProteinName, Label)),
              data = protein_level, size = 1.5, inherit.aes = FALSE) +
    facet_grid(ProteinName ~ Label) +
    theme(legend.position = "bottom") +
    ylab("log-Intensity") +
    xlab("run ID")
}
