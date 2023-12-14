# Basic usage of the package with a breast cancer expression data set

# tools -------

  library(tidyverse)
  library(rlang)
  library(microViz)
  library(trafo)

  library(org.Hs.eg.db)
  library(AnnotationDbi)

  ## for parallelization

  library(furrr)

  select <- dplyr::select
  reduce <- purrr::reduce

# data -------

  data("brca")

  counts <- brca

  ## listing all available genes

  genes <- names(brca)

  genes <- genes[!genes %in% c('sample_id',
                               'patient_id',
                               'timepoint',
                               'metastasis',
                               'histology',
                               'er_status')]

  ## a data frame with log2-transformed variables
  ## increasing by 1 to avoid potential log2(zero)

  log_expression <- counts

  log_expression[genes] <- log_expression[genes] %>%
    map_dfc(~log2(.x + 1))

# Numeric vector stats -------

  ## geometric mean

  Gmean(counts$HERC2)

  counts$HERC2 %>%
    psych::geometric.mean()

  counts[, c("HERC1", "HERC2")] %>%
    as.matrix %>%
    Hmean

  Hmean(counts$HERC2)

  perCI(counts$HERC2)

  counts$HERC2 %>%
    quantile(c(0.025, 0.975))

  bcaCI(counts$HERC2)

  counts$HERC2 %>%
    coxed::bca()

  Gini(counts$HERC2, unbiased = FALSE)

  Gini(counts$HERC2, unbiased = TRUE)

  counts$HERC2 %>%
    DescTools::Gini(unbiased = FALSE)

  freqRatio(counts$TSPAN6)

  percUnique(counts$TSPAN6)

  counts$TSPAN6 %>%
    caret::nzv(saveMetrics = T)

# Column stats and selection of variant genes ---------

  counts[genes] %>%
    colMedians(na.rm = TRUE)

  counts[genes] %>%
    colMins

  counts[genes] %>%
    colMax

  counts[genes] %>%
    colGmeans

  counts[genes] %>%
    colHmeans

  counts[genes] %>%
    colVars

  counts[genes] %>%
    colSDs

  counts[genes] %>%
    colQuantiles(c(0.25, 0.5, 0.75))

  counts[genes] %>%
    colSDs

  counts[genes] %>%
    colGini

  counts[genes] %>%
    colFreqRatios

  counts[genes] %>%
    colPercUniques

  counts[genes] %>%
    colCI(method = 'bca')


  ## selection of variant genes

  colStats <- log_expression[genes] %>%
    distr_stats(freqCut = 9, uniqueCut = 10)

  variant_genes <- colStats %>%
    filter(gini_coef >= 0.1,
           freqRatio < 5) %>%
    .$variable

# Row stats -------

  trans_counts <- counts %>%
    column_to_rownames('sample_id') %>%
    select(all_of(genes)) %>%
    t %>%
    as.data.frame

  microViz::rowMedians(trans_counts)

  rowMins(trans_counts)

  microViz::rowMax(trans_counts)

  rowGmeans(trans_counts)

  rowHmeans(trans_counts)

  rowVars(trans_counts)

  rowSDs(trans_counts)

  rowGini(trans_counts)

  rowQuantiles(trans_counts, c(0, 0.25, 0.5, 0.75, 1))

  rowCI(trans_counts, method = 'percentile')

  rowCI(trans_counts, method = 'bca')

  row_stats(trans_counts)

# Shared elements ---------

  set.seed(12345)

  gene_sets <- 1:10 %>%
    map(function(x) sample(variant_genes[1:1000], size = 300))

  shared_features(gene_sets, m = 6)

# Two-sample test -------

  ## comparing gene expression between ER-positive and ER-negative cancers
  ## two-tailed T test, FDR correction for multiple testing ('BH' stands for
  ## Benjamini-Hochberg)

  er_dge <- log_expression %>%
    filter(!is.na(er_status)) %>%
    test_two_groups(split_fct = 'er_status',
                    type = 't',
                    variables = variant_genes,
                    adj_method = 'BH',
                    .parallel = TRUE)

  ## identification of significantly regulated genes: pFDR < 0.05
  ## and moderate-to-large effect size (d >= 0.5)

  er_significant <- er_dge %>%
    identify_significant(label_variable = 'response',
                         p_variable = 'p_adjusted',
                         regulation_variable = 'effect_size',
                         regulation_level = 0.5)

  ## volcano plot: effect size and significance

  plot_volcano(data = er_dge,
               regulation_variable = 'effect_size',
               p_variable = 'p_adjusted',
               signif_level = 0.05,
               regulation_level = 0.5,
               top_regulated = 10,
               label_variable = 'response',
               x_lab = "Effect size of regulation, Cohen's d, ER- vs ER+",
               plot_title = "Differential gene expression, ER status")

  ## volcano plot: regulation and significance

  plot_volcano(data = er_dge,
               regulation_variable = 'estimate',
               p_variable = 'p_adjusted',
               signif_level = 0.05,
               regulation_level = log2(1.5),
               x_lab = expression("log"[2] * "regulation,  ER- vs ER+"),
               plot_title = "Differential gene expression, ER status")

  ## forest plot for the top-regulated genes

  er_dge %>%
    filter(response %in% unlist(er_significant)) %>%
    plot_top(regulation_variable = 'estimate',
             label_variable = 'response',
             p_variable = 'p_adjusted',
             signif_level = 0.05,
             regulation_level = 0,
             lower_ci_variable = 'lower_ci',
             upper_ci_variable = 'upper_ci',
             top_regulated = 30,
             plot_title = 'Top regulated genes',
             x_lab = expression("log"[2] * "regulation,  ER- vs ER+")) +
    theme(axis.text.y = element_text(face = 'italic'))

  plot_regulated(data = er_dge,
                 regulation_variable = 'effect_size',
                 label_variable = 'response',
                 p_variable = 'p_value',
                 regulation_level = 1,
                 top_regulated = 20,
                 plot_title = 'Top regulated genes, ER status',
                 x_lab = "Effect size, Cohen's d",
                 show_txt = TRUE,
                 txt_size = 2.5,
                 txt_color = 'white') +
    theme(axis.text.y = element_text(face = 'italic'))

# ANOVA -------

  ## comparing tumors of different histologies: ductal, lobular and mixed one

  histo_data <- log_expression %>%
    filter(histology %in% c('IDC', 'MIXED_IDLC', 'ILC')) %>%
    mutate(histology = droplevels(histology))

  histo_dge <- histo_data %>%
    test_anova(split_fct = 'histology',
               variables = variant_genes,
               adj_method = 'BH',
               .parallel = TRUE)

  ## identification of significantly regulated genes: pFDR < 0.05 and
  ## eta-squared >= 0.06 in ANOVA

  histo_significant <- histo_dge$anova %>%
    identify_significant(label_variable = 'response',
                         p_variable = 'p_adjusted',
                         regulation_variable = 'effect_size',
                         signif_level = 0.05,
                         regulation_level = 0.06)

  ## volcano plots

  histo_dge$lm %>%
    filter(level == 'MIXED_IDLC') %>%
    plot_volcano( regulation_variable = 'estimate',
                  p_variable = 'p_adjusted',
                  signif_level = 0.05,
                  regulation_level = log2(1.5),
                  x_lab = expression("log"[2] * "regulation, mixed vs ductal carcinoma"),
                  plot_title = "Differential gene expression, mixed-histology cancers")

  histo_dge$lm %>%
    filter(level == 'ILC') %>%
    plot_volcano( regulation_variable = 'estimate',
                  p_variable = 'p_adjusted',
                  signif_level = 0.05,
                  regulation_level = log2(1.5),
                  x_lab = expression("log"[2] * "regulation, luminal vs ductal carcinoma"),
                  plot_title = "Differential gene expression, luminal cancers")

# Classification of the histology-regulated genes -------

  histo_class <-
    classify(data = histo_data,
             variables = histo_significant,
             split_fct = 'histology')

  histo_class$classification %>%
    group_by(histology) %>%
    top_n(n = 10,
          delta_auc) %>%
    ungroup %>%
    plot_wordcloud(label_variable = 'variable',
                   split_fct = 'histology',
                   size_variable = 'auc',
                   color_variable = 'auc',
                   show.legend = TRUE,
                   shape = 'square',
                   fontface = 'italic',
                   fraction_rotated = 0,
                   nrow = 2)

  heat_map(data = histo_data,
           variables = histo_significant,
           split_fct = 'histology',
           midpoint = 0,
           limits = c(-3, 3),
           oob = scales::squish,
           hide_x_axis_text = TRUE,
           facet = TRUE,
           plot_title = 'Histology-specific genes')

  histo_data <- histo_data %>%
    column_to_rownames('sample_id')

  library(clustTools)

  histo_distances <-
    subset_distance(histo_data,
                    variables = histo_significant,
                    split_fct = 'histology',
                    dist_FUN = calculate_dist,
                    method = 'cosine')

  histo_distances %>% plot

# Modeling with a confounder -------

  ## comparing the T1 and T2 time points
  ## choosing patients with both time points present

  time_data <- log_expression %>%
    select(sample_id,
           patient_id,
           timepoint,
           all_of(genes)) %>%
    filter(timepoint %in% c('T1', 'T2')) %>%
    filter(complete.cases(.)) %>%
    blast(patient_id) %>%
    map_dfr(function(df) if(all(c('T1', 'T2') %in% df$timepoint)) df else NULL)

  time_data <- time_data %>%
    mutate(timepoint = droplevels(timepoint))

  time_dge <- time_data %>%
    test_anova(split_fct = 'timepoint',
               variables = variant_genes,
               confounder = 'patient_id',
               adj_method = 'BH',
               .parallel = TRUE)

  time_signifcant <- time_dge$anova %>%
    identify_significant(label_variable = 'response',
                         p_variable = 'p_value',
                         regulation_variable = 'effect_size',
                         signif_level = 0.05,
                         regulation_level = 0.06)

  heat_map(data = time_data,
           variables = time_signifcant,
           split_fct = 'timepoint')

# GO enrichment --------

  ## A list of Entrez ID vectors

  go_universe <-
    mapIds(org.Hs.eg.db,
           keys = variant_genes,
           keytype = 'SYMBOL',
           column = 'ENTREZID')

  go_universe <- go_universe[!is.na(go_universe)]

  go_input <- list(er = er_dge,
                   histo = histo_dge$anova,
                   time = time_dge$anova) %>%
    map(identify_significant,
        label_variable = 'response',
        p_variable = 'p_value') %>%
    map(mapIds,
        x = org.Hs.eg.db,
        keytype = 'SYMBOL',
        column = 'ENTREZID') %>%
    map(~.x[!is.na(.x)])

  plan('multisession')

  go_enrichment <- go_input %>%
    future_map(GOana,
               universe = go_universe,
               ontology = 'BP',
               .options = furrr_options(seed = TRUE))

  plan('sequential')

  ## significant and top 10 most enriched GOs

  go_enrichment <- go_enrichment %>%
    map(filter, p_adjusted < 0.05)

  top_go <- go_enrichment %>%
    map(top_n,
        n = 10,
        or) %>%
    map(mutate,
        log_inv_p = -log10(p_adjusted))

  go_volcanos <-
    list(data = go_enrichment,
         plot_title = paste('GO enrichment,',
                            c('ER status', 'histology', 'time point'))) %>%
    pmap(plot_volcano,
         regulation_variable = 'or',
         p_variable = 'p_adjusted',
         show_vlines = 'right',
         label_variable = 'term',
         top_regulated = 10,
         label_type = 'text') %>%
    map(~.x +
          scale_x_continuous(trans = 'log2'))

  go_enrichment_bars <-
    plot_regulated(data = go_enrichment$histo,
                   regulation_variable = 'or',
                   p_variable = 'p_adjusted',
                   label_variable = 'term',
                   top_regulated = 20,
                   plot_title = 'Top enriched GOs, histology',
                   x_lab = 'OR, enrichment over the genome',
                   show_txt = TRUE)  +
    scale_y_discrete(labels = function(x) map_chr(x, wrap_text, len = 4))

# Word cloud plots -------

  top_go_words <-
    map2_dfr(top_go, c('ER status', 'histology', 'time point'),
             ~mutate(.x, condition = .y)) %>%
    mutate(condition = factor(condition,
                              c('ER status', 'histology', 'time point')))

  plot_wordcloud(data = top_go_words,
                 split_fct = 'condition',
                 label_variable = 'go_id',
                 color_variable = 'log_inv_p',
                 size_variable = 'or',
                 size_type = 'size',
                 size_range = c(0, 4),
                 wrap = TRUE,
                 cust_theme = theme_minimal(),
                 show.legend = TRUE,
                 color_lab = 'Significance',
                 size_lab = 'Enrichment',
                 shape = 'square',
                 nrow = 2,
                 fraction_rotated = 0.5)




# Bubble plots, top regulated genes -------

  ## top regulated genes in ER- vs ER+ tumor

  top_er_genes <- er_dge %>%
    top_n(n = 30,
          wt = abs(estimate)) %>%
    arrange(estimate) %>%
    .$response

  ## average expression values

  top_er_expression <- log_expression %>%
    select(er_status, all_of(top_er_genes)) %>%
    filter(complete.cases(.)) %>%
    pivot_longer(cols = all_of(top_er_genes),
                 names_to = 'gene',
                 values_to = 'log2_expression') %>%
    mutate(gene = factor(gene, top_er_genes)) %>%
    group_by(er_status, gene) %>%
    summarise(log2_expression = mean(log2_expression)) %>%
    ungroup

  top_er_bubble <- plot_bubble(data = top_er_expression,
                               x_variable = 'er_status',
                               y_variable = 'gene',
                               size_variable = 'log2_expression',
                               color_variable = 'log2_expression',
                               plot_title = 'Top regulated genes, ER status',
                               x_lab = 'ER status',
                               y_lab = 'Gene',
                               color_lab = 'expression',
                               size_lab = 'expression',
                               size_type = 'area') +
    guides(fill = 'legend',
           size = 'legend') +
    theme(axis.text.y = element_text(face = 'italic'))

# Deviations from center --------

  histo_deviations <-
    avg_deviation(data = log_expression,
                  variables = genes,
                  split_fct = 'histology',
                  grand_center = 'mean',
                  split_center = 'mean',
                  ci = 'percentile')

# Aggregation --------

  aggRows(iris[1:4], f = iris$Species, fun = colMeans)

  test_biopsy <- MASS::biopsy

  large_biopsy <- rep(list(test_biopsy), 1000) %>%
    do.call('rbind', .)

  agg_biopsy <- aggRows(as.matrix(test_biopsy[, 2:10]),
                        f = test_biopsy$ID,
                        fun = colMins)

  agg_large_biopsy <- aggRows(as.matrix(large_biopsy[, 2:10]),
                              f = large_biopsy$ID,
                              fun = colMins,
                              na.rm = TRUE,
                              .parallel = FALSE)

  agg_col_biopsy <- aggCols(as.data.frame(agg_large_biopsy[, 1:9]),
                            f = c('cat1', 'cat1',
                                  'cat2', 'cat2',
                                  'cat3', 'cat2',
                                  'cat1', 'cat3',
                                  'cat4'),
                            fun = rowGmeans)

  summary(agg_large_biopsy)
  summary(agg_biopsy)
  summary(test_biopsy)

# END ------
