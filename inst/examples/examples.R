# Basic usage of the package with a breast cancer expression data set

# tools -------

  library(tidyverse)
  library(rlang)
  library(microViz)
  library(trafo)
  library(fastTest)
  library(GOSemSim)

  library(org.Hs.eg.db)
  library(AnnotationDbi)

  ## for parallelization

  library(furrr)
  library(GOSemSim)

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

  Gini(counts$HERC2, unbiased = FALSE)

  Gini(counts$HERC2, unbiased = TRUE)

  counts$HERC2 %>%
    DescTools::Gini(unbiased = FALSE)

  freqRatio(counts$TSPAN6)

  percUnique(counts$TSPAN6)

  counts$TSPAN6 %>%
    caret::nzv(saveMetrics = T)

# Column stats and selection of variant genes ---------

  system.time(counts[genes] %>%
                colMedians(na.rm = TRUE))

  counts[genes] %>%
    colMedians(na.rm = TRUE)

  colMins(counts[genes]) %>% head

  counts[genes] %>% colMax %>% head

  counts[genes] %>% colGmeans %>% head

  counts[genes] %>% colHmeans %>% head

  counts[genes] %>% colVars %>% head

  counts[genes] %>%
    colSDs

  counts[genes] %>%
    colQuantiles(c(0.25, 0.5, 0.75))

  counts[genes] %>%
    colSDs

  counts[genes] %>% colGini %>% head

  counts[genes] %>% colFreqRatios %>% head

  counts[genes] %>% colPercUniques %>% head

  counts[genes] %>%
    colCI(method = 'bca')

  counts[genes] %>%
    colMissing

  ## selection of variant genes

  system.time(distr_stats(log_expression[genes]))

  #system.time(caret::nearZeroVar(log_expression[genes]))

  ## selection of genes with a Gini coefficient cutoff

  colStats <- log_expression[genes] %>%
    distr_stats

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

# Descriptive stats for numeric variables -------

  test_counts <- counts[, c('metastasis', 'histology', genes[1:100])] %>%
    filter(complete.cases(.))

  test_counts[10, genes[5]] <- NA

  colComplete(test_counts)

  fast_num_stats(test_counts,
                 split_fct = 'histology',
                 variables = genes[1:100],
                 format = 'full')

# Visualization of differences in expression of selected genes ---------

  ## differences in median expression between the histology strata

  test_counts %>%
    delta(split_fct = 'histology',
          variables = genes[1:100],
          extreme = FALSE,
          average_fun = colMedians)

  test_counts %>%
    blast(metastasis) %>%
    plot_common_change(variables = genes[1:4],
                       split_fct = 'histology')


# Shared elements ---------

  set.seed(12345)

  gene_sets <- 1:10 %>%
    map(function(x) sample(variant_genes[1:1000], size = 300))

  count_features(gene_sets)

  shared_features(gene_sets, m = 6)

# Two-sample test -------

  ## comparing gene expression between ER-positive and ER-negative cancers
  ## two-tailed T test, FDR correction for multiple testing ('BH' stands for
  ## Benjamini-Hochberg)

  er_dge <-
    f_t_test(log_expression[, variant_genes],
             f = log_expression$er_status,
             adj_method = "BH",
             as_data_frame = TRUE,
             safely = TRUE) %>%
    as_tibble

  ## identification of significantly regulated genes: pFDR < 0.05
  ## and moderate-to-large effect size (d >= 0.5)

  er_significant <- er_dge %>%
    identify_significant(label_variable = 'variable',
                         p_variable = 'p_adjusted',
                         regulation_variable = 'cohen_d',
                         regulation_level = 0.5)

  ## volcano plot: effect size and significance

  plot_volcano(data = er_dge,
               regulation_variable = 'cohen_d',
               p_variable = 'p_adjusted',
               signif_level = 0.05,
               regulation_level = 0.5,
               top_regulated = 10,
               label_variable = 'variable',
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
    filter(variable %in% unlist(er_significant)) %>%
    plot_top(regulation_variable = 'estimate',
             label_variable = 'variable',
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
                 regulation_variable = 'cohen_d',
                 label_variable = 'variable',
                 p_variable = 'p_value',
                 regulation_level = 1,
                 top_regulated = 20,
                 plot_title = 'Top regulated genes, ER status',
                 x_lab = "Effect size, Cohen's d",
                 show_txt = TRUE,
                 txt_size = 2.5,
                 txt_color = 'white') +
    theme(axis.text.y = element_text(face = 'italic'))

# ANOVA with post-hoc T test vs cohort mean -------

  ## comparing tumors of different histologies: ductal, lobular and mixed one

  histo_data <- log_expression %>%
    filter(histology %in% c('IDC', 'MIXED_IDLC', 'ILC')) %>%
    mutate(histology = droplevels(histology))

  ## testing: ANOVA with one-sample post-hoc T test for comparison
  ## of the histologies with the cohort mean

  histo_dge <-
    f_one_anova(histo_data[, variant_genes],
                f = histo_data$histology,
                adj_method = "BH",
                as_data_frame = TRUE,
                safely = TRUE) %>%
    as_tibble

  histo_t_test <- avg_deviation(histo_data,
                                split_fct = "histology",
                                variables = variant_genes,
                                adj_method = "BH")

  ## identification of significantly regulated genes: pFDR < 0.05 and
  ## eta-squared >= 0.06 in ANOVA

  histo_significant <- histo_dge %>%
    identify_significant(label_variable = 'variable',
                         p_variable = 'p_adjusted',
                         regulation_variable = 'etasq',
                         signif_level = 0.05,
                         regulation_level = 0.06)

  ## volcano plots

  histo_volcanos <-
    list(data = histo_t_test %>%
           blast(histology),
         plot_title = paste("Differential gene expression,",
                            levels(histo_t_test$histology),
                            "vs cohort mean")) %>%
    pmap(plot_volcano,
         p_variable = 'p_adjusted',
         regulation_variable = "deviation_center",
         signif_level = 0.05,
         regulation_level = 0.5,
         top_regulated = 20,
         label_variable = 'variable',
         x_lab = expression("log"[2] * "regulation, vs cohort mean"))

# Classification of the histology-regulated genes -------

  histo_class <-
    classify(data = histo_data,
             variables = histo_significant,
             split_fct = 'histology')

  histo_class$classification %>%
    group_by(histology) %>%
    slice_max(delta_auc, n = 10) %>%
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
           plot_title = 'Histology-specific genes') +
    guides(y = guide_axis(n.dodge = 2))

  histo_data <- histo_data %>%
    column_to_rownames('sample_id')

  ## cosine distance between the histologies
  ## in respect to expression of the significantly regulated genes

  library(clustTools)

  histo_distances <-
    subset_distance(histo_data,
                    variables = histo_significant,
                    split_fct = 'histology',
                    dist_FUN = calculate_dist,
                    method = 'cosine')

  histo_distances %>% plot("mean")

# GO enrichment --------

  ## A list of Entrez ID vectors

  go_universe <-
    mapIds(org.Hs.eg.db,
           keys = variant_genes,
           keytype = 'SYMBOL',
           column = 'ENTREZID')

  go_universe <- go_universe[!is.na(go_universe)]

  ## differentially regulated genes: ER status and histology

  go_input <- list(er = er_dge,
                   histo = histo_dge) %>%
    map(identify_significant,
        label_variable = 'variable',
        p_variable = 'p_value') %>%
    map(mapIds,
        x = org.Hs.eg.db,
        keytype = 'SYMBOL',
        column = 'ENTREZID') %>%
    map(~.x[!is.na(.x)])

  ## testing GO enrichment

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
    map(slice_max,
        or,
        n = 10) %>%
    map(mutate,
        log_inv_p = -log10(p_adjusted))

  go_volcanos <-
    list(data = go_enrichment,
         plot_title = paste('GO enrichment,',
                            c('ER status', 'histology'))) %>%
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
    map2_dfr(top_go, c('ER status', 'histology'),
             ~mutate(.x, condition = .y)) %>%
    mutate(condition = factor(condition,
                              c('ER status', 'histology')))

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
    slice_max(abs(estimate),
              n = 30) %>%
    arrange(estimate) %>%
    .$variable

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

# Semantic distances -----

  go_db <- godata(annoDb = org.Hs.eg.db,
                  ont = "BP",
                  computeIC = FALSE)

  gene_sem(go_input$er[1:10],
           semData = go_db,
           measure = 'Wang',
           as_dist = FALSE,
           .parallel = FALSE)

  go_dists <-
    go_sem(go_enrichment$histo$go_id,
           semData = go_db,
           measure = 'Wang',
           as_dist = FALSE,
           .parallel = TRUE)

# Binary data frequencies -----------

  ## samples with at least 100 copies of RNA

  bin_data <- counts[c('er_status', genes[1:5])]

  bin_data[genes[1:5]] <- bin_data[genes[1:5]] %>%
    map_dfc(~ifelse(.x >= 2, 'expressed', 'not_expressed'))

  ## introducing some NAs

  bin_data$TSPAN6[100] <- NA

  count_binary(data = bin_data,
               variables = genes[1:5])

  count_binary(data = bin_data,
               #variables = genes[1:5],
               split_fct = 'er_status')

# Box plot from stats ---------

  stat_tbl <-
    colQuantiles(counts[genes[1:20]],
                 c(0.5, 0.25, 0.75, 0.025, 0.975)) %>%
    set_colnames(c('median', 'perc25', 'perc75', 'perc025', 'perc975')) %>%
    as.data.frame %>%
    rownames_to_column('gene_symbol')

  box_from_stats(stat_tbl,
                 x_variable = 'gene_symbol',
                 fill_variable = 'gene_symbol')

# END ------
