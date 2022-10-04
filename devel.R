# Testing during development

# data ------

  test_data <- tibble::tibble(estimate = rnorm(500, sd = 1) + sample(seq(0, 1, by = 1e-4), 500),
                              gene_name = paste0('gene_', 1:500))

  test_data <- dplyr::mutate(test_data,
                             p_value = 1 - pnorm(abs(estimate),
                                                 sd = 1,
                                                 lower.tail = TRUE),
                             p_value = p_value + sample(seq(0, 0.05, by = 1e-6), 500),
                             lower_conf = estimate + qnorm(0.025),
                             upper_conf = estimate + qnorm(0.975))

  test_expr <- tibble::tibble(group = factor(sample(c('A', 'B'),
                                                    size = 500,
                                                    replace = TRUE)))

  for(i in 1:1000) {

    gene_name <- paste0('gene ', i)

    test_expr <- dplyr::mutate(test_expr,
                               !!gene_name := rnorm(500))

  }

  test_expr$dummy <- 1

# volcano -----

  test_volcano <- plot_volcano(data = test_data,
                               regulation_variable = 'estimate',
                               p_variable = 'p_value',
                               regulation_level = 1,
                               signif_level = 0.05,
                               fill_scale = c('coral3',
                                              'darkolivegreen',
                                              'gray80'),
                               txt_size = 2.5,
                               label_variable = 'gene_name',
                               top_significant = 10,
                               txt_face = 3,
                               plot_title = 'Test volcano',
                               plot_subtitle = 'gene comparison',
                               point_alpha = 0.5,
                               label_type = 'text')

# GSEA -----

  test_sign <- plot_sign(data = test_data,
                         regulation_variable = 'estimate',
                         p_variable = 'p_value',
                         signif_level = 0.05,
                         regulation_level = 1)

# Forest -----

  test_selected_forest <- plot_forest(data = test_data[1:10, ],
                                      regulation_variable = 'estimate',
                                      label_variable = 'gene_name',
                                      p_variable = 'p_value',
                                      signif_level = 0.05,
                                      lower_ci_variable = 'lower_conf',
                                      upper_ci_variable = 'upper_conf',
                                      show_txt = TRUE,
                                      show_ci_txt = TRUE)

  test_forest <- plot_top(data = test_data,
                          regulation_variable = 'estimate',
                          label_variable = 'gene_name',
                          p_variable = 'p_value',
                          signif_level = 0.05,
                          lower_ci_variable = 'lower_conf',
                          upper_ci_variable = 'upper_conf',
                          show_txt = TRUE,
                          show_ci_txt = TRUE,
                          top_regulated = 10)

# Significance bar ------

  test_bar <- plot_signifcant(data = test_data,
                              p_variable = 'p_value',
                              label_variable = 'gene_name',
                              top_significant = 20,
                              signif_level = 0.05)

# testing ------

  microViz:::t_tester(data = test_expr,
                      split_fct = 'group',
                      variable = 'dummy')

  microViz:::wilcox_tester(data = test_expr,
                           split_fct = 'group',
                           variable = 'gene_1',
                           conf.int = TRUE)

  test_comp <- test_two_groups(data = test_expr,
                               split_fct = 'group',
                               variables = paste0('gene_', 1:1000),
                               .parallel = TRUE,
                               type = 'wilcox',
                               conf.int = TRUE)

# ANOVA ------

  microViz:::anova_tester(data = test_expr,
                          split_fct = 'group',
                          variable = 'gene_1')

  microViz:::anova_tester(data = test_expr,
                          split_fct = 'group',
                          variable = 'dummy1')

  test_aov <- test_anova(data = test_expr,
                         split_fct = 'group',
                         variables = paste0('gene ', 1:1000),
                         .parallel = FALSE)

