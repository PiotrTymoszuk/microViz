# Example of a multi-data set analysis

# tools -------

  library(tidyverse)
  library(microViz)

# data ------

  ## stratification of the socioeconomic status (SES) class

  school_data <- MASS::nlschools %>%
    mutate(SES_class = cut(SES,
                           c(-Inf, 25, Inf),
                           c('low', 'high')),
           IQ_sum = IQ + lang,
           IQ_sum_norm = zScores(IQ) + zScores(lang)) %>%
    as_tibble

  ## I'm including only classes with at least 20 students
  ## and with student with both low and high SES

  classes <- school_data %>%
    count(class) %>%
    filter(n >= 20) %>%
    .$class %>%
    as.character

  classes_both <- school_data %>%
    count(class, SES_class)

  classes_both <- split(classes_both, classes_both$class) %>%
    map_dfr(function(x) if(nrow(x) == 2) x else NULL) %>%
    .$class %>%
    unique

  classes <- intersect(classes, classes_both)

  school_data <- school_data %>%
    filter(class %in% classes) %>%
    mutate(class = droplevels(class))

  ## variables of interest

  vars <- c('lang', 'IQ', 'IQ_sum', 'IQ_sum_norm', 'GS')

# Plotting variables of interest as a function of the SES -------

  school_response_test <-
    test_two_groups(data = school_data,
                    split_fct = 'SES_class',
                    variables = vars,
                    type = 't') %>%
    mutate(plot_cap = paste0(effect_size_name,
                             ' = ', signif(effect_size, 2),
                             ', p = ', signif(p_adjusted, 2)))

  school_response_plots <-
    list(x = vars,
         y = c('Language skills',
               'Verbal IQ',
               'Skill sum',
               'Normalized skill sum',
               'Class size'),
         v = school_response_test$plot_cap,
         z = c('test score',
               'points',
               'points',
               'Z-scores',
               'Class size')) %>%
    pmap(function(x, y, v, z) school_data %>%
           ggplot(aes(x = SES_class,
                      y = .data[[x]],
                      fill = SES_class)) +
           geom_boxplot() +
           geom_point(shape = 16,
                      position = position_jitter(0.1, 0.05),
                      size = 2,
                      alpha = 0.25) +
           theme_classic() +
           theme(legend.position = 'none') +
           labs(title = y,
                subtitle = v,
                y = z,
                x = 'Socioeconomic status')) %>%
    set_names(vars)

# Testing for differences in IQ and language skills in each class -------

  class_data <-
    split(school_data, school_data$class)

  ## testing for differences with Mann-Whitney U test
  ## separately for each class

  class_test <- class_data %>%
    map(safely(test_two_groups),
        split_fct = 'SES_class',
        variables = vars,
        type = 't') %>%
    map(~.x$result) %>%
    compact %>%
    map(mutate,
        regulation = ifelse(p_value >= 0.05, 'ns',
                            ifelse(effect_size > 0, 'positive',
                                   ifelse(effect_size < 0,
                                          'negative', 'ns'))),
        regulation = factor(regulation, c('positive', 'negative', 'ns')))

  ## effect sizes and significance in single classes

  estimate_plot <- class_test %>%
    plot_common_estimates(label_variable = 'response',
                          regulation_variable = 'effect_size',
                          regulation_status = 'regulation',
                          plot_title = 'Effect of SES in single classes',
                          x_lab = "Cohen's d") +
    geom_vline(xintercept = 0,
               linetype = 'dashed')

  ## numbers of classes with significant effects and average effect sizes

  effect_plot <- class_test %>%
    plot_common_effect(label_variable = 'response',
                       regulation_variable = 'estimate',
                       regulation_status = 'regulation',
                       plot_title = 'Effect of SES in single classes',
                       x_lab = "difference in score",
                       y_lab = 'Number of classes',
                       plot_subtitle = paste('classes: n =', length(class_test)),
                       cutoff = 10,
                       show_whiskers = FALSE)

# Visualization: multi-class heat map and before-after plots -------

  ## classification of the variables and common heat map

  var_classification <- class_data %>%
    map(classify,
        variables = vars[-5],
        split_fct = 'SES_class')

  var_classification <- var_classification %>%
    map(~.x$classification)

  cm_heat_plot <- class_data %>%
    common_heat_map(variables = vars[-5],
                    split_fct = 'SES_class',
                    normalize = TRUE,
                    norm_center = 'mean',
                    norm_dispersion = 'sem',
                    variable_classification = var_classification)

  ## before - after plots

  common_before_after <- class_data %>%
    plot_common_change(variables = vars[-5],
                       split_fct = 'SES_class',
                       dodge = 0,
                       show_whiskers = TRUE,
                       txt_color = 'black',
                       txt_position = 'right',
                       palette = c('steelblue', 'orangered3'),
                       normalize = FALSE)

# Differences in variables of interest between the SES classes ------

  class_data %>%
    map(delta,
        split_fct = 'SES_class',
        variables = c('lang', 'IQ'),
        average_fun = colMedians)

# END -------
