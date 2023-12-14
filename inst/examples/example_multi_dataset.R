# Example of a multi-dataset analysis

# tools -------

  library(tidyverse)
  library(microViz)

# data ------

  ## stratification of the socioeconomic status (SES) class

  school_data <- MASS::nlschools %>%
    mutate(SES_class = cut(SES,
                           c(-Inf, 25, Inf),
                           c('low', 'high'))) %>%
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

# Plotting variables of interest as a function of the SES -------

  school_response_test <- test_two_groups(data = school_data,
                                          split_fct = 'SES_class',
                                          variables = c('lang', 'IQ'),
                                          type = 't') %>%
    mutate(plot_cap = paste0(effect_size_name,
                             ' = ', signif(effect_size, 2),
                             ', p = ', signif(p_adjusted, 2)))

  school_response_plots <-
    list(x = c('lang', 'IQ'),
         y = c('Language skills', 'Verbal IQ'),
         v = school_response_test$plot_cap,
         z = c('test score', 'points')) %>%
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
    set_names(c('lang', 'IQ'))

# Testing for differences in IQ and language skills in each class -------

  class_data <-
    split(school_data, school_data$class)

  ## testing for differences with Mann-Whitney U test
  ## separately for each class

  class_test <- class_data %>%
    map(safely(test_two_groups),
        split_fct = 'SES_class',
        variables = c('lang', 'IQ'),
        type = 't') %>%
    map(~.x$result) %>%
    compact %>%
    map2_dfr(., names(.),
            ~mutate(.x, class = .y)) %>%
    re_adjust(p_variable = 'p_value', method = 'BH')

  class_test %>%
    draw_estimate_scatter(label_variable = 'response',
                          regulation_variable = 'effect_size',
                          regulation_level = 0.5,
                          p_variable = 'p_value',
                          signif_level = 0.1,
                          central_stat = 'mean',
                          error_stat = 'percentile_ci',
                          plot_title = 'Effect of SES in single classes',
                          x_lab = "Cohen's d") +
    geom_vline(xintercept = 0,
               linetype = 'dashed')
