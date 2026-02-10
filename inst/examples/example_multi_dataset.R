# Example of a multi-data set analysis

# tools -------

  library(tidyverse)
  library(microViz)
  library(fastTest)

# data ------

  ## stratification of the socioeconomic status (SES) class

  school_data <- MASS::nlschools %>%
    mutate(SES_class = cut(SES,
                           c(-Inf, 25, Inf),
                           c("low", "high")),
           IQ_sum = IQ + lang,
           IQ_sum_norm = zScores(IQ) + zScores(lang)) %>%
    as_tibble

  ## I"m including only classes with at least 20 students
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

  vars <- c("lang", "IQ", "IQ_sum", "IQ_sum_norm", "GS")

# Plotting variables of interest as a function of the SES -------

  school_response_test <-
    f_t_test(school_data[, vars],
             f = school_data$SES,
             as_data_frame = TRUE,
             adj_method = "BH") %>%
    mutate(plot_cap = paste0("d = ", signif(cohen_d, 2),
                             ", p = ", signif(p_adjusted, 2))) %>%
    as_tibble

  school_response_plots <-
    list(x = vars,
         y = c("Language skills",
               "Verbal IQ",
               "Skill sum",
               "Normalized skill sum",
               "Class size"),
         v = school_response_test$plot_cap,
         z = c("test score",
               "points",
               "points",
               "Z-scores",
               "Class size")) %>%
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
           theme(legend.position = "none") +
           labs(title = y,
                subtitle = v,
                y = z,
                x = "Socioeconomic status")) %>%
    set_names(vars)

# Testing for differences in IQ and language skills in each class -------

  class_data <-
    split(school_data, school_data$class)

  ## testing for differences with Mann-Whitney U test
  ## separately for each class

  class_test <- class_data %>%
    map(~safely(f_wilcox_test)(.x[, vars],
                       .x[["SES_class"]],
                       as_data_frame = TRUE,
                       adj_method = "BH")) %>%
    map(~.x$result) %>%
    compact %>%
    map(identify_significant,
        label_variable = "variable",
        p_variable = "p_adjusted",
        regulation_variable = "biserial_r",
        return_data = TRUE) %>%
    map(as_tibble)


  ## effect sizes and significance in single classes

  estimate_plot <- class_test %>%
    plot_common_estimates(label_variable = "variable",
                          regulation_variable = "biserial_r",
                          regulation_status = ".regulation",
                          plot_title = "Effect of SES in single classes",
                          x_lab = "Biserial r") +
    geom_vline(xintercept = 0,
               linetype = "dashed")

  ## numbers of classes with significant effects and average effect sizes

  effect_plot <- class_test %>%
    plot_common_effect(label_variable = "variable",
                       regulation_variable = "biserial_r",
                       regulation_status = ".regulation",
                       plot_title = "Effect of SES in single classes",
                       x_lab = "effect size",
                       y_lab = "Number of classes",
                       plot_subtitle = paste("classes: n =", length(class_test)),
                       cutoff = 10,
                       show_whiskers = FALSE)

# Visualization: multi-class heat map and before-after plots -------

  ## classification of the variables and common heat map

  var_classification <- class_data %>%
    map(classify,
        variables = vars[-5],
        split_fct = "SES_class")

  var_classification <- var_classification %>%
    map(~.x$classification)

  cm_heat_plot <- class_data %>%
    common_heat_map(variables = vars[-5],
                    split_fct = "SES_class",
                    normalize = TRUE,
                    norm_center = "mean",
                    norm_dispersion = "sem",
                    variable_classification = var_classification,
                    midpoint = 0,
                    x_lab = "school class")

  ## before - after plots

  common_before_after <- class_data %>%
    plot_common_change(variables = vars[-5],
                       split_fct = "SES_class",
                       dodge = 0,
                       show_whiskers = TRUE,
                       txt_color = "black",
                       txt_position = "right",
                       palette = c("steelblue", "orangered3"),
                       normalize = FALSE,
                       y_lab = "points",
                       x_lab = "social status, SES")

# Differences in variables of interest between the SES classes ------

  class_data %>%
    map(delta,
        split_fct = "SES_class",
        variables = c("lang", "IQ"),
        average_fun = colMedians)

# END -------
