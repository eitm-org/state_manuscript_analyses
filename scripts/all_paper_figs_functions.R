
if(!require(here)){
  install.packages("here")
  library(here)
}

if(!require(knitr)){
  install.packages("knitr")
  library(knitr)
}

if(!require(lqmm)){
  install.packages("lqmm")
  library(lqmm)
}

if(!require(nlme)){
  install.packages("nlme")
  library(nlme)
}

####################################################
#
#        DYNAMIC MODELS VISUALIZATION FUNCTIONS
#
####################################################

wrapper <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}


#substitute x_fix and x_rand with the actual names of these variables
replace_vars <- function(summ_formula) {
  summ_formula <- gsub("x_fix", xvar_fixed, summ_formula)
  summ_formula <- gsub("x_rand", xvar_rand, summ_formula)
  summ_formula <- gsub("y", yvar, summ_formula)
  return(summ_formula)
}

plot_memod <- function(df, char_input) {
  explainer <- var_defs[var_defs$column_me == yvar, "definition"]
  
  mod_call <- paste(yvar, "~", xvar_fixed, "|", xvar_rand)
  
  lme_plot <- ggplot(aes(x = x_fix), data = df) +
    geom_smooth(method = "lm", se = FALSE, aes(y = logrand_fit, color = color_col, group = x_rand), linewidth = .5, alpha = .4) +
    geom_point(aes(y = logy, color = color_col), pch = 16, alpha = .5) +
    labs(caption = wrapper(explainer, width = 45))
  
  test_df <- df %>%
    #for the test, you just want one row for each distinct slope that was made for a draw_id(x_rand)
    distinct(x_rand, rand_m, color_col) %>%
    mutate(m_sign = case_when(rand_m < 0 ~ "negative",
                              TRUE ~ "positive")) 
  
  if (typeof(df[["color_col"]])  == "character" & length(unique(df$color_col)) > 10) {
    lme_plot <- lme_plot + 
      scale_color_viridis_d(guide = "none", end = .9)
    return_list <- list()
  } else if (typeof(df[["color_col"]])  == "character" & length(unique(df$color_col)) <= 10) {
    lme_plot <- lme_plot + 
      scale_color_viridis_d(end = .8) +
      guides(color = guide_legend(title = char_input))
    #compare coefficients between categories
    #stat test
    aov_res <- aov(rand_m ~ color_col, data = test_df)
    aov_summ <- summary(aov_res)
    aov_pval <- aov_summ[[1]]$`Pr(>F)`[[1]]
    #positive/negative in diff groups
    df_summ <- test_df %>%
      group_by(color_col, m_sign) %>%
      mutate(n_m_sign = n()) %>%
      group_by(color_col) %>%
      mutate(n_all = n()) %>%
      filter(m_sign == "positive") %>%
      group_by(color_col, n_all) %>%
      reframe(pos_m_prop = n_m_sign/n_all*100) %>%
      distinct()
    pos_p_prop_sd <- sd(df_summ$pos_m_prop)
    return_list <- list("m_aov_pval" = aov_pval,
                        "m_pos_p_prop_sd" = pos_p_prop_sd)
  } else {
    #color_col is a number (continuous)
    lme_plot <- lme_plot + 
      scale_color_viridis_c(end = .8) +
      guides(color = guide_legend(title = char_input))
    #compare coefficients between categories
    #stat test
    linmod <- lm(rand_m ~ color_col, data = test_df)
    lm_summ <- summary(linmod)
    lm_pval <- lm_summ$coefficients[2,4]
    lm_slope <- lm_summ$coefficients[2,1]
    #not doing positive/negative in diff groups because color_col is continuous
    return_list <- list("m_lm_pval" = lm_pval,
                        "m_lm_slope" = lm_slope)
  }
  
  if (slope < 0.001) {
    graph_slope <- format(slope, scientific = TRUE)
  } else {
    graph_slope <- round(slope, 3)
  }
  
  lme_plot <- lme_plot +
    theme_bw() +
    geom_smooth(aes(y = logfixed_fit), se = FALSE, method = "lm", color = 1, linewidth = 1) +
    xlab(xvar_fixed) +
    ylab(paste("log(", yvar, ")", sep = "")) +
    ggtitle(mod_call) +
    annotate("text", x = min(df$logy), y = max(df$logy), hjust = 0, vjust = 1, label = paste("Pval:", x_pval, "\nSlope:", graph_slope), )
  lme_plot <- cowplot::plot_grid(lme_plot, gridExtra::tableGrob(lme_df, theme = ttheme_minimal()))
  
  return_list <- append(return_list, list("lme_plot" = lme_plot))
  
  return(return_list)
}

get_lmer_vars <- function(lme1) {
  #extract variables from the formula
  lme_form <- lme1$call$fixed
  #get the dependent variable
  #double arrows set the variable globally so I don't have to pass them to the inner function
  yvar <<- all.vars(lme_form)[1]
  #get the fixed independent variable
  xvar_fixed <<- labels(terms(lme_form))
  #get the random independent variable
  #extract random variant call
  rand1 <- as.character(labels(terms(lme1$call$random)))
  #extract the name of the random variable
  xvar_rand <<- substr(rand1, unlist(gregexpr('\\|', rand1))[1] + 2, nchar(rand1))
  
  #was this feature listed as important in laurie's ML analysis?
  imp_feat <- case_when(yvar %in% imp_feats$Feature ~ 1,
                        TRUE ~ 0)
  
  #get output from lmer
  lme_tab <- sjPlot::tab_model(lme1, p.style = "scientific", digits.p = 4)
  lme_df <<- sjtable2df::mtab2df(mtab = lme_tab, n_models = 1, output = "data.frame")
  slope <<- lme1$coefficients$fixed[[2]]
  # as.numeric(lme_df[2, "Estimates"])
  x_pval <<- as.numeric(lme_df[2, "p"])
  marg_cond <- lme_df[11, "Estimates"]
  cond_r2 <- as.numeric(substr(marg_cond, unlist(gregexpr('\\/', marg_cond))[1] + 2, nchar(marg_cond)))
  marg_r2 <- as.numeric(substr(marg_cond, 1, unlist(gregexpr('\\/', marg_cond))[1] - 2))
  # lme_tab_output <<- tableGrob(lme_df, theme = ttheme_minimal())
  #grab the dataframe from the model
  #predict y values based on the model before renaming the variables
  rand_coef <- as.data.frame(lme1$coefficients$random)
  df <- lme1$data
  
  df <- df[,c(yvar, xvar_fixed, xvar_rand, "age", "age_high", "cohort_clinical")] %>%
    mutate(age_high = recode(age_high, 
                             "True" = "age_high",
                             "False" = "age_low"))
  
  names(df) <- c("y", "x_fix", "x_rand", "age", "age_high", "cohort_clinical")
  
  df <- df %>%
    mutate(fixed_b = lme1$coefficients$fixed[[1]],
           fixed_m = lme1$coefficients$fixed[[2]],
           rand_b = rand_coef[x_rand, "draw_id..Intercept."] + fixed_b,
           rand_m = rand_coef[x_rand, "draw_id.draw_month"] + fixed_m,
           rand_fit = rand_b + rand_m*x_fix,
           fixed_fit = fixed_b + fixed_m*x_fix,
           logy = log(y),
           logy = case_when(logy == -Inf ~ min(logy[logy != -Inf]),
                            TRUE ~ logy),
           logrand_fit = log(rand_fit),
           logrand_fit = case_when(is.nan(logrand_fit) ~ min(logrand_fit[!is.nan(logrand_fit)]),
                                   TRUE ~ logrand_fit),
           logfixed_fit = log(fixed_fit),
           logfixed_fit = case_when(logfixed_fit == -Inf ~ min(logfixed_fit[logfixed_fit != -Inf]),
                                    TRUE ~ logfixed_fit))
  #make plot with a line of a different color for each draw_id
  tempdf <- df %>%
    mutate(color_col = x_rand)
  rand_output <- plot_memod(tempdf, "")
  #see models for age, age_high, cohort_clinical
  #age
  tempdf <- df %>%
    mutate(color_col = age)
  age_output <- plot_memod(tempdf, "Age")
  #age_high
  tempdf <- df %>%
    mutate(color_col = age_high)
  ageh_output <- plot_memod(tempdf, "Age, Categorical")
  #cohort_clinical
  tempdf <- df %>%
    mutate(color_col = cohort_clinical)
  clin_output <- plot_memod(tempdf, "Clinical Cohort")
  rm(tempdf)
  
  return_list <- list("yvar" = yvar,
                      "lme_plot" = rand_output$lme_plot,
                      "age_plot" = age_output$lme_plot,
                      "age_m_lm_pval" = age_output$m_lm_pval,
                      "age_m_lm_slope" = age_output$m_lm_slope,
                      "agehigh_plot" = ageh_output$lme_plot,
                      "agehigh_m_aov_pval" = ageh_output$m_aov_pval,
                      "agehigh_m_posneg_comp_sd" = ageh_output$m_pos_p_prop_sd, 
                      "clin_cohort_plot" = clin_output$lme_plot,
                      "clin_cohort_m_aov_pval" = clin_output$m_aov_pval,
                      "clin_cohort_m_posneg_comp_sd" = clin_output$m_pos_p_prop_sd,
                      "slope" = slope,
                      "x_pval" = x_pval,
                      "cond_r2" = cond_r2,
                      "marg_r2" = marg_r2, 
                      "laurie_shortlisted" = imp_feat)
  return(return_list)
}

####################################################
#
#           CONTRIVED SAMPLES FUNCTIONS
#
####################################################
#input: name of the sbsXX_ratio column you're trying to graph
  #things you should have defined in the enviroment:
    #df with columns named "BT474_percent", "replicate", "snv_count", + your input column name
    #gt_lines-- a named list of numeric "ground truth lines", named for the SBS_ratio in the format "SBS[[:digit:]]+"
#output: plot of ratio by sample purity for that sbsXX_ratio
#make sure u use this in a results = "asis" chunk!
plot_sb_ratio <- function(colname, write_output = TRUE) {
  #get data-- metadata + column you're interested in
  df1 <- df[,c("BT474_percent", "replicate", "snv_count", colname)] %>%
    filter(BT474_percent < 100)
  #make snv count table
  snv_table <- df1 %>%
    arrange(BT474_percent)
  #write snv_table to csv
  #make an output_data folder for the coefficient csvs if there isn't one already
  output_path <- file.path("output_data", "fig1", "recovery_figs")

  #rename columns
  names(df1) <- c("BT474_percent", "replicate", "snv_count", "yvar")
  #get SBS name to look for it in the list of ground truth lines (gt_lines)
  sbs_name <- str_remove(colname, "_ratio")
  #make the column name look nice to use it as an axis label in the plot
  nice_colname <- str_replace_all(str_replace_all(str_to_title(str_replace_all(colname, "_", " ")), "Sbs", "SBS"), "Ratio", "Proportion")
  #make plot
  cat("\n")
  cat("\n")
  cat(paste("##", colname))
  cat("\n")
  plot1 <- ggplot(data = df1, aes(y = yvar, x = BT474_percent)) +
    geom_smooth(method = "lm", color = "black", se = FALSE, size = 1) +
    geom_jitter(alpha = .5, size = 3, width = 0.16) +
    theme_classic() +
    ylab(nice_colname) +
    scale_x_continuous(breaks = unique(df1$BT474_percent), name = "Percentage Exogenous DNA", limits = c(0, max(df1$BT474_percent) + .2))
    # ggtitle(paste(nice_colname, "by Sample Purity"))
  print(plot1)
  #add ground truth line
  # if (sbs_name %in% names(gt_lines)) {
  #   gt_lines_df <- as.data.frame(gt_lines)
  #   gt_lines_df["sbs"] <- rownames(gt_lines_df)
  #   gt_lines_df <- gt_lines_df  %>%
  #     mutate(gt = "Ground Truth") %>%
  #     filter(sbs == sbs_name)
  #   plot1 <- plot1 +
  #     geom_hline(data = gt_lines_df, aes(yintercept = gt_lines, color = gt),
  #                linetype = 'dashed',
  #                # color = 'grey',
  #                linewidth = 0.5) +
  #     labs(color = "") +
  #     scale_colour_grey(start = 0.5)
  #     # geom_label(aes(x = 2.5, y =  gt_lines[[sbs_name]], label = "Ground Truth"))
  # }
  
  #write to csv
  if (write_output == TRUE) {
    if(!dir.exists(output_path)) {
      if (!dir.exists(file.path("output_data", "fig1"))) {
        if (!dir.exists(file.path("output_data"))) {
          dir.create(file.path(new_wd, "output_data"))
        }
        dir.create(file.path(new_wd, "output_data", "fig1"))
      }
      dir.create(file.path(new_wd, output_path))
    }
    write_csv(snv_table, file.path(output_path, paste0("snv_table_", colname, ".csv")))
    write_csv(coef_df, file.path(output_path, paste0("coefficient_table_", colname, ".csv")))
    
  }
  #make tables to output
  #make a litle table for model output
  lmod <- lm(yvar ~ BT474_percent, data = df1)
  lmod <- summary(lmod)
  coef_df <- as.data.frame(lmod$coefficients)
  coef_df <- as.data.frame(lmod$coefficients)
  coef_df[["Terms"]] <- row.names(coef_df)
  row.names(coef_df) <- NULL
  coef_df <- coef_df %>% relocate(Terms)
  #make rsq table
  rsq_df <- as.data.frame(c(lmod$r.squared))
  names(rsq_df) <- c("Rsq")
  #print output
  #make sure u use this in a results = "asis" chunk!
  cat("\n")
  pandoc.table(coef_df, caption = "Model Coefficients")
  cat("\n")
  pandoc.table(rsq_df)
  cat("\n")
  ggsave(file.path('plots/fig1', paste0(colname, ".png")), plot1)
  cat("\n\n\\pagebreak\n")
  
  return(plot1)
  
}

####################################################
#
#           MUT-SIGS FUNCTIONS
#
####################################################

find_sbs_others <- function(df, threshold, bound) {
  columns <- colnames(snv_ratios)[-1:-2]
  other <- c()
  
  for (c in columns) {
    col = df[[c]]
    if (sum(col < (sum(col, na.rm=TRUE) * threshold), na.rm=TRUE) > bound * length(col)) {
      other <- append(other, c)
    }
  }
  
  return(other)
}

################################################################
#
#           LQMM MUTSIG vs AGE FUNCTIONS 
#
####################################################

## Get Regression Summary ---------
#' 
#' @param data dataframe for regression
#' @param cols column names for regression
#' @param func_name the regression function
#' @param formula_str independent variable (e.g. '~ age')
#' @param random_formula random formula for mixed effect
#' @param outlier Boolean flag whether to include outliers
#' @param return_type either "summary" or "list": "summary" returns 
#'                    a table of regression stats; "list" returns a list
#'                    of regression objects
#' @examples get_reg_summary(data = agg_snv_global, cols = all_col, 
#'           func_name = "lme", formula_str = "~ age", random_formula <- "~ 1 | draw_id")
#'           
get_reg_summary <- function(data, cols, func_name, formula_str, 
                            random_formula = "~ 1 | draw_id", outlier=TRUE, 
                            return_type = "summary") {
  func <- match.fun(func_name)
  summary_df <- NULL
  table_initialized <- FALSE
  fit_list <- list()
  
  for(var in cols) {
    tryCatch({
      temp_data = data
      formula <- as.formula(paste(var, formula_str))
      # remove outliers
      if (!outlier) {
        Q1 <- quantile(temp_data[[var]], 0.25)
        Q3 <- quantile(temp_data[[var]], 0.75)
        IQR <- Q3 - Q1
        temp_data <- subset(temp_data, temp_data[[var]] > (Q1 - 1.5*IQR) &
                              temp_data[[var]] < (Q3 + 1.5*IQR))
      }
      if(func_name == 'lqmm'){
        fit <- do.call(func, list(formula, random = as.formula(random_formula), 
                                  data = temp_data, group = 'draw_id', tau = 0.5))
      } else if(!is.null(random_formula)){
        fit <- do.call(func, list(formula, random = as.formula(random_formula), 
                                  data = temp_data))
      } else {
        fit <- do.call(func, list(formula, data = temp_data))
      }
      
      fit_list[[var]] <- fit
      # Extract summary statistics
      if(!is.null(random_formula)){
        summary_stats <- summary(fit)$tTable
      } else {
        summary_stats <- summary(fit)$coefficients
      }
      summary_row <- as.data.frame(summary_stats)
      summary_row$Variable <- var  # Add the variable name as a column
      summary_row$Type <- rownames(summary_row)
      row.names(summary_row) <- NULL
      if(!table_initialized) {
        summary_df <- summary_row
        table_initialized <- TRUE
      } else {
        summary_df <- rbind(summary_df, summary_row)
      }
    }, error = function(e) {
      # Skip to the next iteration and report error message
      cat("An error occurred with variable:", var, "Error Message:", e$message, "\n")
    })
  }
  
  if(return_type=="summary"){
    return(summary_df)
  } else {
    return(fit_list)
  }
}

## Anova Tests ---------
#' 
#' @param data dataframe for ANOVA test
#' @param cols column names for ANOVA test
#' @param formula_str1 independent variable (e.g. '~ age')
#' @param formula_str2 independent variable (e.g. '~ age + clinical_cohort')
#' @param random_formula random formula for mixed effect
#' @examples get_anova_summary(data = agg_snv_lme, cols = all_col, formula_str1 = "~age",
#'                             formula_str2 = "~ age + cohort")
get_anova_summary <- function(data, cols, formula_str1, formula_str2,
                              random_formula = "~ 1 | draw_id") {
  table_initialized <- FALSE
  anova_summary <- NULL
  
  for(var in cols) {
    tryCatch({
      form <- as.formula(paste(var, formula_str1))
      form_cc <- as.formula(paste(var, formula_str2))
      
      fit_lme_corr <- do.call(lme, list(fixed = form, random = as.formula(random_formula), data = data, method = "ML"))
      fit_lme_corr_1 <- do.call(lme, list(fixed = form_cc, random = as.formula(random_formula), data = data, method = "ML"))
      
      summary_stats <- anova(fit_lme_corr, fit_lme_corr_1)
      summary_row <- as.data.frame(summary_stats)
      
      summary_row$Variable <- var  # Add the variable name as a column
      row.names(summary_row) <- NULL
      
      if(!table_initialized) {
        anova_summary <- summary_row
        table_initialized <- TRUE
      } else {
        anova_summary <- rbind(anova_summary, summary_row)
      }
    }, error = function(e) {
      cat("An error occurred with variable:", var, "Error Message:", e$message, "\n")
    })
  }
  return(anova_summary)
}

read_agg_data <- function(work_dir) {
  fig6_snvs <- file.path(new_wd, "input_data", "fig6", "joined_snv_clinical.csv")
  agg_snv <- read.csv(file.path(fig6_snvs))
  agg_snv[is.na(agg_snv)] <- 0 # fill all nans with 0
  
  return(agg_snv)
}

lqmm_for_age <- function(work_dir) {
  
  agg_snv <- read_agg_data(work_dir)
  all_col <- names(agg_snv)
  columns_to_remove <- c("draw_age", "draw_month", "draw_id", "cohort") # remove unrelated columns
  all_col <- setdiff(all_col, columns_to_remove)
  
  lqmm_age_summaries <- get_reg_summary(data = agg_snv, random_formula = '~ 1', cols = all_col,
                                        func_name = "lqmm", formula_str = "~ draw_age")
  
  return(lqmm_age_summaries)
}

anova_for_age <- function(work_dir) {
  
  agg_snv <- read_agg_data(work_dir)
  all_col <- names(agg_snv)
  columns_to_remove <- c("draw_age", "draw_month", "draw_id", "cohort") # remove unrelated columns
  all_col <- setdiff(all_col, columns_to_remove)
  
  anova_summary_age <- get_anova_summary(data = agg_snv, cols = all_col, formula_str1 = "~ draw_age", 
                                         formula_str2 = "~ draw_age + cohort")
  
  return(anova_summary_age)
}

# make table and save to input data
make_table <- function(intercepts, slopes, work_dir) {
  lqmm_df <- read_csv(file.path(work_dir, "input_data", "fig6", "lqmm_age_data.csv"))
  anova_df <- read_csv(file.path(work_dir, "input_data", "fig6", "anova_summary_age.csv"))
  
  ## Non-numerical table
  data_table <- data.frame(
    'Variable' = intercepts$Variable,
    'Slope Value' = paste(slopes$Value, "+/-", slopes$`Std. Error`),
    'Intercept Value' = paste(intercepts$Value, "+/-", intercepts$`Std. Error`)
  )
  
  data_table <- merge(data_table, lqmm_df[c("Variable", "lqmm.p")], by = "Variable")
  data_table <- merge(data_table, anova_df[c("Variable", "anova.p")], by = "Variable")
  
  # get significant values based on lqmm and anova
  data_table <- data_table[data_table$'lqmm.p' < 0.05 & data_table$'anova.p' > 0.05, ]
  # remove mutsig variables
  data_table <- data_table[!grepl('mutsig', data_table$Variable, fixed=TRUE), ]
  
  write_csv(data_table, file.path(work_dir, "input_data", "fig6", "age_lqmm_table.csv"))
  
  ## Numerical table
  data_table <- data.frame(
    'Variable' = intercepts$Variable,
    'Slope Value' = slopes$Value,
    'Slope Error' = slopes$`Std. Error`,
    'Intercept Value' = intercepts$Value,
    'Intercept Error' = intercepts$`Std. Error`
  )
  
  data_table <- merge(data_table, lqmm_df[c("Variable", "lqmm.p")], by = "Variable")
  data_table <- merge(data_table, anova_df[c("Variable", "anova.p")], by = "Variable")
  
  # get significant values based on lqmm and anova
  data_table <- data_table[data_table$'lqmm.p' < 0.05 & data_table$'anova.p' > 0.05, ]
  # remove mutsig variables
  data_table <- data_table[!grepl('mutsig', data_table$Variable, fixed=TRUE), ]
  
  write_csv(data_table, file.path(work_dir, "input_data", "fig6", "age_lqmm_numtable.csv"))
}

plot_lin_mods <- function(colname) {
  for_plot <- dynamic_df_grp_wide
  for_plot["yvar"] <- for_plot[colname]
  
  yvar_slope <- paste0(colname, "_slope")
  yvar_int <- paste0(colname, "_intercept")
  
  mean_slope <- mean(unlist(for_plot[yvar_slope]))
  mean_int <- mean(unlist(for_plot[yvar_int]))
  
  mean_lm_df <- data.frame(x = seq(from = min(for_plot$draw_month), to = max(for_plot$draw_month), length.out = 100)) %>%
    mutate(y = mean_int + mean_slope*x)
  
  nice_colname <- str_remove_all(str_replace_all(str_to_title(str_replace_all(colname, "_", " ")), "Sbs", "SBS"), " Count|Adj3|Ratio")
  nice_colname <- paste(nice_colname, "Activity")
  ymax <- max(for_plot$yvar, na.rm = TRUE) + max(for_plot$yvar, na.rm = TRUE)/5

  fig7b_lin <- ggplot() +
    geom_smooth(data = for_plot, aes(y = yvar, x = draw_month, group = draw_id), method = "lm", fill = NA, linewidth = .4, color = "darkgray") + 
    geom_point(data = for_plot, aes(y = yvar, x = draw_month, group = draw_id), size = 1, alpha = .5) + 
    theme(legend.position = "none") +
    theme_classic() +
    geom_smooth(data = mean_lm_df, aes(y = y, x = x),method = "lm", se = FALSE, linewidth = 1.4, color = "cornflowerblue") +
    ylab(nice_colname) +
    xlab("Draw Month") +
    # ggtitle(paste("Dynamic Modelling--", nice_colname)) 
    ylim(c(0, ymax))
  
  return(fig7b_lin)
}
