library(rstatix)
library(readxl)
library(writexl)
library(car)
library(DHARMa)
library(betareg)
library(cowplot)


Fat_body <- read_excel("data/fat_body_fatty_acids.xlsx", sheet = 1)
Liver <- read_excel("data/liver_fatty_acids.xlsx", sheet = 1)

func <- function(phenotype) {

  name <- paste0(deparse(substitute(phenotype)))  
    
  prop <- phenotype %>% mutate_if(is.numeric, list(~ . / 100))
  
  
  tests <- lapply(prop[,c(3:ncol(prop))], function (x) {betareg(x ~ treatment, data = prop)})
  
  summary <- lapply(tests, function (x) {summary(x)})
  anodev <- lapply(tests, function (x) {Anova(x, type = 3)})
  jointtests <- lapply(tests, function (x) {joint_tests(x)})
  
  
  gg <- lapply(3:ncol(phenotype), function(j) {
    cbind(phenotype[j], Treatment = phenotype$treatment)
  })
  names(gg) <- names(tests)
  
  
  sig <- lapply((jointtests), function(i) {
    as.data.frame(i) %>%
      cbind(., Treatment = c("Linuron_H")) %>%
      mutate(.group = case_when(
        p.value < 0.05 & p.value > 0.01 ~ "*",
        p.value < 0.01 & p.value > 0.001 ~ "**",
        p.value < 0.001 & p.value >= 0 ~ "***",
        p.value > 0.05 ~ ""
      ))})
  
  
  graph <- lapply(gg, function(i) {
    
    plot <- ggplot(i, aes(x = Treatment, y = i[, 1])) + 
      geom_boxplot(width = 0.3,
                  position = position_dodge(width = 0.1), 
                   color = 'black', fill = 'gray', 
                   outlier.shape = NA, alpha = 0.7) +
      geom_jitter(position = position_dodge(width = 0.8), 
                  color = "black", size = 3, alpha = 0.7) +
      labs(title = paste0(names(i)[1]), y = paste0("% of total"), x = NULL) +
      theme_minimal() +  # You can change the theme to your preference
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text(data=
                  #letters.gamma[[names(i)[1]]],
                  sig[[names(i)[1]]],
                aes(x=Treatment, y=((max(i[,1],na.rm = TRUE))+(max(i[,1],na.rm = TRUE)*0.1)),label=.group), hjust=1, position = position_dodge(width = 1), angle=0, size=10)
    
    return(plot)
  })
  
  # model.diag.gamma <- (lapply(tests, function(x){simulateResiduals(fittedModel = x, plot = F)}))
  # for(i in 1:length(model.diag.gamma)){plot(model.diag.gamma[[i]], title = names(model.diag.gamma)[i])}
  
  
  summary_table <- phenotype %>% 
    group_by(treatment)%>%
    get_summary_stats(type = "mean_se")
  
  
  summ_ctrl <- subset(summary_table, summary_table$treatment == "Control")
  
  summ_lin <- subset(summary_table, summary_table$treatment == "Linuron_H")
  
  
  stars_sig_table <- function(tests_list, summ_list,
                                 
                                 digits = 2, 
                                 decimal.mark = ".",
                              show_significance = TRUE){
    
    # check arguments
    stopifnot({
      is.numeric(digits)
      digits >= 0
    })
  
    p_values <- lapply(tests_list, function(p) p$p.value) %>% unlist
    p <- data.frame(p.value = p_values, row.names = names(p_values))
    
    summ_list_mean <- summ_list$mean
    
    # transform correlations to specific character format
    Rformatted = formatC(summ_list_mean, format = 'f', digits = digits, decimal.mark = decimal.mark)
    
    # add significance levels if desired
    if (show_significance) {
      # define notions for significance levels; spacing is important.
      stars <- ifelse(is.na(p), "   ", ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "*  ", "   "))))
      Rformatted = paste0(Rformatted, stars)
    }
    
    new_table <- summ_list %>%
      mutate(mean = Rformatted)
  
    return(new_table)
  }
  
  summ_lin_stars <- stars_sig_table(jointtests, summ_lin)
  
  sig_tables <- list(Control = summ_ctrl, Linuron = summ_lin_stars)
  
  combined_df <- do.call(cbind, sig_tables)
  
  df_final <- combined_df[,c(2:5,8:10)]
  
  colnames(df_final) <- c("Variable", "Control n", "Control Mean", "Control SEM", 
                           "Linuron n", "Linuron Mean", "Linuron SEM")
  
  df_final$p.value <- lapply(jointtests, function(p) p$p.value) %>% unlist
  
  vars_to_format <- c("Control Mean", "Control SEM", "Linuron Mean", "Linuron SEM")
  
  df_final <- df_final %>%
    mutate(across(all_of(vars_to_format), ~ if (is.numeric(.)) {
      formatC(.x, format = 'f', digits = 2, decimal.mark = ".")
    } else {
      .x
    })) %>%
    mutate(across(all_of("p.value"), ~ if (is.numeric(.)) {
      formatC(.x, format = 'f', digits = 3, decimal.mark = ".")
    } else {
      .x
    }))
  
  write_xlsx(sig_tables, path = paste0("./tables/", deparse(substitute(phenotype)), "_fatty acids_statistics_table.xlsx"))
  
  assign(paste0(deparse(substitute(phenotype)), "_tests") , tests, envir = .GlobalEnv)
  assign(paste0(deparse(substitute(phenotype)), "_summary") , summary, envir = .GlobalEnv)
  assign(paste0(deparse(substitute(phenotype)), "_anodev") , anodev, envir = .GlobalEnv)
  assign(paste0(deparse(substitute(phenotype)), "_jointtests") , jointtests, envir = .GlobalEnv)
  assign(paste0(deparse(substitute(phenotype)), "_graph") , graph, envir = .GlobalEnv)
  assign(paste0(deparse(substitute(phenotype)), "_result_table") ,df_final, envir = .GlobalEnv)

}

func(Fat_body)
func(Liver)

combined_plot <- plot_grid(plotlist = Fat_body_graph)
combined_plot

lapply(Fat_body_graph, print)
