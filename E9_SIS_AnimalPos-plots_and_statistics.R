#' @title E9 Social Stress Experiment Data Analysis and Visualization
#' @description This script performs data analysis and visualization for the E9 Social Stress experiment.
#' - It loads necessary libraries, sets up directories for saving plots and statistical results,
#'   and defines paths for working and saving directories.
#' - The script then reads lists of control and susceptible animals from CSV files and initializes
#'   tibble objects for storing processed data for social proximity and coil crossing.
#' - It creates tibble objects with consecutive phase rows for further analysis, useful for linear mixed-effects (LME) statistics
#'   and visualizing consecutive phase data.
#' - The script then reads saved tables into a tibble for cage and mice entropy and combines all tibbles into a list.
#' - It creates linear mixed effects models and plots for consecutive phase tibbles and cage change tibbles.
#' - The script also computes rank values for cage change tibbles and generates plots and statistics for cage change tibbles.
#' - The results are saved in the specified directories.
#' - Flags are set to determine whether to save plots and statistics. Directory paths for 
#'   the working environment and saving locations are defined. Modify the paths as needed 
#'   based on the environment (e.g., local machine or server).

message("## create cage change mean plots and statistics ##")

#dynamic plot and statistic for all CC tibbles 

#definitions for dynamic processing
values <- c('SocialProx', 'CoilCrossing', 'MiceEntropy', 'CageEntropy')
phases <- c("active", "inactive")

# Initialisieren einer Liste für die Ergebnisse
#avg hors per cc
avg_plots <- list()
#avg rank in total
rank_plots <- list()

##statistics on rank plots(tobis code, plots again)
allTestResults <- list()
allPlots <- list()
allPosthocResults <- list()

for(i in seq_along(cc_tibbles)){
  value <- values[i]
  color_by <- ifelse(value == 'CageEntropy', 'System', 'Group')
  ifelse(value == 'CageEntropy', color_palette <- c("forestgreen", "goldenrod", "steelblue", "firebrick", "darkorchid"), color_palette <- c("con" = "#1e3791", "res" = "#8aacdb", "sus" = "#f49620"))
  
  message(value)
  
  for(j in seq_along(phases)){
    cc_data <- cc_tibbles[[value]]
    #print(cc_data)
    phase <- phases[j]
    message(phase)
    
    ## PLOTS  ##
    ## AVG PLOT
    avg_phase <- ifelse(phase == "active", 'Avg_act', 'Avg_inact')
    
    #style ribbon
    avg_cc_p <- ggplot(data = cc_data, aes(x = CageChange, y = !!sym(avg_phase), color = !!sym(color_by), group=!!sym(color_by))) + 
      stat_summary(aes(group = !!sym(color_by)), fun = median, geom = "line", linewidth = 1) +
      stat_summary(aes(fill = !!sym(color_by)),
                   fun.min = function(z){quantile(z, 0.25)},
                   fun.max = function(z){quantile(z, 0.75)},
                   fun = median,
                   geom = "ribbon",                              
                   alpha = 0.2, color = NA) +
      labs(title = paste(value, " Cage Change Plot"), subtitle = phase, "phases") +
      scale_fill_manual(values = color_palette) +
      scale_color_manual(values = color_palette) +
      scale_y_continuous(paste(value)) +
      scale_x_discrete(name = "Cage Change")+ 
      facet_grid(Sex~.) +
      theme_minimal(base_size = 14) +  # Increase base font size for readability
      theme(
            panel.grid.major = element_blank(),  # Remove major gridlines
            panel.grid.minor = element_blank(),  # Remove minor gridlines
            plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # Center title
            plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),  # Italic for subtitle
            axis.text.x = element_text(size = 12),  # Adjust text sizes for clarity
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 14, face = "bold"),
            axis.ticks.x = element_line(size = 0.5),
            legend.position = "top",  # Place legend inside the plot area
            panel.background = element_blank())  # Clean panel background
    
    #save plot in plot list
    if(j == 1) avg_plots[[2 * i - 1]] <- avg_cc_p
    else avg_plots[[2 * i]] <- avg_cc_p
    #save plot in file
    if(save_plots) ggsave(filename = paste0(working_directory, "/plots", "/cc_plots/", "avg-cc-",value,"-",phase,".svg"), plot = avg_cc_p, width = 5, height = 5)
    
    # LME-AVG 
    message(paste("cc lme for", value ,", ", avg_phase))
    
    #create dynamic formula
    if(value == 'CageEntropy') {
      formula <- as.formula(paste(avg_phase, "~ System * Sex + CageChange + (1 | Batch)")) ####maybe * CageChange?
    } else {
      formula <- as.formula(paste(avg_phase, "~ Group * Sex + CageChange + (1 | AnimalID)"))
    }
    message(formula)
    
    model <- lmerTest::lmer(formula, data = cc_data)
    #print(summary(model))
    
    # save mixed model results to a file
    mixed_model_results_file <- paste0(working_directory, "/lme/cc-mixed_model_results/", paste0("avg-lme-",value,"-",phase,".txt"))
    model_summary <- capture.output(summary(model))
    if(save_statistics) writeLines(model_summary, con = mixed_model_results_file)
    if(value != 'CageEntropy'){

      ## AVG RANK PLOT
      #define the current rank-value based o the current phase
      rank_phase <- ifelse(phase == "active", 'SystemRank_act', 'SystemRank_inact')
      
      #group data to summarize all Cage changes for each individual
      avg_data <- cc_data %>%
        group_by(Batch, Sex, Group, AnimalID) %>%
        summarise(avg_rank = mean(!!sym(rank_phase)))
      
      #exclude con group from visualisation if  socialProx (here not valuable)
      if(value == 'SocialProx') avg_data <- avg_data %>% filter(Group != "con")
      
      avg_rank_p <- ggplot(data = avg_data, aes(x = Group, y = avg_rank, color = Group, shape = Batch)) +
      geom_jitter(size = 4, alpha = 0.7, width = 0.2, height = 0) +
        scale_shape_manual(values = c(3, 16, 17, 15, 5, 20)) + 
        scale_x_discrete("Animal ID") +
        scale_y_discrete("Rank", limits = c("1", "2", "3", "4")) +
        labs(title = paste("Median rank of each individual"), subtitle = paste(value, "-", phase)) +
        scale_color_manual(values = c("con" = "#1e3791", "res" = "#8aacdb", "sus" = "#f49620")) +
        stat_summary(
          fun.min = function(z) {quantile(z, 0.25)},
          fun.max = function(z) {quantile(z, 0.75)},
          fun = mean,
          color = "black",
          size = 0.8,
          shape = 16) +
        facet_grid(Sex ~ .) +
        theme_minimal(base_size = 14) +  # Increase base font size for readability
        theme(
              panel.grid.major = element_blank(),  # Remove major gridlines
              panel.grid.minor = element_blank(),  # Remove minor gridlines
              plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  # Center title
              plot.subtitle = element_text(hjust = 0.5, size = 14, face = "italic"),  # Italic for subtitle
              axis.text.x = element_text(size = 12),  # Adjust text sizes for clarity
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 14, face = "bold"),
              axis.ticks.x = element_line(size = 0.5),
              legend.position = "top",  # Place legend inside the plot area
              panel.background = element_blank())  # Clean panel background
      
      #save plot in plot list
      if(j == 1) rank_plots[[2 * i - 1]] <- avg_rank_p
      else rank_plots[[2 * i]] <- avg_rank_p
      #save plot in file
      if(save_plots) ggsave(filename = paste0(working_directory, "/plots", "/avg_rank_plots/", "avg_ranks-",value,"-", phase, ".svg"), plot = avg_rank_p, width = 5, height = 5)
      
      #' STATISTICS ON RANK AVG PLOTS
      #' Perform statistical tests and generate plots for rank averages
      #' This function performs statistical tests and generates plots for rank averages based on the specified data and parameters.
      #' @param data The data frame containing the rank averages
      #' @param value The value column name
      #' @param rank The rank column name
      #' @param phase The phase (active or inactive)
      #' @param sex 
      #' @return A list containing the columns testResults, plot, and posthocResults

      for(sex in c("female", "male")){
        
        message("statistics on rank avgs")
        print(sex)
        result <- testAndPlotVariable(avg_data, value, 'avg_rank', phase, sex)#avg_data or cc_data???
        # add result-list(containing the columns testResults, plot, posthocResults) to other fitting list
        if (!is.null(result)) {
          # posthocResults is always NULL for the Wilcoxon test
          if (is.null(result$posthocResults)) {
            allTestResults <- c(allTestResults, list(result$testResults))
          } else {
            allTestResults <- c(allTestResults, list(result$testResults))
            allPosthocResults <- c(allPosthocResults, list(result$posthocResults))
          }
          #add the plot
          allPlots <- c(allPlots, list(result$plot))
        }
      }
    }
  }
}

## PRINT PLOTS IN R ##
# Print all plots individually
for (i in seq_along(avg_plots)) {
  print(avg_plots[[i]])
}
# Arrange plots into grids
gridExtra::grid.arrange(grobs = avg_plots, ncol = 2)   # Average plots
gridExtra::grid.arrange(grobs = rank_plots, ncol = 2)  # Rank plots
gridExtra::grid.arrange(grobs = allPlots, ncol = 4)    # All plots

###################### saving the results #############################################################

# Convert the list of test results to a data frame
allTestResultsDf <- bind_rows(allTestResults)
# Save the test results data frame to a CSV file
# Create the statistics directory if it doesn't exist
statistics_directory <- paste0(working_directory, "/statistics")
if (!dir.exists(statistics_directory)) {
  dir.create(statistics_directory, recursive = TRUE)
}
write.csv(allTestResultsDf, file = paste0(working_directory,"/statistics/test_results.csv"), row.names = FALSE)
# Save the post hoc results to a CSV file
if (!is.null(allPosthocResults) && length(allPosthocResults) > 0) {
  allPosthocResultsDf <- bind_rows(allPosthocResults)
  write.csv(allPosthocResultsDf, file = paste0(working_directory,"/statistics/posthoc_results.csv"), row.names = FALSE)
}