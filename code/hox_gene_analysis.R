library(phytools)
library(rotl)
library(OpenTreeChronograms)
library(ggtree)
library(geiger)
library(nlme)
library(ggpubr)
library(rr2)
library(tidyverse)
library(RRphylo)
library(caper)
library(parallel)
library(PKNCA)
library(phyr)
library(sysfonts)
library(showtext)
library(grid)
library(cowplot)

#### Hox Genes ####
coleoptera_hox = read.table("../data/hox/Coleoptera_HOX_cluster_processed.txt")
diptera_hox = read.table("../data/hox/Diptera_HOX_cluster_processed.txt"); diptera_hox = subset(diptera_hox, V4 != "abd-A")
hymenoptera_hox = read.table("../data/hox/Hymenoptera_HOX_cluster_processed.txt")
lepidoptera_hox = read.table("../data/hox/Lepidoptera_HOX_cluster_processed.txt")
insecta_hox = rbind(coleoptera_hox, diptera_hox, hymenoptera_hox, lepidoptera_hox)

#### N50 Scores to subset HOX data ####
coleoptera_n50 = read.table("../data/hox/Coleoptera_N50_Score.csv", sep = ",")
hymenoptera_n50 = read.table("../data/hox/Hymenoptera_N50_Score.csv", sep = ",")
diptera_n50 = read.table("../data/hox/Diptera_N50_Score.csv", sep = ",")
lepidoptera_n50 = read.table("../data/hox/Lepidoptera_N50_Score.csv", sep = ",")
insecta_n50 = rbind(coleoptera_n50, diptera_n50, hymenoptera_n50, lepidoptera_n50)

#### Merge N50 and Hox data ####
for (df_hox_name in ls(pattern = "_hox$")){
  df_hox = get(df_hox_name)
  colnames(df_hox) = c("Region", "Start", "Stop", "Hox", "Direction", "Score","Accession") # Rename columns
  df_hox = df_hox %>% mutate(Length = Stop - Start) %>% filter(Length >= 150 & Length <= 210 & Score > 70)
  
  df_n50_name = sub("_hox$", "_n50", df_hox_name)
  df_n50 = get(df_n50_name)
  colnames(df_n50) = c("N50", "Accession")
  df_n50 = subset(df_n50, N50 > 10000000)
  df_hox <- left_join(df_hox, df_n50, by = "Accession")
  df_hox = na.omit(df_hox)
  assign(df_hox_name, df_hox)
}
data_frame_list = list(coleoptera_hox, diptera_hox, hymenoptera_hox, lepidoptera_hox, insecta_hox)
mass = read.csv("../data/mass/Insect_Masses_and_Genomes.csv")[,-1]
phylo = read.tree("../data/phylogeny/calibrated_tree.nwk")


#### Hox Gene Analysis #### 
# This code outputs 
# The sections for plotting functions correctly for insecta only, so DO NOT SOURCE THIS DOCUMENT.
# After this, check for which models were significant.
# Be careful that delta AIC < 2 indicates conflicting models, and often the most parsimonious or one with most-normal residual distribution is best.

data_frame_names = c("Coleoptera", "Diptera", "Hymenoptera", "Lepidoptera", "Insecta")
for (i in seq_along(data_frame_names)){
  
  cat("\n #### THE FOLLOWING ANALYSIS IS ON:", data_frame_names[i], "\n")
  
  data_frame_hox = as.data.frame(data_frame_list[i])
  data_frame_mass = mass # keep this
  phylo_tree = phylo # keep this
  
  table <- data_frame_hox %>% # Transform the provided data into better format
    group_by(Accession, Hox) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = Hox, values_from = Count, values_fill = list(Count = 0)) %>%
    ungroup() %>%
    mutate(hoxsums = rowSums(across(-Accession))) # calculate sum of hox genes 
  colnames(table) = sub("-", "_", colnames(table)) # Fix gene naming issues
  
  table$Accession <- gsub("^GCA_|^GCF_|\\.\\d+$", "", table$Accession) # Remove GCA and GCF, as these changed when pulled from NCBI
  data_frame_mass$Accession <- gsub("^GCA_|^GCF_|\\.\\d+$", "", data_frame_mass$Accession) # Remove GCA and GCF, as these changed when pulled from NCBI
  data_frame = merge(data_frame_mass, table, by = "Accession") # Merge data
  
  phylo_tree$tip.label = sub("_", " ", phylo_tree$tip.label) # Fix tip naming issues
  phylo_tree = keep.tip(phylo_tree, data_frame$Species[data_frame$Species %in% phylo_tree$tip.label]) # Subset tree to species in data
  data_frame = data_frame[data_frame$Species %in% phylo_tree$tip.label, ] # Subset data to species present in tree and data
  data_frame = data_frame %>%  select_if(~ n_distinct(.) > 1)
  hox_genes <- names(data_frame)[10:ncol(data_frame)] # Identify HOX genes (columns from 10 onwards)
  hox_genes = sort(hox_genes)
  
  for (hox_gene in hox_genes){
    formula <- as.formula(paste(hox_gene, "~ log10(Mean_Converted_Value) + (1 | Species)")) # PGLMM Linear
    formula_sq <- as.formula(paste(hox_gene, "~ poly(log10(Mean_Converted_Value), 2, raw = TRUE) + (1 | Species)")) # PGLMM Quadratic
    null_formula <- as.formula(paste(hox_gene, "~ 1 + (1 | Species)")) # PGLMM Null
    
    model_fitting_functions <- list( # As data sometimes has poisson or binomial distributions, let's fit 
      pglmmModel_Poisson = function() suppressMessages(pglmm(formula, data = data_frame, family = "poisson", cov_ranef = list(Species = phylo_tree), REML = FALSE)), # Poisson models
      pglmmModel_Poisson_sq = function() suppressMessages(pglmm(formula_sq, data = data_frame, family = "poisson", cov_ranef = list(Species = phylo_tree), REML = FALSE)),
      pglmmModel_Poisson_Null = function() suppressMessages(pglmm(null_formula, data = data_frame, family = "poisson", cov_ranef = list(Species = phylo_tree), REML = FALSE)),
      pglmmModel_Binom = function() suppressMessages(pglmm(formula, data = data_frame, family = "binomial", cov_ranef = list(Species = phylo_tree), REML = FALSE)), # Binomial models
      pglmmModel_Binom_sq = function() suppressMessages(pglmm(formula_sq, data = data_frame, family = "binomial", cov_ranef = list(Species = phylo_tree), REML = FALSE)),
      pglmmModel_Binom_Null = function() suppressMessages(pglmm(null_formula, data = data_frame, family = "binomial", cov_ranef = list(Species = phylo_tree), REML = FALSE))
    )
    
    models <- list()       # To store the fitted models
    model_AICs <- list()   # To store AIC values
    
    print(paste("Fitting models for:", hox_gene))
    
    for (model_name in names(model_fitting_functions)) {  # Loop through each model fitting function
      fit_function <- model_fitting_functions[[model_name]]
      fit_result <- try(fit_function(), silent = TRUE)  # Try fitting the model and catch errors
      
      if (!inherits(fit_result, "try-error")) { # If there is an error then paste the model with the error - most of the time this will be the one with the incompatible data structure
        models[[model_name]] <- fit_result
        model_AICs[[model_name]] <- fit_result$AIC  # Use AIC from the pglmm object
      } else {
        message(paste("Model fitting failed for", model_name))
      }
    }
    
    if (length(models) == 0) { # Check if any models were successfully fitted
      stop(paste("No models were successfully fitted for", hox_gene))
    }
    
    delta_AICs <- list()  # Calculate delta AIC relative to the null models
    
    # Calculate delta AIC for Poisson models
    if (!is.null(model_AICs$pglmmModel_Poisson) && !is.null(model_AICs$pglmmModel_Poisson_Null)) {
      delta_AICs$pglmmModel_Poisson <- model_AICs$pglmmModel_Poisson - model_AICs$pglmmModel_Poisson_Null
    }
    if (!is.null(model_AICs$pglmmModel_Poisson_sq) && !is.null(model_AICs$pglmmModel_Poisson_Null)) {
      delta_AICs$pglmmModel_Poisson_sq <- model_AICs$pglmmModel_Poisson_sq - model_AICs$pglmmModel_Poisson_Null
    }
    
    # Calculate delta AIC for Binomial models
    if (!is.null(model_AICs$pglmmModel_Binom) && !is.null(model_AICs$pglmmModel_Binom_Null)) {
      delta_AICs$pglmmModel_Binom <- model_AICs$pglmmModel_Binom - model_AICs$pglmmModel_Binom_Null
    }
    if (!is.null(model_AICs$pglmmModel_Binom_sq) && !is.null(model_AICs$pglmmModel_Binom_Null)) {
      delta_AICs$pglmmModel_Binom_sq <- model_AICs$pglmmModel_Binom_sq - model_AICs$pglmmModel_Binom_Null
    }
    
    # Identify the model with the lowest AIC
    best_model_name <- names(which.min(unlist(model_AICs)))
    best_model_AIC <- model_AICs[[best_model_name]]
    second_best_model_name <- names(sort(unlist(model_AICs))[2])
    second_best_model_AIC <- model_AICs[[second_best_model_name]]
    
    # Print the best model and its delta AIC relative to the corresponding null model
    if (grepl("_Null$", best_model_name)) {
      print(paste("Best model for", hox_gene, "is:", best_model_name, "with AIC:", round(best_model_AIC,2), "and AIC relative to second best model:", round(second_best_model_AIC - best_model_AIC, 2), "p-value:", round(models[[best_model_name]]$B.pvalue[1], 3)))
    } else {
      delta_AIC <- delta_AICs[[best_model_name]]
      print(paste("Best model for", hox_gene, "is:", best_model_name, "with AIC:", round(best_model_AIC,2), "and delta AIC relative to null model:", round(delta_AIC, 2), "and delta AIC relative to second best model:", round(second_best_model_AIC - best_model_AIC, 2)))
    }
  } 

}


#### THIS IS FOR PLOTTING INSECTA FIGURES FOUND IN THE DISSERTATION - ONLY RUN THIS BELOW IF YOU HAVE RUN THE ABOVE ON INSECTA - THIS IS THE CASE IF YOU RAN EXACTLY THE CODE ABOVE ####
#### All models are fit best to null, apart from Ro, cad, and ftz (though this is an error, as null is actually better fitting)
#### Need to be removed: ShxB, ShxC, ShxD, Antp. Last 4 are actually only in lepidoptera or diptera, so must be removed

for (hox_gene in hox_genes){
  
  pglmmModel_Poisson_Null = pglmm(hox_gene, ~ 1 + (1 | Species), data = data_frame, family = "poisson", cov_ranef = list(Species = phylo_tree), REML = FALSE) # Fit null model

  # Custom latex font
  font_add(family = "lm", regular = "/usr/share/texlive/texmf-dist/fonts/opentype/public/lm/lmroman10-regular.otf")
  font_add(family = "lm-bold", regular = "/usr/share/texmf/fonts/opentype/public/lm/lmroman10-bold.otf")
  showtext_auto()
  
  intercept <- pglmmModel_Poisson_Null$B[1] # Get model intercept
  if (is.na(pglmmModel_Poisson_Null$B[2])){ # If the slope is na (doesnt exit) then make it 0, otherwise set the slope as the second value
    coef_log10_Mean_Converted_Value = 0
  } else {
    coef_log10_Mean_Converted_Value = pglmmModel_Poisson_Null$B[2]
  }
  
  log10_values <- seq(min(log10(data_frame$Mean_Converted_Value)), max(log10(data_frame$Mean_Converted_Value)), length.out = 100) # sequence of values for plotting the line
  new_data <- data.frame(log10_Mean_Converted_Value = log10_values) # make a data frame from them
  new_data$predicted <- exp(intercept + coef_log10_Mean_Converted_Value * new_data$log10_Mean_Converted_Value) # predict model outputs for this sequence
  
  plot_name <- paste0("plot_pglmm_", hox_gene) # set plot names for multiple plots
  
  if (pglmmModel_Poisson_Null$B.pvalue[1] == 0 || pglmmModel_Poisson_Null$B.pvalue[1] < 0.001){ # set p-values from the model outputs
    p_text = "p < 0.001"
  } else {
    p_text = paste0("p = ", round(pglmmModel_Poisson_Null$B.pvalue[1], 3))
  }
  
  hox_gene_label = hox_gene # fix some hox gene naming issues
  if (hox_gene == "hoxsums"){
    hox_gene_label = "Hox Gene Sum"
  } else if (hox_gene == "abd_A"){
    hox_gene_label = "abd-A"
  }
  
  
  
  assign(plot_name, ggplot() + # make the plot for the specific hox gene
           geom_point(data = data_frame, aes(x = log10(Mean_Converted_Value), y = .data[[hox_gene]]), color = "#DC4D10", alpha = 0.6) +
           geom_line(data = new_data, aes(x = log10_Mean_Converted_Value, y = predicted), size = 1.2, color = "black") +
           labs(x = "Log Species Mass", y = paste0(hox_gene)) +
           theme_pubr() +
           theme(
             legend.position = "right",
             text = element_text(family = "lm", size = 14),  # General text font
             axis.title.x = element_blank(),  # X-axis title font
             axis.title.y = element_blank(),  # Y-axis title font
             axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Rotate y-axis tick label
             plot.title = element_text(family = "lm-bold", size = 16),  # Title font
             axis.text = element_text(family = "lm", size = 14)  # Axis text font
           )+
           annotate("text", x = -Inf, y = Inf, 
                    label = p_text, hjust = -0.05, vjust = 5.5, size = 5, fontface = "italic", family = "lm") +
           annotate("text", x = -Inf, y = Inf, 
                    label = "R^2 == 0", hjust = -0.08, vjust = 3.2, size = 5, family = "lm", parse = TRUE) +
           annotate("text", x = -Inf, y = Inf, 
                    label = "PGLMM Poisson", hjust = -0.03, vjust = 2.75, size = 5, family = "lm") +
           annotate("text", x = -Inf, y = Inf, 
                    label = hox_gene_label, hjust = -0.2, vjust = 1.5, size = 5, fontface = "italic", family = "lm-bold")
         
         
  )
  
  print(get(plot_name))
}

remove(plot_pglmm_ShxA, plot_pglmm_ShxB, plot_pglmm_ShxC, plot_pglmm_ShxD, plot_pglmm_Ftz, plot_pglmm_Antp) # remove the plots that didn't work or shouldn't be here - specified above
plot_pglmm_abd_A = plot_pglmm_abd_A + labs(y = "Abd-A") # Fix abd-A label
plot_pglmm_hoxsums = plot_pglmm_hoxsums + labs(y = "Sum of Hox Genes") + theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) # Fix hoxsums plot name

# custom latex fonts
font_add(family = "lm", regular = "/usr/share/texlive/texmf-dist/fonts/opentype/public/lm/lmroman10-regular.otf")
font_add(family = "lm-bold", regular = "/usr/share/texmf/fonts/opentype/public/lm/lmroman10-bold.otf")
showtext_auto()

# Plot Ro outputs
model <- pglmm(Ro ~ log10(Mean_Converted_Value) + (1 | Species), data = data_frame, family = "poisson", cov_ranef = list(Species = phylo_tree), REML = FALSE) 
summary(model)
intercept <- model$B[1]
coef_log10_Mean_Converted_Value <- model$B[2] 
log10_values <- seq(min(log10(data_frame$Mean_Converted_Value)), max(log10(data_frame$Mean_Converted_Value)), length.out = 100)
new_data <- data.frame(log10_Mean_Converted_Value = log10_values)
new_data$predicted_ro <- exp(intercept + coef_log10_Mean_Converted_Value * new_data$log10_Mean_Converted_Value)
plot_pglmm_Ro <- ggplot() +
  geom_point(data = data_frame, aes(x = log10(Mean_Converted_Value), y = Ro), color = "#DC4D10", alpha = 0.6) +
  geom_line(data = new_data, aes(x = log10_Mean_Converted_Value, y = predicted_ro), size = 1.2, color = "black") +
  labs(x = "Log Species Mass", y = "Ro") +
  theme_pubr() +
  theme(
    legend.position = "right",
    text = element_text(family = "lm", size = 14),  # General text font
    axis.title.x = element_blank(),  # X-axis title font
    axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Rotate y-axis tick label
    axis.title.y = element_blank(),
    plot.title = element_text(family = "lm-bold", size = 16),  # Title font
    axis.text = element_text(family = "lm", size = 14)  # Axis text font
  )+
  annotate("text", x = -Inf, y = Inf, 
           label = "p < 0.001", hjust = -0.1, vjust = 5.5, size = 5, fontface = "italic", family = "lm") +
  annotate("text", x = -Inf, y = Inf, 
           label = "R^2 == 0.08", hjust = -0.1, vjust = 3.2, size = 5, family = "lm", parse = TRUE) +
  annotate("text", x = -Inf, y = Inf, 
           label = "PGLMM Poisson", hjust = -0.06, vjust = 2.75, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = Inf, 
           label = "Ro", hjust = -0.2, vjust = 1.5, size = 5, fontface = "italic", family = "lm-bold")
print(plot_pglmm_Ro)

# Plot cad outputs
model <- pglmm(cad ~ log10(Mean_Converted_Value) + (1 | Species), data = data_frame, family = "poisson", cov_ranef = list(Species = phylo_tree), REML = FALSE)
summary(model)
intercept <- model$B[1]
coef_log10_Mean_Converted_Value <- model$B[2]
log10_values <- seq(min(log10(data_frame$Mean_Converted_Value)), max(log10(data_frame$Mean_Converted_Value)), length.out = 100)
new_data <- data.frame(log10_Mean_Converted_Value = log10_values)
new_data$predicted_cad <- exp(intercept + coef_log10_Mean_Converted_Value * new_data$log10_Mean_Converted_Value)
plot_pglmm_cad <- ggplot() +
  geom_point(data = data_frame, aes(x = log10(Mean_Converted_Value), y = cad), color = "#DC4D10", alpha = 0.6) +
  geom_line(data = new_data, aes(x = log10_Mean_Converted_Value, y = predicted_cad), size = 1.2, color = "black") +
  labs(x = "Log Species Mass", y = "cad") +
  theme_pubr() +
  theme(
    legend.position = "right",
    text = element_text(family = "lm", size = 14),  # General text font
    axis.text.x = element_blank(),
    axis.text.y = element_text(angle = 90, vjust = 0.5, hjust = 0.5),  # Rotate y-axis tick label
    axis.title.x = element_blank(),  # X-axis title font
    axis.title.y = element_blank(),
    plot.title = element_text(family = "lm-bold", size = 16),  # Title font
    axis.text = element_text(family = "lm", size = 14)  # Axis text font
  )+
  annotate("text", x = -Inf, y = Inf, 
           label = "p < 0.001", hjust = -0.1, vjust = 5.5, size = 5, fontface = "italic", family = "lm") +
  annotate("text", x = -Inf, y = Inf, 
           label = "R^2 == 0.14", hjust = -0.1, vjust = 3.2, size = 5, family = "lm", parse = TRUE) +
  annotate("text", x = -Inf, y = Inf, 
           label = "PGLMM Poisson", hjust = -0.06, vjust = 2.75, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = Inf, 
           label = "cad", hjust = -0.2, vjust = 1.5, size = 5, fontface = "italic", family = "lm-bold")
print(plot_pglmm_cad)


# Multiple plot axis naming
bottom <- textGrob("Log Species Mass", gp = gpar(fontsize = 19, fontfamily = "lm-bold"))
left <- textGrob("Number of Hox Genes", gp = gpar(fontsize = 19, fontfamily = "lm-bold"), rot = 90)
plot_axis_theme = theme(axis.text.x = element_blank())

# Align the plots
aligned_plots <- align_plots(plot_pglmm_hoxsums+plot_axis_theme, plot_pglmm_abd_A+plot_axis_theme, plot_pglmm_btn+plot_axis_theme, plot_pglmm_cad+plot_axis_theme,
                             plot_pglmm_Dfd+plot_axis_theme, plot_pglmm_eve+plot_axis_theme, plot_pglmm_exex+plot_axis_theme, plot_pglmm_ftz+plot_axis_theme, 
                             plot_pglmm_ind+plot_axis_theme, plot_pglmm_lab+plot_axis_theme, plot_pglmm_Pb+plot_axis_theme, plot_pglmm_Ro+plot_axis_theme, 
                             plot_pglmm_Scr, plot_pglmm_Ubx, plot_pglmm_unpg, plot_pglmm_zen, 
                             align = 'hv')  # Align both horizontally and vertically

combined_plot <- grid.arrange(plot_pglmm_hoxsums+plot_axis_theme, plot_pglmm_abd_A+plot_axis_theme, plot_pglmm_btn+plot_axis_theme, plot_pglmm_cad+plot_axis_theme,
                              plot_pglmm_Dfd+plot_axis_theme, plot_pglmm_eve+plot_axis_theme, plot_pglmm_exex+plot_axis_theme, plot_pglmm_ftz+plot_axis_theme, 
                              plot_pglmm_ind+plot_axis_theme, plot_pglmm_lab+plot_axis_theme, plot_pglmm_Pb+plot_axis_theme, plot_pglmm_Ro+plot_axis_theme, 
                              plot_pglmm_Scr, plot_pglmm_Ubx, plot_pglmm_unpg, plot_pglmm_zen, 
                              ncol = 4, bottom = bottom, left = left)
