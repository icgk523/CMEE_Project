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
library(evobiR)
library(PKNCA)
library(parallel)
library(sysfonts)
library(showtext)
library(cowplot)
source("rsquared_gls_function.R")

calibrated_phylogeny = read.tree("../data/phylogeny/calibrated_tree.nwk"); calibrated_phylogeny$tip.label = sub("_", " ", calibrated_phylogeny$tip.label) # Load entire phylogeny and fix naming
masses = read.csv("../data/mass/Insect_Masses_and_Genomes.csv")[,-1] # Load insect masses without ID column

genomesize = read.csv("../data/genome/ncbi_assembly_information.tsv", sep = "\t"); colnames(genomesize)[1] = "Accession"
masses$Accession <- gsub("^GCA_|^GCF_|\\.\\d+$", "", masses$Accession) # Remove GCA and GCF, as these changed when pulled from NCBI
genomesize$Accession <- gsub("^GCA_|^GCF_|\\.\\d+$", "", genomesize$Accession) # Remove GCA and GCF, as these changed when pulled from NCBI

masses = merge(masses, genomesize, by = "Accession")
masses = distinct(masses)

masses$Log.Assembly.Stats.Total.Sequence.Length = log10(masses$Assembly.Stats.Total.Sequence.Length)
masses$Log_Mean_Converted_Value = log10(masses$Mean_Converted_Value)

masses <- masses %>% group_by(Accession) %>% filter(!is.na(Annotation.BUSCO.Complete) | n() == 1) %>% slice(1) %>% ungroup()
masses = subset(masses, Assembly.Level == "Chromosome" | Assembly.Level == "Complete Genome")

# Custom Fonts for Plots
font_add(family = "lm", regular = "/usr/share/texlive/texmf-dist/fonts/opentype/public/lm/lmroman10-regular.otf")
font_add(family = "lm-bold", regular = "/usr/share/texmf/fonts/opentype/public/lm/lmroman10-bold.otf")
showtext_auto()

# Mass Distributions
log_plot = ggplot(masses) + 
  geom_histogram(aes(x = Log_Mean_Converted_Value, fill = "Log Body Mass"), 
                 color = "black", 
                 alpha = 0.8) +
  theme_classic() +
  labs(x = "Body Mass Values", y = "Number of Species", fill = "Mass Type") +
  facet_wrap(~ Order) +
  scale_fill_manual(values = c("Log Body Mass" = "#DC4D10")) +
  theme(
    legend.position = "top", 
    legend.title = element_blank(),
    text = element_text(family = "lm-bold"),               # Applies font family to all text
    plot.title = element_text(family = "lm-bold"),         # Title text
    plot.subtitle = element_text(family = "lm-bold"),      # Subtitle text
    plot.caption = element_text(family = "lm-bold"),       # Caption text
    axis.title = element_text(family = "lm-bold"),         # Axis titles
    axis.text = element_text(family = "lm-bold"),          # Axis text
    legend.text = element_text(family = "lm-bold"),        # Legend text
    strip.text = element_text(family = "lm-bold")          # Facet text
  )
normal_plot = ggplot(masses) + 
  geom_histogram(aes(x = Mean_Converted_Value, fill = "Body Mass"), 
                 color = "black", 
                 alpha = 0.8) +
  theme_classic() +
  labs(x = "Body Mass Values", y = "Number of Species", fill = "#109FDC") +
  facet_wrap(~ Order) +
  scale_fill_manual(values = c("Body Mass" = "#109FDC")) +
  theme(
    legend.position = "top", 
    legend.title = element_blank(),
    text = element_text(family = "lm-bold"),               # Applies font family to all text
    plot.title = element_text(family = "lm-bold"),         # Title text
    plot.subtitle = element_text(family = "lm-bold"),      # Subtitle text
    plot.caption = element_text(family = "lm-bold"),       # Caption text
    axis.title = element_text(family = "lm-bold"),         # Axis titles
    axis.text = element_text(family = "lm-bold"),          # Axis text
    legend.text = element_text(family = "lm-bold"),        # Legend text
    strip.text = element_text(family = "lm-bold")          # Facet text
  )
combined_plot <- egg::ggarrange(normal_plot, log_plot, ncol = 2)

# Genome Size Distributions
log_plot = ggplot(masses) + 
  geom_histogram(aes(x = Log.Assembly.Stats.Total.Sequence.Length, fill = "Log Genome Size"), 
                 color = "black", 
                 alpha = 0.8) +
  theme_classic() +
  labs(x = "Genome Size (Nt)", y = "Number of Species", fill = "Genome Type") +
  facet_wrap(~ Order) +
  scale_fill_manual(values = c("Log Genome Size" = "#DC4D10")) +
  theme(
    legend.position = "top", 
    legend.title = element_blank(),
    text = element_text(family = "lm-bold"),               # Applies font family to all text
    plot.title = element_text(family = "lm-bold"),         # Title text
    plot.subtitle = element_text(family = "lm-bold"),      # Subtitle text
    plot.caption = element_text(family = "lm-bold"),       # Caption text
    axis.title = element_text(family = "lm-bold"),         # Axis titles
    axis.text = element_text(family = "lm-bold"),          # Axis text
    legend.text = element_text(family = "lm-bold"),        # Legend text
    strip.text = element_text(family = "lm-bold")          # Facet text
  )
normal_plot = ggplot(masses) + 
  geom_histogram(aes(x = Assembly.Stats.Total.Sequence.Length, fill = "Genome Size"), 
                 color = "black", 
                 alpha = 0.8) +
  theme_classic() +
  labs(x = "Body Mass Values", y = "Number of Species", fill = "#109FDC") +
  facet_wrap(~ Order) +
  scale_fill_manual(values = c("Genome Size" = "#109FDC")) +
  theme(
    legend.position = "top", 
    legend.title = element_blank(),
    text = element_text(family = "lm-bold"),               # Applies font family to all text
    plot.title = element_text(family = "lm-bold"),         # Title text
    plot.subtitle = element_text(family = "lm-bold"),      # Subtitle text
    plot.caption = element_text(family = "lm-bold"),       # Caption text
    axis.title = element_text(family = "lm-bold"),         # Axis titles
    axis.text = element_text(family = "lm-bold"),          # Axis text
    legend.text = element_text(family = "lm-bold"),        # Legend text
    strip.text = element_text(family = "lm-bold")          # Facet text
  )
combined_plot <- egg::ggarrange(normal_plot, log_plot, ncol = 1)


########################################################
###### Coleoptera #######
masses_subset = subset(masses, Order == "Coleoptera")
calibrated_phylogeny_subset = keep.tip(calibrated_phylogeny, masses_subset$Species[masses_subset$Species %in% calibrated_phylogeny$tip.label])
masses_subset = masses_subset[masses_subset$Species %in%calibrated_phylogeny_subset$tip.label, ]

formula <- as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ Log_Mean_Converted_Value"))
formula_sq = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)"))
null_formula = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ 1"))

linear_model = lm(formula, data = masses_subset) # Linear model
linear_model_null = lm(null_formula, data = masses_subset) # Linear model
residuals = residuals(linear_model)
names(residuals) = masses_subset$Species
lambda_test <- phylosig(calibrated_phylogeny_subset, residuals, method = "lambda", test = TRUE)

model_fitting_functions <- list(
  pglsModel_Lambda = function() gls(formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Brownian = function() gls(formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins = function() gls(formula, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg = function() gls(formula, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen = function() gls(formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  # pglsModel_Brownian_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  # pglsModel_Lambda_Sq = function() gls(formula_sq, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  # pglsModel_Martins_Sq = function() gls(formula_sq, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  # pglsModel_Blomberg_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  # pglsModel_Grafen_Sq = function() gls(formula_sq, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  pglsModel_Brownian_Null = function() gls(null_formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Null = function() gls(null_formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Null = function() gls(null_formula, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Null = function() gls(null_formula, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Null = function() gls(null_formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML")
)

models <- list()
aic_scores <- list()

for (model_name in names(model_fitting_functions)) {
  fit_function <- model_fitting_functions[[model_name]]
  fit_result <- try(fit_function(), silent = TRUE)
  
  if (!inherits(fit_result, "try-error")) {
    models[[model_name]] <- fit_result
    aic_scores[[model_name]] <- AIC(fit_result)
  } else {
    message(paste("Model fitting failed for", model_name))
  }
}

AIC(linear_model, linear_model_null)
for (model_name in names(aic_scores)) {
  cat(model_name, ":", round(aic_scores[[model_name]], 2), "\n")
}
models[[length(models)+1]] = linear_model
models[[length(models)+1]] = linear_model_null
best_model <- get.best.model(models)

# Find the name of the best model by comparing it to known models
best_model_name <- NULL
for (name in names(models)) {
  if (identical(best_model, models[[name]])) {
    best_model_name <- name
    break
  }
}

range_values <- seq(min(masses_subset$Log_Mean_Converted_Value), max(masses_subset$Log_Mean_Converted_Value), length.out = 100)
new_data <- data.frame(Log_Mean_Converted_Value = range_values)
new_data$Predicted_Values <- predict(best_model, newdata = new_data)

plot_best_model_coleoptera <- ggplot() +
  geom_point(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), alpha = 0.8, color = "#DC4D10") +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Predicted_Values), linewidth = 1.2) +
  geom_smooth(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), method = "lm", linetype = "dashed", se = FALSE, color = "black", linewidth = 1) + 
  labs(x = "Log Species Mass", y = "Log Genome Assembly Length", title = "b) Coleoptera") +
  theme_pubr() +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))
rsquared_coleoptera = rsquared.gls(best_model)[[4]]
p_value_coleoptera = coef(summary(best_model))[[8]]
model_coleoptera = "Brownian"
# Print the plot
print(plot_best_model_coleoptera)



############################################################################################################################################################################################################################
###### Hymenoptera #######
masses_subset = subset(masses, Order == "Hymenoptera")
calibrated_phylogeny_subset = keep.tip(calibrated_phylogeny, masses_subset$Species[masses_subset$Species %in% calibrated_phylogeny$tip.label])
masses_subset = masses_subset[masses_subset$Species %in%calibrated_phylogeny_subset$tip.label, ]

formula <- as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ Log_Mean_Converted_Value"))
formula_sq = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)"))
null_formula = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ 1"))

linear_model = lm(formula, data = masses_subset) # Linear model
linear_model_null = lm(null_formula, data = masses_subset) # Linear model

model_fitting_functions <- list(
  pglsModel_Lambda = function() gls(formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Brownian = function() gls(formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins = function() gls(formula, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg = function() gls(formula, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen = function() gls(formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  pglsModel_Brownian_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Sq = function() gls(formula_sq, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Sq = function() gls(formula_sq, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Sq = function() gls(formula_sq, data = masses_subset, correlation = corGrafen(0.001, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  pglsModel_Brownian_Null = function() gls(null_formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Null = function() gls(null_formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Null = function() gls(null_formula, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Null = function() gls(null_formula, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Null = function() gls(null_formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML")
)

models <- list()
aic_scores <- list()

for (model_name in names(model_fitting_functions)) {
  fit_function <- model_fitting_functions[[model_name]]
  fit_result <- try(fit_function(), silent = TRUE)
  
  if (!inherits(fit_result, "try-error")) {
    models[[model_name]] <- fit_result
    aic_scores[[model_name]] <- AIC(fit_result)
  } else {
    message(paste("Model fitting failed for", model_name))
  }
}

AIC(linear_model, linear_model_null)
for (model_name in names(aic_scores)) {
  cat(model_name, ":", round(aic_scores[[model_name]], 2), "\n")
}
models[[length(models)+1]] = linear_model
models[[length(models)+1]] = linear_model_null

best_model <- get.best.model(models)

# Find the name of the best model by comparing it to known models
best_model_name <- NULL
for (name in names(models)) {
  if (identical(best_model, models[[name]])) {
    best_model_name <- name
    break
  }
}

range_values <- seq(min(masses_subset$Log_Mean_Converted_Value), max(masses_subset$Log_Mean_Converted_Value), length.out = 100)
new_data <- data.frame(Log_Mean_Converted_Value = range_values)
new_data$Predicted_Values <- predict(best_model, newdata = new_data)

plot_best_model_hymenoptera <- ggplot() +
  geom_point(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), alpha = 0.8, color = "#DC4D10") +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Predicted_Values), size = 1.2, linetype = "dashed") +
  # geom_smooth(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) +
  labs(x = "Log Species Mass", y = "Log Genome Assembly Length", title = "e) Hymenoptera") +
  theme_pubr() +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))

# Print the plot
print(plot_best_model_hymenoptera)
rsquared_hymenoptera = 0
p_value_hymenoptera = coef(summary(best_model))[[4]]
model_hymenoptera = ""


############################################################################################################################################################################################################################
###### Diptera #######
masses_subset = subset(masses, Order == "Diptera")
calibrated_phylogeny_subset = keep.tip(calibrated_phylogeny, masses_subset$Species[masses_subset$Species %in% calibrated_phylogeny$tip.label])
masses_subset = masses_subset[masses_subset$Species %in%calibrated_phylogeny_subset$tip.label, ]

formula <- as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ Log_Mean_Converted_Value"))
formula_sq = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)"))
null_formula = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ 1"))

linear_model = lm(formula, data = masses_subset) # Linear model
linear_model_null = lm(null_formula, data = masses_subset) # Linear model

model_fitting_functions <- list(
  pglsModel_Lambda = function() gls(formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Brownian = function() gls(formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins = function() gls(formula, data = masses_subset, correlation = corMartins(0, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg = function() gls(formula, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen = function() gls(formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  # pglsModel_Brownian_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  # pglsModel_Lambda_Sq = function() gls(formula_sq, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  # pglsModel_Martins_Sq = function() gls(formula_sq, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  # pglsModel_Blomberg_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBlomberg(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  # pglsModel_Grafen_Sq = function() gls(formula_sq, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  # 
  pglsModel_Brownian_Null = function() gls(null_formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Null = function() gls(null_formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Null = function() gls(null_formula, data = masses_subset, correlation = corMartins(0, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Null = function() gls(null_formula, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Null = function() gls(null_formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML")
)

models <- list()
aic_scores <- list()

for (model_name in names(model_fitting_functions)) {
  fit_function <- model_fitting_functions[[model_name]]
  fit_result <- try(fit_function(), silent = TRUE)
  
  if (!inherits(fit_result, "try-error")) {
    models[[model_name]] <- fit_result
    aic_scores[[model_name]] <- AIC(fit_result)
  } else {
    message(paste("Model fitting failed for", model_name))
  }
}

AIC(linear_model, linear_model_null)
for (model_name in names(aic_scores)) {
  cat(model_name, ":", round(aic_scores[[model_name]], 2), "\n")
}
models[[length(models)+1]] = linear_model
models[[length(models)+1]] = linear_model_null

best_model <- get.best.model(models)
best_model = gls(formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML")
# Find the name of the best model by comparing it to known models
best_model_name <- NULL
for (name in names(models)) {
  if (identical(best_model, models[[name]])) {
    best_model_name <- name
    break
  }
}

range_values <- seq(min(masses_subset$Log_Mean_Converted_Value), max(masses_subset$Log_Mean_Converted_Value), length.out = 100)
new_data <- data.frame(Log_Mean_Converted_Value = range_values)
new_data$Predicted_Values <- predict(best_model, newdata = new_data)

plot_best_model_diptera <- ggplot() +
  geom_point(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), alpha = 0.8, color = "#DC4D10") +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Predicted_Values), size = 1.2) +
  geom_smooth(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) +
  labs(x = "Log Species Mass", y = "Log Genome Assembly Length", title = "c) Diptera") +
  theme_pubr() +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))

print(plot_best_model_diptera)
rsquared_diptera = rsquared.gls(best_model)[[4]]
p_value_diptera = coef(summary(best_model))[[8]]
model_diptera = "Brownian"

############################################################################################################################################################################################################################
###### Lepidoptera #######
masses_subset = subset(masses, Order == "Lepidoptera")
calibrated_phylogeny_subset = keep.tip(calibrated_phylogeny, masses_subset$Species[masses_subset$Species %in% calibrated_phylogeny$tip.label])
masses_subset = masses_subset[masses_subset$Species %in%calibrated_phylogeny_subset$tip.label, ]

formula <- as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ Log_Mean_Converted_Value"))
formula_sq = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)"))
null_formula = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ 1"))

linear_model = lm(formula, data = masses_subset) # Linear model
linear_model_null = lm(null_formula, data = masses_subset) # Linear model

model_fitting_functions <- list(
  pglsModel_Lambda = function() gls(formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Brownian = function() gls(formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins = function() gls(formula, data = masses_subset, correlation = corMartins(0, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg = function() gls(formula, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen = function() gls(formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  pglsModel_Brownian_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Sq = function() gls(formula_sq, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Sq = function() gls(formula_sq, data = masses_subset, correlation = corMartins(0, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBlomberg(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Sq = function() gls(formula_sq, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  pglsModel_Brownian_Null = function() gls(null_formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Null = function() gls(null_formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Null = function() gls(null_formula, data = masses_subset, correlation = corMartins(0, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Null = function() gls(null_formula, data = masses_subset, correlation = corBlomberg(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Null = function() gls(null_formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML")
)

models <- list()
aic_scores <- list()

for (model_name in names(model_fitting_functions)) {
  fit_function <- model_fitting_functions[[model_name]]
  fit_result <- try(fit_function(), silent = TRUE)
  
  if (!inherits(fit_result, "try-error")) {
    models[[model_name]] <- fit_result
    aic_scores[[model_name]] <- AIC(fit_result)
  } else {
    message(paste("Model fitting failed for", model_name))
  }
}

AIC(linear_model, linear_model_null)
for (model_name in names(aic_scores)) {
  cat(model_name, ":", round(aic_scores[[model_name]], 2), "\n")
}
models[[length(models)+1]] = linear_model
models[[length(models)+1]] = linear_model_null

best_model <- get.best.model(models)

# Find the name of the best model by comparing it to known models
best_model_name <- NULL
for (name in names(models)) {
  if (identical(best_model, models[[name]])) {
    best_model_name <- name
    break
  }
}

range_values <- seq(min(masses_subset$Log_Mean_Converted_Value), max(masses_subset$Log_Mean_Converted_Value), length.out = 100)
new_data <- data.frame(Log_Mean_Converted_Value = range_values)
new_data$Predicted_Values <- predict(best_model, newdata = new_data)

plot_best_model_lepidoptera <- ggplot() +
  geom_point(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), alpha = 0.8, color = "#DC4D10") +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Predicted_Values), size = 1.2) +
  geom_smooth(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) +
  labs(x = "Log Species Mass", y = "Log Genome Assembly Length", title = "f) Lepidoptera") +
  theme_pubr() +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))
print(plot_best_model_lepidoptera)
rsquared_lepidoptera = rsquared.gls(best_model)[[4]]
p_value_lepidoptera = coef(summary(best_model))[[4]]
model_lepidoptera = "Blomberg"

############################################################################################################################################################################################################################
###### Hemiptera #######
masses_subset = subset(masses, Order == "Hemiptera")
calibrated_phylogeny_subset = keep.tip(calibrated_phylogeny, masses_subset$Species[masses_subset$Species %in% calibrated_phylogeny$tip.label])
masses_subset = masses_subset[masses_subset$Species %in%calibrated_phylogeny_subset$tip.label, ]

formula <- as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ Log_Mean_Converted_Value"))
formula_sq = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)"))
null_formula = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ 1"))

linear_model = lm(formula, data = masses_subset) # Linear model
linear_model_null = lm(null_formula, data = masses_subset) # Linear model

model_fitting_functions <- list(
  pglsModel_Lambda = function() gls(formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Brownian = function() gls(formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins = function() gls(formula, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg = function() gls(formula, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen = function() gls(formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  pglsModel_Brownian_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Sq = function() gls(formula_sq, data = masses_subset, correlation = corPagel(0.9, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Sq = function() gls(formula_sq, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBlomberg(0.008, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Sq = function() gls(formula_sq, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  pglsModel_Brownian_Null = function() gls(null_formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Null = function() gls(null_formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Null = function() gls(null_formula, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Null = function() gls(null_formula, data = masses_subset, correlation = corBlomberg(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Null = function() gls(null_formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML")
)

models <- list()
aic_scores <- list()

for (model_name in names(model_fitting_functions)) {
  fit_function <- model_fitting_functions[[model_name]]
  fit_result <- try(fit_function(), silent = TRUE)
  
  if (!inherits(fit_result, "try-error")) {
    models[[model_name]] <- fit_result
    aic_scores[[model_name]] <- AIC(fit_result)
  } else {
    message(paste("Model fitting failed for", model_name))
  }
}

AIC(linear_model, linear_model_null)
for (model_name in names(aic_scores)) {
  cat(model_name, ":", round(aic_scores[[model_name]], 2), "\n")
}
models[[length(models)+1]] = linear_model
models[[length(models)+1]] = linear_model_null

best_model <- get.best.model(models)

# Find the name of the best model by comparing it to known models
best_model_name <- NULL
for (name in names(models)) {
  if (identical(best_model, models[[name]])) {
    best_model_name <- name
    break
  }
}

range_values <- seq(min(masses_subset$Log_Mean_Converted_Value), max(masses_subset$Log_Mean_Converted_Value), length.out = 100)
new_data <- data.frame(Log_Mean_Converted_Value = range_values)
new_data$Predicted_Values <- predict(best_model, newdata = new_data)

plot_best_model_Hemiptera <- ggplot() +
  geom_point(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), alpha = 0.8, color = "#DC4D10") +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Predicted_Values), size = 1.2) +
  geom_smooth(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) +
  labs(x = "Log Species Mass", y = "Log Genome Assembly Length", title = "d) Hemiptera") +
  theme_pubr() +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))
print(plot_best_model_Hemiptera)
rsquared_hemiptera = rsquared.gls(best_model)[[4]]
p_value_hemiptera = coef(summary(best_model))[[8]]
model_hemiptera = "Lambda"
############################################################################################################################################################################################################################
###### Odonata #######
masses_subset = subset(masses, Order == "Odonata")
calibrated_phylogeny_subset = keep.tip(calibrated_phylogeny, masses_subset$Species[masses_subset$Species %in% calibrated_phylogeny$tip.label])
masses_subset = masses_subset[masses_subset$Species %in%calibrated_phylogeny_subset$tip.label, ]

formula <- as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ Log_Mean_Converted_Value"))
formula_sq = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)"))
null_formula = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ 1"))

linear_model = lm(formula, data = masses_subset) # Linear model
linear_model_null = lm(null_formula, data = masses_subset) # Linear model

model_fitting_functions <- list(
  pglsModel_Lambda = function() gls(formula, data = masses_subset, correlation = corPagel(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Brownian = function() gls(formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins = function() gls(formula, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg = function() gls(formula, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen = function() gls(formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  pglsModel_Brownian_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Sq = function() gls(formula_sq, data = masses_subset, correlation = corPagel(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Sq = function() gls(formula_sq, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBlomberg(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Sq = function() gls(formula_sq, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  pglsModel_Brownian_Null = function() gls(null_formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Null = function() gls(null_formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Null = function() gls(null_formula, data = masses_subset, correlation = corMartins(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Null = function() gls(null_formula, data = masses_subset, correlation = corBlomberg(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Null = function() gls(null_formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML")
)

models <- list()
aic_scores <- list()

for (model_name in names(model_fitting_functions)) {
  fit_function <- model_fitting_functions[[model_name]]
  fit_result <- try(fit_function(), silent = TRUE)
  
  if (!inherits(fit_result, "try-error")) {
    models[[model_name]] <- fit_result
    aic_scores[[model_name]] <- AIC(fit_result)
  } else {
    message(paste("Model fitting failed for", model_name))
  }
}

AIC(linear_model, linear_model_null)
for (model_name in names(aic_scores)) {
  cat(model_name, ":", round(aic_scores[[model_name]], 2), "\n")
}
models[[length(models)+1]] = linear_model
models[[length(models)+1]] = linear_model_null

best_model <- get.best.model(models)

# Find the name of the best model by comparing it to known models
best_model_name <- NULL
for (name in names(models)) {
  if (identical(best_model, models[[name]])) {
    best_model_name <- name
    break
  }
}

range_values <- seq(min(masses_subset$Log_Mean_Converted_Value), max(masses_subset$Log_Mean_Converted_Value), length.out = 100)
new_data <- data.frame(Log_Mean_Converted_Value = range_values)
new_data$Predicted_Values <- predict(best_model, newdata = new_data)

plot_best_model_Odonata <- ggplot() +
  geom_point(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), alpha = 0.8, color = "#DC4D10") +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Predicted_Values), size = 1.2) +
  geom_smooth(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) +
  labs(x = "Log Species Mass", y = "Log Genome Assembly Length") +
  theme_pubr() +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))

print(plot_best_model_Odonata)


############################################################################################################################################################################################################################
###### Insecta #######
masses_subset = masses
calibrated_phylogeny_subset = keep.tip(calibrated_phylogeny, masses_subset$Species[masses_subset$Species %in% calibrated_phylogeny$tip.label])
masses_subset = masses_subset[masses_subset$Species %in%calibrated_phylogeny_subset$tip.label, ]

formula <- as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ Log_Mean_Converted_Value"))
formula_sq = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)"))
null_formula = as.formula(paste("Log.Assembly.Stats.Total.Sequence.Length ~ 1"))

linear_model = lm(formula, data = masses_subset) # Linear model
linear_model_null = lm(null_formula, data = masses_subset) # Linear model

model_fitting_functions <- list(
  pglsModel_Lambda = function() gls(formula, data = masses_subset, correlation = corPagel(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Brownian = function() gls(formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins = function() gls(formula, data = masses_subset, correlation = corMartins(0, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg = function() gls(formula, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen = function() gls(formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  pglsModel_Brownian_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Sq = function() gls(formula_sq, data = masses_subset, correlation = corPagel(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Sq = function() gls(formula_sq, data = masses_subset, correlation = corMartins(0, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Sq = function() gls(formula_sq, data = masses_subset, correlation = corBlomberg(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Sq = function() gls(formula_sq, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  
  pglsModel_Brownian_Null = function() gls(null_formula, data = masses_subset, correlation = corBrownian(1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Lambda_Null = function() gls(null_formula, data = masses_subset, correlation = corPagel(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Martins_Null = function() gls(null_formula, data = masses_subset, correlation = corMartins(0, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Blomberg_Null = function() gls(null_formula, data = masses_subset, correlation = corBlomberg(0.1, calibrated_phylogeny_subset, form = ~Species), method = "ML"),
  pglsModel_Grafen_Null = function() gls(null_formula, data = masses_subset, correlation = corGrafen(1, calibrated_phylogeny_subset, form = ~Species), method = "ML")
)

models <- list()
aic_scores <- list()

for (model_name in names(model_fitting_functions)) {
  fit_function <- model_fitting_functions[[model_name]]
  fit_result <- try(fit_function(), silent = TRUE)
  
  if (!inherits(fit_result, "try-error")) {
    models[[model_name]] <- fit_result
    aic_scores[[model_name]] <- AIC(fit_result)
    
    # hist_residuals <- residuals(fit_result)
    # hist(hist_residuals, breaks = 10, main = paste("Residuals Histogram for", model_name),
    #      xlab = "Residuals", col = "blue", border = "black")
  } else {
    message(paste("Model fitting failed for", model_name))
  }
}

AIC(linear_model, linear_model_null)
for (model_name in names(aic_scores)) {
  cat(model_name, ":", round(aic_scores[[model_name]], 2), "\n")
}
models[[length(models)+1]] = linear_model
models[[length(models)+1]] = linear_model_null

best_model <- get.best.model(models)
best_model = gls(formula, data = masses_subset, correlation = corBlomberg(0.2, calibrated_phylogeny_subset, form = ~Species), method = "ML")
# Find the name of the best model by comparing it to known models
best_model_name <- NULL
for (name in names(models)) {
  if (identical(best_model, models[[name]])) {
    best_model_name <- name
    break
  }
}

range_values <- seq(min(masses_subset$Log_Mean_Converted_Value), max(masses_subset$Log_Mean_Converted_Value), length.out = 100)
new_data <- data.frame(Log_Mean_Converted_Value = range_values)
new_data$Predicted_Values <- predict(best_model, newdata = new_data)

plot_best_model_Insecta <- ggplot() +
  geom_point(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), alpha = 0.8, color = "#DC4D10") +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Predicted_Values), size = 1.2) +
  geom_smooth(data = masses_subset, aes(x = Log_Mean_Converted_Value, y = Log.Assembly.Stats.Total.Sequence.Length), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) +
  labs(x = "Log Species Mass", y = "Log Genome Assembly Length", title = "a) Insecta") +
  theme_pubr() +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))

print(plot_best_model_Insecta)
rsquared_insecta = rsquared.gls(best_model)[[4]]
p_value_insecta = coef(summary(best_model))[[4]]
model_insecta = "Blomberg"

############################################################################################################################################################################################################################
###### Plotting #######


font_add(family = "lm", regular = "/usr/share/texlive/texmf-dist/fonts/opentype/public/lm/lmroman10-regular.otf")
font_add(family = "lm-bold", regular = "/usr/share/texmf/fonts/opentype/public/lm/lmroman10-bold.otf")
showtext_auto()



plot_best_model_Insecta_new = plot_best_model_Insecta + scale_x_continuous(breaks = c(-2, 0, 2)) + scale_y_continuous(breaks = c(8.5,9, 9.5))
plot_best_model_coleoptera_new =  plot_best_model_coleoptera + scale_y_continuous(breaks = c(8.4, 8.8, 9.2)) + scale_x_continuous(breaks = c(-1, 0, 1, 2))
plot_best_model_diptera_new = plot_best_model_diptera + scale_x_continuous(breaks = c(0, 1, 2)) + scale_y_continuous(breaks = c(8.2, 8.6, 9))
plot_best_model_Hemiptera_new= plot_best_model_Hemiptera + scale_y_continuous(breaks = c(8.6, 8.8, 9))
plot_best_model_hymenoptera_new = plot_best_model_hymenoptera + scale_x_continuous(breaks = c(0, 1, 2)) + scale_y_continuous(breaks = c(8.4, 8.6, 8.8))
plot_best_model_lepidoptera_new = plot_best_model_lepidoptera + scale_x_continuous(breaks = c(1, 2, 3)) + scale_y_continuous(breaks = c(8.6, 9, 9.4))

annotate_plot <- function(plot, model_name, p_value, rsquared) {
  formatted_rsquared <- ifelse(rsquared == 0, "0", sprintf("%.2f", rsquared))
  formatted_p_value <- ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
  plot +
    annotate("text", x = -Inf, y = Inf, label = model_name, hjust = -0.1, vjust = 2, size = 5, family = "lm-bold") +
    annotate("text", x = -Inf, y = Inf, label = paste("P-value:", formatted_p_value), hjust = -0.1, vjust = 3.5, size = 5, family = "lm") +
    annotate("text", x = -Inf, y = Inf, label = paste("R-squared:", formatted_rsquared), hjust = -0.1, vjust = 5, size = 5, family = "lm")
}


# Annotating each plot with its corresponding values
plot_best_model_Insecta <- annotate_plot(plot_best_model_Insecta_new, model_insecta, p_value_insecta, rsquared_insecta)
plot_best_model_coleoptera <- annotate_plot(plot_best_model_coleoptera_new, model_coleoptera, p_value_coleoptera, rsquared_coleoptera)
plot_best_model_diptera <- annotate_plot(plot_best_model_diptera_new, model_diptera, p_value_diptera, rsquared_diptera)
plot_best_model_Hemiptera <- annotate_plot(plot_best_model_Hemiptera_new, model_hemiptera, p_value_hemiptera, rsquared_hemiptera)
plot_best_model_hymenoptera <- annotate_plot(plot_best_model_hymenoptera_new, model_hymenoptera, p_value_hymenoptera, rsquared_hymenoptera)
plot_best_model_lepidoptera <- annotate_plot(plot_best_model_lepidoptera_new, model_lepidoptera, p_value_lepidoptera, rsquared_lepidoptera)


# Adjust themes for each plot: remove y-axis titles and ensure consistent margins
plot_theme_top = theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),  # Remove y-axis title
  axis.text.x = element_text(size = 20),  # Set x-axis text size to 14
  axis.text.y = element_text(size = 20, angle = 90, hjust = 0.5),  # Set y-axis text size to 14
  plot.title = element_text(family = "lm-bold", size = 19),
  plot.margin = margin(b = 50, l = 10),  # Apply consistent margins to all plots
  legend.position = "none",
  axis.line.x = element_line(color = "black", size = 0.75, lineend = "square"), 
  axis.line.y = element_line(color = "black", size = 1, lineend = "square"), 
  axis.ticks.x = element_line(color = "black", size = 0.5), 
  axis.ticks.y = element_line(color = "black", size = 0.5))
plot_theme_bottom = theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),  # Remove y-axis title
  axis.text.x = element_text(size = 20),  # Set x-axis text size to 14
  axis.text.y = element_text(size = 20, angle = 90, hjust = 0.5),  # Set y-axis text size to 14
  plot.title = element_text(family = "lm-bold", size = 19),
  plot.margin = margin(b = 25, l = 10),  # Apply consistent margins to all plots
  legend.position = "none",
  axis.line.x = element_line(color = "black", size = 0.75, lineend = "square"), 
  axis.line.y = element_line(color = "black", size = 1, lineend = "square"), 
  axis.ticks.x = element_line(color = "black", size = 0.5), 
  axis.ticks.y = element_line(color = "black", size = 0.5))

combined_plot <- egg::ggarrange(plot_best_model_Insecta + plot_theme_top, plot_best_model_coleoptera + plot_theme_top, plot_best_model_diptera+ plot_theme_top, plot_best_model_Hemiptera+ plot_theme_bottom, plot_best_model_hymenoptera+ plot_theme_bottom, plot_best_model_lepidoptera+ plot_theme_bottom, ncol = 3)


plot_grid(combined_plot) +
  annotate("text", label = "Log Species Mass", x = 0.5, y = 0.02, vjust = 1.75, size = 7.5, family = "lm-bold") + 
  annotate("text", label = "Log Genome Size", x = 0.02, y = 0.5, vjust = -2, angle = 90, size = 7.5, family = "lm-bold") +
  theme(plot.margin = margin(l = 30, r = 10, t = 10, b = 20)) + 
  theme(text = element_text(family = "lm"))
