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
library(viridis)
library(extrafont)
library(showtext)
library(piecewiseSEM)
library(gtable)
library(viridis)
library(cowplot)
library(ggplot2)
library(ggplotify)

repeats = read.table("../data/repetitive_elements/Dataframe_ordered_by_tree_final3.txt", header = TRUE) # This is the RE dataset
repeats_names = read.csv("../data/repetitive_elements/outputs.csv", header = TRUE) # These are the genomes we can get data for from the RE dataset
repeats_redmask = read.csv("../data/repetitive_elements/Repeats_Genomes.csv") # These are the genomes we have mass data for
masses = read.csv("../data/mass/Insect_Masses.csv")

repeats = repeats %>% # Merge the three datasets
  mutate(Accession = str_remove(Accession, "\\..*$")) %>%
  inner_join(repeats_names %>%
               mutate(Accession = str_remove(Accession, "\\..*$")), 
             by = "Accession") %>%
  inner_join(masses, by = "Species") %>% filter(TreeOrder != 100) %>%
  mutate(total_TE = LINEs_proportion + SINEs_proportion + LTRs_proportion + DNA_trans_proportion + tandem_repeats_proportion,
         log_total_TE = log10(total_TE))

phylo_tree = read.tree("../data/phylogeny/calibrated_tree.nwk"); phylo_tree$tip.label = sub("_", " ", phylo_tree$tip.label) # Fix tip naming issues
phylo_tree = keep.tip(phylo_tree, repeats$Species[repeats$Species %in% phylo_tree$tip.label]) # Subset tree to species in data
repeats = repeats[repeats$Species %in% phylo_tree$tip.label, ] # Subset data to species present in tree and data

ggplot(repeats) + 
  geom_histogram(aes(x = log_total_TE, fill = "Log Total TE"), 
                 binwidth = 0.1, 
                 color = "black", 
                 alpha = 0.8) +
  geom_histogram(aes(x = total_TE, fill = "Total TE"), 
                 binwidth = 0.05, 
                 color = "black", 
                 alpha = 0.7) +
  theme_classic() +
  labs(x = "TE Values", y = "Number of Species", fill = "TE Type") +
  facet_wrap(~ Order.x) +
  scale_fill_manual(values = c("Log Total TE" = "#DC4D10", "Total TE" = "#109FDC")) +
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


########################################################
###### Insecta #######
repeats_subset = repeats
phylo_tree_subset = keep.tip(phylo_tree, repeats_subset$Species[repeats_subset$Species %in% phylo_tree$tip.label]) # Subset tree to species in data
formula <- as.formula("log_total_TE ~ Log_Mean_Converted_Value") # Linear model
formula_sq <- as.formula("log_total_TE ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)") # Quadratic model
null_formula = as.formula(paste("log_total_TE ~ 1")) # Null model

linear_model <- lm(formula, data = repeats_subset) # Linear model
linear_model_residuals <- residuals(linear_model) # its residuals
names(linear_model_residuals) <- repeats_subset$Species # add species names to residuals

lambda = phylosig(phylo_tree_subset, linear_model_residuals, method = "lambda", test = TRUE); print(lambda) # Test phylogenetic signal

# # If it's proving difficult to find starting values that work, change the model below and test it
# lambdas = seq(0,1,0.01)
# for (lambda_starting in lambdas) { # Loop through each lambda_starting value in lambdas
#   tryCatch({  # Use tryCatch to suppress errors
#     pglsModel_Lambda <- gls(formula_sq, data = repeats_subset, correlation = corBlomberg(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML") # Attempt to create the pgls model
#     print(lambda_starting) # Print the lambda_starting value if successful
#   }, error = function(e) { # Suppress the error and continue
#   })
# }

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = repeats_subset, correlation = corPagel(0.9, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = repeats_subset, correlation = corBrownian(0.9, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = repeats_subset, correlation = corMartins(0.9, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = repeats_subset, correlation = corBlomberg(0.2, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = repeats_subset, correlation = corGrafen(0.9, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = repeats_subset, correlation = corPagel(0.9, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = repeats_subset, correlation = corBrownian(0.9, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = repeats_subset, correlation = corMartins(0.9, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = repeats_subset, correlation = corBlomberg(0.2, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = repeats_subset, correlation = corGrafen(0.9, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = repeats_subset, correlation = corBrownian(0.9, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = repeats_subset, correlation = corPagel(0.9, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = repeats_subset, correlation = corBlomberg(0.2, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = repeats_subset, correlation = corGrafen(0.9, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = repeats_subset, correlation = corMartins(0.9, phylo_tree_subset, form = ~Species), method = "ML")

AIC(linear_model,
    pglsModel_Lambda, pglsModel_Brownian, pglsModel_Martins, pglsModel_Grafen,  pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Martins_sq, pglsModel_Grafen_sq, pglsModel_Blomberg_sq,
    pglsModel_Null_Grafen, pglsModel_Null_Martins, pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Blomberg) # Best model is lambda null after observing the residual distributions

# Define the range of 'Log_Mean_Converted_Value'
intercept <- pglsModel_Null_Pagel$coefficients[1]
insecta = ggplot(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
  geom_point(alpha = 0.8, color = "#DC4D10") +
  geom_hline(yintercept = intercept, color = "black", linewidth = 1) +
  labs(title = "a) Insecta", x = "Log Species Mass", y = "Log Proportion of TE") +
  theme_pubr()+ 
  scale_y_continuous(breaks = c(-1.2, -0.8, -0.4)) 

########################################################
###### Coleoptera #######
repeats_subset = subset(repeats, Order.x == "Coleoptera")
phylo_tree_subset = keep.tip(phylo_tree, repeats_subset$Species[repeats_subset$Species %in% phylo_tree$tip.label]) # Subset tree to species in data
formula <- as.formula("log_total_TE ~ Log_Mean_Converted_Value") # Linear model
formula_sq <- as.formula("log_total_TE ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)") # Quadratic model
null_formula = as.formula(paste("log_total_TE ~ 1")) # Null model

linear_model <- lm(formula, data = repeats_subset) # Linear model
linear_model_residuals <- residuals(linear_model) # its residuals
names(linear_model_residuals) <- repeats_subset$Species # add species names to residuals

lambda = phylosig(phylo_tree_subset, linear_model_residuals, method = "lambda", test = TRUE); print(lambda) # Test phylogenetic signal

# # If it's proving difficult to find starting values that work, change the model below and test it
# lambdas = seq(0,1,0.01)
# for (lambda_starting in lambdas) { # Loop through each lambda_starting value in lambdas
#   tryCatch({  # Use tryCatch to suppress errors
#     pglsModel_Lambda <- gls(formula_sq, data = repeats_subset, correlation = corBlomberg(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML") # Attempt to create the pgls model
#     print(lambda_starting) # Print the lambda_starting value if successful
#   }, error = function(e) { # Suppress the error and continue
#   })
# }

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = repeats_subset, correlation = corPagel(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = repeats_subset, correlation = corBrownian(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = repeats_subset, correlation = corMartins(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = repeats_subset, correlation = corBlomberg(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = repeats_subset, correlation = corGrafen(0.01, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = repeats_subset, correlation = corPagel(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = repeats_subset, correlation = corBrownian(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = repeats_subset, correlation = corMartins(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = repeats_subset, correlation = corBlomberg(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = repeats_subset, correlation = corGrafen(0.01, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = repeats_subset, correlation = corBrownian(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = repeats_subset, correlation = corPagel(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = repeats_subset, correlation = corBlomberg(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = repeats_subset, correlation = corGrafen(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = repeats_subset, correlation = corMartins(0.01, phylo_tree_subset, form = ~Species), method = "ML")

AIC(linear_model,
    pglsModel_Lambda, pglsModel_Brownian, pglsModel_Martins, pglsModel_Grafen,  pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Martins_sq, pglsModel_Grafen_sq, pglsModel_Blomberg_sq,
    pglsModel_Null_Grafen, pglsModel_Null_Martins, pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Blomberg) # Best model is lambda quadratic but it is overfitting to data - best is actually lambda linear

# Define the range of 'Log_Mean_Converted_Value'
range_values <- seq(min(repeats_subset$Log_Mean_Converted_Value), max(repeats_subset$Log_Mean_Converted_Value), length.out = 100) # Make a sequence of data for plotting
predictions <- data.frame( # Make a data frame of predicted values for the data
  Log_Mean_Converted_Value = rep(range_values, 1),
  log_total_TE = pglsModel_Lambda$coefficients[1] +  pglsModel_Lambda$coefficients[2] * range_values,
  Model = rep(c("PGLS Linear"), each = 100)
)

coleoptera = ggplot(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
  geom_point(alpha = 0.8, color = "#DC4D10") +
  geom_line(data = predictions, aes(x = Log_Mean_Converted_Value, y = log_total_TE), linewidth = 1) +
  geom_smooth(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) +
  labs(title = "b) Coleoptera", x = "Log Species Mass", y = "Log Proportion of TE") +
  theme_pubr()


########################################################
###### Diptera #######
repeats_subset = subset(repeats, Order.x == "Diptera")
phylo_tree_subset = keep.tip(phylo_tree, repeats_subset$Species[repeats_subset$Species %in% phylo_tree$tip.label]) # Subset tree to species in data
formula <- as.formula("log_total_TE ~ Log_Mean_Converted_Value") # Linear model
formula_sq <- as.formula("log_total_TE ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)") # Quadratic model
null_formula = as.formula(paste("log_total_TE ~ 1")) # Null model

linear_model <- lm(formula, data = repeats_subset) # Linear model
linear_model_residuals <- residuals(linear_model) # its residuals
names(linear_model_residuals) <- repeats_subset$Species # add species names to residuals

lambda = phylosig(phylo_tree_subset, linear_model_residuals, method = "lambda", test = TRUE); print(lambda) # Test phylogenetic signal

# # If it's proving difficult to find starting values that work, change the model below and test it
# lambdas = seq(0,1,0.01)
# for (lambda_starting in lambdas) { # Loop through each lambda_starting value in lambdas
#   tryCatch({  # Use tryCatch to suppress errors
#     pglsModel_Lambda <- gls(formula_sq, data = repeats_subset, correlation = corBlomberg(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML") # Attempt to create the pgls model
#     print(lambda_starting) # Print the lambda_starting value if successful
#   }, error = function(e) { # Suppress the error and continue
#   })
# }

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = repeats_subset, correlation = corPagel(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = repeats_subset, correlation = corBrownian(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = repeats_subset, correlation = corMartins(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = repeats_subset, correlation = corBlomberg(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = repeats_subset, correlation = corGrafen(1, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = repeats_subset, correlation = corPagel(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = repeats_subset, correlation = corBrownian(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = repeats_subset, correlation = corMartins(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = repeats_subset, correlation = corBlomberg(0.5, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = repeats_subset, correlation = corGrafen(1, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = repeats_subset, correlation = corBrownian(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = repeats_subset, correlation = corPagel(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = repeats_subset, correlation = corBlomberg(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = repeats_subset, correlation = corGrafen(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = repeats_subset, correlation = corMartins(1, phylo_tree_subset, form = ~Species), method = "ML")

AIC(linear_model,
    pglsModel_Lambda, pglsModel_Brownian, pglsModel_Martins, pglsModel_Grafen,  pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Martins_sq, pglsModel_Grafen_sq, pglsModel_Blomberg_sq,
    pglsModel_Null_Grafen, pglsModel_Null_Martins, pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Blomberg) # Best model is brownian null - from residuals best is actually lambda null

# Define the range of 'Log_Mean_Converted_Value'
intercept <- pglsModel_Null_Pagel$coefficients[1]
diptera = ggplot(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
  geom_point(alpha = 0.8, color = "#DC4D10") +
  geom_hline(yintercept = intercept, color = "black", linewidth = 1) +
  labs(title = "c) Diptera", x = "Log Species Mass", y = "Log Proportion of TE") +
  theme_pubr()+ 
  scale_y_continuous(breaks = c(-1.2, -0.8, -0.4)) 

########################################################
###### Hemiptera #######
repeats_subset = subset(repeats, Order.x == "Hemiptera")
phylo_tree_subset = keep.tip(phylo_tree, repeats_subset$Species[repeats_subset$Species %in% phylo_tree$tip.label]) # Subset tree to species in data
formula <- as.formula("log_total_TE ~ Log_Mean_Converted_Value") # Linear model
formula_sq <- as.formula("log_total_TE ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)") # Quadratic model
null_formula = as.formula(paste("log_total_TE ~ 1")) # Null model

linear_model <- lm(formula, data = repeats_subset) # Linear model
linear_model_residuals <- residuals(linear_model) # its residuals
names(linear_model_residuals) <- repeats_subset$Species # add species names to residuals

lambda = phylosig(phylo_tree_subset, linear_model_residuals, method = "lambda", test = TRUE); print(lambda) # Test phylogenetic signal

# # If it's proving difficult to find starting values that work, change the model below and test it
# lambdas = seq(0,1,0.01)
# for (lambda_starting in lambdas) { # Loop through each lambda_starting value in lambdas
#   tryCatch({  # Use tryCatch to suppress errors
#     pglsModel_Lambda <- gls(formula_sq, data = repeats_subset, correlation = corPagel(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML") # Attempt to create the pgls model
#     print(lambda_starting) # Print the lambda_starting value if successful
#   }, error = function(e) { # Suppress the error and continue
#   })
# }

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = repeats_subset, correlation = corPagel(0.04, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = repeats_subset, correlation = corBrownian(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = repeats_subset, correlation = corMartins(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = repeats_subset, correlation = corBlomberg(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = repeats_subset, correlation = corGrafen(0.01, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = repeats_subset, correlation = corPagel(0.04, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = repeats_subset, correlation = corBrownian(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = repeats_subset, correlation = corMartins(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = repeats_subset, correlation = corBlomberg(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = repeats_subset, correlation = corGrafen(0.01, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = repeats_subset, correlation = corBrownian(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = repeats_subset, correlation = corPagel(0.04, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = repeats_subset, correlation = corBlomberg(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = repeats_subset, correlation = corGrafen(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = repeats_subset, correlation = corMartins(0.01, phylo_tree_subset, form = ~Species), method = "ML")

AIC(linear_model,
    pglsModel_Lambda, pglsModel_Brownian, pglsModel_Martins, pglsModel_Grafen,  pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Martins_sq, pglsModel_Grafen_sq, pglsModel_Blomberg_sq,
    pglsModel_Null_Grafen, pglsModel_Null_Martins, pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Blomberg) # Best model is lambda linear - residuals confirm this

# Define the range of 'Log_Mean_Converted_Value'
range_values <- seq(min(repeats_subset$Log_Mean_Converted_Value), max(repeats_subset$Log_Mean_Converted_Value), length.out = 100) # Make a sequence of data for plotting
predictions <- data.frame( # Make a data frame of predicted values for the data
  Log_Mean_Converted_Value = rep(range_values, 1),
  log_total_TE = pglsModel_Lambda$coefficients[1] +  pglsModel_Lambda$coefficients[2] * range_values,
  Model = rep(c("PGLS Linear"), each = 100)
)

hemiptera = ggplot(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
  geom_point(alpha = 0.8, color = "#DC4D10") +
  geom_line(data = predictions, aes(x = Log_Mean_Converted_Value, y = log_total_TE), linewidth = 1) +
  geom_smooth(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) +
  labs(title = "d) Hemiptera", x = "Log Species Mass", y = "Log Proportion of TE") +
  theme_pubr() + 
  scale_y_continuous(breaks = c(-0.9, -0.75, -0.6)) 

########################################################
###### Hymenoptera #######
repeats_subset = subset(repeats, Order.x == "Hymenoptera")
phylo_tree_subset = keep.tip(phylo_tree, repeats_subset$Species[repeats_subset$Species %in% phylo_tree$tip.label]) # Subset tree to species in data
formula <- as.formula("log_total_TE ~ Log_Mean_Converted_Value") # Linear model
formula_sq <- as.formula("log_total_TE ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)") # Quadratic model
null_formula = as.formula(paste("log_total_TE ~ 1")) # Null model

linear_model <- lm(formula, data = repeats_subset) # Linear model
linear_model_residuals <- residuals(linear_model) # its residuals
names(linear_model_residuals) <- repeats_subset$Species # add species names to residuals

lambda = phylosig(phylo_tree_subset, linear_model_residuals, method = "lambda", test = TRUE); print(lambda) # Test phylogenetic signal

# # If it's proving difficult to find starting values that work, change the model below and test it
# lambdas = seq(0,1,0.01)
# for (lambda_starting in lambdas) { # Loop through each lambda_starting value in lambdas
#   tryCatch({  # Use tryCatch to suppress errors
#     pglsModel_Lambda <- gls(formula_sq, data = repeats_subset, correlation = corBlomberg(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML") # Attempt to create the pgls model
#     print(lambda_starting) # Print the lambda_starting value if successful
#   }, error = function(e) { # Suppress the error and continue
#   })
# }

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = repeats_subset, correlation = corPagel(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = repeats_subset, correlation = corBrownian(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = repeats_subset, correlation = corMartins(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = repeats_subset, correlation = corBlomberg(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = repeats_subset, correlation = corGrafen(1, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = repeats_subset, correlation = corPagel(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = repeats_subset, correlation = corBrownian(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = repeats_subset, correlation = corMartins(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = repeats_subset, correlation = corBlomberg(0.5, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = repeats_subset, correlation = corGrafen(1, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = repeats_subset, correlation = corBrownian(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = repeats_subset, correlation = corPagel(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = repeats_subset, correlation = corBlomberg(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = repeats_subset, correlation = corGrafen(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = repeats_subset, correlation = corMartins(1, phylo_tree_subset, form = ~Species), method = "ML")

AIC(linear_model,
    pglsModel_Lambda, pglsModel_Brownian, pglsModel_Martins, pglsModel_Grafen,  pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Martins_sq, pglsModel_Grafen_sq, pglsModel_Blomberg_sq,
    pglsModel_Null_Grafen, pglsModel_Null_Martins, pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Blomberg) # Best model is lambda null

intercept <- pglsModel_Null_Pagel$coefficients[1] 
hymenoptera = ggplot(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
  geom_point(alpha = 0.8, color = "#DC4D10") +
  geom_hline(yintercept = intercept, color = "black", linewidth = 1) +
  labs(title = "e) Hymenoptera", x = "Log Species Mass", y = "Log Proportion of TE") +
  theme_pubr()+ 
  scale_y_continuous(breaks = c(-0.5, -0.9, -1.3))

########################################################
###### Lepidoptera #######
repeats_subset = subset(repeats, Order.x == "Lepidoptera")
phylo_tree_subset = keep.tip(phylo_tree, repeats_subset$Species[repeats_subset$Species %in% phylo_tree$tip.label]) # Subset tree to species in data
formula <- as.formula("log_total_TE ~ Log_Mean_Converted_Value") # Linear model
formula_sq <- as.formula("log_total_TE ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)") # Quadratic model
null_formula = as.formula(paste("log_total_TE ~ 1")) # Null model

linear_model <- lm(formula, data = repeats_subset) # Linear model
linear_model_residuals <- residuals(linear_model) # its residuals
names(linear_model_residuals) <- repeats_subset$Species # add species names to residuals

lambda = phylosig(phylo_tree_subset, linear_model_residuals, method = "lambda", test = TRUE); print(lambda) # Test phylogenetic signal

# # If it's proving difficult to find starting values that work, change the model below and test it
# lambdas = seq(0,1,0.01)
# for (lambda_starting in lambdas) { # Loop through each lambda_starting value in lambdas
#   tryCatch({  # Use tryCatch to suppress errors
#     pglsModel_Lambda <- gls(formula_sq, data = repeats_subset, correlation = corBlomberg(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML") # Attempt to create the pgls model
#     print(lambda_starting) # Print the lambda_starting value if successful
#   }, error = function(e) { # Suppress the error and continue
#   })
# }

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = repeats_subset, correlation = corPagel(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = repeats_subset, correlation = corBrownian(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = repeats_subset, correlation = corMartins(0.01, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = repeats_subset, correlation = corBlomberg(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = repeats_subset, correlation = corGrafen(0.1, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = repeats_subset, correlation = corPagel(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = repeats_subset, correlation = corBrownian(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = repeats_subset, correlation = corMartins(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = repeats_subset, correlation = corBlomberg(0.09, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = repeats_subset, correlation = corGrafen(0.1, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = repeats_subset, correlation = corBrownian(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = repeats_subset, correlation = corPagel(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = repeats_subset, correlation = corBlomberg(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = repeats_subset, correlation = corGrafen(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = repeats_subset, correlation = corMartins(0.01, phylo_tree_subset, form = ~Species), method = "ML")

AIC(linear_model,
    pglsModel_Lambda, pglsModel_Brownian, pglsModel_Martins, pglsModel_Grafen,  pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Martins_sq, pglsModel_Grafen_sq, pglsModel_Blomberg_sq,
    pglsModel_Null_Grafen, pglsModel_Null_Martins, pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Blomberg) # Best model is lambda null

# Define the range of 'Log_Mean_Converted_Value'
intercept <- pglsModel_Null_Pagel$coefficients[1]
lepidoptera = ggplot(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
  geom_point(alpha = 0.8, color = "#DC4D10") +
  geom_hline(yintercept = intercept, color = "black", linewidth = 1) +
  labs(title = "f) Lepidoptera", x = "Log Species Mass", y = "Log Proportion of TE") +
  theme_pubr()+ 
  scale_y_continuous(breaks = c(-1.1, -0.8, -0.5)) 

###########################################
#### Combine plots ####
plot_theme_top <- theme(axis.title.x = element_blank(), 
                        axis.title.y = element_blank(), 
                        plot.title = element_text(family = "lm-bold", size = 19),
                        text = element_text(family = "lm", size = 14), 
                        axis.line.x = element_line(color = "black", size = 1, lineend = "square"), 
                        axis.line.y = element_line(color = "black", size = 1, lineend = "square"), 
                        axis.ticks.x = element_line(color = "black", size = 0.5), 
                        axis.ticks.y = element_line(color = "black", size = 0.5),
                        plot.margin = margin(b = 50, l = 10),
                        axis.text.x = element_text(size = 17),  # Set x-axis text size to 14
                        axis.text.y = element_text(size = 17, angle = 90, hjust = 0.5)  # Set y-axis text size to 14
) 
plot_theme_bottom <- theme(axis.title.x = element_blank(), 
                           axis.title.y = element_blank(), 
                           plot.title = element_text(family = "lm-bold", size = 19),
                           text = element_text(family = "lm", size = 14), 
                           axis.line.x = element_line(color = "black", size = 1, lineend = "square"), 
                           axis.line.y = element_line(color = "black", size = 1, lineend = "square"), 
                           axis.ticks.x = element_line(color = "black", size = 0.5), 
                           axis.ticks.y = element_line(color = "black", size = 0.5),
                           plot.margin = margin(b = 25, l = 10),
                           axis.text.x = element_text(size = 17),  # Set x-axis text size to 14
                           axis.text.y = element_text(size = 17, angle = 90, hjust = 0.5),  # Set y-axis text size to 14
)

font_add(family = "lm", regular = "/usr/share/texlive/texmf-dist/fonts/opentype/public/lm/lmroman10-regular.otf")
font_add(family = "lm-bold", regular = "/usr/share/texmf/fonts/opentype/public/lm/lmroman10-bold.otf")
showtext_auto()

insecta_combined <- insecta + plot_theme_top
coleoptera_combined <- coleoptera + plot_theme_top 
diptera_combined <- diptera + plot_theme_top
hemiptera_combined <- hemiptera + plot_theme_bottom
hymenoptera_combined <- hymenoptera + plot_theme_bottom
lepidoptera_combined <- lepidoptera + plot_theme_bottom

# Combine the plots into a single layout
combined_plot <- egg::ggarrange(insecta_combined, coleoptera_combined, diptera_combined, 
                                hemiptera_combined, hymenoptera_combined, lepidoptera_combined, ncol = 3)

# Combine the plot with the legend and add annotations
final_plot <- plot_grid(combined_plot) +
  annotate("text", label = "Log Species Mass", x = 0.5, y = 0.02, vjust = 1.75, size = 7.5, family = "lm-bold") + 
  annotate("text", label = "Log Proportion of TE per Species", x = 0.02, y = 0.5, vjust = -2, angle = 90, size = 7.5, family = "lm-bold") +
  theme(plot.margin = margin(l = 30, r = 10, t = 10, b = 20)) + 
  theme(text = element_text(family = "lm"))

# Display the final plot
print(final_plot)


