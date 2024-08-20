library(phytools)
library(rotl)
library(OpenTreeChronograms)
library(ggtree)
library(geiger)
library(nlme)
library(ggpubr)
library(dplyr)
library(rr2)
library(grDevices)
library(evobiR)
library(PKNCA)
library(parallel)
library(grid)
library(gridExtra)
library(patchwork)
library(sysfonts)
library(showtext)
library(grid)
library(cowplot)

coleoptera_summary_table = read.csv("../data/toga/summary_table_coleoptera.csv") # Load TOGA output data 
diptera_summary_table = read.csv("../data/toga/summary_table_diptera.csv") # Load TOGA output data 
diptera_summary_table <- subset(diptera_summary_table, !Accession %in% c("GCA_007989325.1", "GCA_030788265.1")) # Remove two problematic species
lepidoptera_summary_table = read.csv("../data/toga/summary_table_lepidoptera.csv") # Load TOGA output data
hymenoptera_summary_table = read.csv("../data/toga/summary_table_hymenoptera.csv") # Load TOGA output data
calibrated_phylogeny = read.tree("../data/phylogeny/calibrated_tree.nwk"); calibrated_phylogeny$tip.label = sub("_", " ", calibrated_phylogeny$tip.label) # Load entire phylogeny and fix naming
masses = read.csv("../data/mass/Insect_Masses_and_Genomes.csv")[,-1] # Load insect masses without ID column
masses$Accession <- gsub("^GCA_|^GCF_|\\.\\d+$", "", masses$Accession) # Remove GCA and GCF, as these changed when pulled from NCBI
masses$Log_Mean_Converted_Value = log10(masses$Mean_Converted_Value) # Create column for log-mass

###### Coleoptera #######
summary_table = coleoptera_summary_table
phylogeny = calibrated_phylogeny
mass = masses
reference_species = "Diorhabda.carinulata"
reference_accession = "GCF_026250575.1"

if (!reference_accession %in% summary_table$Accession){
  new_row <- data.frame(Accession = reference_accession, matrix(0, nrow = 1, ncol = ncol(summary_table) - 1))
  colnames(new_row)[2:ncol(new_row)] <- colnames(summary_table)[2:ncol(summary_table)]
  summary_table <- rbind(summary_table, new_row)
}

summary_table = summary_table %>% # Calculate summaries and combine with masses df
  mutate(
    Accession = gsub("^GCA_|^GCF_|\\.\\d+$", "", Accession),
    Log_one2zero_Ratio = log10(one2zero/(one2one+one2zero+many2one+one2many+many2many)),
    Log_one2zero = log10(one2zero),
    I_and_PI = I + PI,
    Log_I_and_PI = log10(I_and_PI),
    Likely_lost = L + M + PG + PM + UL,
    Log_Likely_lost = log10(Likely_lost),
    Ratio = log10(I_and_PI / (I_and_PI + Likely_lost))) %>%
  inner_join(masses, by = "Accession")

tree = keep.tip(calibrated_phylogeny, summary_table$Species[summary_table$Species %in% calibrated_phylogeny$tip.label]) # Prune phylogeny

distance_matrix = data.frame(cophenetic.phylo(tree)) # Calculate pairwise distance matrix
distance_matrix$Species = rownames(distance_matrix) # Set a column which contains species name
distance_matrix = distance_matrix[, c("Species", reference_species)] # Keep only species names and their distance to reference 
distance_matrix = filter(distance_matrix, distance_matrix[,2] > 0) # Also remove reference distance to itself

rownames(distance_matrix) = NULL; colnames(distance_matrix) = c("Species", "Distance_To_Reference") # Set row and column names
distance_matrix$Log_Distance_To_Reference = log10(distance_matrix$Distance_To_Reference)

summary_table = merge(summary_table, distance_matrix, by = "Species") # Merge working data set and distance matrix data
summary_table = summary_table[!apply(summary_table[, c("Likely_lost", "I_and_PI")], 1, function(x) sum(x) == 0), ] # Remove rows which have failed TOGA runs

summary_table <- summary_table[summary_table$Log_Distance_To_Reference <= quantile(summary_table$Log_Distance_To_Reference, 0.75), ] # Filter out rows that are over the 75th percentile

pruned_tree = keep.tip(tree, summary_table$Species[summary_table$Species %in% tree$tip.label]) 

range_values <- seq(min(summary_table$Log_Mean_Converted_Value), max(summary_table$Log_Mean_Converted_Value), length.out = 100) # Creates a sequence of 100 equally spaced values between the minimum and maximum Log_Mean_Converted_Value
fixed_distances <- mean(summary_table$Log_Distance_To_Reference)

# If it's proving difficult to find starting values that work, change the model below and test it
# lambdas = seq(0,1,0.01)
# for (lambda_starting in lambdas) { # Loop through each lambda_starting value in lambdas
#   tryCatch({  # Use tryCatch to suppress errors
#     pglsModel_Lambda <- gls(formula, data = summary_table, correlation = corPagel(lambda_starting, pruned_tree, form = ~Species), method = "ML") # Attempt to create the pgls model
#     print(lambda_starting) # Print the lambda_starting value if successful
#   }, error = function(e) { # Suppress the error and continue
#   })
# }

### PGLS Models for Log_I_and_PI ###
formula <- as.formula(paste("Log_I_and_PI ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_I_and_PI ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_I_and_PI ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(0.5, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(0.5, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(0.5, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Grafen Quadratic is best

### PGLS Models for Log_one2zero ###
formula <- as.formula(paste("Log_one2zero ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_one2zero ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_one2zero ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(0.5, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(0.5, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(0.5, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Grafen Quadratic is lowest AIC, but Brownian_sq has best residuals distribution - check this yourself

### PGLS Models for Log_one2zero_Ratio ###
formula <- as.formula(paste("Log_one2zero_Ratio ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_one2zero_Ratio ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_one2zero_Ratio ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(0.5, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(0.5, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(0.5, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Brownian Quadratic is lowest AIC

# Plot the best fitting models
pglsModel_Grafen_sq = gls(Log_I_and_PI ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq_one2zero = gls(Log_one2zero ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(Log_one2zero_Ratio ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")

range_values <- seq(min(summary_table$Log_Mean_Converted_Value), max(summary_table$Log_Mean_Converted_Value), length.out = 100)
fixed_distance <- mean(summary_table$Log_Distance_To_Reference)
new_data <- data.frame(
  Log_Mean_Converted_Value = range_values, 
  Log_Distance_To_Reference = fixed_distance
)
new_data$Log_I_and_PI <- predict(pglsModel_Grafen_sq, newdata = new_data)
new_data$Log_one2zero <- predict(pglsModel_Brownian_sq_one2zero, newdata = new_data)
new_data$Log_one2zero_Ratio <- predict(pglsModel_Brownian_sq, newdata = new_data)


distance_label <- "Mean Distance" # Creates labels for the fixed distances
new_data$Distance_Label <- factor(fixed_distance, labels = distance_label)


p_text = paste0("p = ", round(coef(summary(pglsModel_Grafen_sq))[3,4],3))

r_squared_value <- round(rsquared.gls(pglsModel_Grafen_sq)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_i_and_pi_coleoptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  labs(x = "Log Species Mass", y = "Log I/PI Genes") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))+
  annotate("text", x = -Inf, y = Inf, 
           label = p_text, hjust = -0.05, vjust = 5.5, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = Inf, 
           label = r_text, hjust = -0.08, vjust = 3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = Inf, 
           label = "Grafen Q.", hjust = -0.03, vjust = 2.75, size = 5, family = "lm-bold")

p_text = paste0("p = ", round(coef(summary(pglsModel_Brownian_sq_one2zero))[3,4],3))
r_squared_value <- round(rsquared.gls(pglsModel_Brownian_sq_one2zero)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_one2zero_coleoptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_one2zero, linetype = Distance_Label, group = Distance_Label), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Distance Label") +  # Discrete linetype scale for lines
  labs(x = "Log Species Mass", y = "Log One-to-None Genes") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14)) + guides(color = "none", linetype = "none") +
  annotate("text", x = -Inf, y = -Inf, 
           label = p_text, hjust = -0.05, vjust = -2.75, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = -Inf, 
           label = r_text, hjust = -0.08, vjust = -3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = -Inf, 
           label = "Brownian Q.", hjust = -0.03, vjust = -5.6, size = 5, family = "lm-bold")

p_text = paste0("p = ", round(coef(summary(pglsModel_Brownian_sq))[3,4],3))
r_squared_value <- round(rsquared.gls(pglsModel_Brownian_sq)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_one2zero_ratio_coleoptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio, linetype = Distance_Label, group = Distance_Label), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Distance Label") +  # Discrete linetype scale for lines
  labs(x = "Log Species Mass", y = "Log One-to-None Ratio") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))+ guides(color = "none", linetype = "none") +
  annotate("text", x = -Inf, y = -Inf, 
           label = p_text, hjust = -0.05, vjust = -2.75, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = -Inf, 
           label = r_text, hjust = -0.08, vjust = -3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = -Inf, 
           label = "Brownian Q.", hjust = -0.03, vjust = -5.6, size = 5, family = "lm-bold")

############################################################################################################################################################################################################################



###### Diptera #######
summary_table = diptera_summary_table
phylogeny = calibrated_phylogeny
mass = masses
reference_species = "Drosophila.melanogaster"
reference_accession = "GCF_000001215.4"

if (!reference_accession %in% summary_table$Accession){
  new_row <- data.frame(Accession = reference_accession, matrix(0, nrow = 1, ncol = ncol(summary_table) - 1))
  colnames(new_row)[2:ncol(new_row)] <- colnames(summary_table)[2:ncol(summary_table)]
  summary_table <- rbind(summary_table, new_row)
}

summary_table = summary_table %>% # Calculate summaries and combine with masses df
  mutate(
    Accession = gsub("^GCA_|^GCF_|\\.\\d+$", "", Accession),
    Log_one2zero_Ratio = log10(one2zero/(one2one+one2zero+many2one+one2many+many2many)),
    Log_one2zero = log10(one2zero),
    I_and_PI = I + PI,
    Log_I_and_PI = log10(I_and_PI),
    Likely_lost = L + M + PG + PM + UL,
    Log_Likely_lost = log10(Likely_lost),
    Ratio = log10(I_and_PI / (I_and_PI + Likely_lost))) %>%
  inner_join(masses, by = "Accession")

tree = keep.tip(calibrated_phylogeny, summary_table$Species[summary_table$Species %in% calibrated_phylogeny$tip.label]) # Prune phylogeny

distance_matrix = data.frame(cophenetic.phylo(tree)) # Calculate pairwise distance matrix
distance_matrix$Species = rownames(distance_matrix) # Set a column which contains species name
distance_matrix = distance_matrix[, c("Species", reference_species)] # Keep only species names and their distance to reference 
distance_matrix = filter(distance_matrix, distance_matrix[,2] > 0) # Also remove reference distance to itself

rownames(distance_matrix) = NULL; colnames(distance_matrix) = c("Species", "Distance_To_Reference") # Set row and column names
distance_matrix$Log_Distance_To_Reference = log10(distance_matrix$Distance_To_Reference)

summary_table = merge(summary_table, distance_matrix, by = "Species") # Merge working data set and distance matrix data
summary_table = summary_table[!apply(summary_table[, c("Likely_lost", "I_and_PI")], 1, function(x) sum(x) == 0), ] # Remove rows which have failed TOGA runs

summary_table <- summary_table[summary_table$Log_Distance_To_Reference <= quantile(summary_table$Log_Distance_To_Reference, 0.75), ] # Filter out rows that are over the 75th percentile

pruned_tree = keep.tip(tree, summary_table$Species[summary_table$Species %in% tree$tip.label]) 

range_values <- seq(min(summary_table$Log_Mean_Converted_Value), max(summary_table$Log_Mean_Converted_Value), length.out = 100) # Creates a sequence of 100 equally spaced values between the minimum and maximum Log_Mean_Converted_Value
fixed_distances <- c(mean(summary_table$Log_Distance_To_Reference) - sd(summary_table$Log_Distance_To_Reference), # Creates three fixed values: mean - 1 SD, mean, and mean + 1 SD of Log_Distance_To_Reference
                     mean(summary_table$Log_Distance_To_Reference),
                     mean(summary_table$Log_Distance_To_Reference) + sd(summary_table$Log_Distance_To_Reference))

# If it's proving difficult to find starting values that work, change the model below and test it
# lambdas = seq(0,1,0.0001)
# for (lambda_starting in lambdas) { # Loop through each lambda_starting value in lambdas
#   tryCatch({  # Use tryCatch to suppress errors
#     pglsModel_Lambda <- gls(formula_sq, data = summary_table, correlation = corBlomberg(lambda_starting, pruned_tree, form = ~Species), method = "ML") # Attempt to create the pgls model
#     print(lambda_starting) # Print the lambda_starting value if successful
#   }, error = function(e) { # Suppress the error and continue
#   })
# }

### PGLS Models for Log_I_and_PI ###
formula <- as.formula(paste("Log_I_and_PI ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_I_and_PI ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_I_and_PI ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(0.4, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(0.4, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(0.4, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(0.4, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(0.4, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(0.4, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(0.4, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(0.0069, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(0.4, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(0.4, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(0.4, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(0.0069, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(0.4, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(0.4, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Blomberg Quadratic is best but it is completely overfitting to the data, as is Lambda Quadratic - Residuals of Lambda Null model is probably best here

### PGLS Models for Log_one2zero ###
formula <- as.formula(paste("Log_one2zero ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_one2zero ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_one2zero ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(0.01, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(0.009, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(0.01, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(0.01, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Lambda Linear and Null have close AIC's - but residuals of linear plot better

### PGLS Models for Log_one2zero_Ratio ###
formula <- as.formula(paste("Log_one2zero_Ratio ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_one2zero_Ratio ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_one2zero_Ratio ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(0, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(0.5, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(0, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(0.01, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(0.7, pruned_tree, form = ~Species), method = "ML") # high starting value but one of few that fits
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(0, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(0, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(0.01, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(0, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(0.14, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(0.0069, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(0.01, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Lambda Linear is lowest AIC: -406

# Plot the best fitting models
pglsModel_Null_Pagel = gls(Log_I_and_PI  ~ Log_Distance_To_Reference, data = summary_table, correlation = corPagel(0.0001, pruned_tree, form = ~Species), method = "ML")
pglsModel_Lambda_one2zero = gls(Log_one2zero  ~ Log_Mean_Converted_Value + Log_Distance_To_Reference, data = summary_table, correlation = corPagel(0.01, pruned_tree, form = ~Species), method = "ML")
pglsModel_Lambda = gls(Log_one2zero_Ratio  ~ Log_Mean_Converted_Value + Log_Distance_To_Reference, data = summary_table, correlation = corPagel(0, pruned_tree, form = ~Species), method = "ML")

range_values <- seq(min(summary_table$Log_Mean_Converted_Value), max(summary_table$Log_Mean_Converted_Value), length.out = 100)
fixed_distance <- mean(summary_table$Log_Distance_To_Reference)
new_data <- data.frame(
  Log_Mean_Converted_Value = range_values, 
  Log_Distance_To_Reference = fixed_distance
)
new_data$Log_I_and_PI <- predict(pglsModel_Null_Pagel, newdata = new_data)
new_data$Log_one2zero <- predict(pglsModel_Lambda_one2zero, newdata = new_data)
new_data$Log_one2zero_Ratio <- predict(pglsModel_Lambda, newdata = new_data)


distance_label <- "Mean Distance" # Creates labels for the fixed distances
new_data$Distance_Label <- factor(fixed_distance, labels = distance_label)

p_text = paste0("p = ", round(coef(summary(pglsModel_Null_Pagel))[2,4],3))
if (round(coef(summary(pglsModel_Null_Pagel)))[2,4] < 0.001){
  p_text = paste0("p < 0.001")
}

r_squared_value <- round(rsquared.gls(pglsModel_Null_Pagel)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_i_and_pi_diptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  labs(x = "Log Species Mass", y = "Log I/PI Genes") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))+
  annotate("text", x = -Inf, y = Inf, 
           label = p_text, hjust = -0.05, vjust = 5.5, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = Inf, 
           label = r_text, hjust = -0.08, vjust = 3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = Inf, 
           label = "Lambda N.", hjust = -0.03, vjust = 2.75, size = 5, family = "lm-bold")

p_text = paste0("p = ", round(coef(summary(pglsModel_Lambda_one2zero))[2,4],3))
r_squared_value <- round(rsquared.gls(pglsModel_Lambda_one2zero)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_one2zero_diptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_one2zero, linetype = Distance_Label, group = Distance_Label), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Distance Label") +  # Discrete linetype scale for lines
  labs(x = "Log Species Mass", y = "Log One-to-None Genes") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14)) + guides(color = "none", linetype = "none") +
  annotate("text", x = -Inf, y = -Inf, 
           label = p_text, hjust = -0.05, vjust = -2.75, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = -Inf, 
           label = r_text, hjust = -0.08, vjust = -3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = -Inf, 
           label = "Lambda L.", hjust = -0.03, vjust = -5.6, size = 5, family = "lm-bold")

p_text = paste0("p = ", round(coef(summary(pglsModel_Lambda))[2,4],3))
r_squared_value <- round(rsquared.gls(pglsModel_Lambda)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_one2zero_ratio_diptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio, linetype = Distance_Label, group = Distance_Label), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Distance Label") +  # Discrete linetype scale for lines
  labs(x = "Log Species Mass", y = "Log One-to-None Ratio") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))+ guides(color = "none", linetype = "none") +
  annotate("text", x = -Inf, y = -Inf, 
           label = p_text, hjust = -0.05, vjust = -2.75, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = -Inf, 
           label = r_text, hjust = -0.1, vjust = -3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = -Inf, 
           label = "Lambda L.", hjust = -0.03, vjust = -5.6, size = 5, family = "lm-bold")


############################################################################################################################################################################################################################



###### Lepidoptera #######
summary_table = lepidoptera_summary_table
phylogeny = calibrated_phylogeny
mass = masses
reference_species = "Helicoverpa.armigera"
reference_accession = "GCF_030705265.1"

if (!reference_accession %in% summary_table$Accession){
  new_row <- data.frame(Accession = reference_accession, matrix(0, nrow = 1, ncol = ncol(summary_table) - 1))
  colnames(new_row)[2:ncol(new_row)] <- colnames(summary_table)[2:ncol(summary_table)]
  summary_table <- rbind(summary_table, new_row)
}

summary_table = summary_table %>% # Calculate summaries and combine with masses df
  mutate(
    Accession = gsub("^GCA_|^GCF_|\\.\\d+$", "", Accession),
    Log_one2zero_Ratio = log10(one2zero/(one2one+one2zero+many2one+one2many+many2many)),
    Log_one2zero = log10(one2zero),
    I_and_PI = I + PI,
    Log_I_and_PI = log10(I_and_PI),
    Likely_lost = L + M + PG + PM + UL,
    Log_Likely_lost = log10(Likely_lost),
    Ratio = log10(I_and_PI / (I_and_PI + Likely_lost))) %>%
  inner_join(masses, by = "Accession")

tree = keep.tip(calibrated_phylogeny, summary_table$Species[summary_table$Species %in% calibrated_phylogeny$tip.label]) # Prune phylogeny

distance_matrix = data.frame(cophenetic.phylo(tree)) # Calculate pairwise distance matrix
distance_matrix$Species = rownames(distance_matrix) # Set a column which contains species name
distance_matrix = distance_matrix[, c("Species", reference_species)] # Keep only species names and their distance to reference 
distance_matrix = filter(distance_matrix, distance_matrix[,2] > 0) # Also remove reference distance to itself

rownames(distance_matrix) = NULL; colnames(distance_matrix) = c("Species", "Distance_To_Reference") # Set row and column names
distance_matrix$Log_Distance_To_Reference = log10(distance_matrix$Distance_To_Reference)

summary_table = merge(summary_table, distance_matrix, by = "Species") # Merge working data set and distance matrix data
summary_table = summary_table[!apply(summary_table[, c("Likely_lost", "I_and_PI")], 1, function(x) sum(x) == 0), ] # Remove rows which have failed TOGA runs

summary_table <- summary_table[summary_table$Log_Distance_To_Reference <= quantile(summary_table$Log_Distance_To_Reference, 0.75), ] # Filter out rows that are over the 75th percentile

pruned_tree = keep.tip(tree, summary_table$Species[summary_table$Species %in% tree$tip.label]) 

range_values <- seq(min(summary_table$Log_Mean_Converted_Value), max(summary_table$Log_Mean_Converted_Value), length.out = 100) # Creates a sequence of 100 equally spaced values between the minimum and maximum Log_Mean_Converted_Value
fixed_distances <- c(mean(summary_table$Log_Distance_To_Reference) - sd(summary_table$Log_Distance_To_Reference), # Creates three fixed values: mean - 1 SD, mean, and mean + 1 SD of Log_Distance_To_Reference
                     mean(summary_table$Log_Distance_To_Reference),
                     mean(summary_table$Log_Distance_To_Reference) + sd(summary_table$Log_Distance_To_Reference))

# If it's proving difficult to find starting values that work, change the model below and test it
# lambdas = seq(0,1,0.01)
# for (lambda_starting in lambdas) { # Loop through each lambda_starting value in lambdas
#   tryCatch({  # Use tryCatch to suppress errors
#     pglsModel<- gls(formula, data = summary_table, correlation = corMartins(lambda_starting, pruned_tree, form = ~Species), method = "ML") # Attempt to create the pgls model
#     print(lambda_starting) # Print the lambda_starting value if successful
#   }, error = function(e) { # Suppress the error and continue
#   })
# }

### PGLS Models for Log_I_and_PI ###
formula <- as.formula(paste("Log_I_and_PI ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_I_and_PI ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_I_and_PI ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(0, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(0, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(0, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Blomberg Null is best 

### PGLS Models for Log_one2zero ###
formula <- as.formula(paste("Log_one2zero ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_one2zero ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_one2zero ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(0, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(0, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(0, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Blomberg Linear and Null have close AIC's - but residuals of linear plot better

### PGLS Models for Log_one2zero_Ratio ###
formula <- as.formula(paste("Log_one2zero_Ratio ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_one2zero_Ratio ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_one2zero_Ratio ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(0, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML") # high starting value but one of few that fits
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(0, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(0, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Blomberg Linear is lowest AIC but Null fit bests


pglsModel_Null_Blomberg = gls(Log_I_and_PI  ~ Log_Distance_To_Reference, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_one2zero = gls(Log_one2zero  ~ Log_Mean_Converted_Value + Log_Distance_To_Reference, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg_one2zero_ratio = gls(Log_one2zero_Ratio  ~ Log_Distance_To_Reference, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")

range_values <- seq(min(summary_table$Log_Mean_Converted_Value), max(summary_table$Log_Mean_Converted_Value), length.out = 100)
fixed_distance <- mean(summary_table$Log_Distance_To_Reference)
new_data <- data.frame(
  Log_Mean_Converted_Value = range_values, 
  Log_Distance_To_Reference = fixed_distance
)
new_data$Log_I_and_PI <- predict(pglsModel_Null_Blomberg, newdata = new_data)
new_data$Log_one2zero <- predict(pglsModel_Blomberg_one2zero, newdata = new_data)
new_data$Log_one2zero_Ratio <- predict(pglsModel_Null_Blomberg_one2zero_ratio, newdata = new_data)


distance_label <- "Mean Distance" # Creates labels for the fixed distances
new_data$Distance_Label <- factor(fixed_distance, labels = distance_label)

p_text = paste0("p = ", round(coef(summary(pglsModel_Null_Blomberg))[2,4],3))
if (round(coef(summary(pglsModel_Null_Blomberg)))[2,4] < 0.001){
  p_text = paste0("p < 0.001")
}

r_squared_value <- round(rsquared.gls(pglsModel_Null_Blomberg)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_i_and_pi_lepidoptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  labs(x = "Log Species Mass", y = "Log I/PI Genes") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))+
  annotate("text", x = -Inf, y = Inf, 
           label = p_text, hjust = -0.05, vjust = 6.2, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = Inf, 
           label = r_text, hjust = -0.08, vjust = 3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = Inf, 
           label = "Blomberg N.", hjust = -0.03, vjust = 2.75, size = 5, family = "lm-bold")

p_text = paste0("p = ", round(coef(summary(pglsModel_Blomberg_one2zero))[2,4],3))
if (round(coef(summary(pglsModel_Blomberg_one2zero)))[2,4] < 0.001){
  p_text = paste0("p < 0.001")
}
r_squared_value <- round(rsquared.gls(pglsModel_Blomberg_one2zero)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_one2zero_lepidoptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_one2zero, linetype = Distance_Label, group = Distance_Label), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Distance Label") +  # Discrete linetype scale for lines
  labs(x = "Log Species Mass", y = "Log One-to-None Genes") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14)) + guides(color = "none", linetype = "none")+
  annotate("text", x = -Inf, y = -Inf, 
           label = p_text, hjust = -0.05, vjust = -2.75, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = -Inf, 
           label = r_text, hjust = -0.08, vjust = -3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = -Inf, 
           label = "Blomberg L.", hjust = -0.03, vjust = -5.6, size = 5, family = "lm-bold")

p_text = paste0("p = ", round(coef(summary(pglsModel_Null_Blomberg_one2zero_ratio))[2,4],3))
if (round(coef(summary(pglsModel_Null_Blomberg_one2zero_ratio)))[2,4] < 0.001){
  p_text = paste0("p < 0.001")
}
r_squared_value <- round(rsquared.gls(pglsModel_Null_Blomberg_one2zero_ratio)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_one2zero_ratio_lepidoptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio, linetype = Distance_Label, group = Distance_Label), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Distance Label") +  # Discrete linetype scale for lines
  labs(x = "Log Species Mass", y = "Log One-to-None Ratio") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))+ guides(color = "none", linetype = "none")+
  annotate("text", x = -Inf, y = -Inf, 
           label = p_text, hjust = -0.05, vjust = -2.75, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = -Inf, 
           label = r_text, hjust = -0.08, vjust = -3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = -Inf, 
           label = "Blomberg N.", hjust = -0.03, vjust = -5.6, size = 5, family = "lm-bold")


############################################################################################################################################################################################################################
###### Hymenoptera #######
summary_table = hymenoptera_summary_table
phylogeny = calibrated_phylogeny
mass = masses
reference_species = "Apis.cerana"
reference_accession = "GCF_029169275.1"

if (!reference_accession %in% summary_table$Accession){
  new_row <- data.frame(Accession = reference_accession, matrix(0, nrow = 1, ncol = ncol(summary_table) - 1))
  colnames(new_row)[2:ncol(new_row)] <- colnames(summary_table)[2:ncol(summary_table)]
  summary_table <- rbind(summary_table, new_row)
}

summary_table = summary_table %>% # Calculate summaries and combine with masses df
  mutate(
    Accession = gsub("^GCA_|^GCF_|\\.\\d+$", "", Accession),
    Log_one2zero_Ratio = log10(one2zero/(one2one+one2zero+many2one+one2many+many2many)),
    Log_one2zero = log10(one2zero),
    I_and_PI = I + PI,
    Log_I_and_PI = log10(I_and_PI),
    Likely_lost = L + M + PG + PM + UL,
    Log_Likely_lost = log10(Likely_lost),
    Ratio = log10(I_and_PI / (I_and_PI + Likely_lost))) %>%
  inner_join(masses, by = "Accession")

tree = keep.tip(calibrated_phylogeny, summary_table$Species[summary_table$Species %in% calibrated_phylogeny$tip.label]) # Prune phylogeny

distance_matrix = data.frame(cophenetic.phylo(tree)) # Calculate pairwise distance matrix
distance_matrix$Species = rownames(distance_matrix) # Set a column which contains species name
distance_matrix = distance_matrix[, c("Species", reference_species)] # Keep only species names and their distance to reference 
distance_matrix = filter(distance_matrix, distance_matrix[,2] > 0) # Also remove reference distance to itself

rownames(distance_matrix) = NULL; colnames(distance_matrix) = c("Species", "Distance_To_Reference") # Set row and column names
distance_matrix$Log_Distance_To_Reference = log10(distance_matrix$Distance_To_Reference)

summary_table = merge(summary_table, distance_matrix, by = "Species") # Merge working data set and distance matrix data
summary_table = summary_table[!apply(summary_table[, c("Likely_lost", "I_and_PI")], 1, function(x) sum(x) == 0), ] # Remove rows which have failed TOGA runs

summary_table <- summary_table[summary_table$Log_Distance_To_Reference <= quantile(summary_table$Log_Distance_To_Reference, 0.75), ] # Filter out rows that are over the 75th percentile

pruned_tree = keep.tip(tree, summary_table$Species[summary_table$Species %in% tree$tip.label]) 

results <- data.frame(Gene_Variables = character(), Best_Model = character(), P_Value = numeric(), lambda = numeric(), R_squared = numeric(), stringsAsFactors = FALSE)  # Prepare to store results

gene_variables <- c("Log_one2zero_Ratio", "Log_I_and_PI", "Log_one2zero")

range_values <- seq(min(summary_table$Log_Mean_Converted_Value), max(summary_table$Log_Mean_Converted_Value), length.out = 100) # Creates a sequence of 100 equally spaced values between the minimum and maximum Log_Mean_Converted_Value
fixed_distances <- c(mean(summary_table$Log_Distance_To_Reference) - sd(summary_table$Log_Distance_To_Reference), # Creates three fixed values: mean - 1 SD, mean, and mean + 1 SD of Log_Distance_To_Reference
                     mean(summary_table$Log_Distance_To_Reference),
                     mean(summary_table$Log_Distance_To_Reference) + sd(summary_table$Log_Distance_To_Reference))

# If it's proving difficult to find starting values that work, change the model below and test it
# lambdas = seq(0,1,0.01)
# for (lambda_starting in lambdas) { # Loop through each lambda_starting value in lambdas
#   tryCatch({  # Use tryCatch to suppress errors
#     pglsModel<- gls(formula, data = summary_table, correlation = corMartins(lambda_starting, pruned_tree, form = ~Species), method = "ML") # Attempt to create the pgls model
#     print(lambda_starting) # Print the lambda_starting value if successful
#   }, error = function(e) { # Suppress the error and continue
#   })
# }

### PGLS Models for Log_I_and_PI ###
formula <- as.formula(paste("Log_I_and_PI ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_I_and_PI ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_I_and_PI ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(1, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(1, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Brownian Null is best 

### PGLS Models for Log_one2zero ###
formula <- as.formula(paste("Log_one2zero ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_one2zero ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_one2zero ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(0.9, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(0.9, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(0.9, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Very close AICs but Lambda Linear residuals are best

### PGLS Models for Log_one2zero_Ratio ###
formula <- as.formula(paste("Log_one2zero_Ratio ~ Log_Mean_Converted_Value + Log_Distance_To_Reference"))
formula_sq <- as.formula(paste("Log_one2zero_Ratio ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE) + Log_Distance_To_Reference"))
null_formula = as.formula(paste("Log_one2zero_Ratio ~ Log_Distance_To_Reference"))

linear_model = lm(formula, data = summary_table) # Linear model
linear_model_residuals = linear_model$residuals # Extract model residuals

names(linear_model_residuals) = summary_table$Species # Change residual names to species names to match in phylosig
linear_model_residuals = ReorderData(pruned_tree, linear_model_residuals, taxa.names="row names") # Reorder residual names to match phylogeny
lambda_test = phylosig(pruned_tree, linear_model_residuals, method = "lambda", test = TRUE); print(lambda_test[1]) # If lambda is high (near 1) and significant, proceed with PGLS and test for correlation structures

# 0.5 is a good approximate starting value - adjust accordingly to lambda test above
pglsModel_Lambda = gls(formula, data = summary_table, correlation = corPagel(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = summary_table, correlation = corBrownian(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = summary_table, correlation = corMartins(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = summary_table, correlation = corBlomberg(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = summary_table, correlation = corGrafen(0.9, pruned_tree, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = summary_table, correlation = corPagel(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = summary_table, correlation = corBrownian(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = summary_table, correlation = corMartins(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = summary_table, correlation = corBlomberg(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = summary_table, correlation = corGrafen(0.9, pruned_tree, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = summary_table, correlation = corBrownian(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = summary_table, correlation = corPagel(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = summary_table, correlation = corBlomberg(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = summary_table, correlation = corGrafen(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = summary_table, correlation = corMartins(0.9, pruned_tree, form = ~Species), method = "ML")

AIC(linear_model, pglsModel_Lambda, pglsModel_Brownian,  pglsModel_Grafen, pglsModel_Martins, pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Grafen_sq, pglsModel_Martins_sq, pglsModel_Blomberg_sq, 
    pglsModel_Null_Grafen,  pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Martins, pglsModel_Null_Blomberg) # Blomberg Null fits best - check residuals


pglsModel_Null_Brownian = gls(Log_I_and_PI  ~ Log_Distance_To_Reference, data = summary_table, correlation = corBrownian(1, pruned_tree, form = ~Species), method = "ML") # AIC between all models 
pglsModel_Lambda_one2zero = gls(Log_one2zero  ~ Log_Mean_Converted_Value + Log_Distance_To_Reference, data = summary_table, correlation = corPagel(0.9, pruned_tree, form = ~Species), method = "ML")
pglsModel_Null_Blomberg_one2zero_ratio = gls(Log_one2zero_Ratio  ~ Log_Distance_To_Reference, data = summary_table, correlation = corPagel(0.1, pruned_tree, form = ~Species), method = "ML")

range_values <- seq(min(summary_table$Log_Mean_Converted_Value), max(summary_table$Log_Mean_Converted_Value), length.out = 100)
fixed_distance <- mean(summary_table$Log_Distance_To_Reference)
new_data <- data.frame(
  Log_Mean_Converted_Value = range_values, 
  Log_Distance_To_Reference = fixed_distance
)
new_data$Log_I_and_PI <- predict(pglsModel_Null_Brownian, newdata = new_data)
new_data$Log_one2zero <- predict(pglsModel_Lambda_one2zero, newdata = new_data)
new_data$Log_one2zero_Ratio <- predict(pglsModel_Null_Blomberg_one2zero_ratio, newdata = new_data)


distance_label <- "Mean Distance" # Creates labels for the fixed distances
new_data$Distance_Label <- factor(fixed_distance, labels = distance_label)

p_text = paste0("p = ", round(coef(summary(pglsModel_Null_Brownian))[2,4],3))
r_squared_value <- round(rsquared.gls(pglsModel_Null_Brownian)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_i_and_pi_hymenoptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_I_and_PI), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  labs(x = "Log Species Mass", y = "Log I/PI Genes") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))+
  annotate("text", x = -Inf, y = Inf, 
           label = p_text, hjust = -0.05, vjust = 5.5, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = Inf, 
           label = r_text, hjust = -0.08, vjust = 3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = Inf, 
           label = "Brownian N.", hjust = -0.03, vjust = 2.75, size = 5, family = "lm-bold")

p_text = paste0("p = ", round(coef(summary(pglsModel_Lambda_one2zero))[2,4],3))
r_squared_value <- round(rsquared.gls(pglsModel_Lambda_one2zero)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_one2zero_hymenoptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_one2zero, linetype = Distance_Label, group = Distance_Label), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Distance Label") +  # Discrete linetype scale for lines
  labs(x = "Log Species Mass", y = "Log One-to-None Genes") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14)) + guides(color = "none", linetype = "none")+
  annotate("text", x = -Inf, y = -Inf, 
           label = p_text, hjust = -0.05, vjust = -2.75, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = -Inf, 
           label = r_text, hjust = -0.08, vjust = -3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = -Inf, 
           label = "Lambda L.", hjust = -0.03, vjust = -5.6, size = 5, family = "lm-bold")

p_text = paste0("p = ", round(coef(summary(pglsModel_Null_Blomberg_one2zero_ratio))[2,4],3))
r_squared_value <- round(rsquared.gls(pglsModel_Null_Blomberg_one2zero_ratio)[[4]], 2)
r_text <- substitute("R"^2 == r_val, list(r_val = r_squared_value))
plot_log_one2zero_ratio_hymenoptera = ggplot() +
  geom_point(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio, color = Log_Distance_To_Reference), alpha = 1) +
  geom_line(data = new_data, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio, linetype = Distance_Label, group = Distance_Label), size = 1.2) +
  geom_smooth(data = summary_table, aes(x = Log_Mean_Converted_Value, y = Log_one2zero_Ratio), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) + 
  scale_color_viridis_c(option = "turbo", name = "Log Distance to Reference", direction = 1) +  # Continuous color scale for points
  scale_linetype_manual(values = c("solid", "dashed", "dotted"), name = "Distance Label") +  # Discrete linetype scale for lines
  labs(x = "Log Species Mass", y = "Log One-to-None Ratio") +
  theme_pubr() +
  guides(color = guide_colorbar(reverse = TRUE)) +
  theme(legend.position = "right",
        text = element_text(family = "lm", size = 14))+ guides(color = "none", linetype = "none")+
  annotate("text", x = -Inf, y = -Inf, 
           label = p_text, hjust = -0.05, vjust = -2.75, size = 5, fontface = "italic", family = "lm-bold") +
  annotate("text", x = -Inf, y = -Inf, 
           label = r_text, hjust = -0.08, vjust = -3.2, size = 5, family = "lm") +
  annotate("text", x = -Inf, y = -Inf, 
           label = "Blomberg N.", hjust = -0.03, vjust = -5.6, size = 5, family = "lm-bold")


############################################################################################################################################################################################################################
font_add(family = "lm", regular = "/usr/share/texlive/texmf-dist/fonts/opentype/public/lm/lmroman10-regular.otf")
font_add(family = "lm-bold", regular = "/usr/share/texmf/fonts/opentype/public/lm/lmroman10-bold.otf")
showtext_auto()

# Define your plots
plots <- list(
  plot_log_i_and_pi_coleoptera, plot_log_one2zero_coleoptera, plot_log_one2zero_ratio_coleoptera,
  plot_log_i_and_pi_diptera, plot_log_one2zero_diptera, plot_log_one2zero_ratio_diptera,
  plot_log_i_and_pi_hymenoptera, plot_log_one2zero_hymenoptera, plot_log_one2zero_ratio_hymenoptera,
  plot_log_i_and_pi_lepidoptera, plot_log_one2zero_lepidoptera, plot_log_one2zero_ratio_lepidoptera
)

# Adjust themes for each plot: remove y-axis titles and ensure consistent margins
plots <- lapply(plots, function(p) {
  p + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.text.x = element_text(size = 20),  # Set x-axis text size to 14
    axis.text.y = element_text(size = 20, angle = 90, hjust = 0.5),  # Set y-axis text size to 14
    plot.title = element_text(family = "lm-bold"),
    plot.margin = margin(b = 40, l = 18, t = 10),  # Apply consistent margins to all plots
    legend.position = "none",
  ) +
    scale_x_continuous(n.breaks = 3) +
    scale_y_continuous(n.breaks = 3)
})

plots[[3]] <- plots[[3]] +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4)

plots[[8]] <- plots[[8]] +
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(n.breaks = 4)
plots[[10]] <- plots[[10]] +
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(n.breaks = 4)
plots[[12]] <- plots[[12]] +
  scale_x_continuous(n.breaks = 3) +
  scale_y_continuous(n.breaks = 5)
plots[[6]] <- plots[[6]] +
  scale_x_continuous(n.breaks = 4) +
  scale_y_continuous(n.breaks = 4)


plts = c(4:6)
for (i in plts){
  plots[[i]] = plots[[i]] + 
    scale_x_continuous(breaks = c(1, 1.5, 2)) +
    scale_y_continuous(n.breaks = 4)
}

plts = c(7:9)
for (i in plts){
  plots[[i]] = plots[[i]] + 
    scale_x_continuous(breaks = c(1, 1.5, 2)) +
    scale_y_continuous(n.breaks = 4)
}

plts = c(1,2,4,5,7,8,10,11)
for (i in plts){
  plots[[i]] = plots[[i]] + 
    theme(axis.text.x = element_blank())
}



# Create labeled spaces between the columns (for each column)
label_a <- textGrob("a) Coleoptera", y = 0.5, hjust = 0.35, gp = gpar(fontfamily = "lm-bold", fontsize = 21))
label_b <- textGrob("b) Diptera", y = 0.5, hjust = 0.5, just = "center", gp = gpar(fontfamily = "lm-bold", fontsize = 21))
label_c <- textGrob("c) Hymenoptera", y = 0.5, hjust = 0.375,just = "center", gp = gpar(fontfamily = "lm-bold", fontsize = 21))
label_d <- textGrob("d) Lepidoptera", y = 0.5,  hjust = 0.4,just = "center", gp = gpar(fontfamily = "lm-bold", fontsize = 21))

# Arrange labels into a single row
column_labels <- plot_grid(label_a, label_b, label_c, label_d, ncol = 4)

# Create y-axis labels for each row
y_axis_label_row1 <- textGrob("Log I/PI Genes", rot = 90, hjust = 0.32, gp = gpar(fontfamily = "lm-bold", fontsize = 21))
y_axis_label_row2 <- textGrob("Log One-to-None Genes", rot = 90, hjust = 0.42, gp = gpar(fontfamily = "lm-bold", fontsize = 21))
y_axis_label_row3 <- textGrob("Log One-to-None Ratio", rot = 90, hjust = 0.42, gp = gpar(fontfamily = "lm-bold", fontsize = 21))

# Create individual row arrangements with y-axis labels
row1 <- plot_grid(y_axis_label_row1, plot_grid(plots[[1]], plots[[4]], plots[[7]], plots[[10]], ncol = 4), rel_widths = c(0.025, 1))
row2 <- plot_grid(y_axis_label_row2, plot_grid(plots[[2]], plots[[5]], plots[[8]], plots[[11]], ncol = 4), rel_widths = c(0.025, 1))
row3 <- plot_grid(y_axis_label_row3, plot_grid(plots[[3]], plots[[6]], plots[[9]], plots[[12]], ncol = 4), rel_widths = c(0.025, 1))

# Add space between the rows using plot_spacer()
combined_plot <- plot_grid(
  column_labels,  # Column labels at the top
  row1,
  row2,
  row3,
  ncol = 1,
  rel_heights = c(0.1, 1, 1, 1)  # Adjust heights: 0.05 for the spacers, 0.4 for the legend
)

# Add the common x-axis label at the bottom
final_plot <- grid.arrange(
  combined_plot,
  bottom = textGrob("Log Species Mass", hjust = 0.3, vjust = -0.5, gp=gpar(fontfamily = "lm-bold", fontsize = 21)),
  ncol = 1,
  heights = c(2.7, 0.01)
)

# Display the plot
print(final_plot)


