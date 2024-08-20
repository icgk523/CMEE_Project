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

ggplot(repeats, aes(x=log_total_TE)) + 
  geom_histogram(binwidth=0.1, fill = "#eb674a", color = "black") +
  theme_classic() +
  labs(x = "Log Proportion of TE", y = "Number of Species") +
  facet_wrap(~ Order.x)
ggplot(repeats, aes(x=total_TE)) + 
  geom_histogram(binwidth=0.05, fill = "#af674b", color = "black") +
  theme_classic() +
  labs(x = "Log Proportion of TE", y = "Number of Species") +
  facet_wrap(~ Order.x)

# scale_x_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, by = 0.1))
ggplot(repeats, aes(y = log_total_TE, x = Log_Mean_Converted_Value, group = Order.x)) +
  geom_point(aes(color = Order.x)) + 
  theme_pubr() +
  scale_color_viridis(discrete = TRUE) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y = "Log Species Weight", x = "Log Species Genome Size")


repeats_subset = subset(repeats, Order.x == "Lepidoptera")
# repeats_subset = repeats

phylo_tree_subset = keep.tip(phylo_tree, repeats_subset$Species[repeats_subset$Species %in% phylo_tree$tip.label]) # Subset tree to species in data

formula <- as.formula("log_total_TE ~ Log_Mean_Converted_Value")
formula_sq <- as.formula("log_total_TE ~ poly(Log_Mean_Converted_Value, 2, raw = TRUE)")
null_formula = as.formula(paste("log_total_TE ~ 1"))
linear_model_null <- lm(null_formula, data = repeats_subset)
linear_model <- lm(formula, data = repeats_subset)
linear_model_residuals <- residuals(linear_model)
names(linear_model_residuals) <- repeats_subset$Species

# Test phylogenetic signal
lambda = phylosig(phylo_tree_subset, linear_model_residuals, method = "lambda", test = TRUE); print(lambda)
# if (lambda[1] < 0.5){
#   lambda_starting = 0.1
# } else {
#   lambda_starting = 1
# }
lambda_starting = 1
pglsModel_Lambda = gls(formula, data = repeats_subset, correlation = corPagel(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian = gls(formula, data = repeats_subset, correlation = corBrownian(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins= gls(formula, data = repeats_subset, correlation = corMartins(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg = gls(formula, data = repeats_subset, correlation = corBlomberg(0.04, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen = gls(formula, data = repeats_subset, correlation = corGrafen(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Lambda_sq = gls(formula_sq, data = repeats_subset, correlation = corPagel(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Brownian_sq = gls(formula_sq, data = repeats_subset, correlation = corBrownian(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Martins_sq = gls(formula_sq, data = repeats_subset, correlation = corMartins(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Blomberg_sq = gls(formula_sq, data = repeats_subset, correlation = corBlomberg(0.04, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Grafen_sq = gls(formula_sq, data = repeats_subset, correlation = corGrafen(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML")

pglsModel_Null_Brownian = gls(null_formula, data = repeats_subset, correlation = corBrownian(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Pagel = gls(null_formula, data = repeats_subset, correlation = corPagel(1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Blomberg = gls(null_formula, data = repeats_subset, correlation = corBlomberg(0.1, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Grafen = gls(null_formula, data = repeats_subset, correlation = corGrafen(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML")
pglsModel_Null_Martins = gls(null_formula, data = repeats_subset, correlation = corMartins(lambda_starting, phylo_tree_subset, form = ~Species), method = "ML")


models = list(linear_model = linear_model,
              linear_model_null = linear_model_null,
              pglsModel_Lambda = pglsModel_Lambda, 
              pglsModel_Brownian = pglsModel_Brownian,
              pglsModel_Grafen = pglsModel_Grafen, 
              pglsModel_Martins = pglsModel_Martins,
              
              pglsModel_Lambda_sq = pglsModel_Lambda_sq, 
              pglsModel_Brownian_sq = pglsModel_Brownian_sq,
              pglsModel_Grafen_sq = pglsModel_Grafen_sq, 
              pglsModel_Martins_sq = pglsModel_Martins_sq,
              
              pglsModel_Null_Grafen = pglsModel_Null_Grafen, 
              pglsModel_Null_Martins = pglsModel_Null_Martins, 
              pglsModel_Null_Brownian = pglsModel_Null_Brownian, 
              pglsModel_Null_Pagel = pglsModel_Null_Pagel) 



best_model <- get.best.model(models)
summary(best_model)

AIC(linear_model, linear_model_null,
    pglsModel_Lambda, pglsModel_Brownian, pglsModel_Martins, pglsModel_Grafen,  pglsModel_Blomberg, 
    pglsModel_Lambda_sq, pglsModel_Brownian_sq, pglsModel_Martins_sq, pglsModel_Grafen_sq, pglsModel_Blomberg_sq,
    pglsModel_Null_Grafen, pglsModel_Null_Martins, pglsModel_Null_Brownian, pglsModel_Null_Pagel, pglsModel_Null_Blomberg)

# Define the range of 'Log_Mean_Converted_Value'
range_values <- seq(min(repeats_subset$Log_Mean_Converted_Value), max(repeats_subset$Log_Mean_Converted_Value), length.out = 100)

# Define the desired colors manually
viridis_colors <- viridis(4)
model_colors <- c(
  "Simple Linear" = viridis_colors[3],  # Darkest
  "PGLS Quadratic" = viridis_colors[2],    # Second lightest
  "PGLS Linear" = viridis_colors[1],     # Lightest
  "PGLS Null" = "#B12A90FF"
)

# Compute the predicted values for each model
predictions <- data.frame(
  Log_Mean_Converted_Value = rep(range_values, 1),
  log_total_TE = c(
    linear_model$coefficients[1] + linear_model$coefficients[2] * range_values,
    pglsModel_Lambda$coefficients[1] +  pglsModel_Lambda$coefficients[2] * range_values),
  # pglsModel_Lambda_sq$coefficients[1] + pglsModel_Lambda_sq$coefficients[2] * range_values + pglsModel_Lambda_sq$coefficients[3] * (range_values^2)),
  Model = rep(c("Simple Linear", "PGLS Linear"), each = 100)
)
# Extract fitted values and residuals from the null model
intercept <- pglsModel_Null_Pagel$coefficients[1]
null_model_df <- data.frame(
  Log_Mean_Converted_Value = repeats_subset$Log_Mean_Converted_Value,
  log_total_TE = repeats_subset$log_total_TE,
  Fitted = rep(pglsModel_Null_Pagel$coefficients[1], nrow(repeats_subset)),
  Residuals = residuals(pglsModel_Null_Pagel)
)
# # Plot the data using ggplot2
# insecta = ggplot(null_model_df, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
#   geom_point(color = "black") +
#   geom_hline(yintercept = intercept, color = viridis_colors[4], linewidth = 1) +
#   labs(title = "a) Insecta", x = "Log Species Mass", y = "Log Proportion of TE") +
#   theme_pubr()+
#   theme(legend.position = "none")
# save(insecta, file = "../results/insecta_ggplot_object.RData")
# 
# coleoptera = ggplot(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
#   geom_point(color = "black") +
#   geom_line(data = predictions, aes(x = Log_Mean_Converted_Value, y = log_total_TE), linewidth = 1) +
#   geom_smooth(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) +
#   labs(title = "b) Coleoptera", x = "Log Species Mass", y = "Log Proportion of TE") +
#   theme_pubr()+
#   scale_color_manual(values = model_colors)+
#   theme(legend.position = "none")
# save(coleoptera, file = "../results/coleoptera_ggplot_object.RData")

# diptera = ggplot(data = null_model_df, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
#   geom_point(color = "black") +
#   geom_hline(yintercept = intercept, color = viridis_colors[3], linewidth = 1) +
#   labs(title = "c) Diptera", x = "Log Species Mass", y = "Log Proportion of TE") +
#   theme_pubr()+
#   scale_color_manual(values = model_colors)+
#   theme(legend.position = "none")
# save(diptera, file = "../results/diptera_ggplot_object.RData")

# hemiptera = ggplot(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
#     geom_point(color = "black") +
#     geom_line(data = predictions, aes(x = Log_Mean_Converted_Value, y = log_total_TE), linewidth = 1) +
#     geom_smooth(data = repeats_subset, aes(x = Log_Mean_Converted_Value, y = log_total_TE), method = "lm", linetype = "dashed", se = FALSE, color = "black", size = 1) +
#     labs(title = "d) Hemiptera", x = "Log Species Mass", y = "Log Proportion of TE") +
#     theme_pubr()+
#     scale_color_manual(values = model_colors)+
#     theme(legend.position = "none")
# save(hemiptera, file = "../results/hemiptera_ggplot_object.RData")

hymenoptera = ggplot(null_model_df, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
  geom_point(color = "black") +
  geom_hline(yintercept = intercept, color = viridis_colors[4], linewidth = 1) +
  labs(title = "e) Hymenoptera", x = "Log Species Mass", y = "Log Proportion of TE") +
  theme_pubr()+
  theme(legend.position = "none")
# save(hymenoptera, file = "../results/hymenoptera_ggplot_object.RData")

# 
# lepidoptera = ggplot(null_model_df, aes(x = Log_Mean_Converted_Value, y = log_total_TE)) +
#    geom_point(color = "black") +
#    geom_hline(yintercept = intercept, color = viridis_colors[4], linewidth = 1) +
#    labs(title = "f) Lepidoptera", x = "Log Species Mass", y = "Log Proportion of TE") +
#    theme_pubr()+
#    theme(legend.position = "none")
# save(lepidoptera, file = "../results/lepidoptera_ggplot_object.RData")


font_add(family = "lm", regular = "/usr/share/texlive/texmf-dist/fonts/opentype/public/lm/lmroman10-regular.otf")
font_add(family = "lm-bold", regular = "/usr/share/texmf/fonts/opentype/public/lm/lmroman10-bold.otf")
showtext_auto()

load("../results/coleoptera_ggplot_object.RData")
load("../results/lepidoptera_ggplot_object.RData")
load("../results/hymenoptera_ggplot_object.RData")
load("../results/diptera_ggplot_object.RData")
load("../results/hemiptera_ggplot_object.RData")
load("../results/insecta_ggplot_object.RData")

# Function to update plot lines and return the updated gtable
modify_plot_line <- function(plot_gtable, line_color = "black", line_type = "solid") {
  # Identify the panel grob
  panel_index <- which(plot_gtable$layout$name == "panel")
  panel_grob <- plot_gtable$grobs[[panel_index]]
  # Inspect the children of the panel to find any grobs of type 'segments'
  segment_grobs <- which(sapply(panel_grob$children, function(x) inherits(x, "segments")))
  
  # If no segment grobs found, return the original gtable
  if (length(segment_grobs) == 0) {
    message("No segment grobs found for line modification.")
    return(plot_gtable)
  }
  # Loop through segment grobs to modify them
  for (index in segment_grobs) {
    line_grob <- panel_grob$children[[index]]
    # Modify the color and linetype
    line_grob$gp$col <- line_color
    line_grob$gp$lty <- line_type
    # Replace the modified line grob back into the panel
    panel_grob$children[[index]] <- line_grob
  }
  
  # Update the gtable with the modified panel grob
  plot_gtable$grobs[[panel_index]] <- panel_grob
  # Return the modified gtable
  return(plot_gtable)
}
modify_plot_points <- function(plot_gtable, color = "#DC4D10", alpha = 0.6) {
  # Identify the panel grob
  panel_index <- which(plot_gtable$layout$name == "panel")
  panel_grob <- plot_gtable$grobs[[panel_index]]
  # Inspect the children of the panel to find any grobs of type 'points'
  points_grobs <- which(sapply(panel_grob$children, function(x) inherits(x, "points")))
  # If no points grobs found, return the original gtable
  if (length(points_grobs) == 0) {
    message("No points grobs found for modification.")
    return(plot_gtable)
  }
  # Loop through points grobs to modify them
  for (index in points_grobs) {
    point_grob <- panel_grob$children[[index]]
    # Modify the color and alpha
    point_grob$gp$col <- color
    point_grob$gp$alpha <- alpha
    # Replace the modified point grob back into the panel
    panel_grob$children[[index]] <- point_grob
  }
  # Update the gtable with the modified panel grob
  plot_gtable$grobs[[panel_index]] <- panel_grob
  # Return the modified gtable
  return(plot_gtable)
}
plot_theme_top <- theme(axis.title.x = element_blank(), 
                        axis.title.y = element_blank(), 
                        plot.title = element_text(family = "lm-bold", size = 19),
                        text = element_text(family = "lm", size = 14), 
                        axis.line.x = element_line(color = "black", size = 0.5, lineend = "square"), 
                        axis.line.y = element_line(color = "black", size = 0.5, lineend = "square"), 
                        axis.ticks.x = element_line(color = "black", size = 0.5), 
                        axis.ticks.y = element_line(color = "black", size = 0.5),
                        plot.margin = margin(b = 50, l = 10),
                        axis.text.x = element_text(size = 20),  # Set x-axis text size to 14
                        axis.text.y = element_text(size = 20, angle = 90, hjust = 0.5)  # Set y-axis text size to 14
) 
plot_theme_bottom <- theme(axis.title.x = element_blank(), 
                           axis.title.y = element_blank(), 
                           plot.title = element_text(family = "lm-bold", size = 19),
                           text = element_text(family = "lm", size = 14), 
                           axis.line.x = element_line(color = "black", size = 0.5, lineend = "square"), 
                           axis.line.y = element_line(color = "black", size = 0.5, lineend = "square"), 
                           axis.ticks.x = element_line(color = "black", size = 0.5), 
                           axis.ticks.y = element_line(color = "black", size = 0.5),
                           plot.margin = margin(b = 25, l = 10),
                           axis.text.x = element_text(size = 20),  # Set x-axis text size to 14
                           axis.text.y = element_text(size = 20, angle = 90, hjust = 0.5),  # Set y-axis text size to 14
)

insecta_combined <- insecta + plot_theme_top  + scale_y_continuous(breaks = c(-1.2, -0.9, -0.6)) 
coleoptera_combined <- coleoptera + plot_theme_top 
diptera_combined <- diptera + plot_theme_top + scale_y_continuous(breaks = c(-1.2, -0.8, -0.4)) 
hemiptera_combined <- hemiptera + plot_theme_bottom + scale_y_continuous(breaks = c(-0.85, -0.75, -0.65)) 
hymenoptera_combined <- hymenoptera + plot_theme_bottom + scale_y_continuous(breaks = c(-0.5, -0.9, -1.3)) 
lepidoptera_combined <- lepidoptera + plot_theme_bottom+ scale_y_continuous(breaks = c(-1, -0.8, -0.6)) 

insecta_modified_gtable <- modify_plot_line(ggplot_gtable(ggplot_build(insecta_combined)))
insecta_modified_gtable <- modify_plot_points(insecta_modified_gtable, color = "#DC4D10", alpha = 0.75)

diptera_modified_gtable <- modify_plot_line(ggplot_gtable(ggplot_build(diptera_combined)), line_color = "black")
diptera_modified_gtable <- modify_plot_points(diptera_modified_gtable, color = "#DC4D10", alpha = 0.75)

hymenoptera_modified_gtable <- modify_plot_line(ggplot_gtable(ggplot_build(hymenoptera_combined)))
hymenoptera_modified_gtable <- modify_plot_points(hymenoptera_modified_gtable, color = "#DC4D10", alpha = 0.75)

lepidoptera_modified_gtable <- modify_plot_line(ggplot_gtable(ggplot_build(lepidoptera_combined)))
lepidoptera_modified_gtable <- modify_plot_points(lepidoptera_modified_gtable, color = "#DC4D10", alpha = 0.75)

coleoptera_gtable <- modify_plot_points(ggplot_gtable(ggplot_build(coleoptera_combined)), color = "#DC4D10", alpha = 0.75)
hemiptera_gtable <- modify_plot_points(ggplot_gtable(ggplot_build(hemiptera_combined)), color = "#DC4D10", alpha = 0.75)

# Convert the modified gtables back to ggplot objects
insecta_combined <- as.ggplot(insecta_modified_gtable)
coleoptera_combined <- as.ggplot(coleoptera_gtable)  # No modification
diptera_combined <- as.ggplot(diptera_modified_gtable)
hemiptera_combined <- as.ggplot(hemiptera_gtable)  # No modification
hymenoptera_combined <- as.ggplot(hymenoptera_modified_gtable)
lepidoptera_combined <- as.ggplot(lepidoptera_modified_gtable)
insecta_combined <- as.ggplot(insecta_modified_gtable)
coleoptera_combined <- as.ggplot(coleoptera_gtable)  # No modification
diptera_combined <- as.ggplot(diptera_modified_gtable)
hemiptera_combined <- as.ggplot(hemiptera_gtable)  # No modification
hymenoptera_combined <- as.ggplot(hymenoptera_modified_gtable)
lepidoptera_combined <- as.ggplot(lepidoptera_modified_gtable)

# Convert gtable to ggplot objects
insecta_combined <- as.ggplot(insecta_modified_gtable)
coleoptera_combined <- as.ggplot(coleoptera_gtable)  # No modification
diptera_combined <- as.ggplot(diptera_modified_gtable)
hemiptera_combined <- as.ggplot(hemiptera_gtable)  # No modification
hymenoptera_combined <- as.ggplot(hymenoptera_modified_gtable)
lepidoptera_combined <- as.ggplot(lepidoptera_modified_gtable)


# Combine the plots into a single layout
combined_plot <- egg::ggarrange(
  insecta_combined, 
  coleoptera_combined, 
  diptera_combined, 
  hemiptera_combined, 
  hymenoptera_combined, 
  lepidoptera_combined, 
  ncol = 3
)

# Combine the plot with the legend and add annotations
final_plot <- plot_grid(combined_plot) +
  annotate("text", label = "Log Species Mass", x = 0.5, y = 0.02, vjust = 1.75, size = 7.5, family = "lm-bold") + 
  annotate("text", label = "Log Proportion of TE per Species", x = 0.02, y = 0.5, vjust = -2, angle = 90, size = 7.5, family = "lm-bold") +
  theme(plot.margin = margin(l = 30, r = 10, t = 10, b = 20)) + 
  theme(text = element_text(family = "lm"))

# Display the final plot
print(final_plot)


