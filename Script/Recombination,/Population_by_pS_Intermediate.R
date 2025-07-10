# Librairies
library(tidyverse)

# Dossier contenant les fichiers CSV
folder <- "/home/alafitte/Internship/Rapport de stage/Données/Intermediate_recombination"

# Fichiers
files <- file.path(folder, sprintf("resultats_final_%02d.csv", 1:100))

# Constantes
L <- 1000
µ <- 1e-6

# Lecture et combinaison
combined_data <- files %>%
  map_dfr(read_csv, col_types = cols()) %>%
  filter(if_all(c(Pis, branch_length, taille_pop), ~ !is.na(.))) %>%
  filter(taille_pop > 0, Pis > 0, branch_length > 0)

# Agrégation
donnees_par_espece <- combined_data %>%
  group_by(species, taille_pop) %>%
  summarise(
    Pis_total = sum(Pis),
    Sw_sum_total = sum(Sw_sum),
    moyenne_branch_length = mean(branch_length),
    .groups = "drop"
  ) %>%
  filter(moyenne_branch_length > 0) %>%
  mutate(
    N_effectif = (Pis_total * L) / (4 * Sw_sum_total * µ),
    log_branch_length = log10(moyenne_branch_length),
    ecart = N_effectif - taille_pop
  )

# Modèle linéaire
modele_lm <- lm(N_effectif ~ taille_pop, data = donnees_par_espece)
pente_lm <- round(coef(modele_lm)[2], 3)

# R² par rapport à la diagonale y = x
residus_diag <- donnees_par_espece$N_effectif - donnees_par_espece$taille_pop
SSE_diag <- sum(residus_diag^2)
SST_diag <- sum((donnees_par_espece$N_effectif - mean(donnees_par_espece$N_effectif))^2)
r2_diag <- round(1 - SSE_diag / SST_diag, 3)

# Limites pour échelle 1:1
lims <- range(c(donnees_par_espece$taille_pop, donnees_par_espece$N_effectif))

ggplot(donnees_par_espece, aes(x = taille_pop, y = N_effectif)) +
  geom_segment(aes(xend = taille_pop, yend = taille_pop),
               arrow = arrow(length = unit(0.15, "cm")),
               alpha = 0.3, color = "red") +
  geom_point(aes(color = log_branch_length), size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "lightgray") +
  scale_color_gradient(name = "Branch length", low = "blue", high = "yellow") +
  annotate("text", x = lims[1] + 0.05 * diff(lims),
           y = lims[2] - 0.05 * diff(lims),
           label = paste0("R² (vs y = x) = ", r2_diag),
           hjust = 0, size = 4.5) +
  annotate("text", x = lims[1] + 0.05 * diff(lims),
           y = lims[2] - 0.15 * diff(lims),
           label = paste0("Slope = ", pente_lm),
           hjust = 0, size = 4.5) +
  coord_fixed(ratio = 1, xlim = lims, ylim = lims) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 15)  # <---- Double la taille des titres des axes
  ) +
  labs(
    x = "N implemented",
    y = "N estimated"
  )


# Export (facultatif)
ggsave("/home/alafitte/Internship/Rapport de stage/Images/Results/Intermediate_recombination/Report/Picture_7.png",
       width = 8, height = 6, dpi = 300)