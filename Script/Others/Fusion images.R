library(magick)

# --- Charger les images ---
img3 <- image_read("/home/alafitte/Internship/Rapport de stage/Images/Results/Low_recombination/Report/dNdS vs pop_size.png")
img2 <- image_read("/home/alafitte/Internship/Rapport de stage/Images/Results/Free_recombination/Report/dNdS_vs_pop_size.png")
img1 <- image_read("/home/alafitte/Internship/Rapport de stage/Images/Results/Intermediate_recombination/Report/dNdS vs pop_size.png")

# --- Redimensionner les images ---
img1_resized <- image_scale(img1, "500")
img2_resized <- image_scale(img2, "500")
img3_resized <- image_scale(img3, "500")

# --- Créer des étiquettes alignées à gauche ---
title_A <- image_blank(width = 500, height = 50, color = "white") %>%
  image_annotate("A", size = 30, gravity = "northwest", location = "+10+10", color = "black")

title_B <- image_blank(width = 500, height = 50, color = "white") %>%
  image_annotate("B", size = 30, gravity = "northwest", location = "+10+10", color = "black")

title_C <- image_blank(width = 500, height = 50, color = "white") %>%
  image_annotate("C", size = 30, gravity = "northwest", location = "+10+10", color = "black")

# --- Empiler les titres + images verticalement ---
img1_labeled <- image_append(c(title_A, img1_resized), stack = TRUE)
img2_labeled <- image_append(c(title_B, img2_resized), stack = TRUE)
img3_labeled <- image_append(c(title_C, img3_resized), stack = TRUE)

# --- Assembler les 3 ensembles horizontalement ---
final_combined <- image_append(c(img1_labeled, img2_labeled, img3_labeled))

# --- Affichage et sauvegarde ---
print(final_combined)
image_write(final_combined, "/home/alafitte/Images/Captures d’écran/Figure 5.png")

