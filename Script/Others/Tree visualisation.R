# Charger les packages
library(treeio)
library(ggtree)

# Lire l'arbre avec les annotations NHX
arbre <- read.nhx("/home/alafitte/bat_genes_complete_filtered-PhyML_tree.nhx")

# Créer l’arbre avec ggtree
p <- ggtree(arbre, aes(color = log10(as.numeric(W)))) + 
  geom_tiplab(size = 5) +  # Tu peux augmenter la taille du texte si besoin
  scale_color_viridis_c(option = "plasma", name = "Population size") +
  theme_tree2() + 
  ggtitle(" ")

# Sauvegarder en grand format
ggsave("/home/alafitte/Images/Captures d’écran/Figure 2.png", plot = p, width = 37, height = 6, dpi = 300)

# Afficher dans la session si tu veux
print(p)