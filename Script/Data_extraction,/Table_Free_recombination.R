# Load necessary packages for data manipulation and string handling
library(dplyr)
library(stringr)
library(tools)

# === Function to extract species names and branch lengths from a Newick-formatted tree file ===
extract_species_and_branch_lengths <- function(newick_path) {
  # Read the entire Newick tree file as a single string
  newick_str <- readLines(newick_path, warn = FALSE)
  newick_str <- paste(newick_str, collapse = "")
  
  # Use regex to find all matches of "SpeciesName:BranchLength"
  matches <- str_match_all(newick_str, "([A-Za-z]+_[A-Za-z]+):([0-9\\.eE\\-]+)")[[1]]
  
  species <- unique(matches[, 2])
  lengths <- as.numeric(matches[, 3])
  
  df <- data.frame(species = species, branch_length = lengths, stringsAsFactors = FALSE)
  cat("Detected species with their branch lengths:\n")
  print(df)
  return(df)
}

# === Main analysis function ===
analyze_species_from_newick_mod <- function(base_folder, newick_path, output_csv, replicate_number) {
  species_df <- extract_species_and_branch_lengths(newick_path)
  species_list <- species_df$species
  
  final_df <- data.frame(
    species = character(),
    pN_pS = numeric(),
    Pin = numeric(),
    Pis = numeric(),
    Nw_sum = numeric(),
    Sw_sum = numeric(),
    dN_dS = numeric(),
    dN = numeric(),
    dS = numeric(),
    Nw_sum.1 = numeric(),
    Sw_sum.1 = numeric(),
    taille_pop = numeric(),
    branch_length = numeric(),
    stringsAsFactors = FALSE
  )
  
  folder_path <- file.path(base_folder, paste0("sim_", replicate_number))
  
  # If simulation folder does not exist, fill NA for all species
  if (!dir.exists(folder_path)) {
    cat("❌ Simulation folder sim_", replicate_number, " not found. Filling all species with NA.\n")
    final_df <- data.frame(
      species = species_list,
      pN_pS = NA, Pin = NA, Pis = NA,
      Nw_sum = NA, Sw_sum = NA,
      dN_dS = NA, dN = NA, dS = NA,
      Nw_sum.1 = NA, Sw_sum.1 = NA,
      taille_pop = NA,
      branch_length = species_df$branch_length,
      stringsAsFactors = FALSE
    )
    write.csv(final_df, file = output_csv, row.names = FALSE)
    cat("✅ CSV created with NA for sim_", replicate_number, "\n")
    return()
  }
  
  # Process each species
  for (sp in species_list) {
    cat("\n--- Processing species:", sp, "in replicate sim_", replicate_number, "---\n")
    
    suffix <- sprintf("_sim_%d.txt", replicate_number)
    pnps_path <- file.path(folder_path, paste0(sp, "_pNpS", suffix))
    dnds_path <- file.path(folder_path, paste0(sp, "_dNdS", suffix))
    pop_path  <- file.path(folder_path, paste0(sp, "_pop",  suffix))
    
    # If any file is missing, fill NA for this species
    if (!file.exists(pnps_path) | !file.exists(dnds_path) | !file.exists(pop_path)) {
      cat("❌ One or more files missing for", sp, "-> Filling with NA\n")
      final_df <- rbind(final_df, data.frame(
        species = sp,
        pN_pS = NA, Pin = NA, Pis = NA,
        Nw_sum = NA, Sw_sum = NA,
        dN_dS = NA, dN = NA, dS = NA,
        Nw_sum.1 = NA, Sw_sum.1 = NA,
        taille_pop = NA,
        branch_length = species_df$branch_length[species_df$species == sp],
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Read pNpS
    pnps_lines <- readLines(pnps_path, warn = FALSE)
    pnps_lines <- pnps_lines[!grepl("^#", pnps_lines)]
    if (length(pnps_lines) == 0) { cat("⚠️  pNpS file empty\n"); next }
    pnps_values <- as.numeric(strsplit(tail(pnps_lines, 1), "\\s+")[[1]][1:5])
    
    # Read dNdS
    dnds_lines <- readLines(dnds_path, warn = FALSE)
    dnds_lines <- dnds_lines[!grepl("^#", dnds_lines)]
    if (length(dnds_lines) == 0) { cat("⚠️  dNdS file empty\n"); next }
    dnds_values <- as.numeric(strsplit(tail(dnds_lines, 1), "\\s+")[[1]][1:5])
    
    # Read population size
    pop_lines <- readLines(pop_path, warn = FALSE)
    numeric_lines <- pop_lines[!grepl("^#", pop_lines)]
    if (length(numeric_lines) < 3) { cat("⚠️  Not enough lines in pop file\n"); next }
    taille_pop <- as.numeric(strsplit(numeric_lines[3], "\\s+")[[1]][2])
    
    # Append row
    final_df <- rbind(final_df, data.frame(
      species = sp,
      pN_pS = pnps_values[1],
      Pin = pnps_values[2],
      Pis = pnps_values[3],
      Nw_sum = pnps_values[4],
      Sw_sum = pnps_values[5],
      dN_dS = dnds_values[1],
      dN = dnds_values[2],
      dS = dnds_values[3],
      Nw_sum.1 = dnds_values[4],
      Sw_sum.1 = dnds_values[5],
      taille_pop = taille_pop,
      branch_length = species_df$branch_length[species_df$species == sp],
      stringsAsFactors = FALSE
    ))
    
    cat("✅ Done:", sp, "\n")
  }
  
  # Write final CSV
  write.csv(final_df, file = output_csv, row.names = FALSE)
  cat("\n✅ Data exported to:", output_csv, "\n")
}

# === Batch processing for all replicates ===
base_folder_path <- "/home/alafitte/Internship/Results/Recombination/Free_recombination"
newick_tree_file <- "/home/alafitte/bat_genes_complete_filtered-PhyML_tree.nhx"
output_folder <- "/home/alafitte/Internship/Rapport de stage/Données/Free_recombination"

for (i in 1:50) {
  cat("\n===========================\n")
  cat("🔁 Processing replicate:", i, "\n")
  cat("===========================\n")
  
  output_file <- sprintf("%s/resultats_final_%02d.csv", output_folder, i)
  analyze_species_from_newick_mod(base_folder_path, newick_tree_file, output_file, i)
}
