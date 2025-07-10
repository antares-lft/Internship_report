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
  # Pattern: species names like "Genus_species" and numerical branch lengths
  matches <- str_match_all(newick_str, "([A-Za-z]+_[A-Za-z]+):([0-9\\.eE\\-]+)")[[1]]
  
  # Extract unique species names and corresponding branch lengths
  species <- unique(matches[, 2])
  lengths <- as.numeric(matches[, 3])
  
  # Create a dataframe with species and branch length columns
  df <- data.frame(species = species, branch_length = lengths, stringsAsFactors = FALSE)
  
  cat("Detected species with their branch lengths:\n")
  print(df)
  return(df)
}

# === Main analysis function customized for 'Low_recombination/sim_n/' folders ===
analyze_species_from_newick_mod <- function(base_folder, newick_path, output_csv, replicate_number) {
  # Extract species and branch lengths from the Newick tree
  species_df <- extract_species_and_branch_lengths(newick_path)
  species_list <- species_df$species
  
  # Initialize an empty dataframe to collect results for all species
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
  
  # Define path to the replicate folder
  folder_path <- file.path(base_folder, paste0("sim_", replicate_number))
  
  # Loop over each species to read its data files and extract statistics
  for (sp in species_list) {
    cat("\n--- Processing species:", sp, "in replicate sim_", replicate_number, "---\n")
    
    # Build paths to the three expected files for this species in this replicate
    suffix <- sprintf("_sim_%d.txt", replicate_number)
    pnps_path <- file.path(folder_path, paste0(sp, "_pNpS", suffix))
    dnds_path <- file.path(folder_path, paste0(sp, "_dNdS", suffix))
    pop_path  <- file.path(folder_path, paste0(sp, "_pop",  suffix))
    
    # Check if files exist, skip species if any are missing
    if (!file.exists(pnps_path)) { cat("❌ Missing file:", pnps_path, "\n"); next }
    if (!file.exists(dnds_path)) { cat("❌ Missing file:", dnds_path, "\n"); next }
    if (!file.exists(pop_path))  { cat("❌ Missing file:", pop_path,  "\n"); next }
    
    # Read and process pN/pS file (ignore comment lines starting with '#')
    pnps_lines <- readLines(pnps_path, warn = FALSE)
    pnps_lines <- pnps_lines[!grepl("^#", pnps_lines)]
    if (length(pnps_lines) == 0) { cat("⚠️  pNpS file is empty\n"); next }
    # Extract first 5 numeric values from the last line of the file
    pnps_values <- as.numeric(strsplit(tail(pnps_lines, 1), "\\s+")[[1]][1:5])
    
    # Read and process dN/dS file (similarly ignoring comment lines)
    dnds_lines <- readLines(dnds_path, warn = FALSE)
    dnds_lines <- dnds_lines[!grepl("^#", dnds_lines)]
    if (length(dnds_lines) == 0) { cat("⚠️  dNdS file is empty\n"); next }
    dnds_values <- as.numeric(strsplit(tail(dnds_lines, 1), "\\s+")[[1]][1:5])
    
    # Read population size from pop file, expecting at least 3 lines with data
    pop_lines <- readLines(pop_path, warn = FALSE)
    numeric_lines <- pop_lines[!grepl("^#", pop_lines)]
    if (length(numeric_lines) < 3) { cat("⚠️  Not enough lines in pop file\n"); next }
    taille_pop <- as.numeric(strsplit(numeric_lines[3], "\\s+")[[1]][2])
    
    # Retrieve the branch length for this species from the species dataframe
    branch_length <- species_df$branch_length[species_df$species == sp]
    
    # Append the extracted and calculated data to the final dataframe
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
      branch_length = branch_length,
      stringsAsFactors = FALSE
    ))
    
    # Print confirmation and summary for the current species
    cat("✅ pN/pS values:", paste(pnps_values, collapse = ", "), "\n")
    cat("✅ dN/dS values:", paste(dnds_values, collapse = ", "), "\n")
    cat("✅ Population size:", taille_pop, "\n")
    cat("✅ Branch length:", branch_length, "\n")
  }
  
  # Display the complete results table before exporting
  cat("\n=== Final aggregated table ===\n")
  print(final_df)
  
  # Write the results to a CSV file
  write.csv(final_df, file = output_csv, row.names = FALSE)
  cat("\n✅ Data exported to:", output_csv, "\n")
}

# === Loop to generate the 20 CSV output files for replicates 1 to 20 ===
base_folder_path <- "/home/alafitte/Internship/Results/Population_size/N2500/Intermediate_recombination"
newick_tree_file <- "/home/alafitte/bat_genes_complete_filtered-PhyML_tree.nhx"
output_folder <- "/home/alafitte/Internship/Rapport de stage/Données/Population_size/Intermediate_recombination/N2500"

for (i in 1:20) {
  cat("\n===========================\n")
  cat("🔁 Processing replicate:", i, "\n")
  cat("===========================\n")
  
  output_file <- sprintf("%s/resultats_final_%02d.csv", output_folder, i)
  analyze_species_from_newick_mod(base_folder_path, newick_tree_file, output_file, i)
}
