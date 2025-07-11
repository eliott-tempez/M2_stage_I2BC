# For each of the 116 genomes, this script determines the
# closest neighbour genome. The output is saved in neighbour_genomes.csv


library(yaml)
library(ggtree)
library(ape)

config <- yaml.load_file("/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/my_functions/filepaths.yaml")
paths <- config$local_paths
tree_file <- paths$tree_file
output_folder <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/2_initial_data_analysis_figures/"

# Read the tree
tree <- read.tree(tree_file)
# Get the list of genomes
genomes <- tree$tip.label
# Create a data frame to store the neighbour genomes
neighbour_genomes <- data.frame(genome = genomes, neighbour = NA, stringsAsFactors = FALSE)
# For each genome, find the closest neighbour
for (i in seq_along(genomes)) {
  # Get the current genome
  current_genome <- genomes[i]
  # Get the distances to all other genomes
  distances <- cophenetic.phylo(tree)[current_genome, ]
  # Remove the current genome from the distances
  distances <- distances[distances != 0]
  # Find the closest neighbour
  closest_neighbour <- names(which.min(distances))
  # Store the closest neighbour
  neighbour_genomes$neighbour[i] <- closest_neighbour
}
# Save the neighbour genomes to a CSV file
output_file <- file.path(output_folder, "neighbour_genomes.csv")
write.table(neighbour_genomes, output_file, row.names = FALSE, quote = FALSE, sep = "\t")
