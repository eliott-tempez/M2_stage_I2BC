library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)
library(yaml)
library(dplyr)
library(phangorn)
library(viridis)
ggsave <- function(..., bg = "white",
                   width = 1000, height = 1000,
                   units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg,
                  width = width,
                  height = height,
                  units = units,
                  dpi = dpi)
}


config <- yaml.load_file("/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/my_functions/filepaths.yaml")
paths <- config$local_paths



# Functions
# Get the lca node for a given cluster
get_lca_node <- function(tree, tips) {
  tips <- intersect(tips, tree$tip.label)
  # If only one genome in the cluster, return tip point
  if (length(tips) == 1) {
    return(which(tree$tip.label == tips))
  } else if (length(tips) > 1) {
    return(getMRCA(tree, tips))
  } else {
    return(NA)  # No valid tips
  }
}



# Filenames
tree_file <- paths$tree_file
dense_dir <- paths$dense_dir
genomes_list <- paths$genomes_list
fa_dir <- paths$fa_dir
output_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/18_clustering/"
good_candidates <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/good_candidates.txt"
cluster_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/18_clustering/good_candidates_clustering.tsv"
species_file <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/13_plot_dense_results/species.tsv"



# Read genomes list
genomes <- scan(genomes_list, what = "", sep = " ")[1:116]
# Read list of good denovo
good_candidates <- scan(good_candidates, what = "", sep = "\n")
# Get the info for each genome
data <- data.frame(matrix(NA, nrow = length(genomes), ncol = 3))
denovo_genomes <- data.frame(matrix(NA, nrow = 0, ncol = 2))
colnames(denovo_genomes) <- c("denovo", "genome")
rownames(data) <- genomes
colnames(data) <- c("n_denovo", "n_good_candidates", "gc_perc")
for (g in genomes) {
  # Filenames
  dense_file <- paste0(dense_dir, g, "/denovogenes.tsv")
  ## Number of denovo
  denovo <- read.table(dense_file, header = TRUE, sep = "\t")
  denovo_names <- denovo$CDS
  # Add the match denovo/genome to dtf
  for (i in seq_along(denovo_names)) {
    denovo_name <- denovo_names[i]
    denovo_genomes <- rbind(denovo_genomes,
                            data.frame(denovo = denovo_name, genome = g))
  }
  n_denovo <- nrow(denovo)
  # Good candidates
  n_good_candidates <- sum(denovo_names %in% good_candidates)
  ## GC content
  fa_file <- paste0(fa_dir, g, ".fa")
  fasta <- readLines(fa_file)
  # get rid of > lines
  fasta <- fasta[!grepl(">", fasta)]
  # Get the number of nucleotides
  n_nucleotides <- sum(nchar(fasta))
  # Get the GC content
  n_gc <- sum(nchar(gsub("[^GCgc]", "", paste(fasta, collapse = ""))))
  gc <- n_gc / n_nucleotides * 100

  # Populate the data
  data[g, "n_denovo"] <- n_denovo
  data[g, "gc_perc"] <- gc
  data[g, "n_good_candidates"] <- n_good_candidates
}
data$n_bad_candidates <- data$n_denovo - data$n_good_candidates
n_good <- sum(data$n_good_candidates)
n_bad <- sum(data$n_bad_candidates)
# Distinct species
species_data <- read.table(species_file, header = TRUE, sep = "\t")
colnames(species_data) <- c("genome", "genre", "species.ani", "mash", "aai", "rename")
data$species <- species_data[match(rownames(data), species_data$genome), "species.ani"]



# Read the tree
tree <- read.tree(tree_file)

# Extract the clusters
cluster_df <- read.table(cluster_file, header = FALSE, sep = "\t")
colnames(cluster_df) <- c("cluster", "denovo")
# Add the genomes
cluster_df <- left_join(cluster_df, denovo_genomes, by = "denovo")
# Get the lcas
cluster_lcas <- cluster_df %>%
  group_by(cluster) %>%
  summarise(
    genomes_in_cluster = list(unique(genome)),
    .groups = 'drop'
  ) %>%
  mutate(
    lca_node = sapply(genomes_in_cluster, function(tips) get_lca_node(tree, tips))
  ) %>%
  filter(!is.na(lca_node))


# Plot the tree
p <- ggtree(tree, layout = "circular", branch.length = "none")
tree_data <- p$data

# Match LCA nodes to tree coordinates
cluster_lcas_coords <- cluster_lcas %>%
  left_join(tree_data, by = c("lca_node" = "node"))
# Keep only unique lcas and count
cluster_lcas_summary <- cluster_lcas_coords %>%
  group_by(lca_node, x, y) %>%
  summarise(
    count = n(),
    label = first(label),
    .groups = "drop"
  )
cluster_lcas_summary$uni_group <- as.factor("cluster LCAs")



# Add the heatmaps
## Good denovo ##
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "n_good_candidates", drop = FALSE], offset = 0.5,
        width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "#179207",
            name = paste0("Disrupted\nintegrity\n(n = ", n_good, ")"),
            guide = guide_colorbar(order = 1))

## Bad denovo ##
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "n_bad_candidates", drop = FALSE], offset = 1.5,
        width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "white", high = "#4f535a",
            name = paste0("Preserved\nintegrity\n(n =", n_bad, ")"),
            guide = guide_colorbar(order = 2,
                       limits = c(0, max(data$n_bad_candidates))))

## GC content ##
p <- p + new_scale_fill()
p <- gheatmap(p, data[, "gc_perc", drop = FALSE], offset = 2.5,
              width = .05, colnames = FALSE) +
  scale_fill_gradient(low = "#009E73", high = "#6b00b2",
                      name = "GC %",
                      guide = guide_colorbar(order = 3))

# Add the points
p <- p + new_scale_fill()
p <- p + geom_point(data = cluster_lcas_summary,
      aes(x = x, y = y, fill = uni_group),
      inherit.aes = FALSE, shape = 21,
      color = "black", stroke = 0.5, size = 2) +
  scale_fill_manual(values = "#c51515", name = "",
          guide = guide_legend(override.aes = list(size = 4)))

# Adjust the legend position to be below the figure
p <- p + theme(legend.position = "top")

p

#ggsave(paste0(output_dir, "denovo_good_bad.png"))




## Dotplot
# Get number of unique species
species_in_cluster_nb <- cluster_df %>%
  group_by(cluster) %>%
  summarise(n_species_denovo = n_distinct(genome), .groups = "drop")
# Add the lca node
species_in_cluster_nb <- species_in_cluster_nb %>%
  left_join(cluster_lcas, by = "cluster") %>%
  select(cluster, n_species_denovo, lca_node)
# Get the number of species below the node
species_in_cluster_nb <- species_in_cluster_nb %>%
  rowwise() %>%
  mutate(n_species_below = length(Descendants(tree, lca_node, type = "tips")[[1]])) %>%
  ungroup()
# Calculate the coverage
species_in_cluster_nb <- species_in_cluster_nb %>%
  mutate(coverage = n_species_denovo / n_species_below)

# Plot
ggplot(species_in_cluster_nb, aes(x = n_species_denovo, y = coverage)) +
  geom_jitter(fill = "#8261dd", alpha = 0.5, color = "#4c3099", size = 4, width = 0.1, height = 0.02) +
  labs(title = "Number of de novo and coverage for each cluster",
       x = "number of species with a denovo gene",
       y = "Cluster coverage (# denovo / # leaves)") +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(breaks = seq(0, max(species_in_cluster_nb$n_species_denovo), 1)) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

#ggsave(paste0(output_dir, "denovo_coverage.png"))



############# Clusters non orfan #############
# Get the clusters with more than one denovo
clusters_non_orphan <- cluster_df %>%
  group_by(cluster) %>%
  summarise(n_denovo = n(), .groups = "drop") %>%
  filter(n_denovo > 1)
# Filter the data to keep only the non-orphan clusters
clusters_non_orphan_df <- cluster_df %>%
  filter(cluster %in% clusters_non_orphan$cluster) %>%
  select(cluster, denovo, genome)
clusters_non_orphan_lcas <- cluster_lcas %>%
  filter(cluster %in% clusters_non_orphan$cluster) %>%
  select(cluster, lca_node)
clusters_non_orphan_lcas_coords <- cluster_lcas_coords %>%
  filter(cluster %in% clusters_non_orphan$cluster) %>%
  select(lca_node, x, y)

genomes_in_clusters <- data.frame(matrix(NA, nrow = 0, ncol = 3))
colnames(genomes_in_clusters) <- c("cluster", "genome", "has_denovo")


### Individual trees ###
# For each cluster, plot the tree
for (i in seq_len(nrow(clusters_non_orphan))) {
  cluster <- clusters_non_orphan[i, "cluster"][[1]]
  # Get the denovo and genomes
  denovo <- clusters_non_orphan_df %>%
    filter(cluster == !!cluster) %>%
    pull(denovo) %>%
    na.omit()
  genomes <- clusters_non_orphan_df %>%
    filter(cluster == !!cluster) %>%
    pull(genome) %>%
    na.omit()
  # Get the lca node
  lca_node <- clusters_non_orphan_lcas %>%
    filter(cluster == !!cluster) %>%
    pull(lca_node) %>%
    na.omit()
  # Get the coordinates
  coords <- clusters_non_orphan_lcas_coords %>%
    filter(lca_node == !!lca_node)
  # Get the species
  tips_under_lcas <- as.vector(tree_data[tree_data$node %in% Descendants(tree, lca_node, type = "tips")[[1]], "label"])$label
  species <- data.frame(
    species = data[tips_under_lcas, "species"],
    row.names = tips_under_lcas
  )


  # Plot the tree
  p <- ggtree(tree, layout = "circular", branch.length = "none") +
    geom_point(data = coords, aes(x = x, y = y), color = "#000000", 
               size = 5, shape = 21, stroke = 0.5, fill = "#c51515") +
    geom_point(aes(x = x, y = y),
               data = tree_data[tree_data$node %in% Descendants(tree, lca_node, type = "tips")[[1]], ],
               color = "#000000", size = 4, shape = 21,
               stroke = 0.5, fill = "#b5b6b9") +
    geom_point(aes(x = x, y = y),
               data = tree_data[tree_data$label %in% genomes, ],
               color = "#000000", size = 4, shape = 21,
               stroke = 0.5, fill = "#179207")

  p <- p + new_scale_fill()
  p <- gheatmap(p, species[, "species", drop = FALSE], offset = 2,
              width = .1, colnames = FALSE, legend_title = NULL) +
  scale_fill_manual(
    guide = "none",
    values = c(viridis::viridis(length(unique(na.omit(species$species))), option = "turbo"), "white"),
    na.value = "white"
  )

  #ggsave(paste0(output_dir, "clusters/tree_", i, ".png"), plot = p)



  ### Analyse conservation for species under lca with no de novo ###
  genomes_under_lca <- unlist(tree_data[tree_data$node %in% Descendants(tree, lca_node, type = "tips")[[1]], "label"])
  genomes_with_denovo <- unlist(tree_data[tree_data$label %in% genomes, "label"])
  for (focal_genome in genomes_under_lca) {
    has_de_novo <- focal_genome %in% genomes_with_denovo
    if (has_de_novo) {
      denovo_name <- denovo[which(genomes == focal_genome)]
      has_de_novo <- denovo_name
    }
    genomes_in_clusters <- rbind(genomes_in_clusters,
                                 data.frame(cluster = cluster,
                                            genome = focal_genome,
                                            has_denovo = has_de_novo))
  }
}
# Export info about the genomes present in each cluster
write.table(genomes_in_clusters, file = paste0(output_dir, "genomes_in_clusters.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# In order to get the conservation results, execute the analyse_non_orfan_clusters.py script



### Individual trees with conservation results ###
# Read the conservation results
