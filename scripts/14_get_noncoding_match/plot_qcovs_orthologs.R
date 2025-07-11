# Plot the qcov distributions for 1000 random orthologs

library(yaml)
library(ggplot2)

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
dense_dir <- paths$dense_dir
genomes_list <- paths$genomes_list
output_dir <- "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/"


extract_qcovs <- FALSE


if (extract_qcovs) {
  # For each genome combination
  genomes <- scan(genomes_list, what = "", sep = " ")[1:116]
  qcovs <- c()
  for (focal_genome in genomes) {
    for (neighbor_genome in genomes) {
      # Skip if the genomes are the same
      if (focal_genome == neighbor_genome) next
      # Else read the blast file
      blast_file_pattern <- paste0(dense_dir, focal_genome,
                                   "/diamondblast_out/",
                                   neighbor_genome,
                                   "_CDS_BLASTp_*.out")
      blast_files <- Sys.glob(blast_file_pattern)
      blast_file <- blast_files[1]
      data <- read.table(blast_file, header = FALSE, sep = "\t")
      colnames(data) <- c("qseqid", "sseqid", "pident", "length",
                          "mismatch", "gapopen", "qstart",
                          "qend", "sstart", "send", "evalue",
                          "bitscore", "qlen", "qcov")
      # Remove ann non matches
      data <- data[!data$sseqid == "*",]
      # Remove all evals > 1e-5
      data <- data[data$evalue < 1e-5, ]
      # Add the qcovs to the list
      qcovs_data <- data$qcov
      qcovs <- c(qcovs, qcovs_data)
    }
  }
  # Sample 10k qcovs
  print(paste0("Number of qcovs extracted: 10000/", length(qcovs)))
  if (length(qcovs) > 10000) {
    qcovs <- sample(qcovs, 10000)
  }
  # Print qcovs to a file (space separated)
  qcovs_file <- paste0(output_dir, "qcovs_orthologs.txt")
  write(qcovs, file = qcovs_file)
} else {
  # If qcovs are already extracted, read them from the file
  qcovs_file <- paste0(output_dir, "qcovs_orthologs.txt")
  qcovs <- read.table(qcovs_file, header = FALSE)
}


p <- ggplot(qcovs, aes(x = V1)) +
  geom_histogram(fill = "#8d8dd8", color = "black") +
  labs(title = "Distribution of qcovs for orthologs (n = 10000)",
       x = "query coverage (%)", y = "Count") +
  theme_minimal()

# Calculate and print deciles
deciles <- quantile(qcovs$V1, probs = seq(0, 1, 0.1))
print(paste0("70% of qcovs are over ", deciles[4]))
print(paste0("80% of qcovs are over ", deciles[3]))
print(paste0("90% of qcovs are over ", deciles[2]))

# Add them to plot
p <- p + geom_vline(xintercept = deciles[4], linetype = "dashed", color = "red") +
  geom_vline(xintercept = deciles[3], linetype = "dashed", color = "orange") +
  geom_vline(xintercept = deciles[2], linetype = "dashed", color = "green") +
  annotate("text", x = deciles[4] + 5, y = 3800, label = paste0("70%\n(qcov\n", round(deciles[4], 2), "%)"), color = "red") +
  annotate("text", x = deciles[3] + 5, y = 3800, label = paste0("80%\n(qcov\n", round(deciles[3], 2), "%)"), color = "orange") +
  annotate("text", x = deciles[2] + 5, y = 3800, label = paste0("90%\n(qcov\n", round(deciles[2], 2), "%)"), color = "green")

print(p)
ggsave(paste0(output_dir, "qcovs_orthologs_histogram.png"))

