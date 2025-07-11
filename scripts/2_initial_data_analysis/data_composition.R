library(ggplot2)
library(dplyr)
library(tidyr)
ggsave <- function(..., bg = "white", width = 1000, height = 1000, units = "px", dpi = 100) {
  ggplot2::ggsave(..., bg = bg, width = width, height = height, units = units, dpi = dpi)
}

# Import the file
file <- "/home/eliott.tempez/Documents/archaea_data/Thermocomplete_table.csv"
data <- read.csv(file)

# Which columns are filled for at least 95% of all organisms
full_col <- c()
len_max <- dim(data)[1] * .95
for (col in colnames(data)) {
  if (length(which(data[, col] != "")) > len_max) {
    full_col <- c(full_col, col)
  }
}
other_col <- setdiff(colnames(data), full_col)

# Delete non informative columns
data_high <- data[, setdiff(full_col, c("Genome", "Geographic.Region"))]
data_low <- cbind(data[,"Name"], data[, setdiff(other_col, c("Assembly.accession", "acces_EGM", "EGM_name"))])



## Start with columns for which we have a good amount of info
# Genome size
ggplot(data_high, aes(x = Genome_Size_Mb)) +
  geom_histogram(binwidth = 0.1, fill = "#A6C3AA", color = "black") +
  labs(x = "Genome size (Mb)", y = "Count", title = "Distribution of the genome sizes in the dataset") +
  theme_minimal()
# export figure
ggsave("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/initial_data_analysis_figures/genome_size.png")


# GC content
ggplot(data_high, aes(x = GC_content.)) +
  geom_histogram(binwidth = 1, fill = "#C3A6AC", color = "black") +
  labs(x = "GC content (%)", y = "Count", title = "Distribution of the GC content in the dataset") +
  theme_minimal()
# export figure
ggsave("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/initial_data_analysis_figures/GC_content.png")

# Average CDS and intergenic length
# Combine the data into a long format for ggplot
data_long <- data_high %>%
  select(Average_CDS_length, Average_intergenic_length) %>%
  pivot_longer(cols = everything(), names_to = "Type", values_to = "Length")

# Plot histograms with ggplot
ggplot(data_long, aes(x = Length, fill = Type)) +
  geom_histogram(position = "identity", alpha = 0.75, color = "black") +
  scale_fill_manual(values = c("Average_CDS_length" = "#C3A6C3", "Average_intergenic_length" = "#C3C3A6")) +
  labs(x = "Average length (nucl)", y = "Count", title = "Average CDS and intergenic region length in the dataset") +
  theme_minimal() +
  theme(legend.title = element_blank())
# export figure
ggsave("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/initial_data_analysis_figures/dna_legth.png")


# Map of the geographic regions of extraction
world <- map_data("world")
ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "#4e4444", fill = "#e2eef6"
  ) +
  geom_point(
    data = data,
    aes(Geographic_location_longitude, Geographic_location_latitude),
    alpha = 1, color = "#db0000",
    shape = 8, size = 3
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

ggsave("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/initial_data_analysis_figures/map.png", width = 1500, height = 1000)


# Do we have the complete chromosome?
ggplot(data_high, aes(x = tolower(Complete_chromosome))) +
  geom_bar(fill = "#a7c3a6", color = "black") +
  labs(x = "Complete chromosome", y = "Count", title = "Distribution of the presence of a complete chromosome in the dataset") +
  theme_minimal()
  ggsave("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/initial_data_analysis_figures/complete.chr.png")


# Protein coding density
ggplot(data_high, aes(x = X..protein.codingDensity)) +
  geom_histogram(fill = "#c0c3a6", color = "black") +
  labs(x = "Protein coding density", y = "Count", title = "Distribution of the protein coding density in the dataset") +
  theme_minimal()
   ggsave("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/initial_data_analysis_figures/protein_coding_density.png")


# Depth
data_low$depth <- gsub(",", "", data_low$depth)
data_low$depth <- as.numeric(gsub("[a-zA-Z]", "", data_low$depth))
depth_na <- sum(is.na(data_low$depth))

# plot boxplot
ggplot(data_low, aes(y = depth)) +
  geom_boxplot(fill = "#a8a6c3", color = "black", width = 0.9) +
  xlim(-1, 1) +
  labs(y = "Depth (meters)", x = "",
        title = paste0("Collection depth (unknown : ", depth_na, "/122)")) +
  theme_minimal()
  ggsave("/home/eliott.tempez/Documents/M2_Stage_I2BC/results/initial_data_analysis_figures/depth_boxplot.png")
