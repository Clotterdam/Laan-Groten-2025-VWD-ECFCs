#make a pcaplot of ECFC clones
library(openxlsx)
library(readxl)
#read data:
##Read datafile
df.heatmap <- read_excel("I:/Hemostase/Folders Medewerkers TH/Laan, S.N.J/1. Experiments/087-1 full characterization no VWF mutation/2. RNA of ECFC clones (exp 084-1)/heatmap_data.xlsx",.name_repair = "universal")
#READ THE COMBINED EXCEL HERE other info on sample numbers
excel_file <-
  "I:/Hemostase/Folders Medewerkers TH/Laan, S.N.J/1. Experiments/087-1 full characterization no VWF mutation/12. Combined results/Stim_Scratch_Imaging_Plasma v2.xlsx"
merged_data <- read_excel(excel_file, .name_repair = "universal")
#make a database with no NA
#remove all rows with duplicate info (as 1 patient will have multiple clones and multiple conditions in the scratch.)
merged_data <- merged_data[merged_data$Cluster %in% c(1, 2), ]
merged_data <- merged_data[merged_data$Mito == "Medium", ]

# Handle missing values
# Impute missing values with column means
df.heatmap_imputed <- apply(df.heatmap, 2, function(x) ifelse(is.na(x), mean(x, na.rm = TRUE), x))

# Define the desired order of groups
desired_order <- c("Control", "no DDAVP response", "no VWF mutation")  # Add your desired order of groups here
# Reorder the levels of the factor variable according to the desired order
merged_data$patient_group <- factor(merged_data$patient_group, levels = desired_order)
# Extract unique clusters and assign colors
groups <- unique(merged_data$patient_group)
patient_group_colors <- c("Control" = "#F3E660", 
                          "no DDAVP response" = "#C73E4C", 
                          "no VWF mutation" = "#5C136E")

# Map colors to clusters
group_color_map <- setNames(patient_group_colors, groups)
point_colors <- group_color_map[merged_data$patient_group]

# Perform PCA
pca_result <- prcomp(df.heatmap_imputed, scale. = TRUE)

# Identify rows to exclude
rows_to_exclude <- c(8, 10,17,19,21,33, 31 ,16,36)

# Plot PCA
plot(pca_result$x[-rows_to_exclude, 1], pca_result$x[-rows_to_exclude, 2], 
     xlab = "Principal Component 1", ylab = "Principal Component 2",
     main = "PCA Plot",
     col = point_colors[-rows_to_exclude], pch = 16, cex = 2
)
abline(h = 0, v = 0, col = "gray", lty = 2)

# Add row numbers for rows not excluded
text_labels <- as.character((1:nrow(pca_result$x))[-rows_to_exclude])
text(pca_result$x[-rows_to_exclude, 1], pca_result$x[-rows_to_exclude, 2], 
     labels = text_labels, pos = 4, col = "black", cex = 1, )

# Draw lines connecting points to text labels
for (i in 1:length(text_labels)) {
  segments(x0 = pca_result$x[-rows_to_exclude, 1][i], y0 = pca_result$x[-rows_to_exclude, 2][i],
           x1 = pca_result$x[-rows_to_exclude, 1][i], y1 = pca_result$x[-rows_to_exclude, 2][i] - 0.2,
           col = "black")
}

# Check if there are any missing or invalid values in the "Cluster" column
if (any(is.na(merged_data$patient_group)) || any(!merged_data$patient_group %in% clusters)) {
  warning("Missing or invalid values in 'Cluster' column.")
} else {
  # Add legend if no issues found
  legend("topright", legend = levels(factor(merged_data$patient_group)), fill = patient_group_colors, title = "Cluster")
}

# Add percentage of explained variance
explained_var <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
text(0.2, -3, paste("PC1:", explained_var[1], "%"), adj = c(0, 0))
text(-4, 0.85, paste("PC2:", explained_var[2], "%"), adj = c(0, 0))



# Convert PCA result to a data frame to make a better figure.
pca_df <- data.frame(PC1 = pca_result$x[,1], PC2 = pca_result$x[,2], row_number = 1:nrow(pca_result$x))
library(ggrepel)
# Plot PCA using ggplot2
pcaplot <- ggplot(pca_df[-rows_to_exclude,], aes(x = PC1, y = PC2, label = row_number)) +
  geom_point(aes(fill = merged_data$patient_group[-rows_to_exclude]), 
             shape = 21,
             stroke = 1.1,
             size = 5.5) +
  #geom_text(position = position_dodge(width = .5), size = 5) +
  labs(x = "PCA1 (48.98%)", y = "PCA2 (13.14%)", title = "PCA Plot") +
  scale_color_manual(values = patient_group_colors, name = "" ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
  geom_text_repel(size = 5, segment.size = 1, segment.color = "black") +  # Use geom_text_repel for scattered labels
  theme_minimal()+
  scale_color_manual(values = patient_group_colors, name = "Patient group") +
  scale_fill_manual(values = patient_group_colors, name = "Patient group")
pcaplot

ggsave("PCA v2 .png", plot = pcaplot, width = 10, height = 5)  # Adjust width and height as needed

