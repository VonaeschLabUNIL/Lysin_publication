# ===== Load Required Libraries =====
library(tidyverse)
library(readr)
library(ggplot2)
library(gridExtra)
library(ggpubr)

# ===== 1. Define Colors for Each Medium (edit freely) =====
medium_palette <- c(
  "BHI" = "black",  # Reference control
  "AF9" = "#1f78b4",  # dark blue
  "AF10" = "#ff7f00",  # dark orange
  "AF23" = "#33a02c",  # rich green
  "AF45" = "#e31a1c",  # strong red
  "AF79" = "#6a3d9a",  # purple
  "AF94" = "#b15928",  # brown
  "AF111" = "#008080", # dark teal
  "AF236" = "#636363", # charcoal
  "AF316" = "#7b3294", # dark violet
  "K12" = "#ff1493"   # deep pink
)

# ===== 2. Load Metadata and Donor-Species Mapping =====
metadata <- read_tsv("EXP4_Metadata_Table.tsv")
donor_map <- read_tsv("CORRESPOND_DONOR_SPECIES.tsv")
cluster_map <- read_tsv("STRAIN_CLUSTERS_COLORS.tsv")

metadata <- metadata %>%
  mutate(
    WELL = toupper(WELL),
    PLATE_NAME = toupper(PLATE_NAME),
    STRAIN = LABID
  )

# ===== Set Plot Title Mode (choose "DONOR" or "SPECIES_NAME") =====
title_mode <- "SPECIES_NAME"

# ===== 3. Load Plate Files =====
plate_files <- paste0("PLATE_", 1:34, ".tsv")

# ===== 4. Function to Process One Plate =====
process_plate <- function(file_path) {
  plate_name <- tools::file_path_sans_ext(basename(file_path)) %>% toupper()
  message("Processing ", plate_name, "...")
  
  od <- read_tsv(file_path)
  od <- od %>% mutate(across(everything(), ~ str_replace_all(.x, ",", ".")))
  
  od_long <- od %>%
    rename(Time = 1) %>%
    pivot_longer(-Time, names_to = "Well", values_to = "OD") %>%
    mutate(
      Time = as.numeric(Time),
      OD = as.numeric(OD),
      Plate = plate_name,
      Well = toupper(Well)
    )
  
  merged <- od_long %>%
    left_join(metadata, by = c("Plate" = "PLATE_NAME", "Well" = "WELL"))
  
  return(merged)
}

# ===== 5. Process All Plates =====
all_plates <- lapply(plate_files, process_plate)
final_df <- bind_rows(all_plates)

# ===== 6. Improved Blank Subtraction (Plate-aware) =====
blank_means <- final_df %>%
  filter(toupper(STRAIN) == "BLANK") %>%
  group_by(Plate, MEDIUM, Time) %>%
  summarise(blank_mean = mean(OD, na.rm = TRUE), .groups = "drop")

corrected_df <- final_df %>%
  left_join(blank_means, by = c("Plate", "MEDIUM", "Time")) %>%
  mutate(
    OD_corrected = OD - blank_mean,
    OD_corrected = ifelse(OD_corrected < 0, 0, OD_corrected)
  )

# ===== 7. Compute Mean and SD Per Strain/Medium/Time (Excl. BLANK) =====
stats_df <- corrected_df %>%
  filter(toupper(STRAIN) != "BLANK") %>%
  group_by(STRAIN, MEDIUM, Time, TREATMENT) %>%
  summarise(
    mean_OD = mean(OD_corrected, na.rm = TRUE),
    sd_OD = sd(OD_corrected, na.rm = TRUE),
    .groups = "drop"
  )

# ===== 8. EXPORT: Per-Donor Panels Showing All LABID Strains =====
output_folder <- "FOLDER_STREP_IN_SPENT"
dir.create(output_folder, showWarnings = FALSE)

labid_strains <- c("AF9", "AF10", "AF23", "AF45", "AF79",
                   "AF94", "AF111", "AF236", "AF316", "K12")

donor_list <- metadata %>%
  filter(!is.na(TREATMENT), TREATMENT != "BHI", TREATMENT != "EMPTY", MEDIUM != "EMPTY") %>%
  distinct(TREATMENT) %>%
  pull(TREATMENT)

max_y <- max(stats_df$mean_OD + stats_df$sd_OD, na.rm = TRUE) * 1.1

donor_plots <- list()

for (donor in donor_list) {
  df_subset <- stats_df %>%
    filter((TREATMENT == donor | MEDIUM == "BHI") & STRAIN %in% labid_strains)
  
  if (nrow(df_subset) == 0) {
    message("Skipping ", donor, ": no data found.")
    next
  }
  
  label_df <- df_subset %>%
    group_by(STRAIN, MEDIUM) %>%
    filter(Time == max(Time, na.rm = TRUE)) %>%
    ungroup()
  
  species_label <- donor_map %>%
    filter(DONOR == donor) %>%
    pull(SPECIES_NAME) %>%
    unique()
  
  if (title_mode == "SPECIES_NAME" && length(species_label) > 0) {
    message("Using species name for ", donor, ": ", species_label)
    plot_title <- paste("Growth of LABID strains in", species_label)
  } else {
    plot_title <- paste("Growth of LABID strains in donor", donor, "media")
  }
  
  p <- ggplot(df_subset, aes(x = Time, y = mean_OD, color = STRAIN, fill = STRAIN)) +
    #geom_ribbon(aes(ymin = mean_OD - sd_OD, ymax = mean_OD + sd_OD), alpha = 0.2, color = NA) +
    geom_line(size = 0.4, alpha = 0.9) +
    geom_text(data = label_df, aes(label = STRAIN), size = 1.2, hjust = 0, vjust = 1.2, show.legend = FALSE) +
    facet_wrap(~MEDIUM, ncol = 3) +
    scale_color_manual(values = medium_palette) +
    scale_fill_manual(values = medium_palette) +
    coord_cartesian(xlim = c(0, max(stats_df$Time, na.rm = TRUE) + 1), ylim = c(0, max_y)) +
    labs(
      title = plot_title,
      x = "Time (h)",
      y = "Corrected Mean OD600",
      color = "Strain"
    ) +
    theme_classic()
  
  ggsave(
    filename = file.path(output_folder, paste0("AllStrains_in_", donor, ".pdf")),
    plot = p,
    width = 10,
    height = 6
  )
  
  donor_plots[[donor]] <- p
}

# ===== 9. EXPORT: Combined Panel with All Donors in Grid (4 per row) =====
if (length(donor_plots) > 0) {
  combined_panel <- ggpubr::ggarrange(
    plotlist = donor_plots,
    ncol = 4,
    nrow = ceiling(length(donor_plots) / 4),
    common.legend = TRUE,
    legend = "bottom"
  )
  
  ggsave(
    filename = "AllStrains_AllDonors_Combined.pdf",
    plot = combined_panel,
    width = 40,
    height = 6 * ceiling(length(donor_plots) / 4),
    limitsize = FALSE
  )
}



# ===== 10. EXPORT: Per-Donor Panels with Strains Colored by Cluster Colors from File =====
cluster_output_folder <- "FOLDER_CLUSTER_COLOR"
dir.create(cluster_output_folder, showWarnings = FALSE)

strain_colors_clustered <- cluster_map %>%
  filter(LABID %in% unique(stats_df$STRAIN)) %>%
  select(LABID, COLOR) %>%
  deframe()

print("Using the following strain color mapping (by cluster):")
print(strain_colors_clustered)

cluster_plots <- list()

donor_list <- unique(stats_df$TREATMENT)

for (donor in donor_list) {
  df_subset <- stats_df %>%
    filter((TREATMENT == donor | MEDIUM == "BHI") & STRAIN %in% names(strain_colors_clustered))
  
  if (nrow(df_subset) == 0) {
    next
  }
  
  label_df <- df_subset %>%
    group_by(STRAIN, MEDIUM) %>%
    filter(Time == max(Time, na.rm = TRUE)) %>%
    ungroup()
  
  species_label <- donor_map %>%
    filter(DONOR == donor) %>%
    pull(SPECIES_NAME) %>%
    unique()
  
  plot_title <- if (title_mode == "SPECIES_NAME" && length(species_label) > 0) {
    paste("Growth (Cluster-colored) in", species_label)
  } else {
    paste("Growth (Cluster-colored) in donor", donor, "media")
  }
  
  p <- ggplot(df_subset, aes(x = Time, y = mean_OD, color = STRAIN, fill = STRAIN)) +
    #geom_ribbon(aes(ymin = mean_OD - sd_OD, ymax = mean_OD + sd_OD), alpha = 0.2, color = NA) +
    geom_line(size = 0.4, alpha = 0.9) +
    geom_text(data = label_df, aes(label = STRAIN), size = 1.2, hjust = 0, vjust = 1.2, show.legend = FALSE) +
    facet_wrap(~MEDIUM, ncol = 3) +
    scale_color_manual(values = strain_colors_clustered) +
    scale_fill_manual(values = strain_colors_clustered) +
    coord_cartesian(xlim = c(0, max(stats_df$Time, na.rm = TRUE) + 1), ylim = c(0, max(stats_df$mean_OD + stats_df$sd_OD, na.rm = TRUE) * 1.1)) +
    labs(
      title = plot_title,
      x = "Time (h)",
      y = "Corrected Mean OD600",
      color = "Strain"
    ) +
    theme_classic()
  
  ggsave(
    filename = file.path(cluster_output_folder, paste0("ClusterColor_AllStrains_in_", donor, ".pdf")),
    plot = p,
    width = 10,
    height = 6
  )
  
  cluster_plots[[donor]] <- p
}



# ===== 11. EXPORT: SPENT vs BHI Only — Cluster Colors + Strain Names in Facets =====
spent_output_folder <- "ONLY_SPENT_VS_BHI"
dir.create(spent_output_folder, showWarnings = FALSE)

# Create donor → label map: "AFB204: Ruminococcus_lactaris"
donor_strain_label <- metadata %>%
  filter(!is.na(TREATMENT)) %>%
  distinct(DONOR = TREATMENT) %>%
  left_join(donor_map, by = "DONOR") %>%
  mutate(label = paste0(DONOR, ": ", SPECIES_NAME)) %>%
  distinct(DONOR, label) %>%
  deframe()


spent_plots <- list()

for (donor in donor_list) {
  df_subset <- stats_df %>%
    filter((TREATMENT == donor | MEDIUM == "BHI") & STRAIN %in% names(strain_colors_clustered)) %>%
    filter(!str_detect(MEDIUM, "repl|REPL|Repl"))  # exclude replenished media
  
  if (nrow(df_subset) == 0) next
  
  # Assign custom facet labels: "SPENT (AFB204: R. lactaris)" and keep BHI clean
  facet_label <- donor_strain_label[[donor]]
  df_subset <- df_subset %>%
    mutate(
      MEDIUM = ifelse(MEDIUM != "BHI", paste0("SPENT (", facet_label, ")"), "BHI")
    )
  
  label_df <- df_subset %>%
    group_by(STRAIN, MEDIUM) %>%
    filter(Time == max(Time, na.rm = TRUE)) %>%
    ungroup()
  
  species_label <- donor_map %>%
    filter(DONOR == donor) %>%
    pull(SPECIES_NAME) %>%
    unique()
  
  plot_title <- if (title_mode == "SPECIES_NAME" && length(species_label) > 0) {
    paste("Growth in SPENT vs BHI for", species_label)
  } else {
    paste("Growth in SPENT vs BHI for donor", donor)
  }
  
  p <- ggplot(df_subset, aes(x = Time, y = mean_OD, color = STRAIN, fill = STRAIN)) +
    geom_line(size = 0.4, alpha = 0.9) +
    geom_text(data = label_df, aes(label = STRAIN), size = 1.2, hjust = 0, vjust = 1.2, show.legend = FALSE) +
    facet_wrap(~MEDIUM, ncol = 2) +
    scale_color_manual(values = strain_colors_clustered) +
    scale_fill_manual(values = strain_colors_clustered) +
    coord_cartesian(
      xlim = c(0, max(stats_df$Time, na.rm = TRUE) + 1),
      ylim = c(0, max(stats_df$mean_OD + stats_df$sd_OD, na.rm = TRUE) * 1.1)
    ) +
    labs(
      #title = plot_title,
      x = "Time (h)",
      y = "Corrected Mean OD600",
      color = "Strain"
    ) +
    theme_classic()
  
  ggsave(
    filename = file.path(spent_output_folder, paste0("SPENT_vs_BHI_in_", donor, ".pdf")),
    plot = p,
    width = 8,
    height = 5
  )
  
  spent_plots[[donor]] <- p
}

# ===== 12. COMBINED PANEL: SPENT vs BHI All Donors =====
if (length(spent_plots) > 0) {
  combined_panel_spent <- ggpubr::ggarrange(
    plotlist = spent_plots,
    ncol = 6,
    nrow = ceiling(length(spent_plots) / 6),
    common.legend = TRUE,
    legend = "bottom"
  )
  
  ggsave(
    filename = file.path(spent_output_folder, "SPENT_vs_BHI_AllDonors_Combined.pdf"),
    plot = combined_panel_spent,
    width = 60,
    height = 6 * ceiling(length(spent_plots) / 6),
    limitsize = FALSE
  )
}


# ===== 13. COMBINED PANEL: One BHI Reference (Top Left) + SPENT-Only Panels =====

spent_with_bhi_folder <- file.path(spent_output_folder, "SPENT_ONLY_WITH_BHI_TOPLEFT")
dir.create(spent_with_bhi_folder, showWarnings = FALSE)

# ---- Donor to Label Map ----
donor_strain_label <- metadata %>%
  filter(!is.na(TREATMENT)) %>%
  distinct(DONOR = TREATMENT) %>%
  left_join(donor_map, by = "DONOR") %>%
  mutate(label = paste0(DONOR, ": ", SPECIES_NAME)) %>%
  distinct(DONOR, label) %>%
  deframe()

# ---- BHI Panel ----
bhi_plot_df <- stats_df %>%
  filter(
    MEDIUM == "BHI",
    toupper(STRAIN) != "BLANK"
  ) %>%
  mutate(MEDIUM = "BHI")

# Optional: manually assign color if needed
bhi_plot_df <- bhi_plot_df %>%
  filter(STRAIN %in% names(strain_colors_clustered))


p_bhi <- NULL
if (nrow(bhi_plot_df) > 0 && any(!is.na(bhi_plot_df$MEDIUM))) {
  label_df_bhi <- bhi_plot_df %>%
    group_by(STRAIN, MEDIUM) %>%
    filter(Time == max(Time, na.rm = TRUE)) %>%
    ungroup()
  
  p_bhi <- ggplot(bhi_plot_df, aes(x = Time, y = mean_OD, color = STRAIN, fill = STRAIN)) +
    geom_line(size = 0.5, alpha = 0.9) +
    geom_text(data = label_df_bhi, aes(label = STRAIN), size = 1.4, hjust = 0, vjust = 1.2, show.legend = FALSE) +
    facet_wrap(~MEDIUM, ncol = 1) +
    scale_color_manual(values = strain_colors_clustered) +
    scale_fill_manual(values = strain_colors_clustered) +
    coord_cartesian(xlim = c(0, max(stats_df$Time, na.rm = TRUE) + 1), ylim = c(0, 1)) +
    labs(x = "Time (h)", y = "Corrected Mean OD600", color = "Strain") +
    theme_classic() +
    theme(strip.text = element_text(size = 14, face = "bold"))
}

# ---- SPENT-Only Panels ----
spent_only_plots <- list()

for (donor in donor_list) {
  if (donor == "BHI") next
  
  df_subset <- stats_df %>%
    filter(
      TREATMENT == donor,
      STRAIN %in% names(strain_colors_clustered),
      MEDIUM != "BHI",
      !str_detect(MEDIUM, "repl|REPL|Repl")
    )
  
  if (nrow(df_subset) == 0 || all(is.na(df_subset$MEDIUM))) next
  
  facet_label <- donor_strain_label[[donor]]
  if (is.na(facet_label) || is.null(facet_label)) next
  
  df_subset <- df_subset %>%
    mutate(MEDIUM = paste0("SPENT (", facet_label, ")"))
  
  if (all(is.na(df_subset$MEDIUM))) next
  
  label_df <- df_subset %>%
    group_by(STRAIN, MEDIUM) %>%
    filter(Time == max(Time, na.rm = TRUE)) %>%
    ungroup()
  
  p <- ggplot(df_subset, aes(x = Time, y = mean_OD, color = STRAIN, fill = STRAIN)) +
    geom_line(size = 0.4, alpha = 0.9) +
    geom_text(data = label_df, aes(label = STRAIN), size = 1.2, hjust = 0, vjust = 1.2, show.legend = FALSE) +
    facet_wrap(~MEDIUM, ncol = 1) +
    scale_color_manual(values = strain_colors_clustered) +
    scale_fill_manual(values = strain_colors_clustered) +
    coord_cartesian(xlim = c(0, max(stats_df$Time, na.rm = TRUE) + 1), ylim = c(0, 1)) +
    labs(x = "Time (h)", y = "Corrected Mean OD600", color = "Strain") +
    theme_classic() +
    theme(strip.text = element_text(size = 12, face = "bold"))
  
  spent_only_plots[[donor]] <- p
}

# ---- Final Clean-up ----
valid_spent <- Filter(function(p) {
  !is.null(p) && !is.null(p$data) &&
    nrow(p$data) > 0 &&
    "MEDIUM" %in% colnames(p$data) &&
    any(!is.na(p$data$MEDIUM))
}, spent_only_plots)

# ---- Combine and Save ----
if (!is.null(p_bhi) || length(valid_spent) > 0) {
  plot_list <- c(if (!is.null(p_bhi)) list(p_bhi) else list(), valid_spent)
  
  combined_panel <- ggpubr::ggarrange(
    plotlist = plot_list,
    ncol = 6,
    nrow = ceiling(length(plot_list) / 5),
    common.legend = TRUE,
    legend = "bottom"
  )
  
  ggsave(
    filename = file.path(spent_with_bhi_folder, "SPENT_Only_With_BHI_TopLeft.pdf"),
    plot = combined_panel,
    width = 30,
    height = 6 * ceiling(length(plot_list) / 6),
    limitsize = TRUE
  )
  
  # Match original panel size
  custom_width <- 30
  custom_height <- 6 * ceiling(length(plot_list) / 6)
  
  cairo_pdf(
    filename = file.path(spent_with_bhi_folder, "SPENT_Only_With_BHI_TopLeft_FullSize_cairo.pdf"),
    width = custom_width,
    height = custom_height,
    family = "sans"
  )
  print(combined_panel)
  dev.off()
  
  # ---- High-Resolution TIFF Export (Huge Version) ----
  ggsave(
    filename = file.path(spent_with_bhi_folder, "SPENT_Only_With_BHI_TopLeft_Huge_600dpi.tiff"),
    plot = combined_panel,
    width = 30,
    height = 6 * ceiling(length(plot_list) / 6),
    units = "in",
    dpi = 600,
    device = "tiff",
    compression = "lzw"
  )
  
 
} else {
  message("No valid BHI or SPENT plots to display.")
}
