# ============================================================
# Circadian Visualization Utilities
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# ============================================================
# 1. Pathway Circadian Curve
# ============================================================

plot_pathway_circadian <- function(
  input_file,
  gene_set_file,
  output_prefix,
  plot_title = NULL
) {
  circadian_df <- read.csv(input_file, sep = "\t", header = TRUE)

  circadian_df <- circadian_df %>%
    dplyr::select(
      GeneID, mesor, sincoef, coscoef,
      acrophase, amplitude, pvalue, qvalue
    )

  gene_list <- readLines(gene_set_file)
  pathway_df <- circadian_df[circadian_df$GeneID %in% gene_list, ]

  write.csv(
    pathway_df,
    paste0(output_prefix, "_cir.csv"),
    row.names = FALSE
  )

  time_points <- seq(0, 24, by = 0.1)

  fit_circadian_curve <- function(t, sincoef, coscoef) {
    sincoef * sin(2 * pi * t / 24) +
      coscoef * cos(2 * pi * t / 24)
  }

  pathway_long <- do.call(
    rbind,
    lapply(1:nrow(pathway_df), function(i) {
      data.frame(
        time = time_points,
        value = fit_circadian_curve(
          time_points,
          pathway_df$sincoef[i],
          pathway_df$coscoef[i]
        )
      )
    })
  )

  pathway_long$time <- rep(time_points, nrow(pathway_df))

  average_values <- aggregate(value ~ time, data = pathway_long, FUN = mean)

  quantiles_time_points <- pathway_long %>%
    group_by(time) %>%
    summarise(
      q25 = quantile(value, 0.25),
      q75 = quantile(value, 0.75)
    )

  p <- ggplot() +
    geom_line(
      data = average_values,
      aes(x = time, y = value),
      color = "#f79b72",
      size = 2
    ) +
    geom_ribbon(
      data = quantiles_time_points,
      aes(x = time, ymin = q25, ymax = q75),
      fill = "#f79b72",
      alpha = 0.2
    ) +
    labs(
      title = plot_title,
      x = "Time (hours)",
      y = "Expression value"
    ) +
    theme_minimal()

  ggsave(
    paste0(output_prefix, "_circadian_curve.pdf"),
    plot = p,
    width = 8,
    height = 6
  )

  return(p)
}

# ============================================================
# 2. Functional Enrichment Polar Heatmap
# ============================================================

plot_polar_enrichment <- function(
  input_file,
  output_file,
  custom_order,
  plot_title = NULL
) {
  go_df <- read.csv(input_file, stringsAsFactors = FALSE)

  go_long <- go_df %>%
    pivot_longer(
      cols = starts_with("LogP_P"),
      names_to = "Phase",
      values_to = "logP"
    ) %>%
    mutate(
      Phase_num = as.numeric(gsub("LogP_P", "", Phase)),
      scaled_logP = ifelse(logP < 0, -logP, 0),
      Description = factor(Description, levels = rev(custom_order))
    )

  p <- ggplot(go_long, aes(x = Phase_num, y = Description, fill = scaled_logP)) +
    geom_tile(color = "grey60", linewidth = 0.3, height = 0.9) +
    scale_x_continuous(
      breaks = 1:8,
      labels = paste0("CT", seq(0, 21, by = 3), "-", seq(3, 24, by = 3))
    ) +
    scale_fill_gradientn(
      colours = c("white", "#c7e0ed", "#6fafd2", "#327db7", "#134b87", "#053061")
    ) +
    coord_polar(theta = "x", start = -pi / 2, clip = "off") +
    labs(title = plot_title, fill = "-log10(P)") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank()
    )

  ggsave(output_file, p, width = 6, height = 6, dpi = 300)

  return(p)
}

# ============================================================
# 3. Acrophase Radial Plot
# ============================================================

plot_acrophase_radial <- function(
  input_files,
  output_file,
  group_order,
  group_colors
) {
  data_list <- lapply(names(input_files), function(group_name) {
    df <- read.csv(input_files[[group_name]], sep = "\t", header = TRUE)

    df %>%
      dplyr::select(acrophase) %>%
      mutate(Group = group_name)
  })

  combined_data <- do.call(rbind, data_list)
  combined_data$Group <- factor(combined_data$Group, levels = group_order)

  p <- ggplot(combined_data, aes(x = acrophase, fill = Group)) +
    geom_histogram(binwidth = 0.15, alpha = 0.8) +
    scale_fill_manual(values = group_colors) +
    coord_polar() +
    scale_x_continuous(
      breaks = seq(0, 24, by = 3),
      limits = c(0, 24)
    ) +
    theme_minimal() +
    labs(title = "Acrophase Distribution")

  ggsave(output_file, p, width = 6, height = 6)

  return(p)
}

# ============================================================
# 4. Single Gene Circadian Curve with Expression
# ============================================================

plot_gene_circadian_expression <- function(
  circadian_file,
  rna_file,
  gene_name,
  subtype_name,
  subtype_map
) {
  circadian_df <- read.csv(circadian_file, sep = "\t", header = TRUE)

  RNA <- read.csv(rna_file, sep = "\t", header = TRUE)

  RNA_clean <- RNA %>%
    separate(col = 1, into = c("Subtype_full", "ZT"), sep = "\\.") %>%
    mutate(
      ZT = as.numeric(str_remove(ZT, "ZT")),
      Subtype = subtype_map[Subtype_full]
    ) %>%
    filter(!is.na(Subtype))

  RNA_long <- RNA_clean %>%
    pivot_longer(
      cols = -c(Subtype_full, ZT, Subtype),
      names_to = "GeneID",
      values_to = "Expression"
    )

  gene_expr <- RNA_long %>%
    filter(GeneID == gene_name & Subtype == subtype_name)

  gene_fit <- circadian_df[circadian_df$GeneID == gene_name, ]

  time_points <- seq(0, 24, by = 0.1)

  fit_curve <- function(t, sincoef, coscoef, mesor) {
    mesor + sincoef * sin(2 * pi * t / 24) +
      coscoef * cos(2 * pi * t / 24)
  }

  gene_curve <- do.call(
    rbind,
    lapply(1:nrow(gene_fit), function(i) {
      data.frame(
        time = time_points,
        value = fit_curve(
          time_points,
          gene_fit$sincoef[i],
          gene_fit$coscoef[i],
          gene_fit$mesor[i]
        )
      )
    })
  )

  p <- ggplot() +
    geom_point(
      data = gene_expr,
      aes(x = ZT, y = Expression),
      color = "#9071cc"
    ) +
    geom_line(
      data = gene_curve,
      aes(x = time, y = value),
      color = "#9071cc"
    ) +
    theme_minimal() +
    labs(
      title = paste("Circadian expression of", gene_name, "in", subtype_name),
      x = "ZT",
      y = "Expression"
    )

  return(p)
}