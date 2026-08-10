normalize_manuscript_method_label <- function(x) {
  dplyr::recode(
    as.character(x),
    "Fast-Higashi(RWR+cov-conv)" = "Fast-Higashi",
    "Fast-Higashi(RWR+cov-conv)(scaled)" = "Fast-Higashi",
    "Higashi(scaled)" = "Higashi",
    "scHiCluster(scaled)" = "scHiCluster",
    .default = as.character(x)
  )
}

manuscript_method_colors <- c(
  "HiCBZIP-GB(NB)" = "#D1495B",
  "HiCBZIP-N(GS)"  = "#E9C46A",
  "HiCBZIP-N(M)"   = "#F4A261",
  "HiCImpute"      = "#9AA44D",
  "scHiCluster"    = "#4C78A8",
  "Higashi"        = "#7E57C2",
  "Fast-Higashi"   = "#B279C8"
)

manuscript_base_colors <- c(
  "Truth" = "#2A9D8F",
  "Simulated (Raw)" = "#9E9E9E"
)

manuscript_all_colors <- c(manuscript_base_colors, manuscript_method_colors)

manuscript_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey95", color = NA)
    )
}

apply_manuscript_colors <- function(p) {
  present_labels <- ggplot2::ggplot_build(p)$plot$data %>%
    dplyr::select(dplyr::any_of(c("model_label", "method_label", "model", "method"))) %>%
    unlist(use.names = FALSE) %>%
    as.character() %>%
    unique()

  present_labels <- present_labels[present_labels %in% names(manuscript_all_colors)]

  p +
    scale_color_manual(values = manuscript_all_colors, breaks = present_labels, drop = TRUE) +
    scale_fill_manual(values = manuscript_all_colors, breaks = present_labels, drop = TRUE)
}

save_manuscript_plot <- function(plot_obj,
                                 filename,
                                 fig_dir = if (exists("path_here", mode = "function")) path_here("results", "figures") else "results/figures",
                                 width = 12,
                                 height = 6,
                                 dpi = 300,
                                 bg = "white",
                                 save_pdf = F) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  png_path <- file.path(fig_dir, paste0(filename, ".png"))
  ggsave(
    filename = png_path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg
  )

  if (save_pdf) {
    pdf_path <- file.path(fig_dir, paste0(filename, ".pdf"))
    ggsave(
      filename = pdf_path,
      plot = plot_obj,
      width = width,
      height = height,
      bg = bg
    )
  }
}
