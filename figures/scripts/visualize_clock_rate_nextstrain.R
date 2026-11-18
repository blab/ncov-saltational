library(jsonlite)
library(dplyr)
library(ggplot2)

json_file <- '~/Downloads/nextstrain_ncov_open_europe_all-time@2022-05-01.json'

# Fixed cutoffs for the two presentation figures.
early_cutoff <- 2020.75  # approximately 30 September 2020
full_cutoff <- 2021.5

x <- fromJSON(json_file, simplifyVector = FALSE)


get_attr_value <- function(node, attr_name) {
  attrs <- node$node_attrs
  if (is.null(attrs) || is.null(attrs[[attr_name]])) return(NA_real_)
  
  value <- attrs[[attr_name]]
  if (is.list(value) && !is.null(value$value)) value <- value$value
  
  value <- suppressWarnings(as.numeric(value))
  if (length(value) == 0 || is.na(value)) NA_real_ else value[[1]]
}

get_attr_label <- function(node, attr_names) {
  attrs <- node$node_attrs
  if (is.null(attrs)) return(NA_character_)
  
  for (attr_name in attr_names) {
    if (!is.null(attrs[[attr_name]])) {
      value <- attrs[[attr_name]]
      if (is.list(value) && !is.null(value$value)) value <- value$value
      if (length(value) > 0 && !is.null(value[[1]])) {
        value <- as.character(value[[1]])
        if (nzchar(value)) return(value)
      }
    }
  }
  NA_character_
}

node_info <- function(node) {
  data.frame(
    name = if (!is.null(node$name)) node$name else NA_character_,
    num_date = get_attr_value(node, "num_date"),
    div = get_attr_value(node, "div"),
    clade = get_attr_label(node, c("clade", "clade_membership")),
    stringsAsFactors = FALSE
  )
}

# Recursively collect all nodes and parent-child edges, retaining both clade
# labels so that the stem branch entering 20I can be identified. Each edge also
# records the earliest date among its descendant tips. This lets us suppress
# branches in Figure 1 that do not lead to any tip plotted in that figure.
walk_tree <- function(node, parent = NULL) {
  current <- node_info(node)
  children <- node$children
  is_tip <- is.null(children) || length(children) == 0
  current$is_tip <- is_tip
  
  nodes <- list(current)
  edges <- list()
  
  if (is_tip) {
    # These are the same criteria used later when constructing `tips`.
    has_plottable_coordinates <-
      is.finite(current$num_date) &&
      is.finite(current$div) &&
      !is.na(current$clade) &&
      nzchar(current$clade)
    
    descendant_tip_dates <- if (has_plottable_coordinates) {
      current$num_date
    } else {
      numeric(0)
    }
  } else {
    descendant_tip_dates <- numeric(0)
    
    for (child in children) {
      result <- walk_tree(child, current)
      nodes <- c(nodes, result$nodes)
      edges <- c(edges, result$edges)
      descendant_tip_dates <- c(
        descendant_tip_dates,
        result$descendant_tip_dates
      )
    }
  }
  
  if (!is.null(parent)) {
    edges[[length(edges) + 1L]] <- data.frame(
      parent_num_date = parent$num_date,
      parent_div = parent$div,
      parent_clade = parent$clade,
      child_num_date = current$num_date,
      child_div = current$div,
      child_clade = current$clade,
      earliest_plottable_tip_date = if (length(descendant_tip_dates) > 0) {
        min(descendant_tip_dates)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }
  
  list(
    nodes = nodes,
    edges = edges,
    descendant_tip_dates = descendant_tip_dates
  )
}

walked <- walk_tree(x$tree)

nodes <- bind_rows(walked$nodes) %>%
  filter(is.finite(num_date), is.finite(div), num_date <= full_cutoff)

branches <- bind_rows(walked$edges) %>%
  filter(
    is.finite(parent_num_date), is.finite(parent_div),
    is.finite(child_num_date), is.finite(child_div),
    child_num_date <= full_cutoff
  )

# Accept either the full Auspice label or the short clade code.
is_20I <- function(clade) {
  !is.na(clade) & clade %in% c("20I (Alpha, V1)", "20I")
}

branches <- branches %>%
  mutate(
    # TRUE only for the branch crossing from outside 20I into 20I.
    highlight_20I = is_20I(child_clade) & !is_20I(parent_clade),
    plot_group = if_else(highlight_20I, "Alpha", "Other")
  )

tips <- nodes %>%
  filter(is_tip, !is.na(clade), nzchar(clade)) %>%
  mutate(plot_group = if_else(is_20I(clade), "Alpha", "Other"))

# Fit ONCE using only the early data. This exact model is used in both figures.
early_nodes <- nodes %>% filter(num_date <= early_cutoff)

fit <- lm(div ~ num_date, data = early_nodes)
fit_summary <- summary(fit)
clock_rate <- unname(coef(fit)[["num_date"]])
clock_se <- unname(coef(summary(fit))["num_date", "Std. Error"])
clock_p <- unname(coef(summary(fit))["num_date", "Pr(>|t|)"])
r_squared <- unname(fit_summary$r.squared)

# Use one shared set of limits and one shared regression prediction grid so the
# two exported images have the same geometry and can be overlaid in slides.
plot_x_limits <- range(nodes$num_date, finite = TRUE)
plot_x_limits[2] <- full_cutoff
plot_y_limits <- range(nodes$div, finite = TRUE)

# The regression line is shown only over the period used to fit it:
# from the earliest plotted date through September 2020.
prediction_grid <- data.frame(
  num_date = seq(plot_x_limits[1], early_cutoff, length.out = 300)
)
prediction_grid$div <- predict(fit, newdata = prediction_grid)

alpha_values <- c(Alpha = 0.85, Other = 0.35)
color_values <- c(Alpha = "darkcyan", Other = "darkgrey")
linewidth_values <- c(Alpha = 1.8, Other = 0.35)

make_plot <- function(max_date, show_alpha = FALSE) {
  plot_nodes <- nodes %>% filter(num_date <= max_date)
  
  # Alpha is deliberately hidden in Figure 1 so that the second figure can
  # introduce the Alpha clade as the visual transition. The branch and tip
  # data remain available in both cases.
  plot_branches <- branches %>%
    # Keep a branch only when its descendant subtree contains at least one
    # plottable tip at or before the figure cutoff. This prevents Figure 1
    # from showing branches that lead exclusively to omitted tips.
    filter(
      child_num_date <= max_date,
      !is.na(earliest_plottable_tip_date),
      earliest_plottable_tip_date <= max_date
    ) %>%
    # In Figure 2, color every branch whose child node belongs to Alpha.
    # Keep `highlight_20I` separate so only the stem entering Alpha is thick.
    mutate(display_alpha = show_alpha & is_20I(child_clade))
  
  plot_tips <- tips %>%
    filter(num_date <= max_date) %>%
    mutate(
      display_group = if_else(show_alpha & is_20I(clade), "Alpha", "Other"),
      display_group = factor(display_group, levels = c("Alpha", "Other"))
    )
  
  plot_prediction <- prediction_grid %>%
    filter(num_date <= max_date)
  
  p <- ggplot() +
    # Draw ordinary branches first.
    geom_segment(
      data = plot_branches %>% filter(!display_alpha),
      aes(
        x = parent_num_date, y = parent_div,
        xend = child_num_date, yend = child_div
      ),
      color = "darkgrey", alpha = 0.35, linewidth = 0.35,
      na.rm = TRUE
    ) +
    # Draw branches within the Alpha clade in teal, using the ordinary width.
    geom_segment(
      data = plot_branches %>% filter(display_alpha),
      aes(
        x = parent_num_date, y = parent_div,
        xend = child_num_date, yend = child_div
      ),
      color = "darkcyan", alpha = 0.85, linewidth = 0.35,
      na.rm = TRUE
    ) +
    # Draw only the stem entering Alpha on top, with greater thickness.
    geom_segment(
      data = plot_branches %>% filter(show_alpha & highlight_20I),
      aes(
        x = parent_num_date, y = parent_div,
        xend = child_num_date, yend = child_div
      ),
      color = "darkcyan", alpha = 0.85, linewidth = 1.8,
      na.rm = TRUE
    ) +
    geom_point(
      data = plot_tips,
      aes(x = num_date, y = div, color = display_group),
      alpha = 1, size = 1.8
    ) +
    # Show the early regression only through September 2020, in black and dotted.
    geom_line(
      data = plot_prediction,
      aes(x = num_date, y = div),
      color = "black", linewidth = 0.9, linetype = "dashed",
      inherit.aes = FALSE
    ) +
    scale_color_manual(
      name = "Clade",
      values = color_values,
      breaks = c("Alpha", "Other"),
      labels = c("Alpha", "Other"),
      drop = FALSE
    ) +
    scale_x_continuous(
      limits = plot_x_limits,
      breaks = seq(2020, full_cutoff, by = 0.5),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = plot_y_limits,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
      x = "Sampling date",
      y = "Root-to-tip divergence"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
  
  p
}

plot_early <- make_plot(
  early_cutoff,
  show_alpha = FALSE
)

plot_full <- make_plot(
  full_cutoff,
  show_alpha = TRUE
)

png('../figure_model_illustration/root_to_tip_divergence_before_alpha.png',
    height = 4., width = 6.5, res = 350, units = 'in')
plot(plot_early)
dev.off()

png('../figure_model_illustration/root_to_tip_divergence_with_alpha.png',
    height = 4., width = 6.5, res = 350, units = 'in')
plot(plot_full)
dev.off()
