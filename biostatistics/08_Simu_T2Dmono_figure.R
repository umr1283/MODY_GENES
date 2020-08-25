### Environment ====================================================================================
library(here)
project_name <- gsub("(.*)_.*", "\\1", basename(here()))
output_directory <- here("Data", "08_Simu_T2Dmono_figure")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")


### Load packages ==================================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(dplyr)
  library(parallel)
  library(data.table)
  library(scales)
  library(ggtext)
  library(glue)
  library(forcats)
  library(patchwork)
})


### Define theme ===================================================================================
options(
  ggplot2.discrete.colour = scale_colour_viridis_d,
  ggplot2.discrete.fill = scale_fill_viridis_d,
  ggplot2.continuous.colour = scale_colour_viridis_c,
  ggplot2.continuous.fill = scale_fill_viridis_c
)
theme_set(
  theme_light() + 
    theme(
      plot.title.position = "plot", 
      plot.caption.position = "plot",
      plot.title = element_markdown(),
      plot.subtitle = element_markdown(face = "italic", size = rel(0.75)),
      plot.caption = element_markdown(),
      axis.text.y = element_markdown(), 
      strip.text.x = element_markdown(),
      strip.text.y = element_markdown(),
      legend.text = element_markdown(),
      panel.grid.minor = element_blank()
    )
)


### Source functions ===============================================================================


### Analysis =======================================================================================
walk(c("sansRemise", "avecRemise"), function(type_resampling) {
  csv_file <- glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz"))
  if (file.exists(csv_file)) {
    simulation_results <- fread(file = glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz")))
  } else {
    files <- list.files(
      path = here("Data", "08_Simu_T2Dmono", paste0("simu_", type_resampling)), 
      pattern = "_xxr.tsv", 
      recursive = TRUE, 
      full.names = TRUE
    )
    
    simulation_results <- rbindlist(
      mclapply(
        X = files,
        mc.preschedule = FALSE, 
        mc.cores = min(80, detectCores()), 
        FUN = function(x) {
          out <- fread(x)[, -c("qc11_done", "mist_estimate", "mist_statistics")][,
            filename := gsub("_xxr.tsv", "", basename(x))
          ]
        }
      ), 
      fill = TRUE
    )
    
    simulation_results[, 
      c("simu_i", "T2D_def", "n_samples_size", "gene_file") := 
        tstrsplit(filename, split = "_", type.convert = TRUE)
    ]
    
    fwrite(simulation_results, file = glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz")))
  }
  
  simulation_results_long <- melt(
    simulation_results,
    measure.vars = grep("p.value", names(simulation_results)),
    variable.name = "type_p_values", 
    value.name = "p_values"
  )
  
  simulation_results_summary <- simulation_results_long[, 
    .(
      mean_pvalue = mean(p_values, na.rm = TRUE),
      sd_pvalue = sd(p_values, na.rm = TRUE),
      n = sum(!is.na(p_values))
    ), 
    by = c("T2D_def", "n_samples_size", "gene_file", "type_p_values")
  ]
  
  simulation_results_summary[,
    type_p_values := c(
      "p.value.S.pi" = "S(&pi;)", "p.value.S.tau" = "S(&tau;)", "p.value.overall" = "Overall"
    )[type_p_values]
  ]
  
  simulation_results_summary[, 
    col_strip := paste0(
      "FG < ", number(as.numeric(gsub("T2D(.*)", "\\1", T2D_def)), accuracy = 0.1), "<br>",
      "<i style='font-size:6pt;'>(Control)</i>"
    ),
    by = c("T2D_def", "gene_file")
  ]

  simulation_results_summary[, gene_file := paste0("<i>", gene_file, "</i>")]
  
  simulation_results_summary[, 
    linetype := any(mean_pvalue < 0.05 & n_samples_size == max(n_samples_size)), 
    by = c("T2D_def", "gene_file", "type_p_values")
  ]
  
  p <- ggplot(
    data = simulation_results_summary[order(as.numeric(gsub("T2D(.*)", "\\1", T2D_def)))], 
    mapping = aes(x = n_samples_size, y = mean_pvalue)
  ) +
    geom_ribbon(
      mapping = aes(ymin = mean_pvalue - sd_pvalue, ymax = mean_pvalue + sd_pvalue, fill = type_p_values), 
      alpha = 0.05
    ) +
    geom_line(mapping = aes(colour = type_p_values)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "firebrick2", size = 0.25) +
    scale_x_continuous(labels = comma, guide = guide_axis(angle = 25)) +
    scale_y_continuous(
      breaks = function(x) c(setdiff(scales::breaks_extended()(x), 0), 0.05),
      labels = function(x) ifelse(x == 0.05, "<b style='color:firebrick2;'>0.05</b>", x)
    ) +
    scale_colour_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
    scale_fill_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
    facet_grid(
      rows = vars(fct_reorder2(gene_file, n_samples_size, mean_pvalue)), 
      cols = vars(fct_inorder(col_strip))
    ) +
    labs(
      x = "Sample Size", 
      y = "P-value", 
      colour = "P-value", 
      fill = "P-value",
      title = "Rare Variant Simulation Analysis Using MiST",
      caption = glue(
        "In average for {number(simulation_results_summary[, mean(n)], accuracy = 0.1)}",
        " resamplings (up to {number(simulation_results_summary[, max(n)], accuracy = 1)}).<br>",
        "Type 2 diabetes cases: adult with glucose lowering drug or fasting glucose (FG) &ge; 7 mmol/l.<br>",
        "The significance threshold is denoted with a horizontal red dashed line <b style='color:firebrick2;'>---</b>."
      )
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme(
      plot.caption = element_markdown(face = "italic", size = rel(0.5)),
      legend.position = "top"
    )
  
  ggsave(
    filename = glue(file.path(output_directory, "MiST_simulation_{type_resampling}.png")), 
    plot = p, 
    width = 16, 
    height = 24, 
    units = "cm", 
    dpi = 120
  )
  
  p_all <- map(c(FALSE, TRUE), function(x) {
    ggplot(
      data = simulation_results_summary[
        type_p_values == "Overall"
      ][linetype == x][order(as.numeric(gsub("T2D(.*)", "\\1", T2D_def)))], 
      mapping = aes(x = n_samples_size, y = mean_pvalue)
    ) +
      geom_ribbon(
        mapping = aes(
          ymin = mean_pvalue - sd_pvalue, 
          ymax = mean_pvalue + sd_pvalue, 
          fill = fct_reorder2(gene_file, n_samples_size, mean_pvalue)
        ), 
        alpha = 0.05
      ) +
      geom_line(aes(
        colour = fct_reorder2(gene_file, n_samples_size, mean_pvalue),
        linetype = fct_reorder2(gene_file, n_samples_size, mean_pvalue)
      )) +
      geom_hline(yintercept = 0.05, linetype = "dashed", colour = "firebrick2", size = 0.25) +
      scale_x_continuous(labels = comma, guide = guide_axis(angle = 25)) +
      scale_y_continuous(
        breaks = function(x) c(scales::breaks_extended()(x), 0.05),
        labels = function(x) ifelse(x == 0.05, "<b style='color:firebrick2;'>0.05</b>", x)
      ) +
      scale_colour_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
      scale_fill_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
      facet_grid(cols = vars(fct_inorder(col_strip))) +
      labs(
        x = "Sample Size", 
        y = "Overall P-value", 
        colour = NULL, 
        fill = NULL,
        linetype = NULL
      ) +
      coord_cartesian(ylim = c(0, 0.75))
  })
  
  ggsave(
    filename = glue(file.path(output_directory, "MiST_simulation_{type_resampling}_overall.png")), 
    plot = wrap_plots(p_all, nrow = 2, ncol = 1) +
      plot_annotation(
        title = "Rare Variant Simulation Analysis Using MiST",
        subtitle = paste(
          "In <b>A</b> are displayed the genes never showing any significant rare variants cluster association.",
          "In <b>B</b> are displayed the genes showing a significant rare variants cluster association.",
          sep = "<br>"
        ),
        caption = glue(
          "In average for {number(simulation_results_summary[, mean(n)], accuracy = 0.1)}",
          " resamplings (up to {number(simulation_results_summary[, max(n)], accuracy = 1)}).<br>",
          "Type 2 diabetes cases: adult with glucose lowering drug or fasting glucose (FG) &ge; 7 mmol/l.<br>",
          "The significance threshold is denoted with a horizontal red dashed line <b style='color:firebrick2;'>---</b>."
        ), 
        tag_levels = "A", 
        theme = theme(plot.caption = element_markdown(face = "italic", size = rel(0.5)))
      ), 
    width = 16, 
    height = 16, 
    units = "cm", 
    dpi = 120
  )
})



### focus on a gene ================================================================================
walk(c("sansRemise", "avecRemise"), function(type_resampling) {
  
  my_gene <- "GCK"
  
  csv_file <- glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz"))
  if (file.exists(csv_file)) {
    simulation_results <- fread(file = glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz")))
  } else {
    files <- list.files(
      path = here("Data", "08_Simu_T2Dmono", paste0("simu_", type_resampling)), 
      pattern = "_xxr.tsv", 
      recursive = TRUE, 
      full.names = TRUE
    )
    
    simulation_results <- rbindlist(
      mclapply(
        X = files,
        mc.preschedule = FALSE, 
        mc.cores = min(80, detectCores()), 
        FUN = function(x) {
          out <- fread(x)[, -c("qc11_done", "mist_estimate", "mist_statistics")][,
            filename := gsub("_xxr.tsv", "", basename(x))
          ]
        }
      ), 
      fill = TRUE
    )
    
    simulation_results[, 
      c("simu_i", "T2D_def", "n_samples_size", "gene_file") := 
        tstrsplit(filename, split = "_", type.convert = TRUE)
    ]
    
    fwrite(simulation_results, file = glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz")))
  }
  
  simulation_results <- simulation_results[simulation_results$gene_file%in%my_gene, ]
  
  simulation_results_long <- melt(
    simulation_results,
    measure.vars = grep("p.value", names(simulation_results)),
    variable.name = "type_p_values", 
    value.name = "p_values"
  )
  
  simulation_results_summary <- simulation_results_long[, 
    .(
      mean_pvalue = mean(p_values, na.rm = TRUE),
      sd_pvalue = sd(p_values, na.rm = TRUE),
      n = sum(!is.na(p_values))
    ), 
    by = c("T2D_def", "n_samples_size", "gene_file", "type_p_values")
  ]
  
  simulation_results_summary[,
    type_p_values := c(
      "p.value.S.pi" = "S(&pi;)", "p.value.S.tau" = "S(&tau;)", "p.value.overall" = "Overall"
    )[type_p_values]
  ]
  
  simulation_results_summary[, 
    col_strip := paste0(
      "FG < ", number(as.numeric(gsub("T2D(.*)", "\\1", T2D_def)), accuracy = 0.1), "<br>",
      "<i style='font-size:6pt;'>(Control)</i>"
    ),
    by = c("T2D_def", "gene_file")
  ]

  simulation_results_summary[, gene_file := paste0("<i>", gene_file, "</i>")]
  
  simulation_results_summary[, 
    linetype := any(mean_pvalue < 0.05 & n_samples_size == max(n_samples_size)), 
    by = c("T2D_def", "gene_file", "type_p_values")
  ]
  
  p <- ggplot(
    data = simulation_results_summary[order(as.numeric(gsub("T2D(.*)", "\\1", T2D_def)))], 
    mapping = aes(x = n_samples_size, y = mean_pvalue)
  ) +
    geom_ribbon(
      mapping = aes(ymin = mean_pvalue - sd_pvalue, ymax = mean_pvalue + sd_pvalue, fill = type_p_values), 
      alpha = 0.05
    ) +
    geom_line(mapping = aes(colour = type_p_values)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "firebrick2", size = 0.25) +
    scale_x_continuous(labels = comma, guide = guide_axis(angle = 25)) +
    scale_y_continuous(
      breaks = function(x) c(setdiff(scales::breaks_extended()(x), 0), 0.05),
      labels = function(x) ifelse(x == 0.05, "<b style='color:firebrick2;'>0.05</b>", x)
    ) +
    scale_colour_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
    scale_fill_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
    facet_grid(
      cols = vars(fct_inorder(col_strip))
    ) +
    labs(
      x = "Sample Size", 
      y = "P-value", 
      colour = "P-value", 
      fill = "P-value",
      title = glue("Rare Variant Simulation Analysis Using MiST: {my_gene}"),
      caption = glue(
        "In average for {number(simulation_results_summary[, mean(n)], accuracy = 0.1)}",
        " resamplings (up to {number(simulation_results_summary[, max(n)], accuracy = 1)}).<br>",
        "Type 2 diabetes cases: adult with glucose lowering drug or fasting glucose (FG) &ge; 7 mmol/l.<br>",
        "The significance threshold is denoted with a horizontal red dashed line <b style='color:firebrick2;'>---</b>."
      )
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme(
      plot.caption = element_markdown(face = "italic", size = rel(0.5)),
      legend.position = "top"
    )
  ggsave(
    filename = glue(file.path(output_directory, "MiST_simulation_{type_resampling}_{my_gene}.png")), 
    plot = p, 
    width = 16, 
    height = 24/2, 
    units = "cm", 
    dpi = 120
  )
})


### Graphic effect =================================================================================

walk(c("sansRemise", "avecRemise"), function(type_resampling) { 
  
  simulation_results <- fread(file = glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz")))
  
  ## graph effect
  dtaeffect <- simulation_results %>% 
    dplyr::select(
      c("gene_file", "T2D_def", "n_samples_size", "Pi_hat")
    ) %>% 
    group_by(gene_file, T2D_def, n_samples_size) %>% 
    dplyr::summarize(
      Pi_hat_mean = mean(x = Pi_hat, na.rm = TRUE), 
      sd = sd(x = Pi_hat, na.rm = TRUE)
    ) %>% 
    ungroup() %>% 
    mutate(n_samples_size = as.numeric(n_samples_size)) %>% 
    arrange(n_samples_size) %>% 
    mutate(sign = factor(x = Pi_hat_mean>0, levels = c("TRUE", "FALSE"), labels = c("Positive", "Negative")))
  
  effect <- ggplot(data = dtaeffect, mapping = aes(x = n_samples_size, y = Pi_hat_mean)) + 
    geom_point(aes(color = sign)) +
    geom_line(data = dtaeffect, mapping = aes(x = n_samples_size, y = Pi_hat_mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "pink") +
    geom_ribbon(data = dtaeffect, mapping = aes(ymin = Pi_hat_mean - sd, ymax = Pi_hat_mean + sd), alpha = 0.1) +
    scale_x_continuous(labels = comma, guide = guide_axis(angle = 25)) +
    scale_colour_viridis_d(option = "plasma", begin = 0.10, end = 0.90)+ 
    theme(axis.text.x = element_text(angle = 90)) + 
    facet_grid(gene_file ~ T2D_def, scales = "free_y") +
    labs(
      x = "Sample Size", 
      y = expression(paste("Mean  ",  hat(pi), italic("  i.e. `Pi hat` estimates"))),
      caption =  glue(
        "Over 100 resamplings.<br>",
        "Type 2 diabetes cases: adult with glucose lowering drug or fasting glucose (FG) &ge; 7 mmol/l.<br>",
        "The significance threshold is denoted with a horizontal red dashed line <b style='color:pink;'>---</b>."
      )
    ) +
    theme(plot.caption = element_markdown(face = "italic", size = rel(0.5)))
  # effect
  
  ggsave(
    filename = glue(file.path(output_directory, "MiST_effects_{type_resampling}.png")), 
    plot = effect
  )

})



### Analysis only overall et 5.6 / 6.1 et 7   ======================================================
walk(c("sansRemise", "avecRemise"), function(type_resampling) {
  csv_file <- glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz"))
  if (file.exists(csv_file)) {
    simulation_results <- fread(file = glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz")))
  } else {
    files <- list.files(
      path = here("Data", "08_Simu_T2Dmono", paste0("simu_", type_resampling)), 
      pattern = "_xxr.tsv", 
      recursive = TRUE, 
      full.names = TRUE
    )
    
    simulation_results <- rbindlist(
      mclapply(
        X = files,
        mc.preschedule = FALSE, 
        mc.cores = min(80, detectCores()), 
        FUN = function(x) {
          out <- fread(x)[, -c("qc11_done", "mist_estimate", "mist_statistics")][,
            filename := gsub("_xxr.tsv", "", basename(x))
          ]
        }
      ), 
      fill = TRUE
    )
    
    simulation_results[, 
      c("simu_i", "T2D_def", "n_samples_size", "gene_file") := 
        tstrsplit(filename, split = "_", type.convert = TRUE)
    ]
    
    fwrite(simulation_results, file = glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz")))
  }
  
  simulation_results <- simulation_results[simulation_results$T2D_def %in% c("T2D5.6", "T2D6.1", "T2D7"), ] 

  simulation_results_long <- melt(
    simulation_results,
    measure.vars = grep("p.value", names(simulation_results)),
    variable.name = "type_p_values", 
    value.name = "p_values"
  )
  
  simulation_results_long <- simulation_results_long[simulation_results_long$type_p_values %in% c("p.value.overall"), ] 
  
  simulation_results_summary <- simulation_results_long[, 
    .(
      mean_pvalue = mean(p_values, na.rm = TRUE),
      sd_pvalue = sd(p_values, na.rm = TRUE),
      n = sum(!is.na(p_values))
    ), 
    by = c("T2D_def", "n_samples_size", "gene_file", "type_p_values")
  ]
  
  simulation_results_summary[,
    type_p_values := c(
      "p.value.S.pi" = "S(&pi;)", "p.value.S.tau" = "S(&tau;)", "p.value.overall" = "Overall"
    )[type_p_values]
  ]
  
  simulation_results_summary[, 
    col_strip := paste0(
      "FG < ", number(as.numeric(gsub("T2D(.*)", "\\1", T2D_def)), accuracy = 0.1), "<br>",
      "<i style='font-size:6pt;'>(Control)</i>"
    ),
    by = c("T2D_def", "gene_file")
  ]

  simulation_results_summary[, gene_file := paste0("<i>", gene_file, "</i>")]
  
  simulation_results_summary[, 
    linetype := any(mean_pvalue < 0.05 & n_samples_size == max(n_samples_size)), 
    by = c("T2D_def", "gene_file", "type_p_values")
  ]
  
  p <- ggplot(
    data = simulation_results_summary[order(as.numeric(gsub("T2D(.*)", "\\1", T2D_def)))], 
    mapping = aes(x = n_samples_size, y = mean_pvalue)
  ) +
    geom_ribbon(
      mapping = aes(ymin = mean_pvalue - sd_pvalue, ymax = mean_pvalue + sd_pvalue, fill = type_p_values), 
      alpha = 0.05
    ) +
    geom_line(mapping = aes(colour = type_p_values)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", colour = "firebrick2", size = 0.25) +
    scale_x_continuous(labels = comma, guide = guide_axis(angle = 25)) +
    scale_y_continuous(
      breaks = function(x) c(setdiff(scales::breaks_extended()(x), 0), 0.05),
      labels = function(x) ifelse(x == 0.05, "<b style='color:firebrick2;'>0.05</b>", x)
    ) +
    scale_colour_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
    scale_fill_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
    facet_grid(
      rows = vars(fct_reorder2(gene_file, n_samples_size, mean_pvalue)), 
      cols = vars(fct_inorder(col_strip))
    ) +
    labs(
      x = "Sample Size", 
      y = "P-value", 
      colour = "P-value", 
      fill = "P-value",
      title = "Rare Variant Simulation Analysis Using MiST",
      caption = glue(
        "In average for {number(simulation_results_summary[, mean(n)], accuracy = 0.1)}",
        " resamplings (up to {number(simulation_results_summary[, max(n)], accuracy = 1)}).<br>",
        "Type 2 diabetes cases: adult with glucose lowering drug or fasting glucose (FG) &ge; 7 mmol/l.<br>",
        "The significance threshold is denoted with a horizontal red dashed line <b style='color:firebrick2;'>---</b>."
      )
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    theme(
      plot.caption = element_markdown(face = "italic", size = rel(0.5)),
      legend.position = "top"
    )
  
  ggsave(
    filename = glue(file.path(output_directory, "MiST_simulation_{type_resampling}_3defs.png")), 
    plot = p, 
    width = 16, 
    height = 24, 
    units = "cm", 
    dpi = 120
  )
  
  p_all <- map(c(FALSE, TRUE), function(x) {
    ggplot(
      data = simulation_results_summary[
        type_p_values == "Overall"
      ][linetype == x][order(as.numeric(gsub("T2D(.*)", "\\1", T2D_def)))], 
      mapping = aes(x = n_samples_size, y = mean_pvalue)
    ) +
      geom_ribbon(
        mapping = aes(
          ymin = mean_pvalue - sd_pvalue, 
          ymax = mean_pvalue + sd_pvalue, 
          fill = fct_reorder2(gene_file, n_samples_size, mean_pvalue)
        ), 
        alpha = 0.05
      ) +
      geom_line(aes(
        colour = fct_reorder2(gene_file, n_samples_size, mean_pvalue),
        linetype = fct_reorder2(gene_file, n_samples_size, mean_pvalue)
      )) +
      geom_hline(yintercept = 0.05, linetype = "dashed", colour = "firebrick2", size = 0.25) +
      scale_x_continuous(labels = comma, guide = guide_axis(angle = 25)) +
      scale_y_continuous(
        breaks = function(x) c(scales::breaks_extended()(x), 0.05),
        labels = function(x) ifelse(x == 0.05, "<b style='color:firebrick2;'>0.05</b>", x)
      ) +
      scale_colour_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
      scale_fill_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
      facet_grid(cols = vars(fct_inorder(col_strip))) +
      labs(
        x = "Sample Size", 
        y = "Overall P-value", 
        colour = NULL, 
        fill = NULL,
        linetype = NULL
      ) +
      coord_cartesian(ylim = c(0, 0.75))
  })
  
  ggsave(
    filename = glue(file.path(output_directory, "MiST_simulation_{type_resampling}_overall_3defs.png")), 
    plot = wrap_plots(p_all, nrow = 2, ncol = 1) +
      plot_annotation(
        title = "Rare Variant Simulation Analysis Using MiST",
        subtitle = paste(
          "In <b>A</b> are displayed the genes never showing any significant rare variants cluster association.",
          "In <b>B</b> are displayed the genes showing a significant rare variants cluster association.",
          sep = "<br>"
        ),
        caption = glue(
          "In average for {number(simulation_results_summary[, mean(n)], accuracy = 0.1)}",
          " resamplings (up to {number(simulation_results_summary[, max(n)], accuracy = 1)}).<br>",
          "Type 2 diabetes cases: adult with glucose lowering drug or fasting glucose (FG) &ge; 7 mmol/l.<br>",
          "The significance threshold is denoted with a horizontal red dashed line <b style='color:firebrick2;'>---</b>."
        ), 
        tag_levels = "A", 
        theme = theme(plot.caption = element_markdown(face = "italic", size = rel(0.5)))
      ), 
    width = 16, 
    height = 16, 
    units = "cm", 
    dpi = 120
  )
})


### FINAL FIG ======================================================================================
walk(c("sansRemise", "avecRemise"), function(type_resampling) {
  csv_file <- glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz"))
  if (file.exists(csv_file)) {
    simulation_results <- fread(file = glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz")))
  } else {
    files <- list.files(
      path = here("Data", "08_Simu_T2Dmono", paste0("simu_", type_resampling)), 
      pattern = "_xxr.tsv", 
      recursive = TRUE, 
      full.names = TRUE
    )
    
    simulation_results <- rbindlist(
      mclapply(
        X = files,
        mc.preschedule = FALSE, 
        mc.cores = min(80, detectCores()), 
        FUN = function(x) {
          out <- fread(x)[, -c("qc11_done", "mist_estimate", "mist_statistics")][,
            filename := gsub("_xxr.tsv", "", basename(x))
          ]
        }
      ), 
      fill = TRUE
    )
    
    simulation_results[, 
      c("simu_i", "T2D_def", "n_samples_size", "gene_file") := 
        tstrsplit(filename, split = "_", type.convert = TRUE)
    ]
    
    fwrite(simulation_results, file = glue(file.path(output_directory, "MiST_simulation_{type_resampling}.csv.gz")))
  }
  
  simulation_results <- simulation_results[simulation_results$T2D_def %in% c("T2D6.1", "T2D7"), ] 

  simulation_results_long <- melt(
    simulation_results,
    measure.vars = grep("p.value", names(simulation_results)),
    variable.name = "type_p_values", 
    value.name = "p_values"
  )
  
  simulation_results_summary <- simulation_results_long[, 
    .(
      mean_pvalue = mean(p_values, na.rm = TRUE),
      sd_pvalue = sd(p_values, na.rm = TRUE),
      n = sum(!is.na(p_values))
    ), 
    by = c("T2D_def", "n_samples_size", "gene_file", "type_p_values")
  ]
  
  simulation_results_summary[,
    type_p_values := c(
      "p.value.S.pi" = "S(&pi;)", "p.value.S.tau" = "S(&tau;)", "p.value.overall" = "Overall"
    )[type_p_values]
  ]
  
  simulation_results_summary[, 
    col_strip := paste0(
      "FG < ", number(as.numeric(gsub("T2D(.*)", "\\1", T2D_def)), accuracy = 0.1), "<br>",
      "<i style='font-size:6pt;'>(Control)</i>"
    ),
    by = c("T2D_def", "gene_file")
  ]

  simulation_results_summary[, gene_file := paste0("<i>", gene_file, "</i>")]
  
  simulation_results_summary[, 
    linetype := any(mean_pvalue < 0.05 & n_samples_size == max(n_samples_size)), 
    by = c("T2D_def", "gene_file", "type_p_values")
  ]

  
  ## filrer only pvaloverall
  simulation_results_summary <- simulation_results_summary[simulation_results_summary$type_p_values %in% c("Overall"), ] 
  
  
  simulation_results_summary$title <- paste0(
    "Controls (FG < ", gsub("T2D", "", simulation_results_summary$T2D_def), " mmol/l) / Cases with T2D"
  )
  simulation_results_summary$gene_file <- as.factor(simulation_results_summary$gene_file)
  
  tmp <- simulation_results_summary %>% 
    group_by(T2D_def) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(
      plot = map(.x = data, function(df) {
  
      df <- as.data.frame(df)
      
      p_all <- ggplot(
          data = df,
          mapping = aes(x = n_samples_size, y = mean_pvalue, group = gene_file)
        ) +
          geom_ribbon(
            mapping = aes(
              ymin = mean_pvalue - sd_pvalue,
              ymax = mean_pvalue + sd_pvalue,
              fill = gene_file # fct_reorder2(gene_file, n_samples_size, mean_pvalue)
            ),
            alpha = 0.05
          ) +
          geom_line(aes(
            colour = gene_file, #fct_reorder2(gene_file, n_samples_size, mean_pvalue),
            linetype = gene_file #fct_reorder2(gene_file, n_samples_size, mean_pvalue)
          )) +
          geom_hline(yintercept = 0.05, linetype = "dashed", colour = "firebrick2", size = 0.25) +
          scale_x_continuous(labels = comma) +
          scale_y_continuous(
            breaks = function(x) c(scales::breaks_extended()(x), 0.05),
            labels = function(x) ifelse(x == 0.05, "<b style='color:firebrick2;'>0.05</b>", x), 
            expand = expansion(c(0, NA))
          ) +
          scale_colour_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
          scale_fill_viridis_d(option = "plasma", begin = 0.10, end = 0.90) +
          labs(
            x = "Sample Size",
            y = "Overall P-value",
            colour = NULL,
            fill = NULL,
            linetype = NULL
          ) +
          coord_cartesian(ylim = c(0, 0.75)) 
          

      return(p_all)
      }) , 
      save_it = map2(.x = plot, .y = T2D_def, function(plot, T2D_def){
        ggsave(
          filename = glue(file.path(output_directory, "MiST_simulation_{type_resampling}_overall_{T2D_def}.png")),
          plot = plot,
          width = 16,
          height = 12,
          units = "cm",
          dpi = 600
        )

      })
    )


})


