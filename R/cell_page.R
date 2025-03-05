

render_cell_page <- function(output, input, session, cur_cell_type, cell_data, markers, markers_enrich_WormCat.1, expr_data, grid_id) {
  
  all_genes <- lapply(strsplit(rownames(expr_data), " "), function(x) {x[3]})
  cel_genes <- unlist(all_genes[grepl("C. elegans", rownames(expr_data))])
  cbr_genes <- unlist(all_genes[grepl("C. briggsae", rownames(expr_data))])
  
  cell_plot_df = data.frame(gene = unique(unlist(all_genes)))
  rownames(cell_plot_df) <- cell_plot_df$gene
  
  cell_plot_df$ortho <- NA
  cell_plot_df[cell_plot_df$gene %in% cel_genes,]$ortho <- "cel"
  cell_plot_df[cell_plot_df$gene %in% cbr_genes,]$ortho <- "cbr"
  cell_plot_df[cell_plot_df$gene %in% cel_genes & cell_plot_df$gene %in% cbr_genes,]$ortho <- "common"

  cell_plot_df$cel_exp <- NA
  cell_plot_df$cbr_exp <- NA
  
  cell_plot_df[cel_genes,]$cel_exp <- expr_data[paste0("C. elegans ", cel_genes), cur_cell_type]
  cell_plot_df[cbr_genes,]$cbr_exp <- expr_data[paste0("C. briggsae ", cbr_genes), cur_cell_type]
  
  cell_plot_df <- cell_plot_df[cell_plot_df$ortho == "common",]
  
  ### Data processing ###
  WormCat.df <- markers_enrich_WormCat.1[markers_enrich_WormCat.1$cell_type == cur_cell_type,]
  Categories <- levels(factor(markers[["C.elegans"]]$WormCat.1))
  Categories <- Categories[Categories != "Pseudogene"]
  
  WormCat.df <- WormCat.df %>% filter(Category %in% Categories)
  WormCat.df[is.na(WormCat.df$RGS),"RGS"] <- 0
  
  cell_plot_df$cel_marker <- "Not a marker"
  cell_plot_df[cell_plot_df$gene %in% markers[["C.elegans"]][markers[["C.elegans"]]$cell_type == cur_cell_type,]$gene,]$cel_marker <- "C. elegans marker"
  cell_plot_df$cbr_marker <- "Not a marker"
  cell_plot_df[cell_plot_df$gene %in% markers[["C.briggsae"]][markers[["C.briggsae"]]$cell_type == cur_cell_type,]$gene,]$cbr_marker <- "C. briggsae marker"
  
  cell_plot_df$comb_marker <- cell_plot_df$cel_marker
  cell_plot_df$comb_marker <- ifelse(cell_plot_df$cel_marker == "C. elegans marker" & cell_plot_df$cbr_marker == "C. briggsae marker",
                                     "Joint marker",
                                     ifelse(cell_plot_df$cbr_marker == "C. briggsae marker", "C. briggsae marker", cell_plot_df$comb_marker))
  
  exp_max <- max(cell_plot_df$cel_exp, cell_plot_df$cbr_exp)
  
  TPMPlot <- ggplot() +
    geom_point(data = cell_plot_df[cell_plot_df$comb_marker == "Not a marker",],
               aes(x = cel_exp + 1,
                   y = cbr_exp + 1,
                   color = comb_marker),
               alpha = 0.5, size = 1, stroke = 0.3) +
    geom_point(data = cell_plot_df[cell_plot_df$comb_marker != "Not a marker",],
               aes(x = cel_exp + 1,
                   y = cbr_exp + 1,
                   color = comb_marker),
               alpha = 0.5, size = 1, stroke = 0.3) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    scale_x_continuous(name = "C. elegans TPM", limits = c(1, exp_max + 100), trans = "log2") +
    scale_y_continuous(name = "C. briggsae TPM", limits = c(1, exp_max + 100), trans = "log2") +
    scale_color_manual(name = "Markers", breaks = c("C. elegans marker", "C. briggsae marker", "Joint marker", "Not a marker"),
                       labels = c("C. elegans marker", "C. briggsae marker", "Joint marker", "Neither"), values = c("#009E73", "#56B4E9", "#CC79A7", "black")) +
    theme(legend.title= element_blank(),
          legend.position = "bottom",
          rect = element_rect(fill = "transparent"),
          axis.title.x = element_text(color = "#009E73", size = 14),
          axis.title.y = element_text(color = "#56B4E9", size = 14),
          axis.text = element_text(size = 12),
          axis.line = element_line(color="grey80", size=1),
          panel.grid.major = element_line(color="grey80", size=0.25),
          panel.grid.minor = element_line(color="grey80", size=0.05),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"))
  
  WormCatPlot.Fold <- ggplot(WormCat.df, aes(x = Category, y = Fold)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(name = "Worm Category (Tier 1)",
                     breaks = Categories,
                     limits = Categories) +
    scale_y_continuous(name = "Fold Enrich.") +
    theme(legend.title= element_blank(),
          legend.position = "none",
          rect = element_rect(fill = "transparent"),
          axis.line = element_line(color="grey80", size=1),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line("grey80", size = 0.2, linetype = "solid"),
          panel.grid.minor.y = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"))
  
  WormCatPlot.Freq <- ggplot(WormCat.df, aes(x = Category, y = RGS)) +
    geom_bar(stat = "identity") +
    scale_x_discrete(name = "Worm Category (Tier 1)",
                     breaks = Categories,
                     limits = Categories) +
    scale_y_continuous(name = "Count") +
    theme(legend.title= element_blank(),
          legend.position = "top",
          rect = element_rect(fill = "transparent"),
          axis.line = element_line(color="grey80", size=1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line("grey80", size = 0.2, linetype = "solid"),
          panel.grid.minor.y = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  temp_shared_genes <- head(markers[["C.elegans"]] %>%
                              filter(cell_type == cur_cell_type & shared_marker & gene_orthology == "1:1") %>%
                              select("gene", "cel_tpm", "cbr_tpm"), n = 15)
  temp_shared_genes$set <- "shared"
  
  temp_cel_genes <- head(markers[["C.elegans"]] %>%
                           filter(cell_type == cur_cell_type & ! shared_marker & gene_orthology == "1:1") %>%
                           select("gene", "cel_tpm", "cbr_tpm"), n = 5)
  temp_cel_genes$set <- "cel"
  
  temp_cbr_genes <- head(markers[["C.briggsae"]] %>%
                           filter(cell_type == cur_cell_type & ! shared_marker & gene_orthology == "1:1") %>%
                           select("gene", "cel_tpm", "cbr_tpm"), n = 5)
  temp_cbr_genes$set <- "cbr"
  
  temp_gene_set <- rbind(temp_shared_genes, temp_cel_genes, temp_cbr_genes) %>%
    tidyr::pivot_longer(cols = c("cel_tpm", "cbr_tpm"), 
                 names_to = "species", 
                 values_to = "value")
  temp_gene_set$gene <- factor(temp_gene_set$gene, levels = c(temp_shared_genes$gene, temp_cel_genes$gene, temp_cbr_genes$gene))
  
  if(sum(is.na(temp_gene_set$value)) > 0) {
    temp_gene_set[is.na(temp_gene_set$value),]$value <- 0
  }
  
  marker_plot <- ggplot(temp_gene_set, aes(x = gene, y = species, fill = as.numeric(value))) +
    geom_tile() +
    geom_rect(aes(xmin = 0.5, xmax = 0.45 + (1 * 15), ymin = 0.5, ymax = 2.5),
              fill = "transparent", color = "black", size = 1) +
    geom_rect(aes(xmin = 15.55, xmax = 15.45 + (1 * 5), ymin = 0.5, ymax = 2.5),
              fill = "transparent", color = "#009E73", size = 1) +
    geom_rect(aes(xmin = 20.55, xmax = 20.5 + (1 * 5), ymin = 0.5, ymax = 2.5),
              fill = "transparent", color = "#56B4E9", size = 1) +
    scale_fill_viridis_c(name = "TPM", 
                         option="magma",
                         trans = "log2",
                         limits = c(10, max(temp_gene_set$value, na.rm = TRUE) + 1),
                         guide = guide_colorbar(),
                         oob = scales::squish) +
    scale_x_discrete() + 
    scale_y_discrete(labels = c("C. briggsae", "C. elegans")) + 
    # guides(fill = guide_legend(title = "TPM")) +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          #axis.line = element_line(color="grey80", size=1),
          axis.line = element_blank(),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.key.size = unit(0.5, "cm"),  # Adjust the size of the legend key
          legend.key.height = unit(0.3, "cm"))
  
  cellCountPlot <- ggplot() +
    geom_density(data = cell_data, aes(x = cel_cell_count), color = "#009E73", fill = "#009E73", alpha = 0.35) +
    geom_density(data = cell_data, aes(x = cbr_cell_count), color = "#56B4E9", fill = "#56B4E9", alpha = 0.35) +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cel_cell_count"]), color = "#009E73") +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cbr_cell_count"]), color = "#56B4E9") +
    scale_x_continuous(name = "Cell count", trans = "log2") +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  geneCountPlot <- ggplot() +
    geom_density(data = cell_data, aes(x = genes_detected_bootstrap_cel), color = "#009E73", fill = "#009E73", alpha = 0.35) +
    geom_density(data = cell_data, aes(x = genes_detected_bootstrap_cbr), color = "#56B4E9", fill = "#56B4E9", alpha = 0.35) +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "genes_detected_bootstrap_cel"]), color = "#009E73") +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "genes_detected_bootstrap_cbr"]), color = "#56B4E9") +
    scale_x_continuous(name = "Genes detected") +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  umiPlot <- ggplot() +
    geom_density(data = cell_data, aes(x = cel_median_umi), color = "#009E73", fill = "#009E73", alpha = 0.35) +
    geom_density(data = cell_data, aes(x = cbr_median_umi), color = "#56B4E9", fill = "#56B4E9", alpha = 0.35) +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cel_median_umi"]), color = "#009E73") +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cbr_median_umi"]), color = "#56B4E9") +
    scale_x_continuous(name = "Mean num. of UMI", trans = "log2") +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  giniPlot <- ggplot() +
    geom_density(data = cell_data, aes(x = cel_gini_median), color = "#009E73", fill = "#009E73", alpha = 0.35) +
    geom_density(data = cell_data, aes(x = cbr_gini_median), color = "#56B4E9", fill = "#56B4E9", alpha = 0.35) +
    geom_rect(data = cell_data[cur_cell_type,], aes(xmin = cel_gini_lower, xmax = cel_gini_upper, ymin = -Inf, ymax = Inf),
              color = "#00805e", fill = "#009E73", alpha = 0.25) +
    geom_rect(data = cell_data[cur_cell_type,], aes(xmin = cbr_gini_lower, xmax = cbr_gini_upper, ymin = -Inf, ymax = Inf),
              color = "#1d9ae2", fill = "#56B4E9", alpha = 0.25) +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cel_gini_median"]), color = "#009E73") +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cbr_gini_median"]), color = "#56B4E9") +
    scale_x_continuous(name = "Gini coef.") +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  jsdPlot <- ggplot() +
    geom_density(data = cell_data, aes(x = jsd_median)) +
    geom_rect(data = cell_data[cur_cell_type,], aes(xmin = jsd_lower, xmax = jsd_upper, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "jsd_median"]), color = "red") +
    scale_x_continuous(name = "Jensen-Shannon distance") +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  corPlot <- ggplot() +
    geom_density(data = cell_data, aes(x = cor_median)) +
    geom_rect(data = cell_data[cur_cell_type,], aes(xmin = cor_lower, xmax = cor_upper, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cor_median"]), color = "red") +
    scale_x_continuous(name = "Pearson correlation") +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  cosPlot <- ggplot() +
    geom_density(data = cell_data, aes(x = cos_median)) +
    geom_rect(data = cell_data[cur_cell_type,], aes(xmin = cos_lower, xmax = cos_upper, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cos_median"]), color = "red") +
    scale_x_continuous(name = "Cosine distance") +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  countMarkersPlot <- ggplot() +
    geom_density(data = cell_data, aes(x = cel_markers), color = "#009E73", fill = "#009E73", alpha = 0.35) +
    geom_density(data = cell_data, aes(x = cbr_markers), color = "#56B4E9", fill = "#56B4E9", alpha = 0.35) +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cel_markers"]), color = "#009E73") +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cbr_markers"]), color = "#56B4E9") +
    scale_x_continuous(name = "Marker count") +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  privMarkersPlot <- ggplot() +
    geom_density(data = cell_data, aes(x = cel_markers_non_one_to_one/cel_markers), color = "#009E73", fill = "#009E73", alpha = 0.35) +
    geom_density(data = cell_data, aes(x = cbr_markers_non_one_to_one/cbr_markers), color = "#56B4E9", fill = "#56B4E9", alpha = 0.35) +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cel_markers_non_one_to_one"]/
                     cell_data[cur_cell_type, "cel_markers"]), color = "#009E73") +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "cbr_markers_non_one_to_one"]/
                     cell_data[cur_cell_type, "cbr_markers"]), color = "#56B4E9") +
    scale_x_continuous(name = "Ratio of non-1:1 markers") +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  DEGPlot <- ggplot() +
    geom_density(data = cell_data, aes(x = deg)) +
    geom_vline(aes(xintercept = cell_data[cur_cell_type, "deg"]), color = "red") +
    scale_x_continuous(name = "Btwn. species DEG") +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  metric_plots <- plot_grid(cellCountPlot,
                            geneCountPlot,
                            countMarkersPlot,
                            privMarkersPlot,
                            giniPlot,
                            umiPlot,
                            jsdPlot,
                            corPlot,
                            cosPlot,
                            DEGPlot,
                            ncol = 5)
  
  wormcat_plots <- plot_grid(WormCatPlot.Fold, WormCatPlot.Freq,
                             rel_heights = c(0.275, 1), ncol = 1, align = "v", axis = "tblr")
  
  high_level_plots <- plot_grid(TPMPlot, wormcat_plots, rel_widths = c(0.8, 1))
  
  plot_row <- plot_grid(high_level_plots, marker_plot, metric_plots, rel_heights = c(1, 0.35, 1), ncol = 1)
  
  # now add the title
  if(cell_data[cur_cell_type,]$cell_class == "progenitor") {
    title <- ggdraw() + 
      draw_label(
        paste0(cur_cell_type, ": ", cell_data[cur_cell_type,]$div_stage, " Cell stage : ", cell_data[cur_cell_type,]$lineage_group),
        fontface = 'bold', x = 0, hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 7))
  } else {
    title <- ggdraw() + 
      draw_label(
        paste0(cur_cell_type, ": ", cell_data[cur_cell_type,]$cell_class),
        fontface = 'bold', x = 0, hjust = 0) +
      theme(plot.margin = margin(0, 0, 0, 7))
  }
  
  output[[grid_id]] <- renderPlot({
    req(cur_cell_type)  # Ensures cur_cell_type is not NULL
    plot_grid(title, plot_row, ncol = 1, rel_heights = c(0.1, 1))
  }, height = 1000, width = 800, res = 100)
}