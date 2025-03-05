
render_cell_markers <- function(output, input, session, cur_cell_type, cur_species, markers, num_markers = 20, sort_by, grid_id) {
  
  # Define the dynamic height of the plot
  plot_height <- max(400, num_markers * 20) # 20px per marker, minimum 400px
  
  # values_to_use <- ifelse(cur_species == "C.elegans", "cel_tpm_log2fc", "cbr_tpm_log2fc")
  markers$gene <- ifelse(markers$gene == "", markers$briggsae_gene_short_name, markers$gene)
  
  # Pivoting the numeric columns into one with a categorical identifier
  markers <- markers %>%
    tidyr::pivot_longer(
      cols = c(cel_tpm_log2fc, cbr_tpm_log2fc),
      names_to = "species",
      values_to = "log2fc") %>%
    mutate(log2fc_cel_tpm = log2fc * cel_tpm,
           log2fc_cbr_tpm = log2fc * cbr_tpm)
  
  if (sort_by == "Specificity x Exp. (log2FC * log2TPM)") {
    if (cur_species == "C.elegans") {
      temp_markers <- markers %>% filter(species == "cel_tpm_log2fc") %>% arrange(desc(log2fc_cel_tpm))
    } else {
      temp_markers <- markers %>% filter(species == "cbr_tpm_log2fc") %>% arrange(desc(log2fc_cbr_tpm))
    }
  } else if (sort_by == "Specificity (log2FC)") {
    if(cur_species == "C.elegans") {
      temp_markers <- markers %>% filter(species == "cel_tpm_log2fc") %>% arrange(desc(log2fc))
    } else {
      temp_markers <- markers %>% filter(species == "cbr_tpm_log2fc") %>% arrange(desc(log2fc))
    }
  } else if (sort_by == "Specificity (p-value)") {
    if (cur_species == "C.elegans") {
      temp_markers <- markers %>% filter(species == "cel_tpm_log2fc") %>% arrange(p_val_adj.cel)
    } else {
      temp_markers <- markers %>% filter(species == "cbr_tpm_log2fc") %>% arrange(p_val_adj.cbr)
    }
  } else if (sort_by == "Expression") {
    if (cur_species == "C.elegans") {
      temp_markers <- markers %>% filter(species == "cel_tpm_log2fc") %>% arrange(desc(cel_tpm))
    } else {
      temp_markers <- markers %>% filter(species == "cbr_tpm_log2fc") %>% arrange(desc(cbr_tpm))
    }
  }
  
  temp_genes <- temp_markers %>% slice_head(n = num_markers) %>% pull(gene)
  
  markers <- markers %>%
    filter(gene %in% temp_genes) %>%
    arrange(paste0(WormCat.1, WormCat.2, WormCat.3)) %>%
    mutate(factor(species, levels = c("cel_tpm_log2fc", "cbr_tpm_log2fc")))
  
  p <- ggplot(data = markers, aes(x = log2fc, y = gene, fill = species, color = species), alpha = 0.8) +
    geom_bar(stat = "identity", position = position_dodge(), na.rm = TRUE) +
    scale_x_continuous(name = "Log2 Fold-change",
                       limits = c(min(markers$log2fc, na.rm = TRUE), max(markers$log2fc, na.rm = TRUE) * 1.1),
                       oob = scales::squish,
                       expand = c(0, 0)) +
    scale_fill_manual(values = c("cel_tpm_log2fc" = "#009E73", "cbr_tpm_log2fc" = "#56B4E9"),
                      breaks = c("cel_tpm_log2fc", "cbr_tpm_log2fc"),
                      labels = c("*C. elegans*", "*C. briggsae*")) +
    scale_color_manual(values = c("cel_tpm_log2fc" = "#00805e", "cbr_tpm_log2fc" = "#1a8acb"),
                       guide = "none") +
    theme(panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          axis.line.y = element_line(color = "grey80", size = 1),
          axis.title.y = element_blank(),
          axis.text.x = element_markdown(angle = 45, hjust = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.position = "top",
          legend.text = element_markdown(size = 12),
          panel.spacing = unit(0.15, "lines"),
          panel.border = element_rect(linewidth = 1, fill = "transparent", colour = "grey80"),
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.y = element_text(angle = 0, hjust = 0)) +
    facet_grid(WormCat.1 ~ ., scales = "free_y", space = "free_y")
  
  output[[grid_id]] <- renderPlot({
    p
  }, height = plot_height, width = 400)
}