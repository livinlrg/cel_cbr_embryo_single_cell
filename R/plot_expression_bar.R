# plot_expression_bar.R

# ------------------------------------------------------------------------------
# Function to Create the UI Grid of Cards (One Column)
# ------------------------------------------------------------------------------
generate_bar_plot_grid_ui <- function(grid_id, nRows, card_width = "800px", card_height = "400px") {
  if (nRows == 0) {
    return(tags$p("Please select at least one gene."))
  }
  
  row_list <- lapply(seq_len(nRows), function(r) {
    plot_id <- paste0(grid_id, "_plot_r", r)
    
    fluidRow(
      column(
        width = 12,  # Single column layout
        div(
          style = "border: 1px solid #ffffff; margin: 5px; padding: 5px;",
          plotlyOutput(plot_id, height = card_height, width = card_width)
        )
      )
    )
  })
  
  tagList(row_list)
}

# ------------------------------------------------------------------------------
# Function to Generate Plotly Bar Plot
# ------------------------------------------------------------------------------
generate_plotly_bar_plot <- function(gene, expr_data, cell_data, cols, scale_selector, shared_scale, exp_max) {
  gene_expr <- expr_data[gene,]
  
  arrange_by = c("ABp lineage", "MS lineage", "E lineage", "C lineage", "D lineage", "ABal lineage", "ABar lineage",
                 "Ciliated neurons", "Non-ciliated neurons", "Glia and excretory cells", "Muscle cells", "Mesoderm", "Hypodermis and seam", "Intestinal and rectal cells",
                 "Pharyngeal cells", "Germline cells")
  
  df <- data.frame(
    cell = factor(colnames(expr_data), levels = colnames(expr_data)),
    TPM = gene_expr + 1,
    cell_class = factor(cell_data[colnames(expr_data), "division_name"], levels = arrange_by),
    label = paste(
      "Cell:", colnames(expr_data),
      "Cell class:", cell_data[colnames(expr_data), "division_name", drop = TRUE],
      "<br>", gene, " TPM = ", round(gene_expr, 2)
    )
  )
  
  df$cell_color <- cols$cell_class[as.character(df$cell_class)]
  df$cell_bullet <- paste0(
    df$cell, 
    " <span style='color:", df$cell_color, "'>&#9679;</span>"
  )
  
  p <- ggplot(df, aes(x = cell, y = TPM, fill = cell_class, text = label), color = NULL) +
    geom_bar(stat = "identity", width = 1) +
    scale_x_discrete(labels = df$cell_bullet) +
    scale_fill_manual(values = cols$cell_class) +
    theme(panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          axis.line = element_line(color = "grey80", size = 2),
          axis.title.x = element_blank(),
          axis.title.y = if(grepl("C. elegans", gene)) {
            element_text(color = "#009E73", size = 10)
          } else {
            element_text(color = "#56B4E9", size = 10)
          },
          axis.text.x = element_markdown(size = 6, angle = 90, hjust = 1),
          panel.grid.major.y = element_line(color = "grey80", size = 0.25),
          panel.grid.minor.y = element_line(color = "grey80", size = 0.05),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.position = "none",
          plot.margin = margin(t = 10, r = 10, b = -750, l = 10, unit = "pt"))
  
  if (scale_selector == "identity" & shared_scale == "shared") {
    p <- p + scale_y_continuous(name = paste0(gene, " TPM"), limits = c(0, exp_max * 1.1), expand = c(0, 0))
  } else if (scale_selector == "log2" & shared_scale == "shared") {
    p <- p + scale_y_continuous(name = paste0(gene, " log2(TPM)"), limits = c(1, exp_max * 1.1), trans = "log2", expand = c(0, 0))
  } else if (scale_selector == "identity" & shared_scale == "individual") {
    p <- p + scale_y_continuous(name = paste0(gene, " TPM"), limits = c(0, max(df$TPM) * 1.1), expand = c(0, 0))
  } else if (scale_selector == "log2" & shared_scale == "individual") {
    p <- p + scale_y_continuous(name = paste0(gene, " log2(TPM)"), limits = c(1, max(df$TPM) * 1.1), trans = "log2", expand = c(0, 0))
  }
  
  ggplotly(p, tooltip = "text")
}

# ------------------------------------------------------------------------------
# Function to Render the Grid in Shiny Server
# ------------------------------------------------------------------------------
render_expression_bar <- function(output, input, session, expr_data, cell_data, cols, pro_or_term, grid_id = "expression_bar") {
  observeEvent({
    list(input$cel_genes, input$cbr_genes)
  }, {
    all_genes <- c(input$cel_genes, input$cbr_genes)  # Include genes from both species
    nRows <- length(all_genes)
    
    if(nRows > 0) {
      exp_max <- max(expr_data[c(all_genes),], na.rm = TRUE)
    } else {
      exp_max <- Inf
    }
    
    # Generate UI for the grid
    output[[grid_id]] <- renderUI({
      tagList(
        generate_legend_card(cols),
        generate_bar_plot_grid_ui(grid_id, nRows, card_width = ifelse(pro_or_term == "pro", "2800px", "1000px"))
      )
    })
    
    # Generate each plot
    for (r in seq_len(nRows)) {
      local({
        local_r <- r
        plot_id <- paste0(grid_id, "_plot_r", local_r)
        gene <- all_genes[local_r]
        
        output[[plot_id]] <- renderPlotly({
          req(gene)
          generate_legend_card(cols)
          generate_plotly_bar_plot(gene, expr_data, cell_data, cols, input$scale_selector, input$shared_scale, exp_max + 1)
        })
      })
    }
  })
}
