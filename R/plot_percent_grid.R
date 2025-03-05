# plot_percent_grid.R

# ------------------------------------------------------------------------------
# Function to Create the UI Grid of Cards (One Column)
# ------------------------------------------------------------------------------
generate_percent_plot_grid_ui <- function(grid_id, nRows) {
  if (nRows == 0) {
    return(tags$p("Please select at least one gene."))
  }
  
  row_list <- lapply(seq_len(nRows), function(r) {
    plot_id <- paste0(grid_id, "_plot_r", r)
    
    fluidRow(
      column(
        width = 12,  # Single column layout
        div(
          style = "border: 1px solid #ddd; margin: 5px; padding: 5px;",
          shinycssloaders::withSpinner(
          plotlyOutput(plot_id, height = "250px")
          )
        )
      )
    )
  })
  
  tagList(row_list)
}

# ------------------------------------------------------------------------------
# Function to Generate Plotly Scatter Plot (TPM vs. Percent Expressed)
# ------------------------------------------------------------------------------
generate_plotly_percent_plot <- function(gene, expr_data, perc_data, cell_data, cols, scale_selector, shared_scale, exp_max, exp_min) {
  gene_expr <- expr_data[gene,]
  gene_perc <- perc_data[gene,]
  
  df <- data.frame(
    cell = colnames(expr_data),
    TPM = gene_expr + 1,
    PercentExpressed = gene_perc,
    cell_class = cell_data[colnames(expr_data), "cell_class_lineage_group"],
    label = paste(
      "Cell:", colnames(expr_data),
      "<br>", gene, " TPM =", round(gene_expr, 2),
      "<br> Percent Expressed =", round(gene_perc, 2), "%"
    )
  )
  
  if(shared_scale == "shared") {
    plot_limits <- c(max(1, exp_min * 0.9), exp_max * 1.1)
  } else if (shared_scale == "individual") {
    plot_limits = c(max(1, min(df$TPM) * 0.9), max(df$TPM) * 1.1)
  }
  
  p <- ggplot(df, aes(x = PercentExpressed, y = TPM, color = cell_class, text = label)) +
    geom_point(shape = 21) +
    scale_x_continuous(name = "Cell (%) Expressing",
                       limits = c(0, 100)) +
    scale_y_continuous(name = paste0(gene, " TPM"), trans = scale_selector,
                       limits = plot_limits) +
    scale_color_manual(values = cols$cell_class) +
    theme(panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          axis.line = element_line(color = "grey80", size = 1),
          # axis.title.x = element_text(color = "black", size = 10),
          axis.title.y = if(grepl("C. elegans", gene)) {
            element_text(color = "#009E73", size = 10)
            } else {
              element_text(color = "#56B4E9", size = 10)
            },
          panel.grid.major = element_line(color = "grey80", size = 0.25),
          panel.grid.minor = element_line(color = "grey80", size = 0.05),
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.position = "none")
  
  ggplotly(p, tooltip = "text")
}

# ------------------------------------------------------------------------------
# Function to Render the Grid in Shiny Server
# ------------------------------------------------------------------------------
render_percent_grid <- function(output, input, session, expr_data, perc_data, cell_data, cols, pro_or_term, grid_id = "percent_grid") {
  observeEvent({
    list(input$cel_genes, input$cbr_genes)
  }, {
    all_genes <- c(input$cel_genes, input$cbr_genes)
    nRows <- length(all_genes)
    
    if(nRows > 0) {
      exp_max <- max(expr_data[c(all_genes),], na.rm = TRUE)
      exp_min <- min(expr_data[c(all_genes),], na.rm = TRUE)
    } else {
      exp_max <- 0
      exp_min <- 0
    }
    
    # Generate UI for the grid
    output[[grid_id]] <- renderUI({
      tagList(
        generate_legend_card(cols),
        generate_percent_plot_grid_ui(grid_id, nRows)
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
          generate_plotly_percent_plot(gene, expr_data, perc_data, cell_data, cols, input$scale_selector, input$shared_scale, exp_max, exp_min)
        })
      })
    }
  })
}
