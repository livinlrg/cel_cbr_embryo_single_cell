# plot_grid.R

# ------------------------------------------------------------------------------
# Function to Create the UI Grid of Cards
# ------------------------------------------------------------------------------
# Function to Create the UI Grid of Cards with Fixed Size and Scrollable Container
generate_plot_grid_ui <- function(grid_id, nRows, nCols, card_width = "400px", card_height = "400px") {
  if (nCols == 0 || nRows == 0) {
    return(tags$p("Please select at least one gene for each species."))
  }
  
  # Determine the grid dimensions, capping at 10x10
  max_cols <- min(nCols, 10)
  max_rows <- min(nRows, 10)
  
  # Create a scrollable grid container with fixed plot sizes
  scroll_container <- div(
    style = paste0("display: grid; ",
             "grid-template-columns: repeat(", max_cols, ", 400px); ",
             "grid-auto-rows: 400px; ",
             "overflow-x: auto; overflow-y: auto; ",
             "max-height: calc(10 * 400px + 20px); ",
             "max-width: calc(10 * 400px + 20px); ",
             "padding: 10px; border: 1px solid #ddd;"),
    
    # Generate grid items with fixed-size cards
    lapply(seq_len(nRows), function(r) {
      lapply(seq_len(nCols), function(c) {
        plot_id <- paste0(grid_id, "_plot_r", r, "_c", c)
        
        div(
          style = paste0(
            "border: 1px solid #ffffff; margin: 0px; padding: 5px; ",
            "width:", card_width, "; height:", card_height, ";",
            "box-sizing: border-box;"
          ),
          shinycssloaders::withSpinner(
          plotlyOutput(plot_id, height = card_height, width = card_width)
          )
        )
      })
    })
  )
  
  return(scroll_container)
}

# ------------------------------------------------------------------------------
# Function to Generate Plotly Scatter Plot
# ------------------------------------------------------------------------------
generate_plotly_plot <- function(geneX, geneY, expr_data, cell_data, scale_selector, shared_scale, exp_max, exp_min, cols) {
  geneX_expr <- expr_data[geneX, ] + 1
  geneY_expr <- expr_data[geneY, ] + 1
  
  df <- data.frame(
    cell = colnames(expr_data),
    geneX = geneX_expr,
    geneY = geneY_expr,
    cell_class = cell_data[colnames(expr_data), "cell_class_lineage_group"],
    label = paste(
      "Cell:", colnames(expr_data),
      "<br>", geneX, "=", round(geneX_expr, 2),
      "<br>", geneY, "=", round(geneY_expr, 2)
    )
  )
  
  if(shared_scale == "shared") {
    plot_limits <- c(max(1, exp_min * 0.9), exp_max * 1.1)
  } else if (shared_scale == "individual") {
    plot_limits = c(max(1, min(df$geneX, df$geneY) * 0.9), max(df$geneX, df$geneY) * 1.1)
  }
  
  p <- ggplot(df, aes(x = geneX, y = geneY, color = cell_class, text = label)) +
    geom_point(shape = 21) +
    scale_x_continuous(name = paste0(geneX, " TPM"), transform = scale_selector, limits = plot_limits) +
    scale_y_continuous(name = paste0(geneY, " TPM"), transform = scale_selector, limits = plot_limits) +
    scale_color_manual(values = cols$cell_class) +
    theme(panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          axis.line = element_line(color = "grey80", size = 1),
          axis.title.x = element_text(color = "#009E73", size = 10),
          axis.title.y = element_text(color = "#56B4E9", size = 10),
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
render_gene_grid <- function(output, input, session, expr_data, cell_data, cols, pro_or_term, grid_id = "gene_grid") {
  observeEvent({
    list(input$cel_genes, input$cbr_genes)
  }, {
    cel_genes <- input$cel_genes
    cbr_genes <- input$cbr_genes
    
    nCols <- length(cel_genes)  # C. elegans genes
    nRows <- length(cbr_genes)  # C. briggsae genes
    
    if(nCols > 0 | nRows > 0) {
      exp_max <- max(expr_data[c(cel_genes, cbr_genes),], na.rm = TRUE)
      exp_min <- min(expr_data[c(cel_genes, cbr_genes),], na.rm = TRUE)
    } else {
      exp_max <- 0
      exp_min <- 0
    }
    
    # Generate UI for the grid with fixed-size cards
    output[[grid_id]] <- renderUI({
      tagList(
        generate_legend_card(cols),
        generate_plot_grid_ui(grid_id, nRows, nCols, card_width = "400px", card_height = "400px")
      )
    })
    
    # Generate each plot
    for (r in seq_len(nRows)) {
      for (c in seq_len(nCols)) {
        local({
          local_r <- r
          local_c <- c
          plot_id <- paste0(grid_id, "_plot_r", local_r, "_c", local_c)
          
          geneX <- cel_genes[local_c]
          geneY <- cbr_genes[local_r]
          
          output[[plot_id]] <- renderPlotly({
            req(geneX, geneY)
            generate_legend_card(cols)
            generate_plotly_plot(geneX, geneY, expr_data, cell_data, input$scale_selector, input$shared_scale, exp_max, exp_min, cols)
          })
        })
      }
    }
  })
}

