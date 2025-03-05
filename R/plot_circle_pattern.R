# plot_circle_pattern.R

library(shiny)
library(ggplot2)
library(ggraph)
library(dplyr)

# Function to Create a Grid of Circle Pattern Plots (Two Columns)
generate_circle_pattern_grid_ui <- function(grid_id, nPlots, card_width = "400px", card_height = "400px") {
  if (nPlots == 0) {
    return(tags$p("Please select at least one gene."))
  }
  
  # Generate rows with two columns each
  scroll_container <- div(
    style = "display: flex; flex-wrap: wrap; overflow-x: auto; padding: 10px; border: 1px solid #ddd;",
    
    # Generate plots as fixed-size cards within a grid layout
    lapply(seq_len(nPlots), function(plot_index) {
      plot_id <- paste0(grid_id, "_plot_", plot_index)
      div(
        style = paste0("flex: 0 0 calc(50% - 10px); ", 
                       # "display: inline-block; vertical-algin: top; ",
                       "box-sizing: border-box;",
                       "border: 1px solid #ffffff; margin: 5px; padding: 5px;",
                       "width:", card_width, "; height:", card_height, ";"),
        shinycssloaders::withSpinner(
        plotOutput(plot_id, height = card_height, width = card_width)
        )
      )
    })
  )
  
  return(scroll_container)
}

plot_circle_pattern <- function(cur_gene, cur_species, CellTable, lo, expr_data, scale_selector, shared_scale, exp_max) {
  
  plot_species <- ifelse(cur_species == "C.elegans", "C. elegans", "C. briggsae")
  
  temp_exp_df <- data.frame(cell_type = colnames(expr_data),
                            tpm = expr_data[paste0(plot_species, " ", cur_gene),])
  
  CellTable <- left_join(CellTable,
                         temp_exp_df,
                         c("MergedDatasetName" = "cell_type"))
  
  layout <- left_join(lo, CellTable[,c("Lineage", "level", "tpm")],
                      by = c("name" = "Lineage"))
  
  layout <- layout[! is.na(layout$level),]
  layout$level <- as.numeric(layout$level)
  
  if(shared_scale == "shared") {
    plot_limits <- c(8, exp_max)
  } else if (shared_scale == "individual") {
    plot_limits = c(8, max(layout$tpm, na.rm = TRUE))
  }
  
  circle_plot <- ggraph(layout) +
    geom_node_arc_bar(aes(filter = level > 2, fill = tpm + 1),
                      color = "grey20", size = 0.25) +
    coord_fixed() +
    geom_node_text(aes(filter = level < 4 & level > 2,
                       label = name,
                       angle = node_angle(x, y) + 90), size = 3, colour = 'white') +
    annotate(geom = "text", x = 0.025, y = 0.75, hjust = 0.5, size = 9,
             label = paste0("italic('", plot_species, "')"), parse = TRUE,
             color = ifelse(cur_species == "C.elegans", "#009E73", "#56B4E9"),
             fontface = 2) +
    annotate(geom = "text", x = 0.025, y = 0.175, hjust = 0.5, size = 6,
             label = paste0("italic('", cur_gene, "')~' TPM'"), parse = TRUE) +
    scale_fill_viridis(trans = scale_selector, option = "magma",
                       limits = plot_limits, oob = scales::squish) +
    guides(alpha = F, size = F, fill = guide_colorbar(title = "")) +
    theme(legend.position = c(0.49,0.42),
          panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          plot.margin = margin(t = -5, r = -20, b = -5, l = -20, unit = "mm"),
          legend.direction = "horizontal",
          legend.justification = "center",
          legend.key.size = unit(9, 'mm'),
          legend.text = element_text(size = 12))
  
  return(circle_plot)
}

render_circle_pattern_grid <- function(output, input, session, CellTable, lo, expr_data, grid_id = "circle_pattern_grid") {
  
  observeEvent({
    list(input$cel_genes, input$cbr_genes)
  }, {
    all_genes <- c(input$cel_genes, input$cbr_genes)  # Genes from both species
    nPlots <- length(all_genes)
    
    if(nPlots > 0) {
      exp_max <- max(expr_data[c(all_genes),], na.rm = TRUE)
    } else {
      exp_max <- 0
    }
    
    # Generate the UI for the plot grid
    output[[grid_id]] <- renderUI({
      generate_circle_pattern_grid_ui(grid_id, nPlots)
    })
    
    # Render each plot in the grid
    for (i in seq_len(nPlots)) {
      local({
        local_i <- i
        plot_id <- paste0(grid_id, "_plot_", local_i)
        gene <- all_genes[local_i]
        species <- ifelse(grepl("C. elegans", gene), "C.elegans", "C.briggsae")
        
        cur_gene <- strsplit(gene, " ")[[1]][3]
        
        output[[plot_id]] <- renderPlot({
          req(cur_gene)
          plot_circle_pattern(cur_gene, species, CellTable, lo, expr_data, input$scale_selector, input$shared_scale, exp_max)
        })
      })
    }
  })
}