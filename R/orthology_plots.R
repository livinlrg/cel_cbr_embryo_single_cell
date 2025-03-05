

# ------------------------------------------------------------------------------
# Function to Create the UI Grid of Cards (One Column)
# ------------------------------------------------------------------------------
generate_tree_plot_grid_ui <- function(grid_id, length_of_tree) {
  
  row_list <- lapply(1, function(r) {
    plot_id <- paste0(grid_id, "_plot_r", r)
    card_id <- paste0(grid_id, "_card_r", r)
    description_id <- paste0(grid_id, "_description_r", r)  # Add this line
    
    fluidRow(
      column(
        width = 12,
        div(
          class = "tree-plot-container",
          style = "border: 1px solid #ddd; margin: 10px; padding: 10px; background-color: #f9f9f9;",
          uiOutput(card_id),
          uiOutput(description_id),
          plotOutput(plot_id, height = paste0(400 + length_of_tree * 0.17, "px"))
        )
      )
    )
  })
  
  tagList(row_list)
}

# ------------------------------------------------------------------------------
# Function to Render the Grid in Shiny Server
# ------------------------------------------------------------------------------
render_tree_plots <- function(output, input, session, og, og_data, synteny, gene_data, grid_id = "orthology_plot") {
  observeEvent(og(), {
    req(og())
    og_list <- og()
    nRows <- length(og_list)
    
    if (nRows == 0) {
      warning("No orthologous groups (OG) found.")
      return()
    }
    
    length_of_tree <- nchar(og_data[og_list[[1]],"tree"])
    
    species_colors <- c("elegans" = "#009E73",
                        "briggsae" = "#56B4E9",
                        "becei" = "#a6cee3",
                        "zanzibari" = "#b2df8a",
                        "inopinata" = "#fb9a99",
                        "latens" = "#e31a1c",
                        "sulstoni" = "#fdbf6f",
                        "panamensis" = "#ff7f00",
                        "remanei" = "#cab2d6",
                        "tribulationis" = "#6a3d9a",
                        "nigoni" = "#ffff99",
                        "tropicalis" = "#b15928",
                        "brenneri" = "darkgrey")

    # Generate UI for the grid
    output[[grid_id]] <- renderUI({
      tagList(
        generate_tree_legend_card(species_colors),
        generate_tree_plot_grid_ui(grid_id, length_of_tree)
      )
    })
    
    # Generate each plot
    for (r in seq_len(nRows)) {
      local({
        local_r <- r
        plot_id <- paste0(grid_id, "_plot_r", local_r)
        card_id <- paste0(grid_id, "_card_r", local_r)
        description_id <- paste0(grid_id, "_description_r", local_r)
        
        cur_og <- og_list[local_r]
        wb_genes <- strsplit(og_data[cur_og,]$WBgenes, ", ")[[1]]
        
        output[[card_id]] <- renderUI({
          req(cur_og, og_data)
          generate_og_card(cur_og, og_data)
        })
        
        output[[description_id]] <- renderUI({
          req(cur_og, og_data)
          generate_og_description(wb_genes, synteny, gene_data)
        })
        
        output[[plot_id]] <- renderPlot({
          req(cur_og, og_data)
          plot_trees(cur_og, og_data, species_colors)
        })
      })
    }
  })
}

plot_trees <- function(cur_og, og_data, species_colors) {
  tree_text <- og_data[cur_og, "tree"]
  if (is.na(tree_text) || tree_text == "") return(NULL)
  
  cur_tree <- read.tree(text = tree_text)
  
  # Extract species and tip labels
  species_info <- sapply(str_split(cur_tree$tip.label, "_"), function(x) head(x, n = 1))
  gene_labels <- sapply(str_split(cur_tree$tip.label, "_"), function(x) tail(x, n = 1))
  
  # Convert tree labels to gene names
  cur_tree$tip.label <- gene_labels
  
  # Create a data frame with tip metadata
  tip_data <- tibble(
    label = gene_labels,  # Must match tree tip labels
    species = as.factor(species_info),  # Ensure factor for correct coloring
    fontface = ifelse(species_info %in% c("elegans", "briggsae"), "bold", "plain")
  )
  
  gg_tr <- ggtree(cur_tree) %<+% tip_data +  # Attach species metadata
    geom_tiplab(aes(label = label, fontface = fontface), align = TRUE, offset = 0.05) +
    geom_tippoint(aes(fill = species), color = "black", shape = 21, size = 4) +  # Corrected species mapping
    theme_tree() +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = species_colors) +
    theme(plot.margin = margin(10, 150, 10, 10, unit = "pt"),
          legend.position = "none")
  
  return(gg_tr)
}

generate_og_card <- function(cur_og, og_data) {
  og_info <- og_data[cur_og, ]
  
  tags$div(
    class = "og-card",
    style = "border: 1px solid #ddd; margin: 5px; padding: 5px;",
    tags$h4(paste0("Orthologous Group: ", cur_og)),
    tags$p(paste0("Gene names: ", og_info$gene_names, " WBGene names: ", og_info$WBgenes))
  )
}

generate_og_description <- function(wb_genes, synteny, gene_data) {
  synteny_temp <- synteny[synteny$cel_WBG_name %in% wb_genes & synteny$cbr_WBG_name %in% wb_genes,]
  
  description_list <- lapply(seq_len(nrow(synteny_temp)), function(i) {
    cur_synteny <- synteny_temp[i,]
    
    if(cur_synteny$cel_common_name %in% rownames(gene_data)) {
      temp_cbr_name <- gene_data[cur_synteny$cel_common_name, "briggsae_id"]
    } else {
      temp_cbr_name <- "Nope"
    }
    
    if (temp_cbr_name == cur_synteny$cbr_WBG_name) {
      confidence <- ifelse(cur_synteny["one_to_one"] == 1, "High",
                           ifelse(cur_synteny["confident_canonical"] == 1, "Medium",
                                  "Low"))
      
      return(
        tags$div(
          class = "og-card",
          style = "border: 1px solid #ddd; margin: 5px; padding: 5px;",
          tags$p(HTML(paste0("<b>", cur_synteny$cel_common_name, " and ", cur_synteny$cbr_common_name, " are classified as 1:1 orthologs</b>"))),
          tags$p(paste0("Confidence: ", confidence)),
          tags$p(paste0("Syntenic: ", ifelse(cur_synteny$syntenic == "IS_SYNTENIC", "Yes", "No"),
                        " | Mutual best BLAT: ", ifelse(cur_synteny$mutual_best_blat == "1", "Yes", "No"),
                        " | Mutual best Smith-Waterman: ", ifelse(cur_synteny$mutual_best_sw == "1", "Yes", "No"))),
          tags$p(paste0("C. elegans aligned: ", round(cur_synteny$cel_percent_align, 1), "% | C. briggsae aligned: ", round(cur_synteny$cbr_percent_align, 1), "%",
                        " | ", "Percent identity: ", round(cur_synteny$percent_identity, 1), "%"))
        )
      )
    } else {
      return(
        tags$div(
          class = "og-card",
          style = "border: 1px solid #ddd; margin: 5px; padding: 5px;",
          tags$p(HTML(paste0("<b>", cur_synteny$cel_common_name, " and ", cur_synteny$cbr_common_name, " are not 1:1 orthologs</b>"))),
          tags$p(paste0("Syntenic: ", ifelse(cur_synteny$syntenic == "IS_SYNTENIC", "Yes", "No"),
                        " | Mutual best BLAT: ", ifelse(cur_synteny$mutual_best_blat == "1", "Yes", "No"),
                        " | Mutual best Smith-Waterman: ", ifelse(cur_synteny$mutual_best_sw == "1", "Yes", "No"))),
          tags$p(paste0("C. elegans aligned: ", round(cur_synteny$cel_percent_align, 1), "% | C. briggsae aligned: ", round(cur_synteny$cbr_percent_align, 1), "%",
                        " | ", "Percent identity: ", round(cur_synteny$percent_identity, 1), "%"))
        )
      )
    }
  })
  return(tagList(description_list))  # Wrap the entire list in a tagList()
}
