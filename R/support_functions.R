
# ------------------------------------------------------------------------------
# Function to Generate the Legend as a Card
# ------------------------------------------------------------------------------
generate_legend_card <- function(cols) {
  tags$div(
    class = "legend-card",
    style = "display: flex; flex-wrap: wrap; padding: 10px; border: 1px solid #ddd; margin: 5px; background-color: white;",
    lapply(names(cols$cell_class), function(class_label) {
      color_val <- cols$cell_class[[class_label]]
      tags$div(
        style = "display: flex; align-items: center; margin-right: 10px; margin-bottom: 5px;",
        tags$span(
          style = paste0(
            "display: inline-block; width: 15px; height: 15px; ",
            "border-radius: 2px; margin-right: 5px; background-color:", color_val, ";"
          )
        ),
        tags$span(class_label, style = "font-size: 12px;")
      )
    })
  )
}

generate_tree_legend_card <- function(cols) {
  tags$div(
    class = "tree-legend-card",
    style = "display: flex; flex-wrap: wrap; gap: 10px; align-items: center; padding: 10px; border: 1px solid #ddd; margin: 5px; background-color: white; max-width: 100%;",
    lapply(names(cols), function(species_label) {
      color_val <- cols[[species_label]]
      tags$div(
        style = "display: flex; align-items: center; gap: 5px; margin-bottom: 5px;",
        tags$span(
          style = paste0(
            "display: inline-block; width: 15px; height: 15px; ",
            "border-radius: 2px; background-color:", color_val, ";"
          )
        ),
        tags$span(HTML(paste0("<i>C. ", species_label, "</i>")), style = "font-size: 12px; white-space: nowrap;")
      )
    })
  )
}

# Function to perform gene search and selection
search_and_select_gene <- function(input, output, session, cur_gene, gff_list_mrna, expr_data, og_data) {
  req(cur_gene) # Ensure a gene is selected
  
  temp_gene_selection <- strsplit(cur_gene, " ")[[1]][3] # Extract the gene name
  cur_species <- strsplit(cur_gene, " ")[[1]][2]
  
  # determine whether the gene is a cds gene name, cosmid gene name, or wbgene name
  if (grepl("WBGene", temp_gene_selection)) {
    temp_gene_selection <- gff_list_mrna[[cur_species]][gff_list_mrna[[cur_species]]$wbgene == temp_gene_selection, "cds_gene_name"]
    temp_og <- gff_list_mrna[[cur_species]][temp_gene_selection, "OG"]
  } else if (temp_gene_selection %in% gff_list_mrna[[cur_species]]$cds_gene_name) {
    temp_og <- gff_list_mrna[[cur_species]][temp_gene_selection, "OG"]
  } else if (temp_gene_selection %in% gff_list_mrna[[cur_species]]$gene_long_name) {
    temp_gene_selection <- gff_list_mrna[[cur_species]][gff_list_mrna[[cur_species]]$gene_long_name == temp_gene_selection, "cds_gene_name"]
    temp_og <- gff_list_mrna[[cur_species]][temp_gene_selection, "OG"]
  } else if (temp_gene_selection %in% gff_list_mrna[[cur_species]]$gene_short_name) {
    temp_gene_selection <- gff_list_mrna[[cur_species]][gff_list_mrna[[cur_species]]$gene_short_name == temp_gene_selection, "cds_gene_name"]
    temp_og <- gff_list_mrna[[cur_species]][temp_gene_selection, "OG"]
  } else {
    temp_og <- NA
  }
  
  if (!is.null(temp_og) && temp_og %in% rownames(og_data)) {
    temp_og_wb_genes <- strsplit(og_data[temp_og, "WBgenes"], ", ")[[1]]
  } else {
    temp_og_wb_genes <- c()
  }
  
  # Retrieve ortholog genes
  cel_ortho_genes <- unlist(lapply(temp_og_wb_genes, function(x) 
    gff_list_mrna[["elegans"]][grepl(x, gff_list_mrna[["elegans"]]$wbgene),]$cds_gene_name))
  
  cbr_ortho_genes <- unlist(lapply(temp_og_wb_genes, function(x) 
    gff_list_mrna[["briggsae"]][grepl(x, gff_list_mrna[["briggsae"]]$wbgene),]$cds_gene_name))
  
  # Check if either species has more than 10 genes in this orthogroup
  if (length(cel_ortho_genes) > 10 || length(cbr_ortho_genes) > 10) {
    showModal(
      modalDialog(
        title = "Large Orthogroup",
        paste(
          "The selected gene's orthogroup has more than 10 orthologous genes in at least one species.",
          "We are only displaying the first 10 genes per species."
        ),
        easyClose = TRUE
      )
    )
  }
  
  # Limit to at most 10 genes in each species
  cel_ortho_genes <- head(cel_ortho_genes, 10)
  cbr_ortho_genes <- head(cbr_ortho_genes, 10)
  
  # Update Select Inputs
  updateSelectizeInput(session, "cel_genes",
                       selected = paste0("C. elegans ", cel_ortho_genes),
                       choices = rownames(expr_data)[grepl("C. elegans", rownames(expr_data))],
                       server = TRUE
  )
  
  updateSelectizeInput(session, "cbr_genes",
                       selected = paste0("C. briggsae ", cbr_ortho_genes),
                       choices = rownames(expr_data)[grepl("C. briggsae", rownames(expr_data))],
                       server = TRUE
  )
}

search_and_select_cell <- function(input, output, session, cur_cell, CellTable) {
  req(cur_cell)
  
  if (is.na(cur_cell)) {
    return(NULL)
  } else if (cur_cell %in% CellTable$MergedDatasetName) {
    return(cur_cell)
  } else if (cur_cell %in% CellTable$Lineage) {
    return(CellTable[which(CellTable$Lineage == cur_cell), "MergedDatasetName"])
  }
}
