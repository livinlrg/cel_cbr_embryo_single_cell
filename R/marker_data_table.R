

# marker_data_table.R
build_marker_table <- function(input_marker_data, input, output, session) {

  use_tables <- c("gene", "briggsae_gene_short_name",
                  "cell_type",
                  "shared_marker", "gene_orthology",
                  "cel_tpm_log2fc", "cbr_tpm_log2fc", "p_val_adj.cel", "p_val_adj.cbr",
                  "cel_tpm", "cbr_tpm",
                  "cel_max_tpm_term", "cbr_max_tpm_term",
                  "cel_max_tpm_pro", "cbr_max_tpm_pro",
                  "cel_tau_joint", "cbr_tau_joint",
                  "elegans_id", "briggsae_id",
                  "WormCat.1", "WormCat.2", "WormCat.3")

  marker_data <- as.data.table(input_marker_data)[, ..use_tables]

  value_column <- function(maxWidth = 65, color_column, class = NULL, ...) {
    colDef(
      maxWidth = maxWidth,
      align = "center",
      class = paste("cell number", class),
      style = function(value, index) {
        list(
          background = input_marker_data[[color_column]][index],
          color = "#111"
        )
      },
      ...
    )
  }
  
  circle_column <- function(color_column, maxWidth = 55, class = NULL, ...) {
    colDef(
      maxWidth = maxWidth,
      align = "center",
      class = paste("cell number", class),
      html = TRUE,
      cell = function(value, index) {
        color <- input_marker_data[[color_column]][index]
        return(div(class = "circle_value_col", style = list(background = color), value))
      },
      ...
    )
  }
  
  p_val_circle_column <- function(color_column, maxWidth = 55, class = NULL, ...) {
    colDef(
      maxWidth = maxWidth,
      align = "center",
      class = paste("cell number", class),
      html = TRUE,
      cell = function(value, index) {
        color <- input_marker_data[[color_column]][index]
        return(div(class = "circle_value_col",
                   style = list(background = color),
                   sprintf("%.0e", value)))
      },
      ...
    )
  }
  
  marker_circle_column <- function(color_column, maxWidth = 55, class = NULL, ...) {
    colDef(
      maxWidth = maxWidth,
      align = "center",
      class = paste("cell number", class),
      html = TRUE,
      cell = function(value, index) {
        color <- if (value) "#44ab43" else "#ff2700"
        return(div(class = "circle_value_col", style = list(background = color), " "))
      },
      ...
    )
  }

  name_column <- function(maxWidth = 125, minWidth = 95, class = NULL, ...) {
    colDef(
      maxWidth = maxWidth,
      minWidth = minWidth,
      filterable = TRUE,
      align = "center",
      class = paste("cell number", class),
      html = TRUE,
      ...
    )
  }

  reactable(
    marker_data,
    defaultColGroup = colGroup(headerVAlign = "bottom"),
    defaultColDef = colDef(
      vAlign = "center",
      headerVAlign = "bottom",
      class = "cell",
      headerClass = "header"
    ),
    defaultSorted = "cel_tpm_log2fc",
    defaultSortOrder = "desc",
    pagination = TRUE,
    showPageSizeOptions = TRUE,
    highlight = TRUE,
    height = "auto",
    rownames = FALSE,
    columnGroups = list(
      colGroup(name = "Gene name<br>Short", columns = c("gene", "briggsae_gene_short_name"), html = TRUE),
      colGroup(name = "log2 Fold-Change", columns = c("cel_tpm_log2fc", "cbr_tpm_log2fc"), html = TRUE),
      colGroup(name = "Adj. p-value", columns = c("p_val_adj.cel", "p_val_adj.cbr"), html = TRUE),
      colGroup(name = "TPM<br>Cell Type", columns = c("cel_tpm", "cbr_tpm"), html = TRUE),
      colGroup(name = "TPM Max<br>Terminal", columns = c("cel_max_tpm_term", "cbr_max_tpm_term"), html = TRUE),
      colGroup(name = "TPM Max<br>Progenitor", columns = c("cel_max_tpm_pro", "cbr_max_tpm_pro"), html = TRUE),
      colGroup(name = "Tau", columns = c("cel_tau_joint", "cbr_tau_joint"), html = TRUE),
      colGroup(name = "Gene name<br>WBGene", columns = c("elegans_id", "briggsae_id"), html = TRUE)
    ),
    columns = list(gene = name_column(name = "elegans"),
                   briggsae_gene_short_name = name_column(name = "briggsae", class = "border-right"),
                   cell_type = name_column(name = "Cell Type"),
                   shared_marker = marker_circle_column(name = "Shared", filterable = TRUE, maxWidth = 85),
                   gene_orthology = name_column(name = "Ortho.", class = "border-right", maxWidth = 85),
                   cel_tpm_log2fc = circle_column(name = "cel", color_column = "cel_tpm_log2fc_col"),
                   cbr_tpm_log2fc = circle_column(name = "cbr", color_column = "cbr_tpm_log2fc_col"),
                   p_val_adj.cel = p_val_circle_column(name = "cel", color_column = "p_val_adj.cel_col"),
                   p_val_adj.cbr = p_val_circle_column(name = "cbr", color_column = "p_val_adj.cbr_col", class = "border-right"),
                   cel_tpm = value_column(name = "cel", color_column = "cel_max_tpm_term_col"),
                   cbr_tpm = value_column(name = "cbr", color_column = "cbr_max_tpm_term_col"),
                   cel_max_tpm_term = value_column(name = "cel", color_column = "cel_max_tpm_term_col"),
                   cbr_max_tpm_term = value_column(name = "cbr", color_column = "cbr_max_tpm_term_col"),
                   cel_max_tpm_pro = value_column(name = "cel", color_column = "cel_max_tpm_pro_col"),
                   cbr_max_tpm_pro = value_column(name = "cbr", color_column = "cbr_max_tpm_pro_col"),
                   cel_tau_joint = value_column(name = "cel", color_column = "cel_tau_joint_col", class = "border-left"),
                   cbr_tau_joint = value_column(name = "cbr", color_column = "cbr_tau_joint_col"),
                   elegans_id = name_column(name = "elegans", class = "border-left", maxWidth = 175, minWidth = 145),
                   briggsae_id = name_column(name = "briggsae", maxWidth = 175, minWidth = 145),
                   WormCat.1 = name_column(name = "WormCat.1", class = "border-left", maxWidth = 145, minWidth = 125),
                   WormCat.2 = name_column(name = "WormCat.2", maxWidth = 145, minWidth = 125),
                   WormCat.3 = name_column(name = "WormCat.3", maxWidth = 145, minWidth = 125)
    )

    # details = function(index) {
    #   if (index == 3) {
    #     tabsetPanel(
    #       tabPanel("plot", plotOutput("plot")),
    #       tabPanel("subtable", reactable(iris[1:3, 1:2], fullWidth = FALSE))
    #     )
    #   } else if (index == 5) {
    #     paste("Details for row:", index)
    #   }
    # }
  )
}
