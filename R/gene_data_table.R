
# gene_data_table.R
build_gene_table <- function(input_gene_data, input, output, session) {
  
  use_tables <- c("gene", "briggsae_gene_short_name",
                  "filter",
                  "jsd_median_joint", "jsd_median_term", "jsd_median_pro",
                  "cel_tau_median_joint", "cbr_tau_median_joint",
                  "cel_max_tpm_term", "cbr_max_tpm_term",
                  "cel_max_tpm_pro", "cbr_max_tpm_pro",
                  "elegans_id", "briggsae_id",
                  "WormCat.1", "WormCat.2", "WormCat.3")
  
  gene_data <- as.data.table(input_gene_data)[, ..use_tables]
  
  filterInput <- function(values, name, default = "All") {
    unique_values <- unique(values)
    
    options_list <- tagList(
      tags$option(value = "", "All"),  # Reset filter option
      tags$option(value = default, selected = "selected", default),
      tagList(lapply(setdiff(unique_values, default), function(value) {
        tags$option(value = value, value)
      }))
    )
    
    tagList(
      tags$select(
        onchange = sprintf("Reactable.setFilter('gene_table', '%s', event.target.value)", name),
        options_list,
        style = "width: 100%; height: 28px;"
      )
    )
  }
  
  distance_column <- function(color_column, maxWidth = 70, class = NULL, ...) {
    colDef(
      maxWidth = maxWidth,
      align = "center",
      class = paste("cell number", class),
      html = TRUE,
      cell = function(value, index) {
        return(input_gene_data[[color_column]][index])
      },
      ...
    )
  }
  
  value_column <- function(maxWidth = 70, color_column, class = NULL, ...) {
    colDef(
      maxWidth = maxWidth,
      align = "center",
      class = paste("cell number", class),
      style = function(value, index) {
        list(
          background = input_gene_data[[color_column]][index],
          color = "#111"
        )
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
    gene_data,
    defaultColGroup = colGroup(headerVAlign = "bottom"),
    defaultColDef = colDef(
      vAlign = "center",
      headerVAlign = "bottom",
      class = "cell",
      headerClass = "header"
    ),
    defaultSorted = "gene",
    defaultSortOrder = "asc",
    pagination = TRUE,
    showPageSizeOptions = TRUE,
    highlight = TRUE,
    height = "auto",
    rownames = FALSE,
    columnGroups = list(
      colGroup(name = "Gene Name<br>Short", columns = c("gene", "briggsae_gene_short_name"), html = TRUE),
      colGroup(name = "Jensen-Shannon<br>Distance", columns = c("jsd_median_joint", "jsd_median_term", "jsd_median_pro"), html = TRUE),
      colGroup(name = "Tau", columns = c("cel_tau_median_joint", "cbr_tau_median_joint"), html = TRUE),
      colGroup(name = "TPM Max<br>Terminal", columns = c("cel_max_tpm_term", "cbr_max_tpm_term"), html = TRUE),
      colGroup(name = "TPM Max<br>Progenitor", columns = c("cel_max_tpm_pro", "cbr_max_tpm_pro"), html = TRUE),
      colGroup(name = "Gene Name<br>WBGene", columns = c("elegans_id", "briggsae_id"), html = TRUE)
    ),
    columns = list(gene = name_column(name = "elegans"),
                   briggsae_gene_short_name = name_column(name = "briggsae", class = "border-right"),
                   filter = colDef(name = "Filter", filterable = TRUE, align = "center", maxWidth = 65, class = paste("cell number", "border-right"),
                                   filterInput = function(values) {
                                     filterInput(values, "filter")
                                   },
                                   style = function(value) {
                                     list(
                                       background = if (value == "Pass") "#44ab43" else "#ff2700",
                                       color = "#111"
                                     )
                                   },
                                   html = TRUE),
                   jsd_median_joint = distance_column(
                     name = "Joint",
                     color_column = "jsd_median_joint_scaled_col"
                   ),
                   jsd_median_term = distance_column(
                     name = "Term.",
                     color_column = "jsd_median_term_scaled_col"
                   ),
                   jsd_median_pro = distance_column(
                     name = "Pro.",
                     color_column = "jsd_median_pro_scaled_col"
                   ),
                   cel_tau_median_joint = value_column(name = "cel", color_column = "cel_tau_median_joint_col", class = "border-left"),
                   cbr_tau_median_joint = value_column(name = "cbr", color_column = "cbr_tau_median_joint_col"),
                   cel_max_tpm_term = value_column(name = "cel", color_column = "cel_max_tpm_term_col"),
                   cbr_max_tpm_term = value_column(name = "cbr", color_column = "cbr_max_tpm_term_col"),
                   cel_max_tpm_pro = value_column(name = "cel", color_column = "cel_max_tpm_pro_col"),
                   cbr_max_tpm_pro = value_column(name = "cbr", color_column = "cbr_max_tpm_pro_col"),
                   elegans_id = name_column(name = "cel", class = "border-left", maxWidth = 175, minWidth = 145),
                   briggsae_id = name_column(name = "cbr", maxWidth = 175, minWidth = 145),
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
