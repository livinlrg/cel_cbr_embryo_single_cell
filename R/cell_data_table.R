
# cell_data_table.R
build_cell_table <- function(input_cell_data, cols) {
  
  use_tables <- c("cell_type", "cell_class_lineage_group",
                  "jsd_median", "cor_median",
                  "genes_detected_bootstrap_cel", "genes_detected_bootstrap_cbr",
                  "cel_markers", "cbr_markers",
                  "cel_markers_common", "cbr_markers_common",
                  "cel_markers_non_one_to_one", "cbr_markers_non_one_to_one",
                  "cel_cell_count", "cbr_cell_count",
                  "cel_median_umi", "cbr_median_umi")
  
  cell_data <- input_cell_data[,use_tables]
  
  make_color_pal <- function(colors, bias = 1, order = "desc") {
    if (order == "desc") {
      colors <- rev(colors)
    }
    get_color <- colorRamp(colors, bias = bias)
    function(x) rgb(get_color(x), maxColorValue = 255)
  }
  
  distance_color_pal_desc <- make_color_pal(c("#ff2700", "#f8fcf8", "#44ab43"), bias = 1.3, order = "desc")
  distance_color_pal_asce <- make_color_pal(c("#ff2700", "#f8fcf8", "#44ab43"), bias = 1.3, order = "asce")

  value_pct_color <- make_color_pal(c("#ffffff", "#f2fbd2", "#c9ecb4", "#93d3ab", "#35b0ab"), bias = 2, order = "asce")
  
  filter_fct <- function(values, name) {
    
    # Create HTML select element for filtering
    tags$select(
      # Default "All" option
      tags$option(value = "", "All"),  
      
      # Dynamically create options using base R
      lapply(unique(values), function(value) {
        tags$option(value = value, value)
      }),
      
      onchange = glue(
        "Reactable.setFilter(
        'cell_table', 
        '{name}', 
        event.target.value, 
        {{ filterMethod: (rows, id, filterValue) => rows[id] === filterValue }}
      )"
      ),
      
      # Optional styling for the dropdown
      style = "width: 100%; height: 28px;"
    )
  }
  
  
  distance_column <- function(maxWidth = 75, class = NULL, ...) {
    colDef(maxWidth = maxWidth,
           align = "center",
           class = paste("cell number", class),
           ...)
  }
  
  value_column <- function(maxWidth = 85, class = NULL, max_value = 1, ...) {
    colDef(
      maxWidth = maxWidth,
      aggregate = "mean",
      align = "center",
      footer = function(values) { div(round(mean(values), 0))},
      class = paste("cell number", class),
      style = function(value) {
        if (is.numeric(value) && max_value > 0) {
          # Calculate the percentage as a scaled value
          scaled_value <- value / max_value
          list(color = "#111", background = value_pct_color(scaled_value))
        } else {
          # Default style for non-numeric or invalid values
          list(color = "#111", background = "#f8f8f8")
        }
      },
      ...
    )
  }
  
  
  name_column <- function(maxWidth = 150, minWidth = 100, class = NULL, use_filter = TRUE, col_name = NULL, bg_color = FALSE, cols = NULL, ...) {
    colDef(
      maxWidth = maxWidth,
      minWidth = minWidth,
      filterable = TRUE,
      align = "center",
      filterInput = if (use_filter && !is.null(col_name)) {
        function(values) filter_fct(values, col_name)
      } else {
        NULL
      },
      class = paste("cell number", class),
      style = if (bg_color) {
        function(value) {
          if (!is.null(cols) && !is.na(value)) {
            color_value <- cols[[trimws(value)]]
            list(background = color_value, color = "black")
          } else {
            list()
          }
        }
      } else {
        NULL
      },
      ...
    )
  }
  
  reactable(cell_data,
            defaultColGroup = colGroup(headerVAlign = "bottom"),
            defaultColDef = colDef(
              vAlign = "center",
              headerVAlign = "bottom",
              class = "cell",
              headerClass = "header"
            ),
            defaultSorted = "cell_type",
            defaultSortOrder = "desc",
            showPageSizeOptions = TRUE,
            # selection = "single",
            # onClick = "select",
            highlight = TRUE,
            height = "auto",
            rownames = FALSE,
            columnGroups = list(
              colGroup(name = "Distances", columns = c("jsd_median", "cor_median"), html = TRUE),
              colGroup(name = "Genes<br>Detected", columns = c("genes_detected_bootstrap_cel", "genes_detected_bootstrap_cbr"), html = TRUE),
              colGroup(name = "Total<br>Markers", columns = c("cel_markers", "cbr_markers"), html = TRUE),
              colGroup(name = "1:1<br>Markers", columns = c("cel_markers_common", "cbr_markers_common"), html = TRUE),
              colGroup(name = "Non-1:1<br>Markers", columns = c("cel_markers_non_one_to_one", "cbr_markers_non_one_to_one"), html = TRUE),
              colGroup(name = "Cell count", columns = c("cel_cell_count", "cbr_cell_count"), html = TRUE),
              colGroup(name = "Median<br>UMI", columns = c("cel_median_umi", "cbr_median_umi"), html = TRUE)
            ),
            columns = list(cell_type = name_column(name = "Cell type", col_name = "cell_type", maxWidth = 200, minWidth = 100),
                           cell_class_lineage_group = name_column(name = "Cell class", col_name = "cell_class_lineage_group", class = "border-right", bg_color = TRUE, cols = cols,
                                                                  footer = function(values) { div(tags$b("Average: "))}),
                           jsd_median = distance_column(
                             name = "JSD",
                             cell = function(value) {
                               scaled <- (value - min(cell_data$jsd_median)) / (max(cell_data$jsd_median) - min(cell_data$jsd_median))
                               color <- distance_color_pal_desc(scaled)
                               value <- format(round(value, 2), nsmall = 1)
                               div(class = "circle_value_col", style = list(background = color), value)
                             },
                             aggregate = "mean",
                             footer = function(values) { div(round(mean(values), 2))}
                           ),
                           cor_median = distance_column(
                             name = "Corr.",
                             cell = function(value) {
                               scaled <- (value - min(cell_data$cor_median)) / (max(cell_data$cor_median) - min(cell_data$cor_median))
                               color <- distance_color_pal_asce(scaled)
                               value <- format(round(value, 2), nsmall = 1)
                               div(class = "circle_value_col", style = list(background = color), value)
                             },
                             aggregate = "mean",
                             footer = function(values) { div(round(mean(values), 2))},
                           ),
                           genes_detected_bootstrap_cel = value_column(name = "elegans", max_value = max(cell_data[,"genes_detected_bootstrap_cel"]), class = "border-left"),
                           genes_detected_bootstrap_cbr = value_column(name = "briggsae", max_value = max(cell_data[,"genes_detected_bootstrap_cbr"])),
                           cel_markers = value_column(name = "elegans", max_value = max(cell_data[,"cel_markers"])),
                           cbr_markers = value_column(name = "briggsae", max_value = max(cell_data[,"cbr_markers"])),
                           cel_markers_common = value_column(name = "elegans", max_value = max(cell_data[,"cel_markers_common"])),
                           cbr_markers_common = value_column(name = "briggsae", max_value = max(cell_data[,"cbr_markers_common"])),
                           cel_markers_non_one_to_one = value_column(name = "elegans", max_value = max(cell_data[,"cel_markers_non_one_to_one"])),
                           cbr_markers_non_one_to_one = value_column(name = "briggsae", max_value = max(cell_data[,"cbr_markers_non_one_to_one"])),
                           cel_markers_non_one_to_one = value_column(name = "elegans", max_value = max(cell_data[,"cel_markers_non_one_to_one"])),
                           cbr_markers_non_one_to_one = value_column(name = "briggsae", max_value = max(cell_data[,"cbr_markers_non_one_to_one"])),
                           cel_cell_count = value_column(name = "elegans", max_value = max(cell_data[,"cel_cell_count"])),
                           cbr_cell_count = value_column(name = "briggsae", max_value = max(cell_data[,"cbr_cell_count"])),
                           cel_median_umi = value_column(name = "elegans", max_value = max(cell_data[,"cel_median_umi"])),
                           cbr_median_umi = value_column(name = "briggsae", max_value = max(cell_data[,"cbr_median_umi"]))
            ),
            
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
