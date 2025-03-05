# app.R

library(shiny)
library(shinyjs)
library(shinythemes)
library(bslib)
library(ggplot2)
library(plotly)
library(ggraph)
library(igraph)
library(ggtext)
library(reactable)
library(glue)
library(data.table)
library(cowplot)
library(viridis)
library(pheatmap)
library(ape)
library(ggtree)
library(stringr)

# sources all R scripts in the R/ directory
r_source_files <-  c("R/cell_data_table.R", "R/cell_markers_page.R",
                     "R/cell_page.R", "R/gene_data_table.R",
                     "R/marker_data_table.R", "R/orthology_plots.R",
                     "R/plot_circle_pattern.R", "R/plot_expression_bar.R",
                     "R/plot_percent_grid.R", "R/plot_tpm_grid.R",
                     "R/support_functions.R")

lapply(r_source_files, source)

source("R/extra_text.R", local = TRUE)

# ------------------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------------------
expr_data <- readRDS("./data/expr_data.rds")
perc_data <- readRDS("./data/perc_data.rds")

gene_data <- readRDS("./data/gene_data_joint_cell_bg.rds")
cell_data <- readRDS("./data/cell_data_joint_mean_cell_bg.rds")
og_data <- readRDS("./data/og_data.rds")
gff_list_mrna <- readRDS("./data/gff_list_mrna.rds")
markers <- readRDS("./data/joint_markers_cell_bg.rds")

marker_enrich_list <- readRDS("./data/marker_wormcat_enrichment_cell_bg.rds")
markers_enrich_WormCat.1 <- data.frame(marker_enrich_list[["C.elegans"]][["WormCat.1"]])

CellTable <- readRDS("./data/CellTable.rds")
lo <- readRDS("./data/lo.rds")

synteny <- readRDS("./data/synteny_filt.rds")

cell_data$cell_class_lineage_group <- ifelse(is.na(cell_data$lineage_group), cell_data$cell_class, cell_data$lineage_group)
terminal_cell_types <- cell_data[cell_data$cell_class != "progenitor",]$cell_type
progenitor_cell_types <- cell_data[cell_data$cell_class == "progenitor",]$cell_type

all_genes <- unique(c(rownames(expr_data),
               paste0("C. elegans ", gff_list_mrna[["elegans"]]$gene_long_name),
               paste0("C. elegans ", gff_list_mrna[["elegans"]]$wbgene),
               paste0("C. briggsae ", gff_list_mrna[["briggsae"]]$gene_long_name),
               paste0("C. briggsae ", gff_list_mrna[["briggsae"]]$gene_short_name),
               paste0("C. briggsae ", gff_list_mrna[["briggsae"]]$wbgene)))
all_cell_types <- unique(c(rev(cell_data$cell_type), CellTable[! is.na(CellTable$MergedDatasetName),"Lineage"]))

term_cols = list(cell_class_lineage_group = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                              'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                              'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1"))

pro_cols = list(cell_class_lineage_group = c(' ABalpa/ABaraa (Pharynx)' = "#b07aa1", '28 cell or earlier' = "grey30", 'ABpxa/ABarp (Epidermis)' = "#59a14f",
                                             'ABpxp (Neurons)' = "#4e79a7", 'Cx' = "darkred", 'Cxa (Epidermis)' = "#59a14f", 'Cxp/D (Muscle)' = "#76b7b2",
                                             'E (Intestine)' = "#f28e2b", 'MSxa (Pharynx)' = "#b07aa1", 'MSxp (Muscle)' = "#76b7b2", 'Other AB' = "darkred"))

cell_table_cols = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                    'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2",
                    'Mesoderm' = "#9c755f",'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
                    'ABalpa/ABaraa (Pharynx)' = "#b07aa1", '28 cell or earlier' = "grey30", 'ABpxa/ABarp (Epidermis)' = "#59a14f",
                    'ABpxp (Neurons)' = "#4e79a7", 'Cx' = "darkred", 'Cxa (Epidermis)' = "#59a14f", 'Cxp/D (Muscle)' = "#76b7b2",
                    'E (Intestine)' = "#f28e2b", 'MSxa (Pharynx)' = "#b07aa1", 'MSxp (Muscle)' = "#76b7b2", 'Other AB' = "darkred")

bar_term_cols = list(cell_class_lineage_group = c('Muscle cells' = "#117733", 'Ciliated neurons' = "#AA4499", 
                                                  'Germline cells' = "#000000", 'Glia and excretory cells' = "#882255", 'Hypodermis and seam' = "#661100", 'Mesoderm' = "#6699CC",
                                                  'Non-ciliated neurons' = "#DDCC77", 'Pharyngeal cells' = "#117733", 'Intestinal and rectal cells' = "#332288"))

bar_joint_cols = list(cell_class_lineage_group = c('ABal lineage' = "#88CCEE", 'ABar lineage' = "#CC6677", 'ABp lineage' = "#DDCC77", 'Muscle cells' = "#117733",
                                               'C lineage' = "#332288", 'Ciliated neurons' = "#AA4499", 'D lineage' = "#44AA99", 'E lineage' = "#999933",
                                               'Germline cells' = "#000000", 'Glia and excretory cells' = "#882255", 'Hypodermis and seam' = "#661100", 'Mesoderm' = "#6699CC",
                                               'MS lineage' = "#CC6677", 'Non-ciliated neurons' = "#DDCC77", 'Pharyngeal cells' = "#117733", 'Intestinal and rectal cells' = "#332288"))

# ------------------------------------------------------------------------------
# Shiny UI
# ------------------------------------------------------------------------------
ui <- fluidPage(
  shinyjs::useShinyjs(),  # Enable shinyjs for dynamic navigation
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
  ),
  
  div(style="padding: 1px 0px; width: '100%'",
      titlePanel(
        title = "", windowTitle = HTML("C. elegans and C. briggsae Expression Database")
      )
  ),
  
  navbarPage(
    title = div(HTML("<i>C. elegans</i> and <i>C. briggsae</i> Expression Database"),
                # Custom help button on the right side of the navbar
                div(
                  class = "navbar-nav",
                  style = "position: absolute; right: 20px; top: 50%; transform: translateY(-50%);",
                  actionButton(
                    "info_button", "Info",
                    icon = icon("info-circle")
                  )
                )
    ),
    id = "main_tabs",
    theme = shinytheme("spacelab"),
    
    # 1) MAIN TAB -------------------------------------------------------------
    tabPanel("Main",
             fluidPage(
               div(style = "background-color: #F0F8FF; color: #333; 
               border: 1px solid #ddd; border-radius: 8px; padding: 15px;
                             box-shadow: 0px 6px 15px rgba(0, 0, 0, 0.3), 0px 3px 5px rgba(0, 0, 0, 0.2);",
                   div(
                     p("This database contains the processed single-cell expression data from:"),
                     HTML("<b>Lineage-resolved analysis of embryonic gene expression evolution in <i>C. elegans</i> and <i>C. briggsae</i></b><br>
                Christopher R. L. Large, Rupa Khanal, LaDeana Hillier, Chau Huynh, Connor Kubo, Junhyong Kim, Robert Waterston, and John I. Murray<br><br>
                For the raw data, please see this <a href='https://github.com/livinlrg/C.elegans_C.briggsae_Embryo_Single_Cell'>Github repository</a>.")
                   )
               ),
               tags$hr(),
               fluidRow(
                 div(style = "display: flex; gap: 0px; align-items: stretch; width: 100%",
                     column(4, 
                            div(style = "background-color: DeepSkyBlue; color: #333; width: 100%;
               border: 1px solid DeepSkyBlue; border-radius: 8px; padding: 15px;
               box-shadow: 0px 6px 15px rgba(0, 0, 0, 0.3), 0px 3px 5px rgba(0, 0, 0, 0.2); margin-bottom: 20px;
                            flex: 1; display: flex; flex-direction: column; height: 100%",
                                div(
                                  tags$h4(HTML("<b>Is the expression of my gene similar between species?</b><br><br>")),
                                  tags$p(HTML("Search any gene in <i>C. elegans</i> or <i>C. briggsae</i> for its expression patterns.
                           If your gene has other orthologs in either species, they will be selected and plotted.<br>")),
                                ),
                                # Gene Search
                                fluidRow(
                                  column(12, selectizeInput(
                                    inputId = "main_gene_search",
                                    label = "Search Gene:", 
                                    choices = NULL,
                                    multiple = FALSE, 
                                    options = list(maxOptions = 10))),
                                  column(12, selectInput(
                                    inputId = "gene_subtab_choice",
                                    label = "Gene Sub-Tab:",
                                    choices = c("Terminal cell types", "Progenitors", "Orthology"), 
                                    selected = "Terminal cell types")),
                                  div(style = "margin-top: auto; text-align: left; padding: 15px",
                                      actionButton("go_gene", "Go to Gene Tab"))
                                )
                            )
                     ),
                     column(4,
                            div(style = "background-color: DeepSkyBlue; color: #333; width: 100%;
               border: 1px solid DeepSkyBlue; border-radius: 8px; padding: 15px;
               box-shadow: 0px 6px 15px rgba(0, 0, 0, 0.3), 0px 3px 5px rgba(0, 0, 0, 0.2); margin-bottom: 20px;
                            flex: 1; display: flex; flex-direction: column; height: 100%",
                                div(
                                  tags$h4(HTML("<b>What genes are expressed in my cell type?</b><br><br>")),
                                  tags$p(HTML("Search any cell type in the dataset for which genes are expressed
               in it and whether they are conserved between species.")),
                                  # Cell Search
                                  fluidRow(
                                    column(12, selectizeInput(
                                      inputId = "main_cell_search",
                                      label = "Search cell type:",
                                      choices = all_cell_types,
                                      multiple = FALSE,
                                      options = list(maxOptions = 10))),
                                    column(12, selectInput(
                                      inputId = "cell_subtab_choice",
                                      label = "Cell sub-tab:",
                                      choices = c("Cell type page", "Cell type markers"),
                                      selected = "Cell type page")),
                                    div(style = "margin-top: auto; text-align: left; padding: 15px;",
                                        actionButton("go_cell", "Go to Cell Type Tab"))
                                  )
                                )
                                
                            )
                     ),
                     column(4,
                            div(style = "background-color: DeepSkyBlue; color: #333; width: 100%;
               border: 1px solid DeepSkyBlue; border-radius: 8px; padding: 15px;
              box-shadow: 0px 6px 15px rgba(0, 0, 0, 0.3), 0px 3px 5px rgba(0, 0, 0, 0.2); margin-bottom: 20px;
                            flex: 1; display: flex; flex-direction: column; height: 100%",
                                div(tags$h4(HTML("<b>Summary data is available in the table tab including:</b><br><br>")),
                                    tags$p(HTML("&rarr; Cell type similarities between species<br>
                                &rarr; Gene expression pattern similarity between species<br>
                                             &rarr; Lists of 'marker genes' for all cell types")),
                                    
                                    # Table Search
                                    selectInput(
                                      inputId = "table_search",
                                      label = "Tables:",
                                      choices = c("Cell type table", "Gene table", "Marker table"),
                                      multiple = FALSE,
                                      selected = "Cell type table"),
                                    div(style = "margin-top: auto; text-align: left; padding: 15px; width: 100%;",
                                        actionButton("go_table", "Go to Table Tab")
                                    )
                                )
                            )
                     ),
                 ),
               ),
               tags$hr(),
               
               fluidRow(
                 div(style = "display: flex; gap: 0px; align-items: stretch; width: 100%",
                     column(6,
                            div(style = "background-color: #F0F8FF; color: #333; width: 100%;
               border: 1px solid #ddd; border-radius: 8px; padding: 15px;
              box-shadow: 0px 6px 15px rgba(0, 0, 0, 0.3), 0px 3px 5px rgba(0, 0, 0, 0.2); margin-bottom: 20px;
                            flex: 1; display: flex; flex-direction: column; height: 100%",
                                div(tags$h4(HTML("<b>Gene Vignette</b><br>")),
                                    tags$p(HTML("Learn more about how to explore your favorite gene:"))),
                                div(style = "text-align: center;",
                                    tags$img(src = "gene_vignette.png", width = "100%", height = "auto",
                                             style = "max-width: 400px; box-shadow: 0px 4px 10px rgba(0,0,0,0.2);")
                                ),
                                div(style = "margin-top: auto; text-align: left; padding: 15px;",
                                    actionButton("gene_vignette_button", "Gene Vignette")
                                )
                            )
                     ),
                     column(6,
                            div(style = "background-color: #F0F8FF; color: #333; width: 100%;
               border: 1px solid #ddd; border-radius: 8px; padding: 15px;
              box-shadow: 0px 6px 15px rgba(0, 0, 0, 0.3), 0px 3px 5px rgba(0, 0, 0, 0.2); margin-bottom: 20px;
                            flex: 1; display: flex; flex-direction: column; height: 100%",
                                div(tags$h4(HTML("<b>Cell Type Vignette</b><br>")),
                                    tags$p(HTML("Learn more about how to explore your favorite cell type:"))),
                                div(style = "text-align: center;",
                                    tags$img(src = "cell_vignette.png", width = "100%", height = "auto",
                                             style = "max-width: 400px; box-shadow: 0px 4px 10px rgba(0,0,0,0.2);")
                                ),
                                div(style = "margin-top: auto; text-align: left; padding: 15px;",
                                    actionButton("cell_type_vignette_button", "Cell Type Vignette")
                                )
                            )
                     )
                 )
               ),
               tags$hr()
             ),
    ),
    
    # 3) ALL GENES TAB -------------------------------------------------------------
    # Gene Tab (Shared Sidebar for Terminal, Progenitor, Orthology)
    tabPanel("Genes",
             sidebarLayout(
               # Shared Sidebar for Terminal, Progenitor, Orthology
               sidebarPanel(
                 h4(HTML("Select <i>C. elegans</i> genes:")),
                 selectizeInput(
                   inputId = "cel_genes",
                   label = NULL,
                   choices = NULL,
                   multiple = TRUE,
                   options = list(maxItems = 10)
                 ),
                 
                 h4(HTML("Select <i>C. briggsae</i> genes:")),
                 selectizeInput(
                   inputId = "cbr_genes",
                   label = NULL,
                   choices = NULL,
                   multiple = TRUE,
                   options = list(maxItems = 10)
                 ),
                 
                 tags$hr(),
                 h4("Ortholog search:"),
                 selectizeInput(
                   inputId = "gene_search",
                   label = "Search Gene:",
                   choices = NULL,
                   multiple = FALSE,
                   options = list(maxOptions = 10)
                 ),
                 actionButton("search_gene", "Select gene"),
                 tags$hr(),
                 div(style = "display: flex; gap: 20px; align-items: center;",
                     div(
                       radioButtons("scale_selector", "Plotting scale",
                                    choices = list(
                                      "Linear" = "identity",
                                      "Log2" = "log2"
                                    ))
                     ),
                     div( 
                       radioButtons("shared_scale", "Shared scale",
                                    choices = list(
                                      "Individual" = "individual",
                                      "Shared" = "shared"
                                    ))
                     )
                 ),
               ),
               
               # Main Panel with Nested Tabs
               mainPanel(
                 tabsetPanel(id = "gene_tabs",
                             
                             # Terminal cell types
                             tabPanel("Terminal cell types",
                                      tabsetPanel(
                                        id = "terminal_subtabs",
                                        tabPanel("TPM vs. TPM plot", shinycssloaders::withSpinner(uiOutput("cm_terminal_TPM_TPM"))),
                                        tabPanel("TPM vs. percent in plot", shinycssloaders::withSpinner(uiOutput("cm_terminal_TPM_percent"))),
                                        tabPanel("Expression bar plot", shinycssloaders::withSpinner(uiOutput("cm_terminal_expression_bar")))
                                      )
                             ),
                             
                             # Progenitor cell types
                             tabPanel("Progenitors",
                                      tabsetPanel(
                                        id = "progenitor_subtabs",
                                        tabPanel("Circular lineage plot", shinycssloaders::withSpinner(uiOutput("circle_pattern_grid"))),
                                        tabPanel("TPM vs. TPM plot", shinycssloaders::withSpinner(uiOutput("cm_progenitor_TPM_TPM"))),
                                        tabPanel("TPM vs. percent in plot", shinycssloaders::withSpinner(uiOutput("cm_progenitor_TPM_percent"))),
                                        tabPanel("Expression bar plot", shinycssloaders::withSpinner(uiOutput("cm_progenitor_expression_bar")))
                                      )
                             ),
                             
                             # Orthology relationship
                             tabPanel("Orthology",
                                      fluidRow(
                                        column(12, shinycssloaders::withSpinner(uiOutput("orthology_plot")))
                                      )
                             )
                 )
               )
             )
    ),
    
    # 4) CELL TAB -------------------------------------------------------------
    tabPanel("Cell type",
             tabsetPanel(id = "cell_tabs",
                         tabPanel("Cell type page",
                                  br(),
                                  fluidRow(
                                    column(4, selectInput("cell_page_cell_class",
                                                          "Select Cell Class:",
                                                          choices = cell_data %>% arrange(cell_class) %>% pull(cell_class_lineage_group) %>% unique())),
                                    column(4, uiOutput("cell_page_cell_type")),
                                  ),
                                  uiOutput("cell_summary_box"),
                                  div(id = "cell_page_header"),
                                  shinycssloaders::withSpinner(
                                    plotOutput("cell_page", height = "800px", width = "800px")
                                  )
                         ),
                         tabPanel("Cell type markers",
                                  br(),
                                  fluidRow(
                                    column(3, selectInput("cell_marker_cell_class", "Select Cell Class:",
                                                          choices = cell_data %>% arrange(cell_class) %>% pull(cell_class_lineage_group) %>% unique())),
                                    column(3, uiOutput("cell_marker_cell_type")),
                                    column(5, uiOutput("cell_marker_sort_by")),
                                  ),
                                  column(6, style = "min-width: 475px; max-width: 500px; display: flex;",
                                         div(style = "border: 2px solid #eeeeee; border-radius: 8px; margin: 10px; padding: 10px; background-color: #fff;
                                             display: flex; flex-direction: column; flex-grow: 1;",
                                             h4(HTML("<i>C. elegans</i>")),
                                             sliderInput("num_markers_slider_cel", HTML("Top cell type markers for <i>C. elegans:</i>"),
                                                         min = 1, max = 50, value = 20, step = 1),
                                             div(style = "flex-grow: 1;",  # Ensures the plot fills available space
                                                 shinycssloaders::withSpinner(
                                                   plotOutput("cell_marker_cel", height = "auto")
                                                 )
                                             )
                                         )
                                  ),
                                  column(6, style = "min-width: 475px; max-width: 500px;display: flex;",
                                         div(style = "border: 2px solid #eeeeee; border-radius: 8px; margin: 10px; padding: 10px; background-color: #fff;
                                             display: flex; flex-direction: column; flex-grow: 1;",
                                             h4(HTML("<i>C. briggsae</i>")),
                                             sliderInput("num_markers_slider_cbr", HTML("Top cell type markers for <i>C. briggsae:</i>"),
                                                         min = 1, max = 50, value = 20, step = 1),
                                             div(style = "flex-grow: 1;",
                                                 shinycssloaders::withSpinner(
                                                   plotOutput("cell_marker_cbr", height = "auto")
                                                 )
                                             )
                                         )
                                  )
                         )
             )
    ),
    
    # 5) TABLE TAB -------------------------------------------------------------
    tabPanel("Tables",
             tabsetPanel(id = "table_tabs",
                         tabPanel("Cell table",
                                  shinycssloaders::withSpinner(
                                    reactableOutput("cell_table")
                                  )
                         ),
                         tabPanel("Gene table",
                                  shinycssloaders::withSpinner(
                                    reactableOutput("gene_table")
                                  )
                         ),
                         tabPanel("Marker table",
                                  br(),
                                  fluidRow(
                                    column(4, selectInput("marker_species", "Select Species:", choices = c("C. elegans", "C. briggsae"))),
                                    column(4, uiOutput("marker_cell_class")),
                                    column(4, uiOutput("marker_cell_type"))
                                  ),
                                  shinycssloaders::withSpinner(
                                    reactableOutput("marker_table")
                                  )
                         )
             )
    )
  )
)

# ------------------------------------------------------------------------------
# Shiny Server
# ------------------------------------------------------------------------------
server <- function(input, output, session) {
  
  og <- reactiveVal(character())
  
  observeEvent(input$info_button, {
    tab <- input$main_tabs
    
    showModal(modalDialog(
      title = "Information",
      info_text_lookup(tab),
      easyClose = TRUE,
      footer = modalButton("Close"),
      size = "l"
    ))
  })
  
  # gene_subtab_choice
  # terminal_subtabs
  # progenitor_subtabs
  # gene_subtab_choice
  
  # A simple observer that runs once on startup:
  observe({
    # 1) Update main_gene_search
    updateSelectizeInput(
      session, "main_gene_search",
      choices = all_genes,
      server = TRUE
    )
    
    # 1) Update gene_search
    updateSelectizeInput(
      session, "gene_search",
      choices = all_genes,
      server = TRUE
    )
    
    # 2) Update cel_genes (only the C. elegans genes)
    updateSelectizeInput(
      session, "cel_genes",
      choices = rownames(expr_data)[grepl("C. elegans", rownames(expr_data))],
      server = TRUE
    )
    
    # 3) Update cbr_genes (only the C. briggsae genes)
    updateSelectizeInput(
      session, "cbr_genes",
      choices = rownames(expr_data)[grepl("C. briggsae", rownames(expr_data))],
      server = TRUE
    )
  })
  
  observeEvent(input$gene_vignette_button, {
    showModal(modalDialog(
      title = "Gene Vignette",
      gene_vignette_text,
      easyClose = TRUE,
      footer = modalButton("Close"),
      size = "l"
    ))
  })
  
  observeEvent(input$cell_type_vignette_button, {
    showModal(modalDialog(
      title = "Cell Type Vignette",
      cell_type_vignette_text,
      easyClose = TRUE,
      footer = modalButton("Close"),
      size = "l"
    ))
  })
  
  # --- ACTION: Go to Gene Tab (subtab depends on user choice) -----------------
  observeEvent(input$go_gene, {
    updateNavbarPage(session, "main_tabs", selected = "Genes")
    updateTabsetPanel(session, "gene_tabs", selected = input$gene_subtab_choice)
    
    req(input$main_gene_search)
    
    search_and_select_gene(input, output, session, input$main_gene_search, gff_list_mrna, expr_data, og_data)
  })
  
  # --- ACTION: Search gene in gene tab (subtab depends on user choice) -----------------
  # Triggered when the user selects a gene in the ortholog search
  observeEvent(input$search_gene, {
    search_and_select_gene(input, output, session, input$gene_search, gff_list_mrna, expr_data, og_data)
  })
  
  # search and update the orthogroup
  observeEvent(list(input$cel_genes, input$cbr_genes), {
    temp_og <- c()
    if(length(input$cel_genes) > 0) {
      temp_gene_selection <- strsplit(input$cel_genes, " ")[[1]][3]
      temp_og <- c(temp_og, gff_list_mrna[["elegans"]][temp_gene_selection, "OG"])
    }
    
    if(length(input$cbr_genes) > 0) {
      temp_gene_selection <- strsplit(input$cbr_genes, " ")[[1]][3]
      temp_og <- c(temp_og, gff_list_mrna[["briggsae"]][temp_gene_selection, "OG"])
    }
    
    og(unique(temp_og))
  })
  
  # Update the scale selection when the tab changes
  observeEvent(input$terminal_subtabs, {
    if (input$terminal_subtabs == "TPM vs. TPM plot") {
      updateRadioButtons(session, "scale_selector", selected = "log2")
    } else if (input$terminal_subtabs == "TPM vs. percent in plot") {
      updateRadioButtons(session, "scale_selector", selected = "identity")
    } else if (input$terminal_subtabs == "Expression bar plot") {
      updateRadioButtons(session, "scale_selector", selected = "identity")
    }
  })
  
  # Update the scale selection when the tab changes
  observeEvent(input$progenitor_subtabs, {
    if (input$progenitor_subtabs == "Circular lineage plot") {
      updateRadioButtons(session, "scale_selector", selected = "log2")
    } else if (input$progenitor_subtabs == "TPM vs. TPM plot") {
      updateRadioButtons(session, "scale_selector", selected = "log2")
    } else if (input$progenitor_subtabs == "TPM vs. percent in plot") {
      updateRadioButtons(session, "scale_selector", selected = "identity")
    } else if (input$progenitor_subtabs == "Expression bar plot") {
      updateRadioButtons(session, "scale_selector", selected = "identity")
    }
  })
  
  # --- ACTION: Go to Cell Tab (subtab depends on user choice) -----------------
  observeEvent(input$go_cell, {
    updateNavbarPage(session, "main_tabs", selected = "Cell type") # Switch to main Cell tab
    updateTabsetPanel(session, "cell_tabs", selected = input$cell_subtab_choice) # Switch to main Cell tab
    
    # Determine the cell class for the selected cell type
    selected_cell_type <- search_and_select_cell(input, output, session, input$main_cell_search, CellTable)
    selected_cell_class <- cell_data %>% filter(cell_type == selected_cell_type) %>% pull(cell_class_lineage_group) %>% unique()

    # Update the cell type page inputs
    updateSelectInput(session, "cell_page_cell_class", selected = selected_cell_class)
    updateSelectInput(session, "cell_page_cell_type", 
                      choices = cell_data %>% filter(cell_class_lineage_group == selected_cell_class) %>% pull(cell_type),
                      selected = selected_cell_type)
    
    # Update the cell marker page inputs
    updateSelectInput(session, "cell_marker_cell_class", selected = selected_cell_class)
    updateSelectInput(session, "cell_marker_cell_type", 
                      choices = cell_data[cell_data$cell_class_lineage_group == selected_cell_class, "cell_type"],
                      selected = selected_cell_type)
  })
  
  # --- ACTION: Go to Table Tab (subtab depends on user choice) -----------------
  observeEvent(input$go_table, {
    updateNavbarPage(session, "main_tabs", selected = "Tables")
    updateTabsetPanel(session, "table_tabs", selected = input$table_search)
  })
  
  #### Terminal ####
  # TPM vs. TPM plot grid rendering
  render_gene_grid(output, input, session, expr_data[, terminal_cell_types], cell_data, term_cols, "term", grid_id = "cm_terminal_TPM_TPM")
  
  # TPM vs. percent in plot grid rendering
  render_percent_grid(output, input, session, expr_data[, terminal_cell_types], perc_data[, terminal_cell_types], cell_data, term_cols, "term", grid_id = "cm_terminal_TPM_percent")
  
  # Bar plots with just the terminal cell types
  render_expression_bar(output, input, session, expr_data[, terminal_cell_types], cell_data, bar_term_cols, "term", grid_id = "cm_terminal_expression_bar")
  
  #### Progenitor ####
  render_circle_pattern_grid(output, input, session, CellTable, lo, expr_data, grid_id = "circle_pattern_grid")
  
  # TPM vs. TPM plot grid rendering
  render_gene_grid(output, input, session, expr_data[, progenitor_cell_types], cell_data, pro_cols, "pro", grid_id = "cm_progenitor_TPM_TPM")
  
  # TPM vs. percent in plot grid rendering
  render_percent_grid(output, input, session, expr_data[, progenitor_cell_types], perc_data[, progenitor_cell_types], cell_data, pro_cols, "pro", grid_id = "cm_progenitor_TPM_percent")
  
  # Bar plots with all of the cell types
  render_expression_bar(output, input, session, expr_data, cell_data, bar_joint_cols, "pro", grid_id = "cm_progenitor_expression_bar")
  
  #### Orthology ####
  render_tree_plots(output, input, session, og, og_data, synteny, gene_data)
  
  #### Cell ####
  
  # Cell Page
  output$cell_page_cell_class <- renderUI({
    selectInput("cell_page_cell_class", "Select Cell Class:",
                choices = cell_data$cell_class_lineage_group %>% unique() %>% sort(),
                selected = "Ciliated neurons")
  })
  
  # Update item choices based on selected group
  output$cell_page_cell_type <- renderUI({
    req(input$cell_page_cell_class)
    choices_temp <- cell_data[cell_data$cell_class_lineage_group == input$cell_page_cell_class, "cell_type"]
    selectInput("cell_page_cell_type", "Select Cell Type:",
                choices = choices_temp,
                selected = search_and_select_cell(input, output, session, input$main_cell_search, CellTable))
  })
  
  observeEvent(input$cell_page_cell_type, {
    req(input$cell_page_cell_type)
    
    removeUI(selector = "#cell_summary_box", immediate = TRUE)
    
    # Insert a UI box dynamically above render_cell_page
    insertUI(
      selector = "#cell_page_header",
      where = "beforeBegin",
      ui = div(
        id = "cell_summary_box",  # Assign ID for removal on next update
        style = "background-color: #F8F9FA; border: 1px solid #ddd; border-radius: 8px; padding: 15px; 
               margin-bottom: 15px; box-shadow: 0px 4px 10px rgba(0, 0, 0, 0.2);",
        tags$p(paste("You selected:", input$cell_page_cell_type)),
        tags$p(paste("Associated lineages:", paste(CellTable[which(CellTable$MergedDatasetName == input$cell_page_cell_type), "Lineage"], collapse = ", "))),
        tags$p(paste("Terminal cells:", paste(CellTable[which(CellTable$MergedDatasetName == input$cell_page_cell_type & CellTable$Cell != "progenitor"), "Cell"], collapse = ", ")))
      )
    )
    
    render_cell_page(output, input, session, input$cell_page_cell_type, cell_data, markers, markers_enrich_WormCat.1, expr_data, grid_id = "cell_page")
  })
  
  # Cell markers
  output$cell_marker_cell_class <- renderUI({
    selectInput("cell_marker_cell_class", "Select Cell Class:",
                choices = cell_data$cell_class_lineage_group %>% unique() %>% sort(),
                selected = "Ciliated neurons")
  })
  
  output$cell_marker_cell_type <- renderUI({
    req(input$cell_marker_cell_class)
    choices_temp <- cell_data[cell_data$cell_class_lineage_group == input$cell_marker_cell_class, "cell_type"]
    selectInput("cell_marker_cell_type", "Select Cell Type:",
                choices = choices_temp,
                selected = search_and_select_cell(input, output, session, input$main_cell_search, CellTable))
  })
  
  output$cell_marker_sort_by <- renderUI({
    req(input$cell_marker_cell_class)
    selectInput("cell_marker_sort_by", "Sort by:",
                choices = c("Specificity x Exp. (log2FC * log2TPM)", "Specificity (log2FC)", "Specificity (p-value)", "Expression"),
                selected = "Specificity x Exp. (log2FC * log2TPM)")
  })
  
  observeEvent(input$cell_marker_cell_type, {
    req(input$cell_marker_cell_type)
    
    # For C. elegans
    max_markers_cel <- nrow(markers[["C.elegans"]][markers[["C.elegans"]]$cell_type == input$cell_marker_cell_type, ])
    max_markers_cel <- max(1, max_markers_cel) # Ensure at least 1 marker is selectable
    updateSliderInput(session, "num_markers_slider_cel", max = max_markers_cel, value = min(20, max_markers_cel))
    
    # For C. briggsae
    max_markers_cbr <- nrow(markers[["C.briggsae"]][markers[["C.briggsae"]]$cell_type == input$cell_marker_cell_type, ])
    max_markers_cbr <- max(1, max_markers_cbr) # Ensure at least 1 marker is selectable
    updateSliderInput(session, "num_markers_slider_cbr", max = max_markers_cbr, value = min(20, max_markers_cbr))
  })
  
  # Render C. elegans plot
  observeEvent({
    input$marker_plot_species
    input$cell_marker_cell_type
    input$num_markers_slider_cel
    input$cell_marker_sort_by
  }, {
    req(input$cell_marker_cell_type)
    render_cell_markers(output, input, session, 
                        cur_cell_type = input$cell_marker_cell_type, 
                        cur_species = "C.elegans", 
                        markers = markers[["C.elegans"]][markers[["C.elegans"]]$cell_type == input$cell_marker_cell_type,], 
                        num_markers = input$num_markers_slider_cel, 
                        sort_by = input$cell_marker_sort_by,
                        grid_id = "cell_marker_cel")
  })
  
  # Render C. briggsae plot
  observeEvent({
    input$marker_plot_species
    input$cell_marker_cell_type
    input$num_markers_slider_cbr
    input$cell_marker_sort_by
  }, {
    req(input$cell_marker_cell_type)
    render_cell_markers(output, input, session, 
                        cur_cell_type = input$cell_marker_cell_type, 
                        cur_species = "C.briggsae", 
                        markers = markers[["C.briggsae"]][markers[["C.briggsae"]]$cell_type == input$cell_marker_cell_type,], 
                        num_markers = input$num_markers_slider_cbr, 
                        sort_by = input$cell_marker_sort_by,
                        grid_id = "cell_marker_cbr")
  })
  
  #### Tables ####
  
  # Cell Table
  output$cell_table <- renderReactable({
    build_cell_table(cell_data, cols = cell_table_cols)
  })
  
  # Gene table
  output$gene_table <- renderReactable({
    build_gene_table(gene_data, input, output, session)
  })
  
  # Marker Table
  output$marker_cell_class <- renderUI({
    selectInput("marker_cell_class", "Select Cell Class:",
                choices = cell_data %>% arrange(cell_class) %>% pull(cell_class_lineage_group) %>% unique(),
                selected = "Ciliated neurons")
  })
  
  output$marker_cell_type <- renderUI({
    choices_temp <- cell_data[cell_data$cell_class_lineage_group == input$marker_cell_class, "cell_type"]
    selectInput("marker_cell_type", "Select Cell Type:",
                choices = choices_temp,
                selected = "ADE")
  })
  
  output$marker_table <- renderReactable({
    req(input$marker_species, input$marker_cell_class, input$marker_cell_type)
    cur_species <- gsub(" ", "", input$marker_species)
    marker_data <- markers[[cur_species]][markers[[cur_species]]$cell_type == input$marker_cell_type,]
    build_marker_table(marker_data, input, output, session)
  })
  
}

# ------------------------------------------------------------------------------
# Run the application
# ------------------------------------------------------------------------------
# Run the Shiny app
shinyApp(ui = ui, server = server)