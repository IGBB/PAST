library(PAST)
library(shiny)
library(shinythemes)
library(shinyjs) 
library(DT)
library(dplyr)
library(rtracklayer)
library(ggplot2)
library(plotly)

# Define UI for application that draws a histogram
ui <- bootstrapPage(
  useShinyjs(),
  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #95a5a6 !important} .action-button {margin: 5px;}"))
  ),
  navbarPage(
    theme = shinytheme("flatly"), 
    collapsible = TRUE,
    
    # Application title
    "PAST",
    
    # Sidebar with a slider input for number of bins 
    tabPanel("GWAS Input",
             sidebarLayout(
               sidebarPanel(
                 width = 4,
                 h3("Input Type"),
                 tags$p("Select the type of GWAS data you have. \"Single\" has all the data in a single file, \"Two\" has separate associations and effects data with one line per marker in each file, and \"TASSEL\" has separate associations data with one line per markers and effects data with two lines per marker."),
                 selectInput("input_type", NULL,
                             choices = c("single", "two", "TASSEL")),
                 
                 uiOutput("file_header"),
                 uiOutput("file_text"),
                 uiOutput("file_association_upload"),
                 uiOutput("association_error_text"),
                 uiOutput("file_effects_upload"),
                 uiOutput("effects_error_text"),
                 
                 hr(),
                 uiOutput("gwas_column_names")
               ),
               
               # Show a plot of the generated distribution
               mainPanel(
                 DTOutput("gwas_data")
               )
             )
    ),
    
    tabPanel("LD Input",
             sidebarLayout(
               sidebarPanel(
                 width = 4,
                 h3("File"),
                 tags$p("Upload the linkage disequilibrium file. This file can be compressed with BZIP, GZIP, or XZ compression."),
                 fileInput("LD_file", "Linkage Disequilibrium File"),
                 uiOutput("LD_error_text"),
                 hr(),
                 h3("Column Names"),
                 tags$p("Load your data and select the appropriate column for each name."),
                 uiOutput("LD_column_names")
               ),
               
               # Show a plot of the generated distribution
               mainPanel(
                 DTOutput("LD_data")
               )
             )
    ),
    
    tabPanel("Annotations Input",
             sidebarLayout(
               sidebarPanel(
                 width = 4,
                 h3("File"),
                 tags$p("Upload the uncompressed annotation file in GFF3 format."),
                 fileInput("genes_file", "Genes File"),
                 uiOutput("genes_error_text"),
                 hr(),
                 h3("Options"),
                 tags$p("PAST searches for genes within 1kb of SNPs by default, but this number can be changed."),
                 numericInput("window", "Window Size", value = 1000),
                 tags$p(tagList("PAST uses R", tags$sup("2"), " to determine linkage between SNPs.")),
                 sliderInput(
                   "r_squared_cutoff",
                   NULL,
                   min = 0,
                   max = 1,
                   value = 0.8,
                   step = 0.05
                 ),
                 tags$p("PAST filters for a specific type of annotation (genes by default), but you may wish to use another feature."),
                 selectInput("annotation_filter_type", "Annotation Type", selected = "gene",
                             choices = c("gene", "CDS", "mRNA", "exon"))
               ),
               
               # Show a plot of the generated distribution
               mainPanel(
                 DTOutput("genes_data")
               )
             )
    ),
    
    tabPanel("Pathways Input",
             sidebarLayout(
               sidebarPanel(
                 width = 4,
                 h3("File"),
                 tags$p("Upload the pathways file in TSV format. This file can be compressed with BZIP, GZIP, or XZ compression."),
                 fileInput("pathways_file", "Pathways File"),
                 uiOutput("pathways_error_text"),
                 hr(),
                 
                 h3("Gene Cutoff"),
                 tags$p("PAST discards pathways with less than 5 genes by default."),
                 numericInput("gene_cutoff", NULL, value = 5),
                 
                 h3("Sampling"),
                 tags$p("PAST determines pathway signficance by creating random distributions of gene effects and comparing the actual effects to randomly sampled effects."),
                 numericInput("sample", NULL, value = 1000),
               ),
               
               # Show a plot of the generated distribution
               mainPanel(
                 DTOutput("pathways_data")
               )
             )
    ),
    
    tabPanel("Analysis",
             sidebarLayout(
               sidebarPanel(
                 width = 4,
                 h3("Analysis Title"),
                 tags$p("Provide a title for your analysis. This title is used for the downloaded results."),
                 textInput("title", NULL, value = "New Analysis"),
                 
                 h3("Number of Cores"),
                 tags$p("Select the number of cores to be used in run your analysis in parallel."),
                 selectInput("num_cores", NULL, choices=seq(1, parallel::detectCores(TRUE))),
                 
                 h3("Mode"),
                 tags$p("Select the type of analysis you wish to run. \"Increasing\" searches for pathways associated with an increase in the measured trait, and \"decreasing\" searches for pathways associated with a decrease in the measured trait."),
                 selectInput("mode", NULL,
                             choices = c("increasing", "decreasing")),
                 
                 actionButton("begin_genes_analysis", "Begin SNP-gene Assignment"),
                 actionButton("begin_pathways_analysis", "Begin Pathways Analysis"),
                 hr(),
                 h3("Filtering"),
                 tags$p("Choose your method of filtering the results."),
                 selectInput(
                   "filter_type",
                   "Filter Parameter:",
                   choices = c("p-value", "q-value"),
                   selected = "p-value"
                 ),
                 
                 
                 tags$p("Select the value used for filtering results."),
                 # create a slider to specify filtering level
                 sliderInput(
                   "significance_cutoff",
                   "Pathway Signficance Filter",
                   min = 0,
                   max = 1.0,
                   value = 0.05,
                   step = 0.01
                 ),
               ),
               
               # Show a plot of the generated distribution
               mainPanel(
                 tabsetPanel(type = "pills",
                             tabPanel("SNP-gene Assignments", DTOutput("snp_gene_results")),
                             tabPanel("Pathways", 
                                      fluidRow(
                                        column(6,
                                               DTOutput("pathways_results") 
                                        ),
                                        column(6,
                                               uiOutput("pathway_select"),
                                               plotlyOutput("pathway_plot"),
                                               DTOutput("pathway_genes") 
                                        )
                                      )
                             )
                 )
               )
             )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  options(shiny.maxRequestSize = 600 * 1024 ^ 2)
  
  ## GWAS Data
  
  uploaded_files = reactiveValues()
  
  output$file_header <- renderUI({
    if (input$input_type == "single") {
      return(h3("File"))
    } else {
      return(h3("Files"))
    }
  })
  
  output$file_text <- renderUI({
    if (input$input_type == "single") {
      return(tags$p("Upload a file containing a single trait. These files can be compressed with BZIP, GZIP, or XZ compression."))
    } else {
      return(tags$p("Upload files containing a single trait. These files can be compressed with BZIP, GZIP, or XZ compression."))
    }
  })
  
  output$file_association_upload <- renderUI({
    if (input$input_type == "single") {
      return(fileInput("association_file", "GWAS Data File"))
    } else {
      return(fileInput("association_file", "Association File"))
    }
  })
  
  output$file_effects_upload <- renderUI({
    if (input$input_type == "single") {
      return(NULL)
    } else {
      return(fileInput("effects_file", "Effects File"))
    }
  })
  
  observeEvent(input$association_file, {
    uploaded_files$association_file=input$association_file$datapath
    uploaded_files$effects_file=NULL
    reset("effects_file")
  })
  
  observeEvent(input$effects_file, {
    uploaded_files$effects_file=input$effects_file$datapath
  })
  
  observeEvent(input$input_type, {
    uploaded_files$association_file = NULL
    uploaded_files$effects_file = NULL
    output$association_error_text <- renderUI(NULL)
    output$effects_error_text <- renderUI(NULL)
    reset('association_file')
    reset('effects_file')
  })
  
  association_columns = reactive({
    if (is.null(uploaded_files$association_file))
      return(NULL)
    
    tryCatch({
      output$association_error_text <- renderUI(NULL)
      return(names(read.table(uploaded_files$association_file, sep = "\t", header = T, nrows = 2)))
    }, error = function(e) {
      output$association_error_text <- renderUI({
        span(h4(paste0("Error reading file: ", input$association_file$name)),
             p(e$message), 
             style="color:red"
        )
      })
      reset('association_file')
      uploaded_files$association_file = NULL
      return(NULL)
    })
  })
  
  effects_columns = reactive({
    if (is.null(uploaded_files$effects_file))
      return(NULL)
    
    tryCatch({
      output$effects_error_text <- renderUI(NULL)
      return(names(read.table(uploaded_files$effects_file, sep = "\t", header = T, nrows = 2)))
    }, error = function(e) {
      output$effects_error_text <- renderUI({
        span(h4(paste0("Error reading file: ", input$effects_file$name)),
             p(e$message), 
             style="color:red"
        )
      })
      reset('effects_file')
      uploaded_files$effects_file = NULL
    })
    
  })
  
  output$gwas_column_names <- renderUI({
    if (input$input_type == "single") {
      return(
        div(
          h3("Column Names"),
          tags$p("Load your data and select the appropriate column for each name."),
          selectInput("association_trait",
                      tags$b("Trait Name"),
                      selected = NULL,
                      choices = association_columns()),
          selectInput("association_marker",
                      tags$b("Marker Column Name"),
                      selected = NULL,
                      choices = association_columns()),
          selectInput("association_locus",
                      tags$b("Locus/Chromosome Column Name"),
                      choices = association_columns()),
          selectInput("association_site",
                      tags$b("Site/Position Column Name"),
                      choices = association_columns()),
          selectInput("association_p",
                      tags$b("P-value Column Name"),
                      choices = association_columns()),
          selectInput("association_marker_R2",
                      tags$b(tagList("Marker R", tags$sup("2"), " Column Name")),
                      choices = association_columns()),
          selectInput("association_effect",
                      tags$b("Effect Column Name"),
                      choices = association_columns())
        )
      )
    } else {
      return(
        div(
          h3("Association Column Names"),
          tags$p("Load your data and select the appropriate column for each name."),
          selectInput("association_trait",
                      tags$b("Trait Name"),
                      selected = NULL,
                      choices = association_columns()),
          selectInput("association_marker",
                      tags$b("Marker Column Name"),
                      selected = NULL,
                      choices = association_columns()),
          selectInput("association_locus",
                      tags$b("Locus/Chromosome Column Name"),
                      choices = association_columns()),
          selectInput("association_site",
                      tags$b("Site/Position Column Name"),
                      choices = association_columns()),
          selectInput("association_p",
                      tags$b("P-value Column Name"),
                      choices = association_columns()),
          selectInput("association_marker_R2",
                      tags$b(tagList("Marker R", tags$sup("2"), " Column Name")),
                      choices = association_columns()),
          hr(),
          h3("Effects Column Names"),
          tags$p("Load your data and select the appropriate column for each name."),
          selectInput("effects_trait",
                      tags$b("Trait Name"),
                      selected = NULL,
                      choices = effects_columns()),
          selectInput("effects_marker",
                      tags$b("Marker Column Name"),
                      selected = NULL,
                      choices = effects_columns()),
          selectInput("effects_locus",
                      tags$b("Locus/Chromosome Column Name"),
                      choices = effects_columns()),
          selectInput("effects_site",
                      tags$b("Site/Position Column Name"),
                      choices = effects_columns()),
          selectInput("effects_effect",
                      tags$b("Effect Column Name"),
                      choices = effects_columns())
        )
      )
    }
    
  })
  
  get_gwas_data <- reactive({
    update_progress <- function(message = NULL, value = NULL, detail = NULL) {
      progress$set(message = message, value = value, detail = detail)
    }
    if (input$input_type == "single") {
      if (is.null(uploaded_files$association_file))
        return(NULL)
      input_files = c(uploaded_files$association_file)
    } else {
      effects_file <- uploaded_files$effects_file
      if (is.null(uploaded_files$association_file) | is.null(uploaded_files$effects_file))
        return(NULL)
      input_files = c(uploaded_files$association_file, uploaded_files$effects_file)
    }
    
    if (input$input_type == "single") {
      column_set = list(input$association_trait,
                        input$association_marker,
                        input$association_locus,
                        input$association_site,
                        input$association_p,
                        input$association_marker_R2,
                        input$association_effect)
      
      if (sum(vapply(column_set, is.null, TRUE) != 0) | length(column_set %>% unique) != 7)
        return(NULL)
      
      
      tryCatch({
        progress <- shiny::Progress$new(min = 0, max = 1)
        data <- load_GWAS_data(
          input_files,
          input$association_trait,
          input$association_marker,
          input$association_locus,
          input$association_site,
          input$association_p,
          input$association_marker_R2,
          input$association_effect,
          input$input_type,
          update_progress = update_progress)
        progress$close()
        return(data)
      }, error = function(e) {
        showModal(modalDialog(
          title = "Error in parsing GWAS data",
          e$message,
        ))
        return(NULL)
      })
      
      
    } else {
      association_column_set = list(input$association_trait,
                                    input$association_marker,
                                    input$association_locus,
                                    input$association_site,
                                    input$association_p,
                                    input$association_marker_R2)
      
      effects_column_set = list(input$effects_effect,
                                input$effects_trait,
                                input$effects_marker,
                                input$effects_locus,
                                input$effects_site)
      
      if (sum(vapply(association_column_set, is.null, TRUE) != 0) | 
          sum(vapply(effects_column_set, is.null, TRUE) != 0) | 
          length(association_column_set %>% unique) != 6 | 
          length(effects_column_set %>% unique) != 5)
        return(NULL)
      
      tryCatch({
        progress <- shiny::Progress$new(min = 0, max = 1)
        data <- load_GWAS_data(
          input_files,
          input$association_trait,
          input$association_marker,
          input$association_locus,
          input$association_site,
          input$association_p,
          input$association_marker_R2,
          input$effects_effect,
          input$input_type,
          input$effects_trait,
          input$effects_marker,
          input$effects_locus,
          input$effects_site,
          update_progress = update_progress)
        progress$close()
        return(data)
      }, error = function(e) {
        showModal(modalDialog(
          title = "Error in parsing GWAS data",
          e$message,
        ))
        return(NULL)
      })
      
    }
  })
  
  output$gwas_data <- renderDT({
    if(!(is.null(get_gwas_data()))) {
      datatable(get_gwas_data(), rownames = FALSE)
    } else {
      return(NULL)
    }
  })
  
  ## Linkage Disequilibrium Data
  
  observeEvent(input$LD_file, {
    uploaded_files$LD_file=input$LD_file$datapath
  })
  
  LD_columns = reactive({
    if (is.null(uploaded_files$LD_file))
      return(NULL)
    
    tryCatch({
      output$LD_error_text <- renderUI(NULL)
      return(names(read.table(uploaded_files$LD_file, sep = "\t", header = T, nrows = 2)))
    }, error = function(e) {
      message("ERROR in LD file")
      output$LD_error_text <- renderUI({
        span(h4(paste0("Error reading file: ", input$LD_file$name)),
             p(e$message), 
             style="color:red"
        )
      })
      reset('LD_file')
      uploaded_files$LD_file = NULL
      return(NULL)
    })
  })
  
  output$LD_column_names <- renderUI({
    div(
      selectInput("LD_locus_1", 
                  tags$b("Locus1 Column Name"),
                  selected = NULL,
                  choices = LD_columns()),
      selectInput("LD_position_1", 
                  tags$b("Position1 Column Name"),
                  selected = NULL,
                  choices = LD_columns()),
      selectInput("LD_site_1", 
                  tags$b("Site1 Column Name"), 
                  choices = LD_columns()),
      selectInput("LD_position_2", 
                  tags$b("Position2 Column Name"),
                  selected = NULL,
                  choices = LD_columns()),
      selectInput("LD_site_2", 
                  tags$b("Site2 Column Name"), 
                  choices = LD_columns()),
      selectInput("LD_distance", 
                  tags$b("Distance Column Name"), 
                  choices = LD_columns()),
      selectInput("LD_R2", 
                  tags$b("R2 Column Name"), 
                  choices = LD_columns())
    )
  })
  
  get_LD_data <- reactive({
    update_progress <- function(message = NULL, value = NULL, detail = NULL) {
      progress$set(message = message, value = value, detail = detail)
    }
    
    if (is.null(uploaded_files$LD_file))
      return(NULL)
    
    column_set = list(input$LD_locus_1,
                      input$LD_position_1,
                      input$LD_site_1,
                      input$LD_position_2,
                      input$LD_site_2,
                      input$LD_distance,
                      input$LD_R2)
    
    if (sum(vapply(column_set, is.null, TRUE) != 0) | length(column_set %>% unique) != 7)
      return(NULL)
    LD_columns = c(
      input$LD_locus_1,
      input$LD_position_1,
      input$LD_site_1,
      input$LD_position_2,
      input$LD_site_2,
      input$LD_distance,
      input$LD_R2
    )
    LD_columns
    
    tryCatch({
      progress <- shiny::Progress$new(min = 0, max = 1)
      data = load_LD(uploaded_files$LD_file, 
                     LD_columns, 
                     update_progress = update_progress)
      progress$close()
      return(data)
    }, error = function(e) {
      reset('LD_file')
      showModal(modalDialog(
        title = "Error parsing LD file",
        e$message,
      ))
      return(NULL)
    })
    
  })
  
  output$LD_data <- renderDT({
    if (is.null(get_LD_data()))
      return(NULL)
    data = NULL
    for (chromosome in names(get_LD_data())) {
      data = rbind(data, slice_sample(get_LD_data()[[chromosome]], prop = 0.01))
    }
    datatable(data %>% dplyr::arrange(), rownames = FALSE)
  })
  
  ## Genes Data
  
  observeEvent(input$genes_file, {
    uploaded_files$genes_file=input$genes_file$datapath
  })
  
  get_genes_data <- reactive({
    if (is.null(uploaded_files$genes_file))
      return(NULL)
    
    tryCatch({
      output$genes_error_text <- renderUI(NULL)
      return(as.data.frame(rtracklayer::readGFF(uploaded_files$genes_file, 
                                                filter=list(type=input$annotation_filter_type))))
    }, error = function(e) {
      output$genes_error_text <- renderUI({
        span(h4(paste0("Error reading file: ", input$genes_file$name)),
             p(e$message), 
             style="color:red"
        )
      })
      reset('genes_file')
      uploaded_files$genes_file = NULL
      return(NULL)
    })
  })
  
  output$genes_data <- renderDT({
    if(!(is.null(get_genes_data()))) {
      datatable(get_genes_data(), rownames = FALSE)
    } else {
      return(NULL)
    }
  })
  
  ## Pathways Data
  
  observeEvent(input$pathways_file, {
    uploaded_files$pathways_file=input$pathways_file$datapath
  })
  
  get_pathways_data <- reactive({
    
    if (is.null(uploaded_files$pathways_file))
      return(NULL)
    
    tryCatch({
      output$pathways_error_text <- renderUI(NULL)
      data = read.table(uploaded_files$pathways_file, sep = "\t", header = T, quote = "")
    }, error = function(e) {
      output$pathways_error_text <- renderUI({
        span(h4(paste0("Error reading file: ", input$pathways_file$name)),
             p(e$message), 
             style="color:red"
        )
      })
      reset('pathways_file')
      uploaded_files$pathways_file = NULL
      return(NULL)
    })
    
    if (length(data) != 3) {
      output$pathways_error_text <- renderUI({
        span(h4(paste0("Error reading file: ", input$pathways_file$name)),
             p("The pathways file should have three tab-separated columns -- pathway_id, pathway_name, and gene -- with one line for every gene in the pathway."), 
             style="color:red"
        )
      })
      reset('pathways_file')
      uploaded_files$pathways_file = NULL
      return(NULL)
    }
    data  
  })
  
  output$pathways_data <- renderDT({
    if(!(is.null(get_pathways_data()))) {
      datatable(get_pathways_data(), rownames = FALSE)
    } else {
      return(NULL)
    }
  })
  
  ## Analysis
  
  results = reactiveValues()
  
  observe({
    toggleState("begin_genes_analysis", !(is.null(get_genes_data()) | is.null(get_LD_data()) | is.null(get_gwas_data())))
  })
  
  observe({
    toggleState("begin_pathways_analysis", !(is.null(get_pathways_data()) | is.null(results$gene_assignments)))
  })
  
  observeEvent(input$begin_genes_analysis, {
    update_progress <- function(message = NULL, value = NULL, detail = NULL) {
      progress$set(message = message, value = value, detail = detail)
    }
    
    if (is.null(get_genes_data()) | 
        is.null(get_LD_data()) | 
        is.null(get_gwas_data()))
      return(NULL)
    
    disable("begin_genes_analysis")
    disable("begin_pathways_analysis")
    
    progress <- shiny::Progress$new(min = 0, max = 1)
    update_progress(message = "Assigning SNPs to genes", value = 0, "0%")
    results$gene_assignments <- assign_SNPs_to_genes(get_gwas_data(),
                                                     get_LD_data(),
                                                     get_genes_data(),
                                                     c(input$annotation_filter_type),
                                                     input$window,
                                                     input$r_squared_cutoff,
                                                     as.numeric(input$num_cores),
                                                     update_progress = update_progress
    )
    progress$close()
    
    enable("begin_genes_analysis")
    enable("begin_pathways_analysis")
    output$snp_gene_results <- renderDT(datatable(results$gene_assignments, rownames = FALSE))
  })
  
  observeEvent(input$begin_pathways_analysis, {
    update_progress <- function(message = NULL, value = NULL, detail = NULL) {
      progress$set(message = message, value = value, detail = detail)
    }
    if (is.null(get_pathways_data()) | is.null(results$gene_assignments))
      return(NULL)
    progress <- shiny::Progress$new(min = 0, max = 1)
    update_progress(message = "Beginning pathways analysis", value = 0, "0%")
    results$pathways_assignments <- find_pathway_significance(results$gene_assignments,
                                                              get_pathways_data(),
                                                              input$gene_cutoff,
                                                              input$mode,
                                                              input$sample,
                                                              as.numeric(input$num_cores),
                                                              update_progress = update_progress
    )
    progress$close()
    filter_var = ifelse(input$filter_type == "p-value", "pvalue", "qvalue")
    if (is.null(results$pathways_assignments))
      return(NULL)
    if (input$filter_type == "p-value") {
      results$pathways <- dplyr::filter(results$pathways_assignments, pvalue < input$significance_cutoff)
    } else {
      results$pathways <- dplyr::filter(results$pathways_assignments, qvalue < input$significance_cutoff)
    }
  })
  
  observeEvent(c(input$filter_type, input$significance_cutoff), {
    filter_var = ifelse(input$filter_type == "p-value", "pvalue", "qvalue")
    if (is.null(results$pathways_assignments))
      return(NULL)
    if (input$filter_type == "p-value") {
      results$pathways <- dplyr::filter(results$pathways_assignments, pvalue < input$significance_cutoff)
    } else {
      results$pathways <- dplyr::filter(results$pathways_assignments, qvalue < input$significance_cutoff)
    }
  })
  
  output$pathway_select <- renderUI({
    if (is.null(results$pathways))
      return(NULL)
    selectInput("pathway_display", 
                "Select pathway to display", 
                choices = results$pathways$pathway_id %>% unique())
  })
  
  output$pathway_plot <- renderPlotly({
    if (is.null(results$pathways))
      return(NULL)
    
    if (nrow(results$pathways) == 0 | length(input$pathway_display) == 0)
      return(NULL)
    
    pathway_data = results$pathways %>% filter(pathway_id == input$pathway_display)
    title <-
      paste0(unique(as.character(pathway_data$pathway_id)), " - ",
             unique(as.character(pathway_data$pathway_name)))
    intercept <- pathway_data %>%
      dplyr::arrange(desc(.data$running_enrichment_score)) %>%
      dplyr::select(.data$rank)
    intercept <- intercept[, 1][1]
    rugplot <-
      ggplot(pathway_data,
             aes(x = rank,
                 y = running_enrichment_score,
                 label = gene_id)) +
      geom_line(stat = "identity", color = "#fe4365") +
      geom_rug(sides = "t", color = "#fe4365") +
      geom_vline(xintercept = intercept,
                 color = "red",
                 linetype = "longdash") +
      ggtitle(title) +
      geom_text(alpha = 0.0) + 
      labs(x = "Gene Rank", y = "Running Enrichment Score") +
      theme_minimal() + 
      theme(
        axis.text = element_text (color = "#a1a1a1"),
        plot.background = element_rect (color = "#a1a1a1", fill = "#222831"),
        panel.background = element_rect (color = "#a1a1a1", fill = "#222831"),
        plot.title = element_text(color = "#a1a1a1"),
        axis.title = element_text(color = "#a1a1a1"),
        panel.grid.minor = element_blank()
      )
    ggplotly(rugplot)
  })
  
  output$pathway_genes <- renderDT({
    filter_var = ifelse(input$filter_type == "p-value", "pvalue", "qvalue")
    if (is.null(results$pathways))
      return(NULL)
    
    if (nrow(results$pathways) == 0 | length(input$pathway_display) == 0)
      return(NULL)
    
    datatable(results$pathways %>% 
                filter(pathway_id == input$pathway_display) %>%
                mutate(across(c(pvalue, qvalue, running_enrichment_score), round, 3)) %>%
                select(gene_id,
                       rank,
                       running_enrichment_score,
                       as.name(filter_var)),
              rownames = FALSE)
    
  })
  
  output$pathways_results <- renderDT({
    filter_var = ifelse(input$filter_type == "p-value", "pvalue", "qvalue")
    if (is.null(results$pathways))
      return(NULL)
    if (nrow(results$pathways) == 0)
      return(NULL)
    datatable(results$pathways %>%
                dplyr::mutate(across(c(pvalue, qvalue), round, 3)) %>%
                dplyr::select(.data$pathway_id, 
                              .data$pathway_name, 
                              as.name(filter_var)) %>%
                unique(),
              rownames = FALSE)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
