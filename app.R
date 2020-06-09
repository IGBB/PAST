## app.R ##
packages = c("PAST", "DT", "zip", "shiny", "shinydashboard", "gridExtra")
for (package in packages){
  if( !is.element(package, .packages(all.available = TRUE)) ) {
    if (package == "PAST"){
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("PAST")
    } else {
      install.packages(package)
    }
  }
  library(package,character.only = TRUE)
}

library(PAST)
library(zip)
library(dplyr)
library(ggplot2)
library(shiny)
library(shinydashboard)
library(gridExtra)

ui <- dashboardPage(
  dashboardHeader(title = "PAST"),
  dashboardSidebar(
    sidebarMenu(id = "tabs",
                menuItem("Input", tabName = "input", icon = icon("file-upload")),
                menuItem("Results", tabName = "results", icon = icon("table"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML('
      .box {
      overflow-y: hidden;
      }
      .box h1 {
        font-size: 1.2em;
        font-weight: bold;
      }
      .shiny-input-container, .shiny-file-input-progress {
        margin-bottom: 0;
      }
      .shiny-notification {
        height: 100%;
        width: 100%;
        position:fixed;
        top: 0;
        left: 0;
        font-size: 250%;
        text-align: center;
        background-color: rgba(0, 0, 0, 0.9);
        color: #FFFFFF;
      }
      .shiny-notification-content {
        margin-top: 50vh;
        margin-left: auto;
        margin-right: auto;
        width: 75%;
      }
      .shiny-notification-close {
      display: none;
      }
      .h1, .h2, .h3, h1, h2, h3 {
      margin-top: 0px !important;
      }
      .form-group {
      margin: 10px 0;
      }
    '))),
    tabItems(
      # First tab content
      tabItem(tabName = "input",
              fluidRow(
                box(
                  title = "Analysis Options", status = "primary", solidHeader = TRUE, width = 3, height = 300,
                  style = "height:270px; overflow-y: scroll;",
                  h1("Analysis Title"),
                  tags$p("Provide a title for your analysis. This title is used for the downloaded results."),
                  textInput("title", NULL, value = "New Analysis"),
                  
                  h1("B73 Genome Version"),
                  selectInput("genome_version",
                              label = NULL,
                              choices = list("RefGen_v2","RefGen_v3","RefGen_v4")),
                  
                  # number of cores
                  # h1("Number of Cores"),
                  # tags$p("Select the number of cores to be used in run your analysis in parallel."),
                  # numericInput("num_cores", NULL, value = 1),
                  
                  # analysis type
                  h1("Mode"),
                  tags$p("Select the type of analysis you wish to run. \"Increasing\" searches for pathways associated with an increase in the measured trait, and \"decreasing\" searches for pathways associated with a decrease in the measured trait."),
                  selectInput("mode", NULL,
                              choices = c("increasing", "decreasing"))
                ),
                box(
                  title = "GWAS", status = "primary", solidHeader = TRUE, width = 3, height = 300,
                  style = "height:270px; overflow-y: scroll;",
                  h1("Files"),
                  tags$p("Upload files containing a single trait. These files can be compressed with BZIP, GZIP, or XZ compression."),
                  fileInput("association_file", "Association File"),
                  fileInput("effects_file", "Effects File"),
                  h1("Association Column Names"),
                  tags$p("The column names below must match the column names of your data."),
                  textInput("association_trait",
                            "Trait Column Name",
                            value = "Trait"),
                  textInput("association_marker",
                            "Marker Column Name",
                            value = "Marker"),
                  textInput("association_locus",
                            "Locus Column Name",
                            value = "Locus"),
                  textInput("association_site",
                            "Site Column Name",
                            value = "Site"),
                  textInput("association_p",
                            "P-value Column Name",
                            value = "p"),
                  textInput("association_marker_R2",
                            tagList("Marker R", tags$sup("2"), " Column Name"),
                            value = "marker_R2"),
                  
                  h1("Effects Column Names"),
                  tags$p("The column names below must match the column names of your data."),
                  textInput("effects_trait",
                            "Trait Column Name",
                            value = "Trait"),
                  textInput("effects_marker",
                            "Marker Column Name",
                            value = "Marker"),
                  textInput("effects_locus",
                            "Locus Column Name",
                            value = "Locus"),
                  textInput("effects_site",
                            "Site Column Name",
                            value = "Site"),
                  textInput("effects_effect",
                            "Effect Column Name",
                            value = "Effect")
                ),
                box(
                  title = "Linkage Disequilibrium", status = "primary", solidHeader = TRUE, width = 3, height = 300,
                  style = "height:270px; overflow-y: scroll;",
                  tags$p("Upload the linkage disequilibrium file. This file can be compressed with BZIP, GZIP, or XZ compression."),
                  fileInput("LD_file", "Linkage Disequilibrium File"),
                  h1("LD Column Names"),
                  tags$p("The column names below must match the column names of your data."),
                  textInput("LD_locus_1",
                            "Locus1 Column Name",
                            value = "Locus1"),
                  textInput("LD_position_1",
                            "Position1 Column Name",
                            value = "Position1"),
                  textInput("LD_site_1",
                            "Locus Column Name",
                            value = "Site1"),
                  textInput("LD_position_2",
                            "Site Column Name",
                            value = "Position2"),
                  textInput("LD_site_2",
                            "Site Column Name",
                            value = "Site2"),
                  textInput("LD_distance",
                            "Site Column Name",
                            value = "Dist_bp"),
                  textInput("LD_R2",
                            "Effect Column Name",
                            value = "R.2")
                ),
                box(
                  title = "Gene Assignment", status = "primary", solidHeader = TRUE, width = 3, height = 300,
                  style = "height:270px; overflow-y: scroll;",
                  h1("Gene Assignment Options"),

                  # window size to search for genes
                  tags$p("PAST searches for genes within 1kb of SNPs by default, but this number can be changed."),
                  numericInput("window", "Window Size", value = 1000),
                  
                  # r^2 for LD
                  tags$p(tagList("PAST uses R", tags$sup("2"), " to determine linkage between SNPs.")),
                  sliderInput(
                    "r_squared_cutoff",
                    NULL,
                    min = 0,
                    max = 1,
                    value = 0.8,
                    step = 0.05
                  )
                ),
                box(
                  title = "Pathway Discovery", status = "primary", solidHeader = TRUE, width = 3, height = 300,
                  style = "height:270px; overflow-y: scroll;",
                  # minimum number of genes in a pathway
                  h1("Gene Cutoff"),
                  tags$p("PAST discards pathways with less than 5 genes by default."),
                  numericInput("gene_cutoff", NULL, value = 5),
                  
                  # number of times to sample effects for generating null distribution
                  h1("Effects"),
                  tags$p("PAST determines pathway signficance by creating 1000 random distributions of gene effects and comparing the actual effects to randomly sampled effects."),
                  numericInput("sample", NULL, value = 1000)
                ),
                box(
                  title = "Begin Analysis", status = "primary", solidHeader = TRUE, width = 3, height = 300,
                  style = "height:250px; overflow-y: scroll;",
                  tags$p("When all parameters have been set, click the button below to begin the analysis. The results of the analysis can be viewed in the Results tab."),
                  br(),
                  actionButton("start_button", "Begin Analysis")
                ),
                box(
                  title = "Load Saved Analysis", status = "primary", solidHeader = TRUE, width = 3, height = 300,
                  style = "height:250px; overflow-y: scroll;",
                  tags$p("Upload previously downloaded results The downloaded file should be unzipped, and the full pathways should be uploaded."),
                  fileInput("results_file", "Results File"),
                  tags$p("Once results have been uploaded, click the button below to view the saved analysis in the Results tab."),
                  br(),
                  actionButton("plot_button", "View Saved Results")
                ),
                box(
                  title = "About PAST", status = "primary", solidHeader = TRUE, width = 3, height = 300,
                  style = "height:270px; overflow-y: scroll;",
                  h1("Please cite us if you use PAST."),
                  tags$p(tagList("Thrash A, Tang JD, DeOrnellis M, Peterson DG, Warburton ML (2020). “PAST: The Pathway Association Studies Tool to Infer Biological Meaning from GWAS Datasets.” Plants, 9(1), 58.",
                                 tags$a(href="https://doi.org/10.3390/plants9010058", "https://doi.org/10.3390/plants9010058")
                  )),
                  h1("Report an Issue"),
                  tags$p(tags$a("https://github.com/IGBB/PAST/issues")),
                  h1("Full Documentation"),
                  tags$p(tags$a("https://igbb.github.io/PAST/"))
                )
              )
      ),
      
      # Second tab content
      tabItem(tabName = "results",
              
              column(width = 4,
                     fluidRow(
                       box(
                         title = "Filtering", status = "primary", solidHeader = TRUE, width=12, height = 300,
                         style = "height:270px;overflow-y: scroll;",
                         # create a selector for filtering type
                         "Choose your method of filtering the results.",
                         selectInput(
                           "filter_type",
                           "Filter Parameter:",
                           choices = c("p-value", "q-value"),
                           selected = "p-value"
                         ),
                         
                         
                         "Select the value used for filtering results.",
                         # create a slider to specify filtering level
                         sliderInput(
                           "significance_cutoff",
                           "Pathway Signficance Filter",
                           min = 0,
                           max = 0.3,
                           value = 0.05,
                           step = 0.01
                         )
                       )
                     ),
                     # create the table view
                     fluidRow(
                       box(
                         title = textOutput("box_title_table"), width=12,
                         status = "primary",
                         solidHeader = TRUE,
                         DT::dataTableOutput("pathways"),
                         style = "height:calc(82vh - 320px);overflow-y: scroll;"
                       )
                     )
              ),
              
              # create the graph view
              box(
                title = textOutput("box_title_plot"),
                status = "primary",
                color = "red",
                solidHeader = TRUE,
                width = 8,
                plotOutput("plots", height = "auto"),
                style = "height:82vh; overflow-y: scroll;"
              ),
              
              # add the download button
              uiOutput("download_button")
              
      )
    )
  )
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 600 * 1024 ^ 2)
  # output$title <- renderText({
  #   input$title
  # })
  
  analysis_type <- reactiveVal(NULL)
  output$download_button <- renderUI({
    req(analysis_type()) # requires readfile to be non-null before running code below
    downloadButton("download_data",
                   "Download Results",
                   style = "margin-left: 15px")  })
  
  # set the table box title based on whether we are doing a new analysis
  # or loading saved data
  output$box_title_table <- renderText({
    if (is.null(analysis_type())) {
      return("Pathways")
    }
    if (analysis_type() == "new") {
      return(paste0(input$mode, " Pathways"))
    }
    else if (analysis_type() == "saved"){
      return(paste0(input$results_file$name, " Pathways"))
    }
  })
  
  # set the plot box title based on whether we are doing a new analysis
  # or loading saved data
  output$box_title_plot <- renderText({
    if (is.null(analysis_type())) {
      return("Plots")
    }
    if (analysis_type() == "new") {
      return(paste0(input$mode, " plots"))
    }
    else if (analysis_type() == "new"){
      return(paste0(input$results_file$name, " plots"))
    } else {
      return("plots")
    }
  })
  
  # get GWAS data
  get_gwas_data <- reactive({
    association_file <- input$association_file
    effects_file <- input$effects_file
    association_columns = c(
      input$association_trait,
      input$association_marker,
      input$association_locus,
      input$association_site,
      input$association_p,
      input$association_marker_R2
    )
    
    effects_columns = c(
      input$effects_trait,
      input$effects_marker,
      input$effects_locus,
      input$effects_site,
      input$effects_effect
    )
    
    if (is.null(association_file) | is.null(effects_file))
      return(NULL)
    
    data <- try(load_GWAS_data(
      association_file$datapath,
      effects_file$datapath,
      association_columns,
      effects_columns
    ))
    data
  })
  
  # load the LD data
  get_LD <- reactive({
    LD_file <- input$LD_file
    if (is.null(LD_file))
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
    data <- try({load_LD(LD_file$datapath, LD_columns)})
    data
    
  })
  
  # assign SNPs to genes
  assign_genes <- reactive({
    all_data <- get_gwas_data()
    LD <- get_LD()
    if (input$genome_version == "RefGen_v2") {
      gene_file_path <- "./data/Gene_model_files/Refgen_v2_gene_models.gff"
    } else if (input$genome_version == "RefGen_v3") {
      gene_file_path <- "./data/Gene_model_files/Refgen_v3_gene_models.gff"
    } else {
      gene_file_path <- "./data/Gene_model_files/Refgen_v4_gene_models.gff"
    }
    
    if (is.null(gene_file_path) | is.null(LD))
      return(NULL)
    window <- input$window
    r_squared_cutoff <- input$r_squared_cutoff
    # num_cores <- input$num_cores
    num_cores <- 8
    filter_type = "gene"
    data <- try({assign_SNPs_to_genes(all_data,
                                      LD,
                                      gene_file_path,
                                      c(filter_type),
                                      window,
                                      r_squared_cutoff,
                                      num_cores)
    })
    data
  })
  
  find_pathways <- reactive({
    genes <- assign_genes()
    if (input$genome_version == "RefGen_v2") {
      pathway_file_path <- "./data/Pathway_files/Refgen_v2_pathways.txt"
    } else if (input$genome_version == "RefGen_v3") {
      pathway_file_path <- "./data/Pathway_files/Refgen_v3_pathways.txt"
    } else {
      pathway_file_path <- "./data/Pathway_files/Refgen_v4_pathways.txt"
    }
    if (is.null(pathway_file_path) | is.null(genes))
      return(NULL)
    gene_cutoff <- input$gene_cutoff
    mode <- input$mode
    sample <- input$sample
    # num_cores <- input$num_cores
    num_cores <- 8
    
    data <- try({find_pathway_significance(genes,
                                           pathway_file_path,
                                           gene_cutoff,
                                           mode,
                                           sample,
                                           num_cores)
    })
    data
    
    
  })
  
  observeEvent(input$start_button, {
    analysis_type("new")
    withProgress(message = 'Beginning analysis', value = 0,detail="0%", {
      
      # load GWAS
      gwas_data <- get_gwas_data()
      if (is.null(gwas_data)) {
        showModal(modalDialog(
          title = "Error",
          "Something went wrong loading your GWAS data. Please make sure your files exist and that the uploads were allowed to complete.",
        ))
        return(NULL);
      } else if (class(gwas_data) == "try-error"){
        showModal(modalDialog(
          title = "Error",
          "Something went wrong loading your GWAS data. Please check that the column names match the names provided in the GWAS options.",
          gwas_data,
          br(),
          tags$b("Association File Column Names: "),
          read.table(input$association_file$datapath, nrows = 1, colClasses = "character", header = F),
          br(),
          tags$b("Effects File Column Names: "),
          read.table(input$effects_file$datapath, nrows = 1, colClasses = "character", header = F),
        ))
        return(NULL);
      }
      incProgress(1/4,detail = paste0(25,"%"), message = "Loaded GWAS data")
      
      # load LD
      LD_data <- get_LD()
      if (is.null(LD_data)) {
        showModal(modalDialog(
          title = "Error",
          "Something went wrong loading your LD data. Please make sure your file exists and that the upload was allowed to complete.",
        ))
        return(NULL);
      } else if (class(LD_data) == "try-error"){
        showModal(modalDialog(
          title = "Error",
          "Something went wrong loading your LD data. Please check that the column names match the names provided in the GWAS options.",
          LD_data,
          br(),
          tags$b("LD File Column Names: "),
          read.table(input$LD_file$datapath, nrows = 1, colClasses = "character", header = F)
        ))
        return(NULL);
      }
      incProgress(1/4,detail = paste0(50,"%"), message = "Loaded LD data")
      
      # assign genes to SNPs
      genes <- assign_genes()
      if (is.null(genes)) {
        showModal(modalDialog(
          title = "Error",
          "Something went wrong loading your annotation data. Please make sure your file exists and that the upload was allowed to complete.",
        ))
        return(NULL);
      } else if (class(genes) == "try-error"){
        showModal(modalDialog(
          title = "Error",
          "Something went wrong while assigning SNPs to genes. Please make sure your annotations are correct and that the first column matches the chromosome/locus of GWAS and LD data.",
          genes,
          br(),
          tags$b("Annotations File Column Names: "),
          read.table(gene_file_path, nrows = 1, colClasses = "character", header = F)
        ))
        return(NULL);
      }
      incProgress(1/4,detail = paste0(75,"%"), message = "Assigned SNPs to genes")
      
      # find pathways
      pathways <- find_pathways()
      if (is.null(pathways)) {
        showModal(modalDialog(
          title = "Error",
          "Something went wrong loading your pathways data. Please make sure your file exists and that the upload was allowed to complete.",
        ))
        return(NULL);
      } else if (class(pathways) == "try-error"){
        showModal(modalDialog(
          title = "Error",
          "Something went wrong assigning genes to pathways. Please make sure that the gene names in your annotations file match the gene names in your pathways file.",
          pathways,
          br(),
          tags$b("Pathways File Column Names: "),
          read.table(pathway_file_path, nrows = 1, colClasses = "character", header = F)
        ))
        return(NULL);
      }
      incProgress(1/4,detail = paste0(100,"%"), message = "Discovered significant pathways")
      newtab <- switch(input$tabs, "input" = "results")
      updateTabItems(session, "tabs", newtab)
      
    })
  })
  
  # load pathway sifnificance data from file
  load_pathways_from_file <- reactive({
    pathway_file <- input$results_file
    if (is.null(pathway_file))
      return(NULL)
    data <- try({read.table(pathway_file$datapath, header = TRUE, sep = "\t")})
    if (class(data) != "try-error"){
      data
    } else {
      NULL
    }
  })
  
  observeEvent(input$plot_button, {
    analysis_type("saved")
    # find pathways
    pathways <- load_pathways_from_file()
    if(is.null(pathways)) {
      showModal(modalDialog(
        title = "Error",
        "Something went wrong loading your data. Please make sure your pathways file is correct."
      ))
      return(NULL);
    }      
    newtab <- switch(input$tabs, "input" = "results")
    updateTabItems(session, "tabs", newtab)
  })
  
  # render pathway table
  output$pathways <- DT::renderDataTable({
    # determine source of data based on analysis type
    if (analysis_type() == "new") {
      pathways <- find_pathways()
    } else {
      pathways <- load_pathways_from_file()
    }
    
    # get filter and cut-off
    filter_type <- input$filter_type
    significance_cutoff <- input$significance_cutoff
    if (is.null(pathways))
      return(NULL)
    
    # get the pathway, p-value, and q-value
    pathways <- pathways %>%
      dplyr::select(.data$pathway_id, .data$pvalue, .data$qvalue) %>%
      unique()
    
    # filter pathways
    if (filter_type == "p-value") {
      pathways <- dplyr::filter(pathways, pvalue < significance_cutoff)
    } else {
      pathways <- dplyr::filter(pathways, qvalue < significance_cutoff)
    }
    
    # render data table without paging and no rownames
    DT::datatable(pathways,
                  rownames = FALSE,
                  options = list(paging = FALSE))
  })
  
  # render plots
  output$plots <- renderPlot({
    # determine source of data based on analysis type
    if (is.null(analysis_type())) {
      return(NULL)
    }
    if (analysis_type() == "new") {
      pathways <- find_pathways()
    } else {
      pathways <- load_pathways_from_file()
    }
    
    # get filter and cut-off
    filter_type <- input$filter_type
    significance_cutoff <- input$significance_cutoff
    if (is.null(pathways))
      return(NULL)
    # filter pathways
    rugplots_data <-
      pathways %>% dplyr::arrange(.data$pathway_number)
    if (filter_type == "p-value") {
      rugplots_data <- dplyr::filter(rugplots_data,
                                     pvalue < significance_cutoff)
    } else {
      rugplots_data <- dplyr::filter(rugplots_data,
                                     qvalue < significance_cutoff)
    }
    
    # split rugplots data into a list of pathways
    rugplots_split <-
      split(rugplots_data, rugplots_data$pathway_number)
    
    # set up a list to store plot instructions
    plots_list <- list()
    
    # for each pathway, draw rugplot
    for (rank in names(rugplots_split)) {
      # get data
      temp_data <- rugplots_split[[rank]]
      
      # title is "PWY-ID - Pathway Description"
      title <- paste0(unique(as.character(temp_data$pathway_id)),
                      " - ",
                      unique(as.character(temp_data$pathway_name)))
      
      # intercept should be at rank of maximum enrichment score
      intercept <- temp_data %>%
        dplyr::arrange(desc(.data$running_enrichment_score)) %>%
        dplyr::select(.data$rank)
      intercept <- intercept[, 1][1]
      
      # set up the rugplot
      rugplot <-
        ggplot(temp_data, aes(x = rank, y = running_enrichment_score)) +
        geom_line(stat = "identity") +
        geom_rug(sides = "t", position = "jitter") +
        geom_vline(xintercept = intercept,
                   color = "black",
                   linetype = "longdash") +
        ggtitle(title) +
        labs(x = "Gene Rank", y = "Running Enrichment Score") +
        scale_x_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000)) +
        theme(
          axis.text = element_text (color = "black"),
          panel.background = element_rect (color = "black", fill = "pink")
        )
      
      # store the rugplot in the list
      plots_list[[rank]] <- rugplot
    }
    
    # draw a blank if no rugplots were plotted, otherwise, draw rugplots in grid
    if (nrow(rugplots_data) == 0) {
      ggplot() + geom_blank()
    } else {
      columns <- floor(sqrt(length(plots_list)))
      do.call("grid.arrange", c(plots_list, ncol = columns))
    }
  } # complicated way to set the height based on the number of columns
  , height = reactive({
    # get number of pathways being plotted
    if (analysis_type() == "new") {
      pathways <- find_pathways()
    } else {
      pathways <- load_pathways_from_file()
    }
    
    # filter
    filter_type <- input$filter_type
    significance_cutoff <- input$significance_cutoff
    
    # return 100 if there are no pathways to keep the function from breaking
    if (is.null(pathways))
      return(100)
    
    # get rugplots data for column calculation
    rugplots_data <-
      pathways %>% dplyr::arrange(.data$pathway_number)
    if (filter_type == "p-value") {
      rugplots_data <- dplyr::filter(rugplots_data,
                                     pvalue < significance_cutoff)
    } else {
      rugplots_data <- dplyr::filter(rugplots_data,
                                     qvalue < significance_cutoff)
    }
    
    # return 100 if there are no pathways to keep the function from breaking
    if (nrow(rugplots_data) == 0)
      return(100)
    
    # split rugplots data into a list of pathways
    rugplots_split <-
      split(rugplots_data, rugplots_data$pathway_number)
    
    # get number of columns
    columns <- floor(sqrt(length(names(rugplots_split))))
    
    # return height
    return(length(names(rugplots_split)) / columns * 200)
    
  }))
  
  # download handler
  output$download_data <- downloadHandler(
    # make filename = analysis_title.zip
    filename = function() {
      paste(input$title, "zip", sep = ".")
    },
    # set up content
    content = function(filename) {
      # get filter type and cutoff
      filter_type <- input$filter_type
      significance_cutoff <- input$significance_cutoff
      
      # set up empty vector of files and move to tempdir()
      fs <- c()
      setwd(tempdir())
      
      # get pathway significance data based on analysis type
      # R Shiny will either calculate this if something has changed or use
      # what's already on the screen in most cases
      if (analysis_type() == "new") {
        pathways <- find_pathways()
      } else {
        pathways <- load_pathways_from_file()
      }
      
      # write full pathways file
      write.table(pathways, paste0(input$title, ".pathways.tsv"), sep = "\t")
      
      # filter pathways
      if (filter_type == "p-value") {
        filtered_pathways <- dplyr::filter(pathways,
                                           pvalue < significance_cutoff)
      } else {
        filtered_pathways <- dplyr::filter(pathways,
                                           qvalue < significance_cutoff)
      }
      
      # write filtered pathways file
      write.table(filtered_pathways,
                  paste0(input$title, ".pathways.filtered.tsv"),
                  sep = "\t")
      
      # store pathways files
      fs <-
        c(
          paste0(input$title, ".pathways.tsv"),
          paste0(input$title, ".pathways.filtered.tsv")
        )
      
      # break until we have pathways information
      if (is.null(pathways))
        return(NULL)
      
      # set up rugplots by arranged by pathway_number and filtering
      rugplots_data <-
        pathways %>% dplyr::arrange(.data$pathway_number)
      if (filter_type == "p-value") {
        rugplots_data <- dplyr::filter(rugplots_data,
                                       pvalue < significance_cutoff)
      } else {
        rugplots_data <- dplyr::filter(rugplots_data,
                                       qvalue < significance_cutoff)
      }
      
      # split based on pathway number
      rugplots_split <-
        split(rugplots_data, rugplots_data$pathway_number)
      
      # for each pathway, draw rugplot
      for (rank in names(rugplots_split)) {
        # get data
        temp_data <- rugplots_split[[rank]]
        
        # title is "PWY-ID - Pathway Description"
        title <- paste0(unique(as.character(temp_data$pathway_id)),
                        " - ",
                        unique(as.character(temp_data$pathway_name)))
        
        # intercept should be at rank of maximum enrichment score
        intercept <-
          temp_data %>% dplyr::arrange(desc(.data$running_enrichment_score)) %>%
          dplyr::select(rank)
        intercept <- intercept[, 1][1]
        
        # set up the rugplot
        rugplot <-
          ggplot(temp_data, aes(x = rank, y = running_enrichment_score)) +
          geom_line(stat = "identity") +
          geom_rug(sides = "t", position = "jitter") +
          geom_vline(xintercept = intercept,
                     color = "black",
                     linetype = "longdash") +
          ggtitle(title) +
          labs(x = "Gene Rank", y = "Running Enrichment Score") +
          scale_x_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000)) +
          theme(
            axis.text = element_text (color = "black"),
            panel.background = element_rect (color = "black", fill = "pink")
          )
        
        # set up the output path for the rugplot
        path <- paste0(input$title,
                       ".",
                       unique(as.character(temp_data$pathway_id)),
                       ".png")
        
        # add the rugplot to the files to be zipped and save it in tempdir()
        fs <- c(fs, path)
        ggsave(path, rugplot)
      }
      
      # zip the file using the name and files specified earlier
      zipr(zipfile = filename, files = fs)
    },
    contentType = "application/zip"
  )
}


shinyApp(ui, server)