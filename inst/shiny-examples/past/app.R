ui <- dashboardPage(
  title = "PAST",
  shinydashboard::dashboardHeader(title = textOutput("title")),
  shinydashboard::dashboardSidebar(
    sidebarMenu(
      textInput("title", "Analysis Title", value = "New Analysis"),
      selectInput(
        "type",
        "Analysis Type:",
        choices = c("new", "saved"),
        selected = "new"
      ),
      menuItemOutput("menuitem"),
      menuItem(
        "Plot",
        sliderInput(
          "significance_cutoff",
          "Pathway Signficance Filter",
          min = 0,
          max = 0.2,
          value = 0.05,
          step = 0.01
        )
      )
    )
  ),
  shinydashboard::dashboardBody(fluidRow(
    box(
      title = textOutput("box_title_table"),
      width = 4,
      status = "primary",
      solidHeader = TRUE,
      DT::dataTableOutput("pathways"),
      style = "height:82vh; overflow-y: scroll;"
    ),
    box(
      title = textOutput("box_title_plot"),
      status = "primary",
      color = "red",
      solidHeader = TRUE,
      width = 8,
      plotOutput("plots", height = "auto"),
      style = "height:82vh; overflow-y: scroll;"
    ),
    downloadButton("download_data", "Download Results",
                   style = "margin-left: 15px")
  ))
)

server <- function(input, output) {
  options(shiny.maxRequestSize = 600 * 1024 ^ 2)
  output$title <- renderText({
    input$title
  })
  output$menuitem <- renderMenu({
    if (input$type == "new")
      return(
        menuItem(
          "Parameters",
          tabName =  "live",
          fileInput("association_file", "Association File"),
          fileInput("effects_file", "Effects File"),
          fileInput("LD_file", "Linkage Disequilibrium File"),
          fileInput("genes_file", "Genes File"),
          fileInput("pathway_file", "Pathways File"),
          numericInput("num_cores", "Number of Cores", value = 4),
          selectInput("mode", "Mode:",
                      choices = c("increasing", "decreasing")),
          menuItem(
            "Advanced",
            numericInput("window", "Window Size", value = 1000),
            sliderInput(
              "r_squared_cutoff",
              "R^2 cutoff:",
              min = 0,
              max = 1,
              value = 0.8,
              step = 0.05
            ),
            numericInput("gene_cutoff", "Gene Cutoff", value = 5),
            numericInput("sample", "Effects", value = 1000)
          )
        )
      )
    else
      (return (menuItem("Saved",
                        fileInput("load_file", "Data"))))
  })
  output$download_data <- downloadHandler(
    filename = function() {
      print(input$title)
      paste(input$title, ".zip", sep = "")
    },
    content = function(filename) {
      significance_cutoff <- input$significance_cutoff
      fs <- c()
      setwd(tempdir())
      if (input$type == "new") {
        pathways <- find_pathways()
      } else {
        pathways <- load_pathways_from_file()
      }
      write.table(pathways, paste0(input$title, ".pathways.tsv"), sep = "\t")
      write.table(
        pathways %>% dplyr::filter(data$pvalue < significance_cutoff),
        paste0(input$title, ".pathways.filtered.tsv"),
        sep = "\t"
      )
      significance_cutoff <- input$significance_cutoff
      fs <-
        c(
          paste0(input$title, ".pathways.tsv"),
          paste0(input$title, ".pathways.filtered.tsv")
        )
      if (is.null(pathways))
        return(NULL)
      rugplots_data <-
        pathways %>% dplyr::arrange(.data$pathway_number) %>%
        dplyr::filter(.data$pvalue < significance_cutoff)
      rugplots_split <- split(rugplots_data, rugplots_data$pathway_number)
      for (rank in names(rugplots_split)) {
        temp_data <- rugplots_split[[rank]]
        title <-
          paste0(unique(as.character(temp_data$pathway_id)),
                 " - ",
                 unique(as.character(temp_data$pathway_name)))
        print(title)
        intercept <-
          temp_data %>% dplyr::arrange(desc(.data$running_enrichment_score)) %>%
          dplyr::select(rank)
        intercept <- intercept[, 1][1]
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
        path <-
          paste0(input$title,
                 ".",
                 unique(as.character(temp_data$pathway_id)),
                 ".png")
        print(path)
        fs <- c(fs, path)
        ggsave(path, rugplot)
      }
      zip(zipfile = filename, files = fs)
    },
    contentType = "application/zip"
  )
  output$box_title_table <- renderText({
    if (input$type == "new") {
      return(paste0(input$mode, " pathways"))
    }
    else {
      return(paste0(input$load_file$name, " pathways"))
    }
  })
  output$box_title_plot <- renderText({
    if (input$type == "new") {
      return(paste0(input$mode, " plots"))
    }
    else {
      return(paste0(input$load_file$name, " plots"))
    }
  })
  get_all_data <- reactive({
    association_file <- input$association_file
    effects_file <- input$effects_file
    if (is.null(association_file) | is.null(effects_file))
      return(NULL)
    print("Merging data")
    merge_data(association_file$datapath, effects_file$datapath)
  })
  LD <- reactive({
    LD_file <- input$LD_file
    if (is.null(LD_file))
      return(NULL)
    print("Parsing LD")
    parse_LD(LD_file$datapath)
  })
  genes <- reactive({
    all_data <- get_all_data()
    LD <- LD()
    genes_file <- input$genes_file
    if (is.null(genes_file) | is.null(LD))
      return(NULL)
    window <- input$window
    r_squared_cutoff <- input$r_squared_cutoff
    num_cores <- input$num_cores
    print("Finding genes")
    parse_SNP(all_data,
              LD,
              genes_file$datapath,
              window,
              r_squared_cutoff,
              num_cores)
  })
  find_pathways <- reactive({
    genes <- genes()
    pathway_file <- input$pathway_file
    if (is.null(pathway_file) | is.null(genes))
      return(NULL)
    gene_cutoff <- input$gene_cutoff
    mode <- input$mode
    sample <- input$sample
    num_cores <- input$num_cores
    print("Finding pathways")
    analyze_pathways(genes,
                     pathway_file$datapath,
                     gene_cutoff,
                     mode,
                     sample,
                     num_cores)
  })
  load_pathways_from_file <- reactive({
    pathway_file <- input$load_file
    if (is.null(pathway_file))
      return(NULL)
    read.table(pathway_file$datapath, header = TRUE, sep = "\t")
  })
  output$pathways <- DT::renderDataTable({
    if (input$type == "new") {
      pathways <- find_pathways()
    } else {
      pathways <- load_pathways_from_file()
    }
    significance_cutoff <- input$significance_cutoff
    if (is.null(pathways))
      return(NULL)
    pathways <- pathways %>%
      dplyr::select(.data$pathway_id, .data$pvalue, .data$qvalue) %>%
      unique() %>% dplyr::filter(.data$pvalue < significance_cutoff)
    datatable(pathways,
              rownames = FALSE,
              options = list(paging = FALSE))
  })
  output$plots <- renderPlot({
    if (input$type == "new") {
      pathways <- find_pathways()
    } else {
      pathways <- load_pathways_from_file()
    }
    significance_cutoff <- input$significance_cutoff
    if (is.null(pathways))
      return(NULL)
    rugplots_data <- pathways %>%
      dplyr::arrange(.data$pathway_number) %>%
      dplyr::filter(.data$pvalue < significance_cutoff)
    rugplots_split <-
      split(rugplots_data, rugplots_data$pathway_number)
    p <- list()
    for (rank in names(rugplots_split)) {
      temp_data <- rugplots_split[[rank]]
      title <-
        paste0(unique(as.character(temp_data$pathway_id)),
               " - ",
               unique(as.character(temp_data$pathway_name)))
      intercept <- temp_data %>%
        dplyr::arrange(desc(.data$running_enrichment_score)) %>%
        dplyr::select(.data$rank)
      intercept <- intercept[, 1][1]
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
      p[[rank]] <- rugplot
    }
    if (nrow(rugplots_data) == 0) {
      ggplot() + geom_blank()
    } else {
      columns <- floor(sqrt(length(p)))
      do.call("grid.arrange", c(p, ncol = columns))
    }
  }
  , height = reactive({
    columns <- floor(sqrt(length(p)))
    if (input$type == "new") {
      pathways <- find_pathways()
    } else {
      pathways <- load_pathways_from_file()
    }
    significance_cutoff <-
      input$significance_cutoff
    if (is.null(pathways))
      return(100)
    pathways <-
      pathways %>%
      dplyr::select(.data$pathway_id, .data$pvalue, .data$qvalue) %>%
      unique() %>%
      dplyr::filter(.data$pvalue < significance_cutoff)
    if (nrow(pathways) == 0)
      return(100)
    return(nrow(pathways) / columns * 200)
  }))
}
shinyApp(ui = ui, server = server)
