#' Interactive visualization of DE results
#'
#' @param results List of result data frames
#' @param dimred List of dimension reduction results
#' @param genemodels Gene models
#' @param bwFiles Vector of bigWig files
#' @param bwCond Vector of conditions for bigWig files
#' @param appTitle App title
#'
#' @import shiny ggplot2 ggrepel dplyr tidyr Gviz
#' @export
#'
iResViewer <- function(results = list(), dimred = list(), genemodels = NULL,
                       bwFiles = NULL, bwCond = NULL, appTitle = "iResViewer") {
  options(ucscChromosomeNames = FALSE)
  options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)

  p_layout <- function() {
    shinydashboard::dashboardPage(
      skin = "blue",

      shinydashboard::dashboardHeader(title = appTitle,
                                      titleWidth = nchar(appTitle) * 20),

      shinydashboard::dashboardSidebar(
        textInput(inputId = "sel.gene", label = "Selected gene")
      ),

      shinydashboard::dashboardBody(fluidRow(
        do.call(shinydashboard::tabBox,
                c(
                  width = 12,
                  lapply(names(dimred), function(w)
                    tabPanel(w, selectInput(paste0(w, "_pcx"), "x-axis component",
                                            choices = 1:7,
                                            selected = 1),
                             selectInput(paste0(w, "_pcy"), "y-axis component",
                                         choices = 1:7,
                                         selected = 2),
                             uiOutput(paste0(w, "_dimred_ui")))),
                  lapply(names(results), function(w)
                    tabPanel(w, DT::dataTableOutput(paste0(w, "_restable")),
                             p(class = "text-center",
                               downloadButton(paste0(w, "_restable_download"), "Download")))),
                  list(tabPanel("Coverage plot",
                                uiOutput("coverage.plot.ui")))

                ),
        )
      ))
    )
  }

  server_function <- function(input, output, session) {
    options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)

    ## Generate result tables
    Map(function(nm) {
      output[[paste0(nm, "_restable")]] <- DT::renderDataTable(
        results[[nm]],
        filter = "top",
        rownames = FALSE,
        options = list(scrollX = TRUE))
    }, names(results))

    ## Generate download buttons for result tables
    Map(function(nm) {
      output[[paste0(nm, "_restable_download")]] <- downloadHandler(
        paste0(nm, "_results_filtered.csv"),
        content = function(file) {
          tmp <- input[[paste0(nm, "_restable_rows_all")]]
          write.csv(results[[nm]][tmp, , drop = FALSE], file)
        })
    }, names(results))

    ## Generate coverage plot
    output$coverage.plot <- renderPlot({
      if (input$sel.gene == "")
        return(NULL)
      else
        plotTracks(mygene = input$sel.gene, genemodels = genemodels,
                   genemodels2 = NULL,
                   gtf_file = NULL, rnaseq_datafiles = bwFiles,
                   rnaseq_condition = bwCond, show_chr = NULL,
                   min_coord = NULL, max_coord = NULL,
                   pdf_filename = NULL, pdf_width = 7, pdf_height = 7)
    })

    output$coverage.plot.ui <- renderUI({
      plotOutput("coverage.plot", height = "800px")
    })

    ## Generate dimension reduction plots
    Map(function(nm) {
      output[[paste0(nm, "_dimred")]] <- renderPlot(
        ggplot(dimred[[nm]], aes_string(x = grep(paste0(input[[paste0(nm, "_pcx")]], "$"),
                                                 colnames(dimred[[nm]]), value = TRUE),
                                        y = grep(paste0(input[[paste0(nm, "_pcy")]], "$"),
                                                 colnames(dimred[[nm]]), value = TRUE))) +
          geom_point() + coord_fixed())
    }, names(dimred))

    Map(function(nm) {
      output[[paste0(nm, "_dimred_ui")]] <- renderUI(
        plotOutput(paste0(nm, "_dimred"), height = "600px"))
    }, names(dimred))

    # output$pca.plot <- renderPlot({
    #   if (!is.null(input$pcx) & !is.null(input$pcy)) {
    #     ggplot(res$pca, aes_string(x = input$pcx, y = input$pcy, color = "trttime",
    #                                label = "sample_id")) +
    #       geom_point(size = 7) + geom_label_repel() + theme_bw() +
    #       xlab(paste0(input$pcx, res$pcavar[input$pcx])) +
    #       ylab(paste0(input$pcy, res$pcavar[input$pcy])) +
    #       coord_fixed() + theme(axis.text = element_text(size = 14),
    #                             axis.title = element_text(size = 16)) +
    #       scale_colour_manual(values = c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D", "black"))
    #   } else {
    #     NULL
    #   }
    # })




  }

  shinyApp(ui = p_layout, server = server_function)
}
