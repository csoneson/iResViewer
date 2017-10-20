#' Interactive visualization of differential gene expression results
#'
#' @param wideResults Named list of result data frames. These will be displayed
#'   as searchable tables, each in its own tab.
#' @param longResults Named list of result data frames in "long" form (one row
#'   per gene/contrast combination). These will be used to generate volcano
#'   plots. Each data frame should have at least four columns: gene (the gene
#'   ID), logFC (log-fold change), FDR (adjusted p-value) and mlog10PValue
#'   (-log10(nominal p-value)).
#' @param dimReds Named list of dimension reduction results. Each element of the
#'   list must have at least three columns: two columns with coordinates in the
#'   low-dimensional space and one column with sample IDs or group labels.
#' @param geneModels A GRanges object with gene models, typically generated from
#'   a gtf file.
#' @param geneInfo Data frame with gene annotation information. Should have at
#'   least two columns: gene (the gene ID) and symbol (the gene symbol).
#' @param bwFiles Named vector with paths to bigWig files. These will be used to
#'   generate coverage plots.
#' @param bwCond Named vector corresponding to \code{bwFiles}, giving the group
#'   label for each bigWig file. These will be used to color the coverage plots
#'   by group.
#' @param abundances Named list with data frames in "long" format (one row per
#'   gene/sample combination), containing abundance estimates. These will be
#'   used to illustrate the abundance pattern for selected genes. Each data
#'   frame must have at least four columns: sample (the sample ID), gene (the
#'   gene ID), value (the abundance) and group (a sample group label, used to
#'   order and color the points in the plot).
#' @param appTitle App title
#' @param ... Additional arguments (currently not used)
#'
#' @author Charlotte Soneson
#'
#' @import shiny ggplot2
#' @importFrom utils write.csv
#' @importFrom dplyr %>% filter arrange
#' @importFrom ggrepel geom_label_repel
#' @importFrom shinydashboard dashboardPage dashboardHeader dashboardSidebar
#'   dashboardBody tabBox
#' @export
#'
iResViewer <- function(wideResults = list(), longResults = list(),
                       dimReds = list(), geneModels = NULL, geneInfo = NULL,
                       bwFiles = NULL, bwCond = NULL, abundances = list(),
                       appTitle = "iResViewer", ...) {
  options(ucscChromosomeNames = FALSE)
  options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)

  pLayout <- function() {
    shinydashboard::dashboardPage(
      skin = "purple",

      shinydashboard::dashboardHeader(title = appTitle,
                                      titleWidth = nchar(appTitle) * 20),

      shinydashboard::dashboardSidebar(
        shiny::textInput(inputId = "sel.gene", label = "Selected gene")
      ),

      shinydashboard::dashboardBody(shiny::fluidRow(
        do.call(shinydashboard::tabBox,
                c(
                  width = 12,

                  ## ======================================================== ##
                  ## Tabs with dimension reduction results
                  ## ======================================================== ##
                  lapply(names(dimReds), function(w)
                    shiny::tabPanel(
                      w,
                      shiny::fluidRow(
                        shiny::column(
                          4, shiny::selectInput(paste0(w, "_pcx"), "x-axis",
                                                choices = colnames(dimReds[[w]]),
                                                selected = colnames(dimReds[[w]])[2])),
                        shiny::column(
                          4, shiny::selectInput(paste0(w, "_pcy"), "y-axis",
                                                choices = colnames(dimReds[[w]]),
                                                selected = colnames(dimReds[[w]])[3])),
                        shiny::column(
                          4, shiny::selectInput(paste0(w, "_pccol"), "color by",
                                                choices = colnames(dimReds[[w]]),
                                                selected = colnames(dimReds[[w]])[ncol(dimReds[[w]])]))
                      ),
                      shiny::fluidRow(
                        shiny::column(
                          4, shiny::selectInput(paste0(w, "_pclab"), "Labels",
                                                choices = colnames(dimReds[[w]]),
                                                selected = colnames(dimReds[[w]])[1])),
                        shiny::column(
                          4, shiny::checkboxInput(paste0(w, "_dopclab"), "Show labels",
                                                  value = FALSE))
                      ),
                      shiny::uiOutput(paste0(w, "_dimred_ui")))
                  ),

                  ## ======================================================== ##
                  ## Tabs with result tables
                  ## ======================================================== ##
                  lapply(names(wideResults), function(w)
                    shiny::tabPanel(
                      w,
                      DT::dataTableOutput(paste0(w, "_restable")),
                      shiny::p(class = "text-center",
                               shiny::downloadButton(paste0(w, "_restable_download"),
                                                     "Download")))
                    ),

                  ## ======================================================== ##
                  ## Tabs with abundances
                  ## ======================================================== ##
                  lapply(names(abundances), function(w)
                    shiny::tabPanel(
                      w,
                      shiny::div(style = "position:relative",
                                 shiny::uiOutput(paste0(w, "_abundance_ui")),
                                 shiny::uiOutput(paste0(w, "_abundance_hover_info"))))
                  ),

                  ## ======================================================== ##
                  ## Tabs with volcano plots
                  ## ======================================================== ##
                  lapply(names(longResults), function(w)
                    shiny::tabPanel(
                      paste0(w, " volcano"),
                      shiny::div(style = "position:relative",
                                 shiny::uiOutput(paste0(w, "_volcano_ui")),
                                 shiny::uiOutput(paste0(w, "_volcano_hover_info"))))
                  ),

                  ## ======================================================== ##
                  ## Tab with coverage plot
                  ## ======================================================== ##
                  list(shiny::tabPanel(
                    "Coverage plot",
                    shiny::uiOutput("coverage.plot.ui"))
                  )

                )
        )
      ))
    )
  }

  server_function <- function(input, output, session) {
    options(ucscChromosomeNames = FALSE, envir = .GlobalEnv)

    ## ====================================================================== ##
    ## Generate result tables
    ## ====================================================================== ##
    Map(function(nm) {
      output[[paste0(nm, "_restable")]] <-
        DT::renderDataTable(
          wideResults[[nm]],
          filter = "top",
          rownames = FALSE,
          options = list(scrollX = TRUE))
    }, names(wideResults))

    ## ====================================================================== ##
    ## Generate download buttons for result tables
    ## ====================================================================== ##
    Map(function(nm) {
      output[[paste0(nm, "_restable_download")]] <-
        shiny::downloadHandler(
          paste0(nm, "_results_filtered.csv"),
          content = function(file) {
            tmp <- input[[paste0(nm, "_restable_rows_all")]]
            utils::write.csv(wideResults[[nm]][tmp, , drop = FALSE], file)
          })
    }, names(wideResults))

    ## ====================================================================== ##
    ## Generate coverage plot
    ## ====================================================================== ##
    output$coverage.plot <- shiny::renderPlot({
      if (input$sel.gene == "")
        return(NULL)
      else
        plotGvizTracks(showGene = gsub("\\.[0-9]+$", "", input$sel.gene),
                       geneModels = geneModels, geneModels2 = NULL,
                       gtfFile = NULL, bwFiles = bwFiles,
                       bwCond = bwCond, showChr = NULL,
                       minCoord = NULL, maxCoord = NULL,
                       pdfFilename = NULL, pdfWidth = 7, pdfHeight = 7)
    })

    output$coverage.plot.ui <- shiny::renderUI({
      shiny::plotOutput("coverage.plot", height = "800px")
    })

    ## ====================================================================== ##
    ## Generate dimension reduction plots
    ## ====================================================================== ##
    Map(function(nm) {
      output[[paste0(nm, "_dimred")]] <- shiny::renderPlot({
        p <- ggplot(dimReds[[nm]], aes_string(x = input[[paste0(nm, "_pcx")]],
                                              y = input[[paste0(nm, "_pcy")]],
                                              color = input[[paste0(nm, "_pccol")]],
                                              label = input[[paste0(nm, "_pclab")]])) +
          geom_point(size = 5) + coord_fixed() + theme_bw() +
          theme(axis.text = element_text(size = 14),
                axis.title = element_text(size = 16))
        if (input[[paste0(nm, "_dopclab")]]) p <- p + ggrepel::geom_label_repel()
        p
      })
    }, names(dimReds))

    Map(function(nm) {
      output[[paste0(nm, "_dimred_ui")]] <- shiny::renderUI(
        shiny::plotOutput(paste0(nm, "_dimred"), height = "600px"))
    }, names(dimReds))

    ## ====================================================================== ##
    ## Abundance plots
    ## ====================================================================== ##
    Map(function(nm) {
      output[[paste0(nm, "_abundance_ui")]] <- shiny::renderUI({
        shiny::plotOutput(paste0(nm, "_abundance"), hover = paste0(nm, "_abundance_hover"),
                          width = "90%", height = "500px")
      })
    }, names(abundances))

    Map(function(nm) {
      output[[paste0(nm, "_abundance")]] <- shiny::renderPlot({
        if (input$sel.gene == "") return(NULL)
        else {
          if (!(tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) %in%
                c(gsub("\\.[0-9]+$", "", tolower(geneInfo$gene)),
                  tolower(geneInfo$symbol)))) return(NULL)
          id <- (geneInfo %>% filter(gsub("\\.[0-9]+$", "",
                                          tolower(gene)) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) |
                                       tolower(symbol) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene))))$gene
          if (is.null(id)) return(NULL)
          else {
            df <- abundances[[nm]] %>% dplyr::filter(gene %in% id) %>%
              dplyr::arrange(group, gene) %>%
              dplyr::mutate(sample = factor(sample, levels = unique(sample)))
            ggplot(df, aes(x = sample, y = value, group = gene, col = group)) +
              geom_line(col = "black") + geom_point(size = 3) +
              theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
                                 legend.position = "bottom",
                                 axis.text.y = element_text(size = 14),
                                 axis.title.y = element_text(size = 16)) +
              guides(fill = guide_legend(nrow = 2, byrow = TRUE)) + xlab("") +
              ylab("abundance")
          }
        }
      })
    }, names(abundances))

    Map(function(nm) {
      output[[paste0(nm, "_abundance_hover_info")]] <- shiny::renderUI({
        if (input$sel.gene == "") return(NULL)
        else {
          if (!(tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) %in%
                c(gsub("\\.[0-9]+$", "", tolower(geneInfo$gene)),
                  tolower(geneInfo$symbol)))) return(NULL)
          id <- (geneInfo %>% filter(gsub("\\.[0-9]+$", "",
                                          tolower(gene)) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) |
                                       tolower(symbol) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene))))$gene
          if (is.null(id)) return(NULL)
          else {
            hover <- input[[paste0(nm, "_abundance_hover")]]
            df <- abundances[[nm]] %>% dplyr::filter(gene %in% id) %>%
              dplyr::arrange(group, gene) %>%
              dplyr::mutate(sample = factor(sample, levels = unique(sample)))
            res2 <- shiny::nearPoints(df, hover, threshold = 5, maxpoints = 1)
            res2$genename <- geneInfo$symbol[match(res2$gene, geneInfo$gene)]
            left_pct <- (hover$x - hover$domain$left)/(hover$domain$right - hover$domain$left)
            top_pct <- (hover$domain$top - hover$y)/(hover$domain$top - hover$domain$bottom)
            left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
            top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
            style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                            "left:", left_px + 2, "px; top:", top_px + 2, "px;")
            shiny::wellPanel(
              style = style,
              shiny::p(shiny::HTML(paste0("<b> Sample: </b>", res2$sample, "<br/>",
                                          "<b> ID: </b>", res2$gene, "<br/>",
                                          "<b> Gene: </b>", res2$genename, "<br/>",
                                          "<b> abundance: </b>", round(res2$value, 3), "<br/>")))
            )
          }
        }
      })
    }, names(abundances))

    ## ====================================================================== ##
    ## Volcano plots
    ## ====================================================================== ##
    Map(function(nm) {
      output[[paste0(nm, "_volcano_ui")]] <- shiny::renderUI({
        shiny::plotOutput(paste0(nm, "_volcano"), hover = paste0(nm, "_volcano_hover"),
                          click = paste0(nm, "_volcano_click"),
                          width = "90%", height = "500px")
      })
    }, names(longResults))

    Map(function(nm) {
      pp <- shiny::reactive(
        shiny::nearPoints(longResults[[nm]], input[[paste0(nm, "_volcano_click")]], threshold = 5,
                          maxpoints = 1, panelvar1 = "contrast"))
      shiny::observeEvent(input[[paste0(nm, "_volcano_click")]], {
        shiny::updateTextInput(session, "sel.gene", value = gsub("\\.[0-9]+$", "", pp()$gene))
      })
    }, names(longResults))

    Map(function(nm) {
      output[[paste0(nm, "_volcano")]] <- shiny::renderPlot({
        p <- ggplot(longResults[[nm]], aes(x = logFC, y = mlog10PValue,
                                           col = as.character(FDR <= 0.05))) +
          theme_bw() + ylab("-log10(PValue)") + facet_grid(~contrast) +
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                             name = "FDR <= 0.05") +
          theme(axis.text = element_text(size = 14),
                axis.title = element_text(size = 16),
                strip.text = element_text(size = 15))
        if (input$sel.gene == "") return(p + geom_point(size = 1))
        else {
          if (!(tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) %in%
                c(gsub("\\.[0-9]+$", "", tolower(longResults[[nm]]$gene)),
                  tolower(longResults[[nm]]$symbol))))
            return(p + geom_point(size = 1))
          id <- (longResults[[nm]] %>%
                   dplyr::filter(gsub("\\.[0-9]+$", "",
                                      tolower(gene)) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene)) |
                                   tolower(symbol) == tolower(gsub("\\.[0-9]+$", "", input$sel.gene))))$gene
          if (is.null(id)) return(p + geom_point(size = 1))
          p + geom_point(size = 1, alpha = 0.1) +
            geom_point(data = longResults[[nm]] %>% dplyr::filter(gene %in% id), size = 4, pch = 21,
                       col = "yellow", aes(fill = as.character(FDR <= 0.05))) +
            guides(fill = FALSE) +
            scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                              name = "FDR <= 0.05")
        }
      })
    }, names(longResults))

    Map(function(nm) {
      output[[paste0(nm, "_volcano_hover_info")]] <- shiny::renderUI({
        vhover <- input[[paste0(nm, "_volcano_hover")]]
        res1 <- shiny::nearPoints(longResults[[nm]], vhover, threshold = 5,
                           maxpoints = 1, panelvar1 = "contrast")
        left_pct <- (vhover$x - vhover$domain$left)/(vhover$domain$right - vhover$domain$left)
        top_pct <- (vhover$domain$top - vhover$y)/(vhover$domain$top - vhover$domain$bottom)
        left_px <- vhover$range$left + left_pct * (vhover$range$right - vhover$range$left)
        top_px <- vhover$range$top + top_pct * (vhover$range$bottom - vhover$range$top)
        style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                        "left:", left_px + 2, "px; top:", top_px + 2, "px;")
        shiny::wellPanel(
          style = style,
          shiny::p(shiny::HTML(paste0("<b> ID: </b>", res1$gene, "<br/>",
                                      "<b> Gene: </b>", res1$symbol, "<br/>",
                                      "<b> logFC: </b>", round(res1$logFC, 3), "<br/>",
                                      "<b> PValue: </b>", signif(res1$PValue, 3), "<br/>",
                                      "<b> FDR: </b>", signif(res1$FDR, 3), "<br/>")))
        )
      })
    }, names(longResults))


  }

  shinyApp(ui = pLayout, server = server_function)
}
