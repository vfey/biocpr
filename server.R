## Required libraries ##

library(shiny)
library(shinyBS)
library(plyr)
library(RColorBrewer)
library(readR)
library(medseqr)
library(heatmapGen2)
library(genefilter)
library(DT)
library(PAMhm)
library(coreheat)
library(digest)
options(shiny.maxRequestSize = 40 * 1024 ^ 2)#increase maximum upload size for fileInput() to 40 MB
#library(shinysky)


## Shiny server logic ##

shinyServer(function(input, output, session) {
  message("\n@@ Starting...\n")
  
  shinyjs::hide('doBiomart')
  shinyjs::hide('imgSize')
  shinyjs::hide('chooseSelType')
  shinyjs::hide('geneSelUI')
  shinyjs::hide('nSurrGenes')
  shinyjs::hide('addRect4genes')
  shinyjs::hide('addStars')
  shinyjs::hide('nvarUI')
  shinyjs::hide('textSize')
  shinyjs::hide('toggleAdv')
  shinyjs::hide('plotTitle')
  shinyjs::hide('sortPlotData')
  shinyjs::hide('toPlot')
  shinyjs::hide('toggleBMset')
  
  raw.data <- reactiveValues()
  plot.data <- reactiveValues()
  plot.data2 <- reactiveValues()
  plot.data3 <- reactiveValues()
  plot.data4 <- reactiveValues()
  adjPlot.data <- reactiveValues()
  cormat0.data <- reactiveValues()
  cormat1.data <- reactiveValues()
  cormat.data <- reactiveValues()
  filt.values.hash <- reactiveValues()
  dat.object.hash <- reactiveValues()
  hideall <- reactiveValues(v = FALSE)
  
  output$idColSelect <- renderUI({
    if (is.null(input$uploadData))
      return()
    ext <- sub(".+(\\.[a-z]{3,4}$)", "\\1", input$uploadData$name)
    if (!file.exists(paste(input$uploadData$datapath, ext, sep = ""))) {
      file.rename(input$uploadData$datapath,
                  paste(input$uploadData$datapath, ext, sep = ""))
    }
    cat("Reading raw data...\n")
    ds <-
      read.to.list(file.path(dirname(input$uploadData$datapath), dir(dirname(
        input$uploadData$datapath
      ))), stringsAsFactors = FALSE)
    ds <- ds[[1]]
    raw.data[["ds"]] <- ds
    id.choices <- c("<row names>", colnames(ds))
    cat("done\n")
    output$nvarUI <- renderUI({
      if (is.null(input$uploadData))
        return()
      if (is.null(input$IDcolName))
        return(NULL)
      if (is.null(input$columnSel))
        return(NULL)
      #									if (input$chooseSelType=="By gene symbol (using pre-selected number of genes)") return(NULL)
      cat("Creating gene # slider...\n")
      sliderInput('nvar',
                  "Number of Genes With Highest Variance",
                  50,
                  nrow(ds),
                  100,
                  50)
    })
    
    tags$html(
      tags$hr(),
      selectInput(
        'IDcolName',
        tags$label("Name of ID column", style = "font-size: 14px;"),
        id.choices,
        id.choices[2],
        FALSE
      )
    )
  })
  
  output$columnSelect <- renderUI({
    if (is.null(input$uploadData))
      return()
    if (is.null(input$IDcolName))
      return()
    ds <- reactiveValuesToList(raw.data)[[1]]
    cat("Generating column selectizer...")
    if (input$IDcolName == "<row names>") {
      cat("<row names>\n")
      dat.cols <- c(1:ncol(ds))
    } else {
      cat("<ID column>\n")
      id.col <- which(colnames(ds) %in% input$IDcolName)
      dat.cols <- c(1:ncol(ds))[-id.col]
    }
    cat("done\n")
    shinyjs::show('doBiomart')
    shinyjs::show('toggleBMset')
    updateButton(session,
                 'sortPlotData',
                 "Sort Plot Data",
                 value = FALSE,
                 style = "info")
    selectizeInput(
      'columnSel',
      label = tags$label("Available Data Columns", style = "font-size: 14px;"),
      choices = names(ds)[dat.cols],
      selected = names(ds)[dat.cols],
      multiple = TRUE
    )
  })
  
  output$inputDT <- renderUI({
    if (is.null(input$uploadData))
      return(NULL)
    if (is.null(input$IDcolName))
      return(NULL)
    if (is.null(input$columnSel))
      return(NULL)
    ds <- reactiveValuesToList(raw.data)[[1]]
    selCols <- which(colnames(ds) %in% input$columnSel)
    if (input$IDcolName == "<row names>" &&
        ncol(ds[selCols]) > length(selCols))
      return(NULL)
    cat("Generating upload overview...\n")
    if (nrow(ds) > 12000) {
      createAlert(
        session,
        anchorId = 'alert_anchor1',
        alertId = "upload_nrow_alert",
        title = "NOTE!",
        content = paste(
          "The data set is very large (",
          nrow(ds),
          " rows). Processing may be slow!",
          sep = ""
        ),
        style = "info",
        append = FALSE
      )
    } else {
      closeAlert(session, 'upload_nrow_alert')
    }
    cat(" (Generating Plot Data...")
    if (input$IDcolName == "<row names>") {
      ds <- data.frame(ds[selCols], stringsAsFactors = FALSE)
    } else {
      id.col <- which(colnames(ds) %in% input$IDcolName)
      ds <-
        as.data.frame(ds[, c(id.col, selCols)], stringsAsFactors = FALSE)
    }
    cat("done)\n")
    if (length(selCols) < 2) {
      createAlert(
        session,
        anchorId = 'alert_anchor1',
        alertId = "hm_ncol_alert",
        title = "Warning!",
        content = "At least 2 data columns are needed for a heatmap! No plot data are generated.",
        style = "warning",
        append = FALSE
      )
    }
    else {
      closeAlert(session, 'hm_ncol_alert')
    }
    plot.data[["ds"]] <- ds
    cat(" Rendering input data...")
    closeAlert(session, 'upload_nrow_alert')
    shiny::withProgress({
      output$uploadTable <-
        DT::renderDataTable(ds, options = list(pageLength = 20))
    }, message = "Rendering data table...", detail = "(Input data)", value =
      0.5)
    shinyjs::show('doBiomart')
    cat("done\n")
    dataTabStyle <-
      "max-height: 1024px; max-width: 99%; border-bottom-left-radius: 1.9%; font-size: 90%; "
    tags$div(DT::dataTableOutput('uploadTable'), style = dataTabStyle)
  })
  
  observeEvent(input$doBiomart, {
    if (is.null(input$uploadData))
      return()
    if (is.null(input$IDcolName))
      return(NULL)
    ds <- reactiveValuesToList(raw.data)
    if (length(ds) == 0)
      return()
    ds <- isolate(ds[[1]])
    cat("Annotating data...\n ")
    createAlert(
      session,
      anchorId = 'alert_anchor1',
      alertId = "biomart_query_alert1",
      title = "NOTE!",
      content = paste(
        "Querying biomart for HGNC Gene Symbols. Depending on your network bandwidth and server workload processing may be slow!\n",
        "You will be automatically redirected to the 'Plot Data' tab.",
        sep = ""
      ),
      style = "info",
      append = FALSE
    )
    createAlert(
      session,
      anchorId = 'alert_anchor2',
      alertId = "biomart_query_alert2",
      title = "NOTE!",
      content = paste(
        "Querying biomart for HGNC Gene Symbols. Depending on your network bandwidth and server workload processing may be slow!",
        sep = ""
      ),
      style = "info",
      append = FALSE
    )
    if (input$IDcolName == "<row names>") {
      id.col <- "row.names"
      hgnc.col <- "ID"
    } else {
      id.col <- hgnc.col <- input$IDcolName
    }
    shiny::withProgress({
      ds <-
        try(convert.bm(
          ds,
          id.col,
          host = input$biomHost,
          biom.filter = input$biomFilt,
          biom.attributes = c("ensembl_gene_id", "hgnc_symbol"),
          rm.dups = FALSE
        ))
    }, message = "Getting information from Biomart...", detail = "This may take a short while...", value =
      0.5)
    if (is(ds, "try-error")) {
      aaid <- switch(
        input$mainNavbarPage,
        "Data Table" = 'alert_anchor1',
        "Plot Data" = 'alert_anchor2',
        "Correlation Heatmap" = 'alert_anchor3'
      )
      aid <- switch(
        input$mainNavbarPage,
        "Data Table" = "bm_fail1",
        "Plot Data" = "bm_fail2",
        "Correlation Heatmap" = "bm_fail3"
      )
      createAlert(
        session,
        anchorId = aaid,
        alertId = aid,
        title = "Warning!",
        content = "The biomart query failed!!'.",
        style = "warning",
        append = FALSE
      )
      return(NULL)
    }
    closeAlert(session, 'biomart_query_alert1')
    closeAlert(session, 'biomart_query_alert2')
    rownames(ds) <- make.unique(ds[, "ensembl_gene_id"])
    ds <- ds[,-which(colnames(ds) %in% "ensembl_gene_id")]
    colnames(ds)[which(colnames(ds) == "hgnc_symbol")] <- hgnc.col
    cat("done\n")
    if (input$mainNavbarPage == "Data Table") {
      updateNavbarPage(session, "mainNavbarPage", "Plot Data")
    }
    #						shinyjs::hide('doBiomart')
    #						shinyjs::hide('toggleBMset')
    #						shinyjs::hide('biomHost')
    #						shinyjs::hide('biomFilt')
    cat("updatting button\n")
    updateButton(session,
                 'sortPlotData',
                 "Unsort Plot Data",
                 style = "success",
                 value = TRUE)
    updateCheckboxInput(session, 'addRect4genes', value = FALSE)
    plot.data[["ds"]] <- ds
  })
  
  observe({
    if (is.null(input$uploadData))
      return()
    if (is.null(input$IDcolName))
      return(NULL)
    if (is.null(input$columnSel))
      return(NULL)
    ds <- reactiveValuesToList(plot.data)
    if (length(ds) == 0)
      return()
    cat("Plot data filtering...")
    ds <- isolate(ds[["ds"]])
    selCols <- which(colnames(ds) %in% input$columnSel)
    if (input$IDcolName == "<row names>") {
      cat("using <row names>...")
      if ("ID" %in% colnames(ds)) {
        cat("and ID column...")
        id.col <- which(colnames(ds) %in% "ID")
        ds <-
          as.data.frame(ds[, c(id.col, selCols)], stringsAsFactors = FALSE)
      } else {
        ds <- data.frame(ds[selCols], stringsAsFactors = FALSE)
      }
    } else {
      id.col <- which(colnames(ds) %in% input$IDcolName)
      ds <-
        as.data.frame(ds[, c(id.col, selCols)], stringsAsFactors = FALSE)
    }
    nvar <- input$nvar
    if (is.null(nvar)) {
      nvar <- 100
    }
    if ((input$IDcolName != "<row names>" &&
         ncol(ds) < 3) || (input$IDcolName == "<row names>" &&
                           ncol(ds) < 2)) {
      d <- ds[, ncol(ds)]
      z <-
        try((d - median(d, na.rm = TRUE)) / (mad(d - median(d, na.rm = TRUE))), silent =
              TRUE)
      if (is(z, "try-error"))
        return(NULL)
      select <-
        order(abs(z), decreasing = TRUE)[seq_len(min(nvar, length(z)))]
    } else {
      if (input$IDcolName == "<row names>") {
        if ("ID" %in% colnames(ds)) {
          Pvars <- try(genefilter::rowVars(ds[, 2:ncol(ds)]), silent = TRUE)
          if (is(Pvars, "try-error"))
            return(NULL)
        } else {
          Pvars <- try(genefilter::rowVars(ds[, 1:ncol(ds)]), silent = TRUE)
          if (is(Pvars, "try-error"))
            return(NULL)
        }
      } else {
        Pvars <- try(genefilter::rowVars(ds[, 2:ncol(ds)]), silent = TRUE)
        if (is(Pvars, "try-error"))
          return(NULL)
      }
      select <-
        order(Pvars, decreasing = TRUE)[seq_len(min(nvar, length(Pvars)))]
    }
    ds.filt <- ds[select, , drop = FALSE]
    plot.data2[["ds2"]] <- ds.filt
    cat("done\n")
  })
  
  observeEvent(input$toPlot, {
    updateNavbarPage(session, "mainNavbarPage", "Correlation Heatmap")
  })
  
  output$geneSelUI <- renderUI({
    if (is.null(input$uploadData))
      return()
    if (is.null(input$IDcolName))
      return()
    if (input$chooseSelType == "By number of genes")
      return()
    ds <- reactiveValuesToList(plot.data2)
    if (length(ds) == 0)
      return()
    ds.filt <- isolate(ds[["ds2"]])
    if (input$IDcolName == "<row names>") {
      if ("ID" %in% colnames(ds.filt)) {
        choices <- ds.filt$ID
      } else {
        choices <- rownames(ds.filt)
      }
    } else {
      choices <- ds.filt[[input$IDcolName]]
    }
    sel <- sort(choices[!choices %in% ""])[1]
    selectInput(
      'geneSel',
      "Select gene to display",
      choices = sort(choices),
      selected = sel,
      TRUE
    )
  })
  
  observe({
    if (is.null(input$uploadData))
      return()
    if (is.null(input$IDcolName))
      return()
    ds <- reactiveValuesToList(plot.data2)
    if (length(ds) == 0)
      return()
    ds.filt <- ds[["ds2"]]
    if (input$chooseSelType == "By number of genes") {
      cat("  Filtering by number of genes...\n")
      updateCheckboxInput(session, 'addRect4genes', value = FALSE)
      select <- 1:nrow(ds.filt)
    } else if (input$chooseSelType == "By gene symbol (using pre-selected number of genes)") {
      cat("  Filtering by gene symbol...\n")
      syms <- input$geneSel
      if (input$IDcolName == "<row names>") {
        if ("ID" %in% colnames(ds.filt)) {
          dat <- as.matrix(ds.filt[,-which(colnames(ds.filt) %in% "ID")])
          lab <- ds.filt["ID"]
        } else {
          dat <- as.matrix(ds.filt)
          lab <-
            data.frame(ID = rownames(ds.filt),
                       stringsAsFactors = FALSE)
          rownames(lab) <- rownames(ds.filt)
        }
      } else {
        dat <-
          as.matrix(ds.filt[,-which(colnames(ds.filt) %in% input$IDcolName)])
        lab <- ds.filt[input$IDcolName]
      }
      # Create a Progress object
      progress <- shiny::Progress$new()
      progress$set(message = "Generating correlation matrix for filtering...", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      updateProgress <- function(value = NULL,
                                 detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = detail)
      }
      cat("    Calling 'cormap2' for gene selection...\n")
      cm0 <- isolate(reactiveValuesToList(cormat0.data))
      if (length(cm0) != 0) {
        cat("  Existing correlation values...\n")
        cm0 <- cm0[[1]]
        if (nrow(dat) != nrow(cm0)) {
          cat("    Input data has changed. Recalculating correlation matrix...\n")
          if (nrow(dat) > 8000) {
            createAlert(
              session,
              anchorId = 'alert_anchor2',
              alertId = "cormat_nrow_alert",
              title = "WARNING!",
              content = paste(
                "Calculating correlation map for",
                nrow(dat),
                "genes. Processing may be slow!"
              ),
              style = "warning",
              append = FALSE
            )
          }
          cormat <-
            cormap2(
              dat,
              biomart = FALSE,
              doPlot = FALSE,
              updateProgress = updateProgress
            )
          cormat0.data[["cm0"]] <- cormat
        } else {
          cat("   Using existing matrix...\n")
          cormat <- cm0
        }
      } else {
        cat("    Calculating correlation matrix for the first time...\n")
        if (nrow(dat) > 8000) {
          createAlert(
            session,
            anchorId = 'alert_anchor2',
            alertId = "cormat_nrow_alert",
            title = "WARNING!",
            content = paste(
              "Calculating correlation map for",
              nrow(dat),
              "genes. Processing may be slow!"
            ),
            style = "warning",
            append = FALSE
          )
        }
        cormat <-
          cormap2(
            dat,
            biomart = FALSE,
            doPlot = FALSE,
            updateProgress = updateProgress
          )
        cormat0.data[["cm0"]] <- cormat
      }
      closeAlert(session, 'cormat_nrow_alert')
      cat("  Starting ID selection...\n")
      cormat <- cormat[nrow(cormat):1,]
      rn <- rownames(cormat)
      symSel <- plyr::llply(syms, function(x) {
        id <- rownames(lab)[lab[, 1] %in% x]
        plyr::llply(id, function(y) {
          rns <- which(rn %in% y)
          if (length(syms) > 2 && input$nSurrGenes == 0) {
            rns
          } else {
            if (rns < input$nSurrGenes) {
              rns1 <- 1:(rns - 1)
            } else {
              rns1 <- (rns - input$nSurrGenes):(rns - 1)
            }
            if (rns > (length(rn) - input$nSurrGenes)) {
              rns2 <- (rns + 1):length(rn)
            } else {
              rns2 <- (rns + 1):(rns + input$nSurrGenes)
            }
            c(rns1, rns, rns2)
          }
        })
      })
      symSel <- unique(unlist(symSel))
      select <- which(rownames(ds.filt) %in% rn[symSel])
    }
    cat("   Setting filtered plot data...\n")
    plot.data3[["ds3"]] <- ds.filt[select, , drop = FALSE]
  })
  
  observe({
    if (is.null(input$uploadData))
      return()
    if (is.null(input$IDcolName))
      return(NULL)
    ds <- reactiveValuesToList(plot.data3)
    if (length(ds) == 0)
      return()
    cat("  Triggering sorting...\n")
    ds.filt <- isolate(ds[["ds3"]])
    if (input$sortPlotData) {
      cat("    Sorting...\n")
      plot.data4[["ds4"]] <-
        ds.filt[order(rownames(ds.filt)), , drop = FALSE]
      updateButton(session, 'sortPlotData', "Unsort Plot Data", style =
                     "success")
    } else {
      cat("    Not sorting\n")
      plot.data4[["ds4"]] <- ds.filt
      updateButton(session, 'sortPlotData', "Sort Plot Data", style = "info")
    }
  })
  
  output$plotDT <- renderUI({
    shinyjs::hide('nvarUI')
    shinyjs::hide('geneSelUI')
    shinyjs::hide('nSurrGenes')
    shinyjs::hide('sortPlotData')
    shinyjs::hide('chooseSelType')
    shinyjs::hide('toPlot')
    if (is.null(input$uploadData))
      return()
    if (is.null(input$IDcolName))
      return(NULL)
    ds <- reactiveValuesToList(plot.data4)
    if (length(ds) == 0)
      return()
    ds.filt <- isolate(ds[["ds4"]])
    if (ncol(ds.filt) < 2)
      return(NULL)
    cat("Rendering plot data...")
    shinyjs::hide('noplotdata')
    if (input$IDcolName != "<row names>" && ncol(ds.filt) < 3) {
      createAlert(
        session,
        anchorId = 'alert_anchor2',
        alertId = "hm_ncol_alert2",
        title = "Warning!",
        content = "At least 2 data columns are needed for a heatmap!",
        style = "warning",
        append = FALSE
      )
    }
    else {
      closeAlert(session, 'hm_ncol_alert2')
    }
    shiny::withProgress({
      output$plotTable <-
        renderDataTable(ds.filt, options = list(pageLength = 20))
    }, message = "Rendering data table...", detail = "(Plot data)", value =
      0.5)
    shinyjs::show('nvarUI')
    shinyjs::show('geneSelUI')
    shinyjs::show('nSurrGenes')
    tabStyle <-
      "max-height: 1024px; max-width: 99%; border-bottom-left-radius: 1.9%; font-size: 90%; "
    cat("done\n\n")
    shinyjs::show('sortPlotData')
    shinyjs::show('chooseSelType')
    shinyjs::show('toPlot')
    tags$div(dataTableOutput('plotTable'), style = tabStyle)
  })
  
  
  
  observeEvent(input$show1, {
    showModal(
      modalDialog(
        HTML(
          '<img src="https://gitlab.utu.fi/dhajam/corplot/blob/master/Markdown/FAQ.md">'
        ),
        title = "Example expression file",
        easyClose = TRUE,
        footer = NULL
      )
    )
  })
  
  
  
  observe({
    if (is.null(input$uploadData))
      return()
    if (input$mainNavbarPage != "Correlation Heatmap")
      return(NULL)
    ds <- reactiveValuesToList(plot.data4)
    if (length(ds) == 0)
      return()
    cat("Getting plot adjustments...\n")
    ds.filt <- isolate(ds[["ds4"]])
    if (input$IDcolName == "<row names>") {
      cat("  IDs in <row names>\n")
      cat("   using <row names>...\n")
      if ("ID" %in% colnames(ds.filt)) {
        cat("  ...and ID column...\n")
        dat <-
          as.matrix(ds.filt[,-which(colnames(ds.filt) %in% "ID"), drop = FALSE])
      } else {
        dat <- as.matrix(ds.filt)
      }
    } else {
      cat("  IDs in ID column\n")
      dat <-
        as.matrix(ds.filt[,-which(colnames(ds.filt) %in% input$IDcolName), drop =
                            FALSE])
    }
    if (ncol(dat) < 2) {
      cat("  break1\n\n")
      hideall[[1]] <- as.logical(TRUE)
      return(NULL)
    }
    if (nrow(dat) == 0) {
      cat("  break2\n\n")
      hideall[[1]] <- as.logical(TRUE)
      return(NULL)
    }
    adj.l <- plotAdjust(dat)
    cex <- round(adj.l$r.cex * 2.1, 1)
    lw <- adj.l$labelwidth * 8 ^ log10(cex)
    adjPlot.data[["pdf.w"]] <- adj.l$pdf.width
    adjPlot.data[["cex"]] <- cex
    adjPlot.data[["lw"]] <- lw
    cat("done\n")
    updateSliderInput(session, 'textSize', value = cex)
  })
  
  observeEvent(input$textSize, {
    if (input$mainNavbarPage != "Correlation Heatmap")
      return(NULL)
    cex <- isolate(reactiveValuesToList(adjPlot.data))[["cex"]]
    if (input$textSize == cex)
      return(NULL)
    cat("Adjusting plot width...\n")
    adj.l.lw <- isolate(reactiveValuesToList(adjPlot.data))[["lw"]]
    labelwidth <- (adj.l.lw) * 8 ^ log10(input$textSize)
    adjPlot.data[["lw2"]] <- labelwidth
  })
  
  output$main_plot_ui <- renderUI({
    shinyjs::hide('toggleAdv')
    shinyjs::hide('imgSize')
    shinyjs::hide('textSize')
    shinyjs::hide('plotTitle')
    shinyjs::hide('chooseSelType')
    shinyjs::hide('addStars')
    shinyjs::hide('addRect4genes')
    shinyjs::hide('downlButID')
    if (is.null(input$uploadData))
      return()
    ds <- reactiveValuesToList(plot.data4)
    if (length(ds) == 0)
      return()
    if (is.null(input$nvar))
      return()
    shinyjs::hide('nohmdata')
    cat("\nPlotting...\n")
    ds.filt <- isolate(ds[["ds4"]])
    if (input$IDcolName == "<row names>") {
      cat("using <row names>...\n")
      if ("ID" %in% colnames(ds.filt)) {
        cat("  ...and ID column...\n")
        dat <-
          as.matrix(ds.filt[,-which(colnames(ds.filt) %in% "ID"), drop = FALSE])
        lab <- ds.filt["ID"]
      } else {
        dat <- as.matrix(ds.filt)
        lab <-
          data.frame(ID = rownames(ds.filt), stringsAsFactors = FALSE)
        rownames(lab) <- rownames(ds.filt)
      }
    } else {
      dat <-
        as.matrix(ds.filt[,-which(colnames(ds.filt) %in% input$IDcolName), drop =
                            FALSE])
      lab <- ds.filt[input$IDcolName]
    }
    if (ncol(dat) < 2) {
      createAlert(
        session,
        anchorId = 'alert_anchor3',
        alertId = "hm_ncol_alert3",
        title = "Warning!",
        content = "At least 2 data columns are needed for a heatmap!",
        style = "warning",
        append = FALSE
      )
      return(NULL)
    } else {
      closeAlert(session, 'hm_ncol_alert3')
    }
    
    # dat hash
    cat("  Generating hash sum of dat object...")
    dat.hash <- isolate(reactiveValuesToList(dat.object.hash))
    dat.hash1 <- digest::sha1(dat)
    cat("done\n")
    if (length(dat.hash) == 0) {
      cat ("   No saved dat hash sum. Saving...\n")
      dat.object.hash[["dh1"]] <- dat.hash1
    } else {
      dat.hash <- isolate(dat.hash[[1]])
    }
    
    # Plot adjustments
    cat("  Loading plot adjustments...\n")
    if (length(reactiveValuesToList(adjPlot.data)[["lw2"]]) == 0) {
      cat("    Using standard adjustments...\n")
      labelwidth <-
        isolate(reactiveValuesToList(adjPlot.data))[["lw"]]
    } else {
      cat("    Using modified adjustments...\n")
      labelwidth <-
        isolate(reactiveValuesToList(adjPlot.data))[["lw2"]]
    }
    width <- height <- 1024 * input$imgSize / 100
    if (isolate(input$imgSize) == 100) {
      hmStyle <- "overflow-y: hidden; overflow-x: hidden; "
    } else {
      hmStyle <- "overflow: scroll; "
    }
    
    # Filters
    corThr <- input$corThr
    cutThr <- input$cutThr
    if (!input$doCorFilt) {
      corThr <- NULL
    }
    if (!input$doCutFilt) {
      cutThr <- NULL
    }
    cat("  Generating filters list...\n")
    filt.hash <- isolate(reactiveValuesToList(filt.values.hash))
    filt.l1 <- list(
      doClust = input$doClust,
      naFrac = input$naFrac,
      corThr = corThr,
      corMar = input$corMar,
      cutThr = cutThr,
      cut.size = input$cutSize
    )
    filt.hash1 <- digest::sha1(filt.l1)
    if (length(filt.hash) == 0) {
      cat ("   No saved filters. Saving...\n")
      filt.values.hash[["fvh"]] <- filt.hash1
    } else {
      filt.hash <- isolate(filt.hash[[1]])
    }
    cat("    Calling 'cormap2' for heatmap...\n     First without plotting...\n")
    # Create a Progress object
    progress0 <- shiny::Progress$new()
    progress0$set(message = "Generating heatmap...", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress0$close())
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress0$getValue()
        value <- value + (progress0$getMax() - value) / 5
      }
      progress0$set(value = value, detail = detail)
    }
    
    # Correlation matrix
    cm1 <- isolate(reactiveValuesToList(cormat1.data))
    if (length(cm1) != 0) {
      cat("  Existing plot correlation values...\n")
      cm1 <- cm1[[1]]
      if (dat.hash != dat.hash1 || filt.hash != filt.hash1) {
        if (dat.hash != dat.hash1) {
          cat("    Plot data has changed.\n")
        }
        if (filt.hash != filt.hash1) {
          cat("    Plot filters were modified.\n")
        }
        cat("     Recalculating plot data correlation matrix and saving changes...\n")
        if (nrow(dat) > 8000) {
          createAlert(
            session,
            anchorId = 'alert_anchor2',
            alertId = "cormat_nrow_alert",
            title = "WARNING!",
            content = paste(
              "Calculating correlation map for",
              nrow(dat),
              "genes. Processing may be slow!"
            ),
            style = "warning",
            append = FALSE
          )
        }
        cormat <-
          cormap2(
            dat,
            biomart = FALSE,
            cluster_correlations = filt.l1$doClust,
            minfrac = filt.l1$naFrac,
            cor.thr = filt.l1$corThr,
            cor.mar = filt.l1$corMar,
            cut.thr = filt.l1$cutThr,
            cut.size = filt.l1$cutSize,
            doPlot = FALSE,
            updateProgress = updateProgress
          )
        cormat1.data[["cm1"]] <- cormat
        dat.object.hash[["dh1"]] <- dat.hash1
        filt.values.hash[["fvh"]] <- filt.hash1
      } else {
        cat("   Using existing plot data correlation matrix...\n")
        cormat <- cm1
      }
    } else {
      cat("    Calculating plot data correlation matrix for the first time...\n")
      if (nrow(dat) > 8000) {
        createAlert(
          session,
          anchorId = 'alert_anchor2',
          alertId = "cormat_nrow_alert",
          title = "WARNING!",
          content = paste(
            "Calculating correlation map for",
            nrow(dat),
            "genes. Processing may be slow!"
          ),
          style = "warning",
          append = FALSE
        )
      }
      cormat <-
        cormap2(
          dat,
          biomart = FALSE,
          cluster_correlations = filt.l1$doClust,
          minfrac = filt.l1$naFrac,
          cor.thr = filt.l1$corThr,
          cor.mar = filt.l1$corMar,
          cut.thr = filt.l1$cutThr,
          cut.size = filt.l1$cutSize,
          doPlot = FALSE,
          updateProgress = updateProgress
        )
      cat("    Saving correlation matrix...")
      cormat1.data[["cm1"]] <- cormat
      cat("done\n")
    }
    closeAlert(session, 'cormat_nrow_alert')
    output$main_plot <- renderPlot({
      cat("  Starting plot...\n")
      # Create a Progress object
      progress <- shiny::Progress$new()
      progress$set(message = "Generating heatmap...", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())
      updateProgress <- function(value = NULL,
                                 detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / 5
        }
        progress$set(value = value, detail = detail)
      }
      cat("    Calling 'cormap2' for heatmap...\n     Now with plotting...\n")
      genes2highl <- NULL
      if (input$addRect4genes &&
          input$chooseSelType == "By gene symbol (using pre-selected number of genes)") {
        cat("   -> Highlighting genes:\n    ")
        genes2highl <- input$geneSel
        print(genes2highl)
      }
      cormat <- cormat[nrow(cormat):1,]
      cm <-
        cormap2(
          dat,
          cormat = cormat,
          lab = lab,
          biomart = FALSE,
          main = input$plotTitle,
          cluster_correlations = FALSE,
          cex = isolate(input$textSize),
          labelheight = labelwidth,
          labelwidth = labelwidth,
          add.sig = input$addStars,
          genes2highl = genes2highl,
          updateProgress = updateProgress
        )
      cat("  done\n\n")
      cormat.data[["cm"]] <- cm
      #									if (is(cm, "try-error")) return(NULL)
    }, width, height)
    shinyjs::show('toggleAdv')
    shinyjs::show('imgSize')
    shinyjs::show('textSize')
    shinyjs::show('plotTitle')
    shinyjs::show('chooseSelType')
    shinyjs::show('addStars')
    shinyjs::show('downlButID')
    if (isolate(input$chooseSelType) == "By gene symbol (using pre-selected number of genes)") {
      shinyjs::show('addRect4genes')
    } else if (isolate(input$chooseSelType) == "By number of genes") {
      shinyjs::hide('addRect4genes')
    }
    tags$div(style = hmStyle, plotOutput('main_plot', height = 1024))
  })
  
  output$cormatUI <- renderUI({
    if (is.null(input$uploadData))
      return()
    if (input$mainNavbarPage != "Correlation Matrix")
      return(NULL)
    cm <- isolate(reactiveValuesToList(cormat.data))
    if (length(cm) == 0)
      return()
    shinyjs::hide('nocmdata')
    cat("Rendering correlation matrix data table...\n")
    cm <- isolate(cm[["cm"]])
    tabStyle <-
      "max-height: 1024px; max-width: 99%; border-bottom-left-radius: 1.9%; font-size: 90%; overflow-x: scroll; "
    shiny::withProgress({
      output$cormatDT <-
        renderDataTable(
          format(
            round(cm, 2),
            width = 12,
            digits = 1,
            scientific = FALSE
          ),
          options = list(pageLength = 20)
        )
    }, message = "Rendering data table...", detail = "(Correlation matrix)", value =
      0.5)
    tags$div(dataTableOutput('cormatDT'), style = tabStyle)
  })
  
  output$downloadPlot <-
    downloadHandler(
      filename = "correlationHeatmap.pdf",
      content = function(file) {
        if (is.null(input$uploadData))
          return()
        if (input$mainNavbarPage != "Correlation Heatmap")
          return(NULL)
        ds <-
          reactiveValuesToList(plot.data4)
        if (length(ds) == 0)
          return()
        cat("DOWNLOADING PLOT\n")
        ds.filt <-
          isolate(ds[["ds4"]])
        if (input$IDcolName == "<row names>") {
          cat("  using <row names>...\n")
          if ("ID" %in% colnames(ds.filt)) {
            cat("  ...and ID column...\n")
            dat <-
              as.matrix(ds.filt[,-which(colnames(ds.filt) %in% "ID")])
          } else {
            dat <- as.matrix(ds.filt)
          }
        } else {
          dat <-
            as.matrix(ds.filt[,-which(colnames(ds.filt) %in% input$IDcolName)])
        }
        genes2highl <- NULL
        if (input$addRect4genes &&
            input$chooseSelType == "By gene symbol (using pre-selected number of genes)") {
          genes2highl <- input$geneSel
        }
        cm <-
          isolate(reactiveValuesToList(cormat.data))
        if (length(cm) == 0)
          return()
        cm <-
          isolate(cm[["cm"]][nrow(cm[["cm"]]):1,])
        cat("  Loading plot adjustments...\n")
        if (length(reactiveValuesToList(adjPlot.data)[["lw2"]]) ==
            0) {
          cat("    Using standard adjustments...\n")
          labelwidth <-
            isolate(reactiveValuesToList(adjPlot.data))[["lw"]]
        } else {
          cat("    Using modified adjustments...\n")
          labelwidth <-
            isolate(reactiveValuesToList(adjPlot.data))[["lw2"]]
        }
        corThr <-
          isolate(input$corThr)
        cutThr <-
          isolate(input$cutThr)
        if (!isolate(input$doCorFilt)) {
          corThr <- NULL
        }
        if (!isolate(input$doCutFilt)) {
          cutThr <- NULL
        }
        pdf.width <-
          min(18, isolate(reactiveValuesToList(adjPlot.data)[["pdf.w"]]))
        shiny::withProgress({
          pdf(file, width = pdf.width, height = pdf.width)
          cat("    Calling 'cormap2' for download...\n")
          invisible(
            cormap2(
              dat,
              cormat = cm,
              lab = NULL,
              biomart = FALSE,
              main = isolate(input$plotTitle),
              cluster_correlations = FALSE,
              cex = isolate(input$textSize),
              minfrac =
                isolate(input$naFrac),
              cor.thr = corThr,
              cor.mar = isolate(input$corMar),
              cut.thr = cutThr,
              cut.size = isolate(input$cutSize),
              labelheight = labelwidth,
              labelwidth =
                labelwidth,
              add.sig = isolate(input$addStars),
              genes2highl = genes2highl
            )
          )
          cat("  done\n\n")
          dev.off()
        }, message = "Creating PDF for download...", detail =
          "(Correlation Heatmap)", value = 0.5)
      },
      contentType = NULL
    )
  
  observe({
    cat("Triggering control selector...\n")
    hide <- reactiveValuesToList(hideall)[[1]]
    if (hide) {
      shinyjs::hide('chooseSelType')
      shinyjs::hide('geneSelUI')
      shinyjs::hide('nSurrGenes')
      shinyjs::hide('sortPlotData')
      shinyjs::hide('toPlot')
      shinyjs::hide('nvarUI')
      shinyjs::hide('IDcolName')
      shinyjs::hide('columnSel')
      shinyjs::hide('addStars')
      shinyjs::hide('toggleAdv')
      shinyjs::hide('imgSize')
      shinyjs::hide('textSize')
      shinyjs::hide('plotTitle')
    }
    if (input$mainNavbarPage == "Correlation Heatmap" &&
        !is.null(input$uploadData) && !is.null(input$nvar)) {
      #							if (input$sortPlotData && input$doBiomart!=0) {
      #								shinyjs::hide('doBiomart')
      #								shinyjs::hide('toggleBMset')
      #								shinyjs::hide('biomHost')
      #								shinyjs::hide('biomFilt')
      #							}
      shinyjs::show('chooseSelType')
      shinyjs::show('geneSelUI')
      shinyjs::show('nSurrGenes')
      shinyjs::show('addStars')
      shinyjs::show('sortPlotData')
      shinyjs::show('toggleAdv')
      shinyjs::show('imgSize')
      shinyjs::show('nvarUI')
      shinyjs::show('textSize')
      shinyjs::show('plotTitle')
      #							shinyjs::hide('toPlot')
    }
    if (input$mainNavbarPage == "Plot Data" &&
        !is.null(input$uploadData)) {
      shinyjs::show('chooseSelType')
      shinyjs::show('geneSelUI')
      shinyjs::show('nSurrGenes')
      shinyjs::show('sortPlotData')
      shinyjs::show('toPlot')
      shinyjs::show('nvarUI')
      shinyjs::show('IDcolName')
      shinyjs::show('columnSel')
      shinyjs::hide('addStars')
      #							if (input$sortPlotData && input$doBiomart!=0) {
      #								shinyjs::hide('doBiomart')
      #								shinyjs::hide('toggleBMset')
      #								shinyjs::hide('biomHost')
      #								shinyjs::hide('biomFilt')
      #							}
    }
    if (input$mainNavbarPage == "Data Table" &&
        !is.null(input$uploadData)) {
      shinyjs::show('IDcolName')
      shinyjs::show('columnSel')
      #							shinyjs::hide('chooseSelType')
      #							shinyjs::hide('geneSelUI')
      #							shinyjs::hide('addStars')
      #							shinyjs::hide('toPlot')
      shinyjs::show('toggleBMset')
    }
  })
})
