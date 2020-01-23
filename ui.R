library(shiny)
library(shinyBS)
library(shinythemes)
library(knitr)
library(rmarkdown)

navbarPage(
  title= div(h4("Correlation Heatmaps",style = "margin-top: 0px;"),
             img(src = "ty_logo.png", height = "50px",
                 style = "position: relative; top: -43px; right: -1000px;")),
  windowTitle = "Correlation Heatmaps",
  theme = shinytheme("sandstone"),
  fluid = TRUE,
  tabPanel("Data Table",
           shinyjs::useShinyjs(),
           sidebarLayout(
             sidebarPanel(
               fileInput('uploadData', 'Choose Data File'),
               actionLink(inputId = "showEx", label = "Example file"),
               uiOutput('idColSelect'),
               uiOutput('columnSelect'),
               shinyBS::bsButton('doBiomart', "Get Gene Symbols", style="primary", size="small", type="action"),
               #shinyBS::bsButton('toggleBMset', "Advanced Biomart Controls", style="primary", size="extra-small", type="toggle", value=FALSE),
               width=2
             ),
             mainPanel(
               bsAlert('alert_anchor1'),
               uiOutput("inputDT"),
               width=10
             )
           )
  ),
  tabPanel("Plot Data",
           sidebarLayout(
             sidebarPanel(
               tags$div(id='noplotdata', tags$label("No data.")),
               radioButtons('chooseSelType', "Selection Criteria", c("By number of genes", "By gene symbol (using pre-selected number of genes)"), "By number of genes"),
               uiOutput('nvarUI'),
               uiOutput('geneSelUI'),
               sliderInput('nSurrGenes', tags$label("Numer of genes surrounding selection", style="font-size: 14px;"), 0, 50, 10, 1),
               shinyBS::bsButton('sortPlotData', "Sort Plot Data", style="info", size="small", type="toggle", value=FALSE),
               shinyBS::bsButton('toPlot', "Open Heatmap", style="primary", size="small", type="action"),
               conditionalPanel("input.toggleBMset", {
                 tags$html(
                   textInput('biomHost', "Biomart Host URL", "www.ensembl.org", placeholder="www.ensembl.org"),
                   textInput('biomFilt', "Biomart Search Filter", "ensembl_gene_id", placeholder="ensembl_gene_id")
                 )
               }),
               width=2
             ),
             mainPanel(
               bsAlert('alert_anchor2'),
               uiOutput("plotDT"),
               width=10
             )
           )
  ),
  tabPanel("Correlation Heatmap",
           sidebarLayout(
             sidebarPanel(
               tags$div(id='nohmdata', tags$label("No data.")),
               checkboxInput('addStars', "Add Significance Stars", FALSE),
               checkboxInput('addRect4genes', "Highlight Selected Genes", FALSE),
               sliderInput('imgSize', tags$label("Image Size (%)", style="font-size: 14px;"), 100, 800, 100, 20),
               sliderInput('textSize', "Gene Label Size", 0.1, 3, 0.8, 0.1),
               textInput('plotTitle', "Enter Plot Title", "", placeholder=""),
               tags$div(id='downlButID', style="text-align: center; ", downloadButton('downloadPlot', tags$label("Download as PDF", style="font-size: 14px;"))),
               tags$br(),
               shinyBS::bsButton('toggleAdv', "Show advanced controls", style="primary", size="small", type="toggle", value=FALSE),
               conditionalPanel("input.toggleAdv", {
                 tags$html(
                   checkboxInput('doClust', "Cluster and Filter Correlations", TRUE),
                   sliderInput('naFrac', "Allowed Fraction of NAs per Row (%)", 0.1, 1, 0.1, 0.05),
                   checkboxInput('doCorFilt', "Filter by Correlation Value", FALSE),
                   conditionalPanel("input.doCorFilt", {
                     tags$html(
                       sliderInput('corThr', "Correlation Filter Threshold", -1, 1, 0.5, 0.05),
                       sliderInput('corMar', "Correlation Filter Margin per Row (%)", 0.001, 1, 0.05, 0.005)
                     )
                   }),
                   checkboxInput('doCutFilt', "Filter by Cutting the Dendrogram", FALSE),
                   conditionalPanel("input.doCutFilt", {
                     tags$html(
                       sliderInput('cutThr', "Threshold for Tree Cutting", 0, 2, 1, 0.1),
                       sliderInput('cutSize', "Number of Genes on a Tree Branch to be Considered a Cluster", 1, 100, 1, 1)
                     )
                   })
                 )
               })
               , width=2
             ),
             mainPanel(
               bsAlert('alert_anchor3'),
               uiOutput('main_plot_ui'),
               width=10)
           )
  ),
  tabPanel("Correlation Matrix",
           tags$div(id='nocmdata', tags$label("No data.")),
           uiOutput('cormatUI')
  ),
  navbarMenu("Help",
             icon = icon("info"),
             tabPanel("About",
                      h4("About"),
                      hr(),
                      includeMarkdown("Markdown/README.md")),
             tabPanel("FAQ",
                      h4("Frequently Asked Quesitons"),
                      hr(),
                      includeMarkdown("Markdown/FAQ.md")),
             tabPanel("Session Information",
                      h4("R session information"),
                      hr(),
                      includeMarkdown("Markdown/RSessionInfo.Rmd"))
             ),
  id='mainNavbarPage', footer=list(tags$hr(), tags$table(style="width: 100%; ", tags$tbody(tags$tr(tags$th(style="text-align: center; ", tags$label(style="font-family: Verdana; font-size: 12pt; font-weight: normal; color: grey; ", "This app has been created and is maintained by the Institute of Biomedicine, University of Turku."))))))
)
