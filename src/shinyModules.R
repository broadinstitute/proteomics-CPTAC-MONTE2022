###############################################################################
## Filename: shinyModules.R
## Created: July 6, 2022
## Author: Stephanie Vartany
## Purpose: Data viewer for "MONTE enables serial immunopeptidome, ubiquitylome,
## proteome, phosphoproteome, acetylome analyses of sample-limited tissues"
## This file contains functions for the user interface and servers
## Note: code is heavily adapted from Karsten Krug's previous shiny apps
###############################################################################

## User interface for each tab in the data viewer
viewerTabUI <- function(id, label = "Viewer Tab", params) {

  ns <- NS(id) # namespace function
  
  tagList(
    ## UI for heatmap
    h3("Heatmap"),
    fluidRow(column(3, wellPanel(
     
     ## text input
     selectizeInput(ns('genes'), 
                    label=paste('Enter your genes of interest (max. ',
                                GENEMAX,')', sep=''), 
                    selected=params$genes.start,
                    choices = NULL,
                    multiple=T),
     
     br(),
     
     ## inputs to customize heatmap
     fluidRow(
       column(6, radioButtons(ns('zscore'), 
                              label='Z-score', 
                              choices=c('row', 'none'), 
                              selected='row')),
       
       column(6, 
              radioButtons(ns('PTMsites'), 
                           label='PTM sites', 
                           choices=c('most variable', 'all'), 
                           selected='most variable'))
       ), # end fluidRow
     
     fluidRow(
       column(6, textInput(ns('min.val'), 
                           label='min', 
                           value=-2, 
                           width='80%')),
       column(6, textInput(ns('max.val'), 
                           label='max', 
                           value=2, 
                           width='80%'))
     ), #end fluidRow
     
     ## inputs for sorting
     fluidRow(
       column(12, selectizeInput(ns('sort.after'), 'Sort by', 
                                 choices= params$annotations,
                                 selected=params$annotations.start, 
                                 multiple=FALSE))  
     ), #end fluidRow
     
     br(),
     
     # ## update plot
     # fluidRow(column(12,
     #                 actionButton(ns("update"), "Update Heatmap",
     #                 icon("paper-plane"), 
     #                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))),
     # 
     # br(),
     
     ## download buttons
     fluidRow(column(6, downloadButton(ns('downloadHM'), 'Download PDF')),
              column(6, downloadButton(ns('downloadTab'), 'Download Excel'))),
     
     br(),
     
     ## Instructions
     p(strong("Getting started")),
     p("Enter your gene names of interest (official gene symbols, e.g. EGFR) into the text field. You can enter up to 20 genes."),
     br(),
     p(strong("Dataset")),
     HTML(params$data_description),
     
    ) #end wellPanel
    ), #end column
    
    ## Heatmap
    column(9, plotOutput(ns("hm"),
                         height = 'auto',
                         width = 'auto'),
    ) # end column
    ), #end fluidRow
    
    
    br(),
    
    
################ COMMENTED OUT FOR NOW
    # ## UI for scatter plot
    # h3("Scatter Plot"),
    # fluidRow(column(3, wellPanel(
    #   
    #   ## pick gene for x and y axis
    #   fluidRow(
    #   column(6, 
    #          ## text input
    #          selectizeInput(ns('scatter.gene.x'), 
    #                         label="Enter x-axis gene", 
    #                         choices = NULL,
    #                         multiple=F)
    #          ),
    #   column(6, 
    #          ## text input
    #          selectizeInput(ns('scatter.gene.y'), 
    #                         label="Enter y-axis gene", 
    #                         choices = NULL,
    #                         multiple=F)
    #   )), #end fluidRow and column
    #  
    #  ## pick ome for x and y axis
    #  fluidRow(
    #    column(6, selectizeInput(ns('scatter.x'), 
    #                             label="x-axis", 
    #                             choices=NULL,
    #                             multiple=F)),
    #    column(6, selectizeInput(ns('scatter.y'), 
    #                             label="y-axis", 
    #                             choices=NULL,
    #                             multiple=F),)
    #  ), #end fluidRow
    #  
    #  # fluidRow(column(12,
    #  #                 actionButton(ns("scatter.update"), "Update Plot",
    #  #                 icon("paper-plane"), 
    #  #                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))),
    #  
    #  ## Instructions
    #  br(),
    #  p(strong("Getting started")),
    #  p("Enter your gene name of interest (official gene symbol, e.g. PIK3CA) into the text field. Then, pick two -omes to create a 2D scatter plot"),
    #  
    # ) #end wellPanel
    # ), #end column
    # 
    # ## scatter plot
    # column(5, plotOutput(ns("scatter"))),
    # column(4, "There will be a table with HLA sequences here")
    # ) #end fluidRow

################

  ) #end tagList
    
}

## server for each tab in the data viewer
viewerTabServer <- function(id, table, params) {
  moduleServer(
    id,
    ## module function
    function(input, output, session) {
      
      # update selectizeInput in ui
      updateSelectizeInput(session,
                           inputId = "genes",
                           choices = params$genes.all,
                           selected=params$genes.start,
                           server=T)
      
      ## get heatmap parameters
      HM.params <- reactive({list(genes.char = input$genes,
                           zscore = input$zscore,
                           PTMsites = input$PTMsites,
                           min.val = as.numeric(input$min.val),
                           max.val = as.numeric(input$max.val),
                           sort.after = input$sort.after,
                           id = id)})
      
      ## generate complex heatmap
      HM.out <- reactive({
          validate(
            need(HM.params()$genes.char, "Input genes to see results"),
            need(HM.params()$min.val < HM.params()$max.val, "Input valid min and max")
          )
          myComplexHeatmap(table, HM.params())
        })
      
      ## render complex heatmap
      observe({
        output$hm <- renderPlot({
          validate(
            need(HM.params()$genes.char, "Input genes to see results"),
            need(HM.params()$min.val < HM.params()$max.val, "Input valid min and max"),
          )
          draw(HM.out()$HM, annotation_legend_side='bottom')
        },
        height = dynamicHeightHM(nrow(HM.out()$Table)))
      })
      
      ## download HM pdf
      output$downloadHM <- downloadHandler(
        filename = function() {paste(FILENAMESTRING,
                                     '-',
                                     gsub(' |\\:','-', Sys.time()),
                                     '.pdf',
                                     sep='')},
        content = function(file) {
          pdf(file = file,
              width = 1400/72,
              height = (dynamicHeightHM(nrow(HM.out()$Table))+48)/72)
          draw(HM.out()$HM, annotation_legend_side='bottom')
          dev.off()
        }
      )
      
      ## download excel
      output$downloadTab <- downloadHandler(
        filename = function(){paste( FILENAMESTRING,
                                     '-',
                                     gsub(' |\\:','-', Sys.time()),
                                     '.xlsx',
                                     sep='')},
        content = function(file) {
          tab = HM.out()$Table
          WriteXLS('tab',
                   ExcelFileName=file,
                   SheetNames=FILENAMESTRING,
                   FreezeCol=4,
                   row.names=T)
        }
      )
      
########### COMMENTED OUT FOR NOW 
      # ## scatter plot
      # updateSelectizeInput(session,
      #                      inputId = "scatter.gene.x",
      #                      choices = params$genes.all,
      #                      selected=params$genes.start[1],
      #                      server=T)
      # updateSelectizeInput(session,
      #                      inputId = "scatter.gene.y",
      #                      choices = params$genes.all,
      #                      selected=params$genes.start[1],
      #                      server=T)
      # observe({
      #   updateSelectizeInput(session,
      #                        inputId = "scatter.x",
      #                        choices = getGeneRowLabels(input$scatter.gene.x, table),
      #                        server = T)
      # })
      # observe({
      #   updateSelectizeInput(session,
      #                        inputId = "scatter.y",
      #                        choices = getGeneRowLabels(input$scatter.gene.y, table),
      #                        server = T)
      # })
      # 
      # scatter.params <- reactive({list(gene.x = input$scatter.gene.x,
      #                                  gene.y = input$scatter.gene.y,
      #                                  x = input$scatter.x,
      #                                  y = input$scatter.y)})
      # observe({
      #   output$scatter <- renderPlot({
      #     validate(
      #       need(scatter.params()$gene.x, "Nothing to show"),
      #       need(scatter.params()$gene.y, "Nothing to show"),
      #       need(scatter.params()$x, "Nothing to show"),
      #       need(scatter.params()$y, "Nothing to show")
      #     )
      #     myScatterPlot(table, scatter.params())
      #   })})
##########################################
      
      
    }
  )
}

## user interface for HLA table in MONTE tab
HLATableUI <- function(id, label = "HLA Table", params) {
  ns <- NS(id) # namespace function
  
  tagList(
    ## UI for heatmap
    h3("HLA Table"),
    fluidRow(column(3, wellPanel(
      
      ## text input
      selectizeInput(ns('hla.gene'), 
                     label="Enter your gene of interest", 
                     selected=params$genes.start[1],
                     choices = NULL,
                     multiple=F),
      
      br(),
      
      ## download button
      fluidRow(column(12, downloadButton(ns('downloadHLA'), 'Download HLA Tables'))),
      
      br(),
      
      ## Instructions
      p(strong("Getting started")),
      p("Enter your gene name of interest (official gene symbols, e.g. PIK3CA) to view associated HLA class 1 and class 2 sequences."),
      br()
    ) #end wellPanel
    ), #end column
    
    column(9,
           fluidRow(fluidRow(h4("HLA Class 1 Sequences")), 
                    fluidRow(tableOutput(ns("hla.cls1")))),
           fluidRow(fluidRow(h4("HLA Class 2 Sequences")), 
                    fluidRow(tableOutput(ns("hla.cls2")))),
           style='padding-left:25px; padding-right:5px; padding-top:0px; padding-bottom:0px'
    ) # end column
    ) # end fluidRow
  )
}

## server for HLA table
HLATableServer <- function(id, hla.table, params) {
  moduleServer(
    id,
    ## module function
    function(input, output, session) {
      # update selectizeInput in ui
      updateSelectizeInput(session,
                           inputId = "hla.gene",
                           choices = params$genes.all,
                           selected=params$genes.start[1],
                           server=T)
      
      HLA.out <- reactive({makeHLATables(input$hla.gene, hla.table)})
      
      output$hla.cls1 <- renderTable({
        validate(
          need(nrow(HLA.out()$cls1) > 0, 
               "No HLA class 1 sequences found for inputted gene")
        )
        HLA.out()$cls1
      },
      striped = T,
      bordered = T)
      
      output$hla.cls2 <- renderTable({
        validate(
          need(nrow(HLA.out()$cls2) > 0, 
               "No HLA class 2 sequences found for inputted gene")
        )
        HLA.out()$cls2
      },
                                     striped = T,
                                     bordered = T)
      
      ## download class 1
      output$downloadHLA <- downloadHandler(
        filename = function(){paste( FILENAMESTRING,
                                     '-HLA-',
                                     input$hla.gene,
                                     '-',
                                     gsub(' |\\:','-', Sys.time()),
                                     '.xlsx',
                                     sep='')},
        content = function(file) {
          cls1.table <- HLA.out()$cls1
          cls2.table <- HLA.out()$cls2
          
          # convert \U2713 or ' ' to 0/1
          cls1.table[cls1.table == '\U2713'] <- '1'
          cls1.table[cls1.table == ' '] <- '0'
          cls2.table[cls2.table == '\U2713'] <- '1'
          cls2.table[cls2.table == ' '] <- '0'
          
          
          WriteXLS(c('cls1.table', 'cls2.table'),
                   ExcelFileName = file,
                   SheetNames = c(paste(input$hla.gene, 'Class 1 Sequences', sep= ' '),
                                  paste(input$hla.gene, 'Class 2 Sequences', sep= ' ')),
                   row.names = F)
        }
      )
    })
}