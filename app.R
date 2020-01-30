#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
## Get installed packages
#
# installed.packages()[, "Package"]
pkg <- installed.packages()[, "Package"]


# Packages used
usedpkg <-
  c(
    "MALDIquant",
    "MALDIquantForeign",
    "xlsx",
    "stringr",
    "shiny",
    "shinyWidgets",
    "plotly",
    "DT"
    
  )

#See if packages is installed an load package

#Client side
# for (pack in usedpkg) {
#   if (!(pack %in% pkg)) {
#     install.packages(pack)
#   }
#   library(pack, character.only = T)
# }



#Shinyio side
library("MALDIquant")
library("MALDIquantForeign")
library("xlsx")
library("stringr")
library("shiny")
library("shinyWidgets")
library("plotly")
library("DT")
#mzR



# Define UI for application that draws a histogram
ui <- fluidPage(
  
  setBackgroundImage(src = "https://cdn.shopify.com/s/files/1/0068/8216/4806/articles/article_6_Health_potential_of_whey_1080x1620.jpg?v=1563879886"),
  
  tags$head(
    tags$style(HTML('li>a{background-color:white}')),
    tags$style(HTML('form{background-color:steelblue !important; border-color:rgb(60,110,160) !important}')),
    tags$style(HTML('form>div>label{color:white !important; font-size:20px}')),
    tags$style(HTML('h1{font-family: Impact, Charcoal, sans-serif;}')),
    tags$style(HTML('.shiny-datatable-output{background-color:rgb(221,221,221) !important;}')),
    tags$style(HTML('div.dataTables_wrapper{background-color:white !important;}'))
  ),
  
  # Application title
    tags$div(tags$h1(
    tags$img(src = "https://www.howest.be/sites/default/files/styles/width_500/public/howest-university-of-applied-sciences-logo.png?itok=scypq5Z0",
             width = "15%", style = "z-index:5"),
    "MS-PEP data analyzer"
  )),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      
      fileInput(inputId = "data", label = "Raw data:"),
      
      #Test sequence: ARTRARRPESKATNATLDPRSFLLRNPNDK
      textInput("seq", "Sequence:", value = "", placeholder = "Enter sequence here"),
      
      actionButton("Submit", "Submit", value = 0),
      
      tags$br(),
      tags$br(),
      
      dropdown(label = "Advanced Settings",
        sliderInput(
          "Relintco",
          "Relative Intensity Cutoff:",
          min = 0,
          max = 100,
          value = 20
        ),
        sliderInput(
          "errco",
          "Maximum error between theoretical and practical:",
          min = 0,
          max = 1,
          value = 0.3,
          step = 0.01
        ),
        checkboxInput("SNM", "Show not matching sequences"),
        checkboxInput("Z", "Add z-value of 2"),
        checkboxInput(
          "proteolitic",
          "Is proteolytically cleaved? (adds OH to b type and subtracts H from y type fragments)",
          value = T
        ),
        #numericInput("Nterm","N-term modification",0), #[WIP]
        numericInput("Cterm", "C-term modification", -1),
        radioButtons("RelInt", label = "Relative intensity reference", choices = c("Highest peak", "Selected peak"),selected ="Highest peak")
        
      ),
      tags$br(),
      downloadButton("Exp", "Export to excel"),
      tags$br(),
      tags$br(),
      htmlOutput("<h1>selPeak</h1>"),
      dataTableOutput("Selected")
      
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        type = "tabs",
        tabPanel(title = "MS-Plot", plotOutput("peakPlot")),
        tabPanel("Theoretical fragments", dataTableOutput("TFrag")),
        tabPanel("Practical fragments", DT::dataTableOutput("PFrag")),
        tabPanel("Visualization", plotlyOutput("plly"))
      )
    )
  ),
  tags$div(
    tags$a(
      "Howest bioinformatics",
      href = "https://www.howest.be/en/programmes/advanced-bachelor/bioinformatics",
      style = "color:white",
      target = "_blank"
    ),
    "|" ,
    tags$a("bit@bio-informatica.be", href = "mailto:bit@bio-informatica.be", style =
             "color:white"),
    style = "position: fixed;
    left: 0;
    bottom: 0;
    width: 100%;
    height: 2em;
    background-color: darkgray;
    color: white;
    text-align: center;
    align-items: center"
  )
)
# Define server logic required to draw a histogram
server <- function(input, output, session) {
  #Extra functions needed
  source("mzR-functions.R")
  
  #Z value selector, changes on input
  Zval <- reactiveVal(1)
  reactive({
    #If ticked: both z == 1 and z == 2, else only 1
    if (input$Z == T) {
      Zval <- 1:2
    }
    else {
      Zval <- 1
    }
  })
  
  # "Object" including info about file submitted
  File <-
    eventReactive(input$Submit, {
      list(
        name = input$data$name,
        location = input$data$datapath,
        basefile = strsplit(input$data$name , ".txt")[1]
      )
    })
  # First spectogram of file (always since only txt import is implemented)
  pepspec1 <- reactive(importTxt(File()$location)[[1]])
  
  # peaks from pepspec1, Signal to noise 10
  pepspec1_peaks_30_10 <- reactive(detectPeaks(
    pepspec1(),
    method = "SuperSmoother",
    halfWindowSize = 30,
    SNR = 10
  ))
  
  # m/z and intensity dataframe of input file. including relative intensity
  peaks.df.reac <- reactive({
    sub_spec1 <- pepspec1()#[1:20000]
    sub_spec_peaks_30_10 <-
      detectPeaks(
        sub_spec1,
        method = "SuperSmoother",
        halfWindowSize = 30,
        SNR = 10
      )
    inter <- data.frame(
      #Specifying packages because other package with mz function
      MALDIquant::mz(sub_spec_peaks_30_10),
      MALDIquant::intensity(sub_spec_peaks_30_10)
    )
    colnames(inter) <- c("m/z", "int")
    
    ## relative intensity
    inter$rel.int <-
      if (input$RelInt == "Highest peak") {
        inter$int / max(inter$int) * 100
      } else {
        if (!is.null(input$rows[1])){
          inter$int / as.numeric(input$rows[3]) * 100
          
        }else {
          inter$int / max(inter$int) * 100
        }
      }
      
    
    inter
  })
  
  #Input sequence (same as input$seq)
  mypeptide <- reactive(input$seq)
  
  # Function to calculate theorethical fragments
  theor_frag_calc <- function(mypeptide) {
    theor_frag.df <-
      tail(
        calculateFragments(
          paste0("A", mypeptide),
          type = c("y"),
          # default
          z = 1:2,
          # default
          modifications = c(Cterm = input$Cterm, Nterm = input$Nterm),
          neutralLoss = NULL
        ),
        2
      )
    theor_frag.df$ion <- ""
    theor_frag.df$type <- ""
    theor_frag.df <-
      rbind(theor_frag.df,
            calculateFragments(
              mypeptide,
              type = c("b", "y"),
              # default
              z = 1:2,
              # default
              modifications = c(Cterm = input$Cterm, Nterm = input$Nterm)
            ))
    ## Filter only b type
    # theor_frag.df[theor_frag.df$type=="b",]
    ## Filter only y type
    # theor_frag.df[theor_frag.df$type=="y",]
    
    ## Add internal fragments
    for (n in (nchar(mypeptide)):2) {
      theor_frag.df <-
        rbind(theor_frag.df,
              calculateFragments(
                substr(mypeptide, 1, n),
                type = c("b", "y"),
                z = 1:2
              ))
    }
    if (input$proteolitic) {
      theor_frag.df[which(theor_frag.df$type == "b"), "mz"] <-
        theor_frag.df[which(theor_frag.df$type == "b"), "mz"] + (15.994915 + 1.007825)
      theor_frag.df[which(theor_frag.df$type == "y"), "mz"] <-
        theor_frag.df[which(theor_frag.df$type == "y"), "mz"] - (1.007825)
    }
    theor_frag.df
  }
  
  #Actually calculate the fragments and only keep unique fragments
  theor_frag.df.react <-
    reactive(unique(theor_frag_calc(mypeptide())))
  
  # expanded dataframe of peaks.df
  peaks.df2.react <- reactive({
    peaks.df2 <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("m/z", "int", "rel.int"))
    for (i in 1:nrow(peaks.df.reac())) {
      current <- peaks.df.reac()[i, "m/z"]
      #High default large error
      seq <- ""
      modification <- ""
      found <- FALSE
      
      for (j in which(sapply(theor_frag.df.react()$mz,FUN = function(x) abs(current - x)) < 1)) {
        #for (j in 1:nrow(theor_frag.df)) {
        #if (abs(current-theor_frag.df[j,"mz"])<1) {
        found <- TRUE
        
        seq <- theor_frag.df.react()[j, "seq"]
        zval <- theor_frag.df.react()[j, "z"]
        # * in ion column stands for ammonia fragment
        # _ in ion column stands for water fragment
        # Can't have both at same time
        if (grepl("\\*", theor_frag.df.react()[j, "ion"])) {
          modification <- "[1xAmonia]"
        } else if (grepl("\\_", theor_frag.df.react()[j, "ion"])) {
          modification <- "[1xWater]"
        }
        
        peaks.df2[nrow(peaks.df2) + 1, "m/z"] <-
          current
        peaks.df2[nrow(peaks.df2), "int"] <-
          peaks.df.reac()[i, "int"]
        peaks.df2[nrow(peaks.df2), "rel.int"] <-
          peaks.df.reac()[i, "rel.int"]
        peaks.df2[nrow(peaks.df2), "sequence"] <- seq
        peaks.df2[nrow(peaks.df2), "Z"] <- zval
        peaks.df2[nrow(peaks.df2), "modification"] <-
          modification
        peaks.df2[nrow(peaks.df2), "error"] <-
          round(abs(current - theor_frag.df.react()[j, "mz"]), 6)
        #}
        
        
      }
      if (!found) {
        peaks.df2[nrow(peaks.df2) + 1, "m/z"] <- peaks.df.reac()[i, "m/z"]
        peaks.df2[nrow(peaks.df2), "int"] <-
          peaks.df.reac()[i, "int"]
        peaks.df2[nrow(peaks.df2), "rel.int"] <-
          peaks.df.reac()[i, "rel.int"]
        peaks.df2[nrow(peaks.df2), "sequence"] <- ""
        peaks.df2[nrow(peaks.df2), "Z"] <- ""
        peaks.df2[nrow(peaks.df2), "modification"] <-
          ""
        peaks.df2[nrow(peaks.df2), "error"] <- ""
      }
    }
    peaks.df2 <-
      peaks.df2[order(peaks.df2$error, decreasing = F),]
      unique(peaks.df2[order(peaks.df2$rel.int, decreasing = T),])
  }
  )
  
  # Filter for export
  peaks.df.filtered.react <- reactive({
    inter <- peaks.df2.react()[which(
      peaks.df2.react()$sequence != "" &
        peaks.df2.react()$rel.int >= input$Relintco &
        peaks.df2.react()$error <= input$errco &
        peaks.df2.react()$Z %in% (if (input$Z) {
          c("", 1, 2)
        } else{
          c("", 1)
        })
    ),]
    # Define color gradient
    cols <- heat.colors(1000)[1000:1]
    inter$color <- cols[round(inter$rel.int)]
    
    for (n in 1:nrow(inter)) {
      inter[n, "start"] <-
        str_locate(mypeptide(), inter[n, "sequence"])[1]
      inter[n, "stop"] <-
        str_locate(mypeptide(), inter[n, "sequence"])[2]
    }
    inter
  }
  
  )
  
  #This includes peaks without match
  peaks.df3.react <- reactive({
    peaks.df3 <- peaks.df2.react()
    peaks.df3$sequence <-
      ifelse(peaks.df3$error > input$errco, "", peaks.df3$sequence)
    peaks.df3$modification <-
      ifelse(peaks.df3$error > input$errco,
             "",
             peaks.df3$modification)
    peaks.df3$Z <-
      ifelse(peaks.df3$error > input$errco, "", peaks.df3$Z)
    peaks.df3$error <-
      ifelse(peaks.df3$error > input$errco,
             paste0(peaks.df3$error, "*"),
             peaks.df3$error)
    peaks.df3
  })
  
  #When pressing submit, following will be executed
  
  observeEvent(input$Submit, {
    withProgress(message = 'Calculating',
                 value = 0,
                 detail = "Reading data",
                 {
                   start <- proc.time()
                   if (is.null(File()) | input$Submit == 0) {
                     return(NULL)
                   }
                   
                   incProgress(0.2, detail = "Detecting Peaks")
                   
                   
                   
                   peaks.df <- peaks.df.reac()
                   
                   
                   
                   ## Show mass values of peaks
                   # MALDIquant::mass(pepspec1_peaks_30_10)
                   
                   ## Show intensity values of peaks
                   # MALDIquant::intensity(pepspec1_peaks_30_10)
                   
                   incProgress(0.2, detail = "Calculating theorethical fragments")
                   
                   
                   
                   theor_frag.df <- theor_frag.df.react()
                   incProgress(0.3, detail = "Matching theoretical with practical")
                   
                   
                   peaks.df2<-reactive(peaks.df2.react())
                   print(proc.time() - start)
                   
                   
                   #incProgress(0.2, detail = "Generating Excel")
                   
                   incProgress(0.3, detail = "Finishing up")
                 })
  })
  
  
  
  #MS-plot output
  
  output$peakPlot <- renderPlot({
    if (input$Submit != 0) {
      ## Create plots
      par(mfrow = c(1, 1))
      # labels 5 highest peaks
      top5 <-
        MALDIquant::intensity(pepspec1_peaks_30_10()) %in% sort(MALDIquant::intensity(pepspec1_peaks_30_10()),
                                                                decreasing = TRUE)[1:5]
      ylim_max = max(MALDIquant::intensity(pepspec1_peaks_30_10())) + max(MALDIquant::intensity(pepspec1_peaks_30_10())) *
        0.15
      
      xlim_range <- range(MALDIquant::mass(pepspec1()))
      plot(
        pepspec1(),
        main = paste0("Peak detection on ", File()$name),
        sub = "halfWinSize = 20, SNR = 2",
        xlim = xlim_range,
        ylim = c(0, ylim_max)
      )
      points(pepspec1_peaks_30_10())
      noise <- MALDIquant::estimateNoise(pepspec1())
      lines(noise, col = "red")
      lines(noise[, 1], noise[, 2] * 2, col = "green")
      labelPeaks(pepspec1_peaks_30_10(), index = top5)
      
    }
  })
  
  #Output theoretical fragments
  output$TFrag <- renderDataTable({
    if(mypeptide() != "") {
      theor_frag.df <- theor_frag.df.react()
      if (input$Submit != 0) {
        theor_frag.df[which(theor_frag.df$z %in% (if (input$Z) {
          c("", 1, 2)
        } else{
          c("", 1)
        })),]
      }
    }
  })
  
  #Output practical fragments
  output$PFrag <- DT::renderDataTable({
    if (exists("peaks.df2.react") & input$Submit != 0) {
      peaks.df2 <- peaks.df2.react()
      test <<- peaks.df2
      peaks.df3 <- peaks.df3.react()
      if (input$SNM) {
        DT::datatable(unique(peaks.df3[which(peaks.df3$rel.int >= input$Relintco &
                                 peaks.df3$Z %in% (if (input$Z) {
                                   c("", 1, 2)
                                 } else{
                                   c("", 1)
                                 })),]),callback = JS("table.on('click.dt', 'tr', function() {
            var row_=table.row(this).data();
            var data = [row_];
           Shiny.onInputChange('rows',data );
    });"))
        
      }
      else {
        DT::datatable(peaks.df2[which(
          peaks.df2$seq != "" &
            peaks.df2$rel.int >= input$Relintco &
            peaks.df2$error <= input$errco &
            peaks.df2$Z %in% (if (input$Z) {
              c("", 1, 2)
            } else{
              c("", 1)
            })
        ),],selection = "single",callback = JS("table.on('click.dt', 'tr', function() {
            var row_=table.row(this).data();
            var data = [row_];
           Shiny.onInputChange('rows',data );
    });"))
      }
    }
    else {
      peaks.df2 <<- NULL
    }
  },)
  
  output$Selected <- renderDataTable({
    if ((!is.null(input$rows)) && (input$RelInt != "Highest peak")){
      peaks.df2.react()[which(rownames(peaks.df2.react()) ==as.numeric(input$rows[1])),1:2]
    }
    else {
      NULL
    }
    }, 
      filter="none",
      selection="none",
      autoHideNavigation=TRUE,
      options=list(searching=F,ordering=F,paging=F,lengthChange=F,info=F)
    )
  output$selPeak <- renderUI({
    if ((!is.null(input$rows)) && (input$RelInt != "Highest peak")){
      "Reference peak for 100% intensity:"
    }
    else {
      NULL
    }
  }
  )
  
  # observeEvent(input$rows, {
  #   print(input$rows[1])
  #   print(Sys.time())
  #   
  # })
  
  #Output plotly graph
  
  output$plly <- renderPlotly({
    if (exists("peaks.df2.react") & input$Submit != 0) {
      peaks.df.filtered <- peaks.df.filtered.react()
      
      letters <- as.list(unlist(strsplit(mypeptide(), "")))
      
      plot_ly(peaks.df.filtered) %>%
        add_segments(
          x = ~ start,
          xend = ~ stop,
          y = ~ sequence,
          yend = ~ sequence,
          showlegend = FALSE
        ) %>%
        layout(
          title = "",
          xaxis = list(
            title = "position of sequence",
            range = c(0, nchar(mypeptide())),
            ticktext = letters,
            tickvals = 1:30
          ),
          margin = list(l = 65)
        )
    }
    else {
      plotly_empty()
    }
  })
  
  # Generate excel file, save as temp.xlsx
  
  generate <- function() {
    peaks.df2 <- peaks.df2.react()
    peaks.df3 <- peaks.df3.react()
    basefile <- File()$basefile
    theor_frag.df <- theor_frag.df.react()
    wb <- createWorkbook(type = "xlsx")
    
    ## Create a new sheet in the workbook
    sheet <- createSheet(wb, sheetName = "Peak data")
    addDataFrame(if (input$SNM) {
      unique(peaks.df3[which(peaks.df3$rel.int >= input$Relintco &
                               peaks.df3$Z %in% (if (input$Z) {
                                 c("", 1, 2)
                               } else{
                                 c("", 1)
                               })),])
      
    }
    else {
      peaks.df2[which(
        peaks.df2$seq != "" &
          peaks.df2$rel.int >= input$Relintco &
          peaks.df2$error <= input$errco &
          peaks.df2$Z %in% (if (input$Z) {
            c("", 1, 2)
          } else{
            c("", 1)
          })
      ),]
    },
    sheet,
    col.names = TRUE,
    row.names = TRUE)
    
    ## Change column width
    setColumnWidth(sheet,
                   colIndex = 5,
                   colWidth = nchar(mypeptide()))
    setColumnWidth(sheet, colIndex = 6, colWidth = 12)
    
    ## Create a new sheet in the workbook
    sheet <- createSheet(wb, sheetName = "Visualisation")
    
    ## Style
    CELL_STYLE <- CellStyle(wb) +
      Font(wb, isBold = TRUE) +
      Alignment(wrapText = TRUE, horizontal = "ALIGN_CENTER")
    
    ## Add first row with peptide sequence
    # Sequence vector as first row
    pepseq <- unlist(str_split(mypeptide(), pattern = ""))
    pepseq.df <- data.frame()
    pepseq.df <- rbind(pepseq)
    row1 <- createRow(sheet, rowIndex = 1)
    for (c in 1:length(pepseq)) {
      row1_cell <- createCell(row1, colIndex = c)
      setCellValue(row1_cell[[1, 1]], pepseq[c])
      setCellStyle(row1_cell[[1, 1]], CELL_STYLE)
    }
    row1_cell <- createCell(row1, colIndex = c + 1)
    setCellValue(row1_cell[[1, 1]], "m/z")
    row1_cell <- createCell(row1, colIndex = c + 2)
    setCellValue(row1_cell[[1, 1]], "Sequence")
    row1_cell <- createCell(row1, colIndex = c + 3)
    setCellValue(row1_cell[[1, 1]], "Error")
    row1_cell <- createCell(row1, colIndex = c + 4)
    setCellValue(row1_cell[[1, 1]], "Relative intensity")
    
    setColumnWidth(sheet,
                   colIndex = c + 2,
                   colWidth = nchar(mypeptide()))
    setColumnWidth(sheet, colIndex = 6, colWidth = 12)
    ## Add fragment data
    
    # Filter on matched sequence
    peaks.df.filtered <- peaks.df.filtered.react()
    
    # Define color gradient
    cols <- heat.colors(1000)[1000:1]
    
    #str_locate(mypeptide,peaks.df.filtered$sequence)
    
    for (n in 1:nrow(peaks.df.filtered)) {
      # Get start and stop indexes
      positions <-
        str_locate(mypeptide(), peaks.df.filtered[n, "sequence"])
      
      # Make a fragment based on these indexes
      fragment <-
        c(
          rep(0, positions[1] - 1),
          rep(1, positions[2] - positions[1] + 1),
          rep(0, nchar(mypeptide()) - positions[2])
        )
      # print(fragment)
      
      # Define color based on intensity
      color <- cols[round(peaks.df.filtered[n, "rel.int"]/max(peaks.df.filtered$rel.int) * 1000)]
      
      # Style if no match in fragment
      CELL_STYLE <- CellStyle(wb) +
        Font(wb, isBold = FALSE) +
        Alignment(wrapText = TRUE, horizontal = "ALIGN_CENTER")
      
      # Style if match in fragment
      CELL_STYLE_1 <- CellStyle(wb) +
        Font(wb, isBold = FALSE) +
        Alignment(wrapText = TRUE, horizontal = "ALIGN_CENTER") +
        Fill(
          foregroundColor = color,
          backgroundColor = color,
          pattern = "SOLID_FOREGROUND"
        )
      
      #Create new row
      row <- createRow(sheet, rowIndex = n + 1)
      for (c in 1:length(pepseq)) {
        row_cell <- createCell(row, colIndex = c)
        if (fragment[c] == 1) {
          setCellStyle(row_cell[[1, 1]], CELL_STYLE_1)
        } else {
          setCellStyle(row_cell[[1, 1]], CELL_STYLE)
        }
      }
      row_cell <- createCell(row, colIndex = c + 1)
      setCellValue(row_cell[[1, 1]], peaks.df.filtered[n, "m/z"])
      row_cell <- createCell(row, colIndex = c + 2)
      setCellValue(row_cell[[1, 1]], peaks.df.filtered[n, "sequence"])
      row_cell <- createCell(row, colIndex = c + 3)
      setCellValue(row_cell[[1, 1]], peaks.df.filtered[n, "error"])
      row_cell <- createCell(row, colIndex = c + 4)
      setCellValue(row_cell[[1, 1]], peaks.df.filtered[n, "rel.int"])
    }
    
    ## Change column width
    setColumnWidth(sheet, colIndex = c(1:ncol(pepseq.df)), colWidth = 4)
    
    sheet <- createSheet(wb, sheetName = "All comb")
    addDataFrame(theor_frag.df[which(theor_frag.df$z %in% (if (input$Z) {
      c("", 1, 2)
    } else{
      c("", 1)
    })),], sheet)
    
    ## Save the workbook to a file
    saveWorkbook(wb, paste0("temp.xlsx"))
  }
  
  #When exitting, remove temp.xlsx
  
  session$onSessionEnded(function() {
    if (file.exists("temp.xlsx")) {
      file.remove(paste0("temp.xlsx"))
    }
    
  })
  
  # If press on download
  output$Exp <- downloadHandler(
    filename = function() {
      paste0("r-report-", File()$basefile, ".xlsx")
    },
    
    content = function(file) {
      generate()
      file.copy(paste0("temp.xlsx"), file)
    }
  )
}
# Run the application
shinyApp(ui = ui, server = server)
