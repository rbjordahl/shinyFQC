library(shiny)


shinyUI(fluidPage(
  titlePanel("Quality Control for Flow Cytometry Data"),
  
  fluidRow(
    column(3,
           fileInput('fcsFiles', strong('Choose fcs file(s):'), multiple = TRUE, 
                     accept = c('text/fcs', '.fcs')),
           
           actionButton("goButton", "Submit!"),
           hr(),
           downloadButton('downloadMarkers', 'Download markers table'),
           br(),
           downloadButton('downloadFCS', 'Download new FCS files'),
           hr(),
           
           ## sample limits: 50
           uiOutput("sample_select"),
           
           lapply(1:50, function(i) {
             uiOutput(paste0('timeSlider', i))
           }),
           
           hr(),
           uiOutput("marker_select"),
           
           hr(),
           div(style = "margin-top: 30px; width: 200px; ", HTML("Developed by")),
           div(style = "margin-top: 10px; ", 
               HTML("<img style='width: 150px;' src='http://archild.sign.a-star.edu.sg/images/logo.png'>"))
    ),
    column(9,
           tabsetPanel(type = "pills",
                       
                       tabPanel("Cell numbers", 
                                hr(),
                                htmlOutput("cnplot")),
                                          
                       tabPanel("Time flow", fluidPage(
                         hr(),
                         htmlOutput("tfplot"),
                         hr(),
                         
                         fluidRow(
                           column(4, offset = 1,
                                  numericInput("tf_binSize", "Bin size:", value = NA)
                           ),
                           column(4, 
                                  numericInput("tf_varCut", "Variation cut:", value = 1)
                           ) 
                         ),
                         
                         textOutput("tf_text")
                       )),
                       
                       tabPanel("Timeline", fluidPage(
                         hr(),
                         uiOutput("tl_sample_choose"),
                         
                         hr(),
                         h4("Timeline plot:"),
                         htmlOutput("tlplot"),
                         
                         hr(),
                         h4("Expression table plot:"),
                         plotOutput("tabplot"),
                         
                         hr(),
                         fluidRow(
                           column(4, offset = 1,
                                  numericInput("tl_binSize", "Bin size:", value = NA)
                           ),
                           column(4,
                                  numericInput("tl_varCut", "Variation cut:", value = 1)
                           ) 
                         ),
                         textOutput("tl_text")
                       )),
                       
                       tabPanel("Margin Events", 
                                fluidPage(
                                  hr(),
                                  fluidRow(
                                    column(4, offset = 1,
                                           selectInput("side", "Select checking side:",
                                                       choices = c("both", "upper", "lower"), 
                                                       selected = "both")
                                    ),
                                    column(4, 
                                           numericInput("tol", "Tolerance value:", -.Machine$double.eps)
                                    ) 
                                  ),
                                  hr(),
                                  plotOutput("meplot"),
                                  hr(),
                                  
                                  h4("Time Distribution of Margin Events:"),
                                  fluidRow(
                                    column(5, offset = 1, uiOutput("me_sample_choose") ),
                                    column(3, uiOutput("me_channel_choose") )
                                  ),
                                  
                                  hr(),
                                  htmlOutput("upMeTimePlot"),
                                  htmlOutput("lowMeTimePlot"),
                                  br(),
                                  br(),
                                  textOutput("meTime_text")
                                  
                                )),
                       
                       tabPanel("QA score", fluidPage(
                         hr(),
                         uiOutput("score_sample_choose"),
                         
                         h4("QA score summary:"),
                         plotOutput("scoreplot"),
                         fluidRow(
                           column(width = 3,
                                  numericInput("score_nrbins", "Row bins:", value = 100)
                           ),
                           column(width = 3, offset = 1, 
                                  numericInput("scoreThres", "Threshold score:", value = 3, step = 0.1)
                           ),
                           column(width = 4, offset = 1, 
                                  uiOutput("sort_ID")
                           ) 
                         )
                       )),
                       
                       tabPanel("Summary",
                                fluidPage(
                                  verticalLayout(
                                    hr(),
                                    h4("Cell Number check:"),
                                    htmlOutput("cntable"),
                                    hr(),
                                    h4("Margin Events check:"),
                                    htmlOutput("metable"),
                                    hr(),
                                    h4("Time flow check:"),
                                    plotOutput("s_tfplot"),
                                    hr(),
                                    h4("Timeline check:"),
                                    plotOutput("s_tlplot")
                                  )
                                )),
                       
                       tabPanel("Gating Heirarchy", 
                                            hr(),
                                            plotOutput("gating")),
                       tabPanel("Data Plots", 
                                hr(),
                                plotOutput("scatter"))
                       
           )
    )
  )
))