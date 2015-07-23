library(openCyto)
library(data.table)
library(ggplot2)
library(plyr)
library(reshape2)
library(flowIncubator)
require(gridExtra)
require(lme4)
require(contrast)
require(epiR)
require(multcomp)
require(MASS)
require(scales)
library(data.table)

## max data size
options(shiny.maxRequestSize=1024^3) 

shinyServer(function(input, output, session) {
  
  ##-------------------------------------------------------------------------
  
  ## load flowset data
  set <- reactive({
    if (input$goButton == 0)
      return()
    isolate({fcsFiles <- input$fcsFiles
             if (is.null(fcsFiles))
               return(NULL)
             set <- read.flowSet(fcsFiles$datapath)
             sampleNames(set) <- fcsFiles$name})
    
    return(set)
  })
 ##_________________________ 
 
 TEMPLATE <-"C:/R Testing/Shiny/flowApp/gates.csv"
 gating_template<-gatingTemplate(TEMPLATE)
 
 Auto <- reactive({
    if(is.null(set()))
      return(NULL)
    GatingSet(set)
  })
 
    reactive({
   gating(gating_template, Auto)
      })
  
 ##_____________________________ 
  
  ## channel names
  channels <- reactive({
    if(is.null(set()))
      return(NULL)
    pd <- set()[[1]]@parameters@data
    channels <- pd$name
    return(channels)
  })
  
  ## marker selection
  markerNames <- reactive({
    if(is.null(set()))
      return(NULL)
    pd <- set()[[1]]@parameters@data
    markers <- paste("<", pd$name, ">:", pd$desc, sep = "")
    return(markers)
  })
  
  output$marker_select <- renderUI({
    if(is.null(markerNames())){
      return(NULL)
    }else{
      checkboxGroupInput('paras', strong('Select markers:'), 
                         markerNames(), selected = markerNames())
    }   
  })
  
  channelSelect <- reactive({
    if(is.null(markerNames()))
      return(NULL)
    setdiff(channels()[match(input$paras, markerNames())], c('Time', "time"))
  })
  
  ## sample selection 
  output$sample_select <- renderUI({
    if(is.null(set())){
      return(NULL)
    }else{
      checkboxGroupInput('samples', strong('Select samples:'), 
                         sampleNames(set()), selected = sampleNames(set()))
    }   
  })
  
  ##Automated Gating
  #________________________________
  
  
   


  
  #__________________________________
  ## time cut UI
  lapply(1:50, function(i) {
    output[[paste0('timeSlider', i)]] <- renderUI({
      if(is.null(set()))
        return(NULL)
      if (i <= length(set())){
        x <- set()[[i]]
        time <- findTimeChannel(x)
        mint <- min(exprs(x)[, time])
        maxt <- max(exprs(x)[, time])
        sliderInput(paste0('timeCut', i), strong(paste0('Time cut for sample ', i," :")),
                    min = mint, max = maxt, value = c(mint, maxt), step = 1)
      }
    })
  })
  
  ## time threshold get from time cut UI
  timeThres <- reactive({
    if(is.null(set()))
      return(NULL)
    timeThres <- list()
    for (i in 1:length(set())){
      sample <- sampleNames(set())[i]
      timeRange <- input[[paste0('timeCut', i)]]
      timeThres[[sample]] <- timeRange
    }
    return(timeThres)
  })
  
  ## time filtered flowset
  fset <- reactive({
    if(is.null(set()))
      return(NULL)
    flowList <- list()
    for (i in 1:length(set())){
      y <- set()[[i]]
      name <- sampleNames(set())[i]
      time <- findTimeChannel(y)
      params <- parameters(y)
      keyval <- keyword(y)
      sub_exprs <- exprs(y)
      timeRange <- timeThres()[[name]]
      minTime <- timeRange[1]
      maxTime <- timeRange[2]
      okCellid <- sub_exprs[,time] >= minTime & sub_exprs[,time] <= maxTime
      sub_exprs <- sub_exprs[okCellid, ]
      flowList[[name]] <- flowFrame(exprs = sub_exprs, parameters = params, description=keyval)
    }
    flowSet <- as(flowList, "flowSet")
  })
  
  
  ## reactive QC analysis
  fsScore <- reactive({
    if(is.null(set()))
      return(NULL)
    flowQA_firstScoreInit(set())})
  
  fcn <- reactive({
    if(is.null(set()))
      return(NULL)
    flowQA_cellnum(set())})
  
  fsm <- reactive({
    if(is.null(fset()))
      return(NULL)
    flowQA_marginevents(fset(), side = input$side, tol = input$tol)})
  
  fsm_score <- reactive({
    if(is.null(set()))
      return(NULL)
    marginEventsScore(set(), parms = channelSelect(),
                      side = input$side, tol = input$tol)
  })
  
  fstf <- reactive({
    if(is.null(set()))
      return(NULL)
    flowQA_timeflow(set(), binSize = input$tf_binSize)})
  
  fstf_score <- reactive({
    if(is.null(fstf()))
      return(NULL)
    timeFlowScore(set(), binSize = fstf()$binSize, 
                  varCut = input$tf_varCut)})
  
  fstf_gvis <- reactive({
    if(is.null(set()))
      return(NULL)
    flowS <- set()
    binSize <- fstf()$binSize
    mint <- min(fsApply(flowS, function(x){
      time <- findTimeChannel(x)
      min(exprs(x)[, time])}))
    maxt <- max(fsApply(flowS, function(x){
      time <- findTimeChannel(x)
      max(exprs(x)[, time])}))
    nrBins <- floor(max(fsApply(flowS, nrow, use.exprs = TRUE)) / binSize)
    tbins <- seq(mint, maxt, len=nrBins + 1)    # time bins
    
    counts <- fsApply(flowS, function(x){
      time <- findTimeChannel(x)
      xx <- sort(exprs(x)[, time])   # time channel
      tbCounts <- hist(xx, tbins, plot = FALSE)$counts  # number of events per time bin
    })
    tcord <- as.data.frame(t(counts))
    tcord$time <- tbins[-1]
    return(tcord)
  })
  
  fstl <- reactive({
    if(is.null(set()))
      return(NULL)
    flowQA_timeline(set(), binSize = input$tl_binSize,
                    varCut= input$tl_varCut)})
  
  qaScore <- reactive({
    if(is.null(fstl()))
      return(NULL)
    fsScore <- flowQA_scoreUpdate(fsScore(), fsm_score())
    fsScore <- flowQA_scoreUpdate(fsScore, fstf_score())
    fsScore <- flowQA_scoreUpdate(fsScore, fstl()$tlScore)
    qaScore <- flowQA_scoreSummary(fsScore)
    return(qaScore) })
  
  ## render cell number plot
  output$cnplot <- renderGvis({
    if(is.null(set()))
      return(NULL)
    cnframe <- fcn()
    data <- as.data.frame(cnframe[match(input$samples, cnframe$sampleName), ] )
    gvisColumnChart(data, 
                    xvar = "sampleName", yvar = "cellNumber",
                    options=list(title="Cell Number for Each Sample",
                                 width = 900, height = 500, vAxis="{minValue: 0}")
    )
  })
  
  ##__________________________________________________
  
  output$gating <- renderPlot({
    if(is.null(set()))
      return(NULL)
    plot(gating_template, fontsize=32)
  })
  
 
 plot1 <- reactive({
   if(is.null(Auto()))
     return("Error?")
   plotGate(Auto,"lymph",arrange=FALSE, smooth=FALSE)
 })
  
 output$scatter <- renderPlot({
    plot1
 })

  
  ##_______________________________________________________
  output$cntable <- renderGvis({
    if(is.null(set()))
      return(NULL)
    cnframe <- fcn()
    data <- as.data.frame(cnframe[match(input$samples, cnframe$sampleName), ] )
    gvisTable(data)
  })
  
  ## render margin events plot        
  output$meplot <- renderPlot({
    if(is.null(set()))
      return(NULL)
    perc <- fsm()
    perc <- as.matrix(perc[channelSelect(), input$samples])
    col.regions=colorRampPalette(c("white",  "darkblue"))(256)
    print(levelplot(perc*100, scales = list(x = list(rot = 45), y = list(rot = 45)),
                    xlab="", ylab="", main="Percentage of margin events",
                    col.regions=col.regions))
  })
  
  output$metable <- renderGvis({
    if(is.null(set()))
      return(NULL)
    perc <- fsm()
    perc <- perc[channelSelect(), input$samples]
    gvisTable(as.data.frame(t(perc)))
  })
  
  output$me_sample_choose <- renderUI({
    if(is.null(set())){
      return(NULL)
    }else{
      selectInput('me_sample_choose', 'Choose a sample:', 
                  choices = input$samples, width = "100%")
    }   
  })
  
  output$me_channel_choose <- renderUI({
    if(is.null(set())){
      return(NULL)
    }else{
      selectInput('me_channel_choose', 'Choose a Channel:', 
                  choices = channelSelect(), width = "100%")
    }   
  })
  
  meTimeData <- reactive({
    if(is.null(set()))
      return(NULL)
    fcs <- fset()[[input$me_sample_choose]]
    exp <- fcs@exprs
    para <- pData(fcs@parameters)
    ranges <- range(fcs)
    time <- findTimeChannel(fcs)
    channel <- input$me_channel_choose
    tc <- as.data.frame(exp[ ,c(time, channel)])
    tc_range <- ranges[ ,channel]
    
    tc_neg <- tc[tc[,2] <= tc_range[1] - input$tol, ] 
    tc_pos <- tc[tc[,2] >= tc_range[2] + input$tol, ] 
    colnames(tc_pos) <- c("Time", "Upper Margin Events")
    colnames(tc_neg) <- c("Time", "Lower Margin Events")
    res <- merge(tc_pos, tc_neg, all = TRUE)
  })
  
  output$upMeTimePlot <- renderGvis({
    if(is.null(meTimeData()))
      return(NULL)
    upMeData <- as.data.frame(meTimeData())[,c(1,2)]
    gvisScatterChart(upMeData, options = list(pointSize = 0.5))
  })
  
  output$lowMeTimePlot <- renderGvis({
    if(is.null(meTimeData()))
      return(NULL)
    lowMeData <- as.data.frame(meTimeData())[,c(1,3)]
    if(nrow(lowMeData) > 2000){
      lowMeData <- lowMeData[sample(1:nrow(lowMeData), 2000), ]
    }
    gvisScatterChart(lowMeData, options = list(pointSize = 0.5)) 
  })
  
  output$meTime_text <- renderText({
    channelRange <- range(fset()[[input$me_sample_choose]])[ ,input$me_channel_choose]
    lower_thres <- round(channelRange[1] - input$tol, digits = 2)
    higher_thres <- round(channelRange[2] + input$tol, digits = 2)
    paste0("Lower threshold: ", lower_thres, ";    Upper threshold:", higher_thres)
  })
  
  ## render time flow plot
  output$s_tfplot <- renderPlot({
    if(is.null(set()))
      return(NULL)
    timeFlowData <- fstf()$timeFlowData[input$samples]
    vcut <- input$tf_varCut
    timeFlowPlot(timeFlowData, timeThres(), vcut) 
  })
  
  output$tfplot <- renderGvis({
    if(is.null(set()))
      return(NULL)
    tcord <- fstf_gvis()
    time <- tcord$time
    for (i in 1:length(set())){
      sample <- sampleNames(set())[i]
      timeRange <- timeThres()[[sample]]
      minTime <- timeRange[1]
      maxTime <- timeRange[2]
      badCellid <- time < minTime | time > maxTime
      tcord[[sample]][badCellid] <- NA
    }
    tcord <- subset(tcord, select = c(input$samples, "time"))
    gvisLineChart(tcord, xvar = "time", options=list(title="Cell Flow vs. Time") ) 
  })
  
  output$tf_text <- renderText({
    paste0("Bin size: ", fstf()$binSize)
  })
  
  
  ## render time line plot
  output$tl_sample_choose <- renderUI({
    if(is.null(fstl()$timeLineData)){
      return(NULL)
    }else{
      selectInput('tl_sample_choose', 'Choose a sample:', 
                  choices = input$samples, width = "100%")
    }   
  })
  
  
  output$s_tlplot <- renderPlot({
    if(is.null(fstl()$timeLineData))
      return(NULL)
    timeLineData <- fstl()$timeLineData[input$samples]
    timeLinePlot(lapply(timeLineData, function(x) x[[1]]), timeThres(),
                 channels = channelSelect())  
  })
  
  
  output$tlplot <- renderGvis({
    if(is.null(fstl()$timeLineData))
      return(NULL)
    timeLineData <- fstl()$timeLineData
    x <- timeLineData[[input$tl_sample_choose]]$res
    range <- timeThres()[[input$tl_sample_choose]]
    time <- x[ ,1]
    badCellid <- time < range[1] | time > range[2]
    x[badCellid,-1] <- NA
    x <- x[, c("time", channelSelect())]
    gvisLineChart(as.data.frame(x), xvar = "time", 
                  options=list(title="Timeline plot") )      
  })
  
  output$tl_text <- renderText({
    paste0("Bin size: ", fstl()$binSize)
  })
  
  output$tabplot <- renderPlot({
    if(is.null(fset()))
      return(NULL)
    x <- fset()
    flowQA_tabplot(x, input$tl_sample_choose, channels = channelSelect(),
                   binSize = fstl()$binSize)    
  })
  
  ## render score plot
  output$score_sample_choose <- renderUI({
    if(is.null(qaScore())){
      return(NULL)
    }else{
      selectInput('score_sample_choose', 'Choose a sample:', 
                  choices = input$samples, width = "100%")
    }   
  })
  
  output$sort_ID <- renderUI({
    if(is.null(qaScore())){
      return(NULL)
    }else{ 
      selectInput('sort_id', 'Sort column:', 
                  choices = colnames(qaScore()[[1]]) )
    }   
  })
  
  output$scoreplot <- renderPlot({
    if(is.null(qaScore()))
      return(NULL)
    sample <- input$score_sample_choose
    sid <- match(input$sort_id, colnames(qaScore()[[sample]]))
    
    flowQA_scoreplot(qaScore(), input$score_sample_choose, timeThres(), 
                     scoreThres = input$scoreThres,
                     sortID = sid, nBins = input$score_nrbins)    
  })
  
  ## data download
  output$downloadMarkers <- downloadHandler(
    filename = function() { 
      "markers.txt"
    },
    content = function(file) {
      write.table(channelSelect(), file, quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
    }
  )
  
  output$downloadFCS <- downloadHandler(
    filename = function(){
      paste0("new_fcs_files", ".tar")
    },
    content = function(file){
      fsep = .Platform$file.sep
      tempdir = paste0(tempdir(), fsep, gsub("\\D", "_", Sys.time()))
      data <- fset()[input$samples]
      if(is.null(data)){
        return(NULL)
      }
      write.flowSet(data, tempdir)
      tar(tarfile = file, files = tempdir)
    }
  )
  
  
})


