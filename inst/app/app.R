
library(shiny)
library(ranger)
library(kmer)
library(ape)
library(DT)
library(plotly)
options(shiny.maxRequestSize = 1000*1024^2)
library(rintrojs)
library(shinythemes)
# library(magrittr)
# library(dplyr)
library(shinyjs)
#for mac, without this Rstudio crashes
Sys.setenv(LIBGL_ALWAYS_SOFTWARE=1)
# library(tibble)
library(RColorBrewer)
library(parallel)

# library(shinylogs)

ui <- fluidPage(
 theme = shinytheme("simplex"),
  introjsUI(), # must include in UI  
  introBox(titlePanel(div("Covidex",
            img(style="width: 50px", src = "covidex_coronavirusazul2.png"))),div("Enabled by data from",tags$a(img(style="width: 50px", src = "GISAID.png"),href="https://www.gisaid.org/")),
           data.step = 1,
           data.intro = "Covidex is an ultra fast and accurate subtyping tool of hcov-SARS-CoV-2 genomes. The classification is performed using 3 machine learning models (random forest) of a k-mer database"
           ),
  
 

  sidebarLayout(
    
    sidebarPanel(
     actionButton("btn","Help"),
    
     introBox( fileInput("file", label = h3("Query file (multi-fasta format)"),
                accept = c(".text",".fasta",".fas",".fasta")),
                data.step = 2,
                data.intro = "Please select the fasta file with the sequences to subtype (max size 10 Mb if running in server)"
                 ),
     textInput("user", "User name (optional)", "anonymous"),
    
     introBox(  actionButton("go", "RUN"), 
                data.step = 3,
                data.intro = "When ready press RUN to begin the analysis"
     ),
     strong(textOutput("text")),
     strong("Questions?"), 
     tags$a(actionButton(inputId = "email1", label = "Contact Admin",
                         icon = icon("envelope", lib = "font-awesome")),
            href="mailto:cacciabue.marco@inta.gob.ar;aguilera.pablo@inta.gob.ar"),
     br(), 
     strong("or visit"),
     tags$a(img(style="width: 100px", src = "sourceforge.png"),href="https://sourceforge.net/projects/covidex/"),
     br(), 
     strong("or by Twitter"),
     tags$a(img(style="width: 25px", src = "twitter.png"),href="https://twitter.com/marcocacciabue"),
     div(class = "footer",
         includeHTML("www/footer.html")),
    
  

),
       
 
    
   mainPanel(
     introBox( downloadButton('report',"Generate report (ENG)"),
               data.step = 5,
               data.intro="Press this button to generate a report with classification data."),
     
      tabsetPanel(
       tabPanel(title="Results",   introBox( DT::dataTableOutput("table"),
                data.step = 4,
                data.intro="In the Results tab, 
                all uploaded sequences with the classification results will appear with an associated confidence value.")),
       tabPanel("Classification model information",DT::dataTableOutput("table2")))
    )

  ),
fluidRow(
  column(12, align="center",
         strong("GISAID data provided on this website are subject to GISAID's"),
         tags$a(href="https://www.gisaid.org/DAA/",strong("Terms and Conditions")),
         
  )
)
)
server <- shinyServer(function(input, output, session) {
  # track_usage(storage_mode = store_json(path = "logs/"))
  introjs(session)
  hintjs(session, options = list("hintButtonLabel"="Done!"))
  observeEvent(input$btn,
               introjs(session))
  observeEvent(input$go, {
    if(length(input$file)==0){
      showModal(modalDialog(
        title = "Important message", easyClose = TRUE,
        "Please load the fasta file first and then press RUN.
Also, remember that the file must NOT exceed 10 MB in size (about 350 complete sequences). If you need to process larger files consider downloading Covidex from the Sourceforge repository.
"
      ))}})

  model1 <- readRDS(paste("models/",sep="","1.rds"))
  model2 <- readRDS(paste("models/",sep="","2.rds"))
  model3 <- readRDS(paste("models/",sep="","3.rds"))
  model1_data<-data.frame(Model=model1$info,date=model1$date,trees=model1$num.trees,Oob= round(model1$prediction.error,4))
  model2_data<-data.frame(Model=model2$info,date=model2$date,trees=model2$num.trees,Oob= round(model2$prediction.error,4))
  model3_data<-data.frame(Model=model3$info,date=model3$date,trees=model3$num.trees,Oob= round(model3$prediction.error,4))
  
  model_data<-rbind(model1_data,model2_data,model3_data)
  

  data_reactive<- eventReactive(input$go,{
    req(input$file)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "processing data", value = 0)
    progress$inc(0.2, detail = paste("Reading fasta"))
    


    #Read input file in fasta format
    query <- read.FASTA(input$file[1,4], type = "DNA")
    
    model<-model1
    if (length(query)>2000){
      # progress$inc(0.4, detail = paste("Sample larger than 1000 sequences: runing in multithread mode"))
      #get number of available cores in the machine to see how much to split the data.
      progress$inc(0.4, detail = paste("Counting kmers in multithread mode"))
      cores<-detectCores()-1
      
      step<-ceiling(length(query)/cores)
      last<-step
      first<-1
      query_split<-0
      #split query sequences
      for (j in 1:cores) {
        #check to get the correct number of sequences at the end  
        if ((first+step)<length(query)){
          query_split[j]<-list(query[first:last])
          first<-first+step
          last<-last+step}else{
            query_split[j]<-list(query[first:length(query)])
            
          }
        
      }
      #create core cluster
      cl <- makeCluster(detectCores()-1)
      # progress$inc(0.6, detail = paste("Kmer counting in multithread mode"))
      
      #Use parallel function to send each split dataset to a different thread.
      query_count_split<-parLapply(cl,query_split,count_parallel,model$kmer)
      stopCluster(cl)
      progress$inc(0.8, detail = paste("merging outputs"))
      
      #Combine the results from the different threads in a single dataset
      query_count<-as.numeric()
      for (i in 1:length(query_count_split))  {
        query_count<-rbind(query_count, query_count_split[[i]])
      }
      
      
      #normalize Kmer counts acording to kmer size and sequence length
      genome_length<-0
      n_length<-0
      for(H in 1:length(query_count[,1])){
        
        k<-query[H]
        k<-as.matrix(k)
        query_count[H,]<- query_count[H,]*model$kmer/(length(k))
        genome_length[H]<-length(k)
        n_length[H]<-round(100*base.freq(k,all = TRUE)[15],2)
      }}
    else{
      
      #Calculate k-mer counts from query sequences
      # progress$inc(0.4, detail = paste("Runing in Single-thread mode"))
      progress$inc(0.4, detail = paste("Counting kmers in single-thread mode"))
      query_count<-kcount(query , k=model$kmer)
      genome_length<-0
      n_length<-0
      for(i in 1:length(query_count[,1])){
        
        k<-query[i]
        k<-as.matrix(k)
        query_count[i,]<- query_count[i,]*model$kmer/(length(k))
        genome_length[i]<-length(k)
        n_length[i]<-round(100*base.freq(k,all = TRUE)[15],2)
      }
      
    }
    
    progress$inc(0.9, detail = paste("Predicting"))
    
    
    calling<-predict(model,query_count)
    #Run the predict method from de Ranger package, retaining the classification result from each tree in the model (to calculate a probability value for each classification)
    calling_all<-predict(model,query_count,predict.all = TRUE)
    probability <- rep(0, length(calling_all$predictions[,1]))
    
    for (i in 1:length(calling_all$predictions[,1])) {
      #extract predictions for each query sample in temp vector,
      #count the number of correct predictions and divide by number of trees to get a probability.
      temp<-calling_all$predictions[i,]
      probability[i] <- sum(temp==which(model$forest$levels==calling$predictions[i]))/model$num.trees
      
    }
   
    
    model<-model2
    calling2<-predict(model,query_count)
    #Run the predict method from de Ranger package, retaining the classification result from each tree in the model (to calculate a probability value for each classification)
    calling_all2<-predict(model,query_count,predict.all = TRUE)
    probability2 <- rep(0, length(calling_all2$predictions[,1]))
    
    for (i in 1:length(calling_all2$predictions[,1])) { #
      #extract predictions for each query sample in temp vector,
      #count the number of correct predictions and divide by number of trees to get a probability.
      temp2<-calling_all2$predictions[i,]
      probability2[i] <- sum(temp2==which(model$forest$levels==calling2$predictions[i]))/model$num.trees
      
    }
    
    model<-model3
    calling3<-predict(model,query_count)
    #Run the predict method from de Ranger package, retaining the classification result from each tree in the model (to calculate a probability value for each classification)
    calling_all3<-predict(model,query_count,predict.all = TRUE)
    probability3 <- rep(0, length(calling_all3$predictions[,1]))
    
    for (i in 1:length(calling_all3$predictions[,1])) {       #count the number of correct predictions and divide by number of trees to get a probability.
      temp3<-calling_all3$predictions[i,]
      probability3[i] <- sum(temp3==which(model$forest$levels==calling3$predictions[i]))/model$num.trees
      
    }
    
    
    
    
    
    # QC<-as.character(probability>input$QC)
    n_QC<-(n_length<1)
    Length_QC<-(genome_length>29500)
    VOC<-(calling2$prediction=="B.1.351")|(calling2$prediction=="B.1.1.7")|(calling2$prediction=="P.1")|(calling2$prediction=="B.1.427")|(calling2$prediction=="B.1.429")
    VOI<-(calling2$prediction=="P.2")|(calling2$prediction=="B.1.525")|(calling2$prediction=="B.1.526")
    FLAG<-(VOI|VOC)
    data_out <- data.frame(label= row.names(query_count),Nextstrain=calling$prediction, probability_N=probability,Rambaut=calling2$prediction,probability_R=probability2,GISAID=calling3$prediction,probability_G=probability3,FLAG=FLAG,VOC=VOC,VOI=VOI,Length=genome_length,Length_QC=Length_QC,N=n_length,N_QC=n_QC)
    rm(query_count,calling,calling_all,calling2,calling_all2)

  list(message="Done!",data_out=data_out)
                })
  output$text <- renderText({
    data_reactive()$message
  })
  
  table<-reactive({ 
    
    data_out<-data_reactive()$data_out
    
    data_out
    
    })
  
  output$table <- DT::renderDataTable({

    col<-brewer.pal(5,"Blues")
    col2<-brewer.pal(5,"Reds")
    
    table<-table()
    datatable(table,selection = 'single',
              options = list(
                columnDefs = list(list(targets = c(8,9,10,12,14), visible = FALSE)),
                dom = 'Bfrtip',
                lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
                pageLength = 15,
                buttons =
                  list('copy', 'print', list(
                    extend = 'collection',
                    buttons = list(
                      list(extend = 'csv', filename = paste(input$file[1,1],"covidex_results"),sep=""),
                      list(extend = 'excel', filename = paste(input$file[1,1],"covidex_results"),sep=""),
                      list(extend = 'pdf', filename = paste(input$file[1,1],"covidex_results"),sep="")),
                    text = 'Download'
                  )
                  )
              ))%>% formatStyle("Rambaut","FLAG",
                                backgroundColor = styleEqual(c(0, 1), c(col[1], col[3])))%>% formatStyle("Length","Length_QC",
                                                                                                         backgroundColor = styleEqual(c(0, 1), c(col2[3], col2[1])))%>% formatStyle("N","N_QC",
                                                                                                                                                                                    backgroundColor = styleEqual(c(0, 1), c(col2[3], col2[1]))) 
  })
  # 


# output$table <- DT::renderDataTable({
#   req(input$go)
#   req(input$file)
#   data_reactive()$data_out
# 
#   
# })
model_table<-reactive({
  model_data
})

output$table2 <- DT::renderDataTable({
  datatable(model_table())
  
  
})

output$report <- downloadHandler(
  # For PDF output, change this to "report.pdf"
  filename = "report.html",
  content = function(file) {
    # Copy the report file to a temporary directory before processing it, in
    # case we don't have write permissions to the current working dir (which
    # can happen when deployed).
    tempReport <- file.path(tempdir(), "report.rmd")
    file.copy("report.rmd", tempReport, overwrite = TRUE)
    
    # Set up parameters to pass to Rmd document
    
    # Knit the document, passing in the `params` list, and eval it in a
    # child of the global environment (this isolates the code in the document
    # from the code in this app).
    rmarkdown::render(tempReport, output_file = file)
    
  }
)
version<-reactive({
  version<-2.1
  return(version)})

count_parallel<-function(x,kmer) ({ 
  
  
  library(kmer)
  kcount(x,k= kmer)})
##uncomment this to end R session when the windows closes
# session$onSessionEnded(function() { 
#  stopApp()
#     q("no") 
#  })
 
})



shinyApp(ui, server)
