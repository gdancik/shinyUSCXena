#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(UCSCXenaTools)
library(DT)
library(dplyr)
library(edgeR)

SAMPLES_TO_KEEP <- '01A'
DEFAULT_GENE <- 'TP53'

MyXenaData <- dplyr::filter(XenaData, grepl('GDC TCGA', XenaCohorts)) %>%
              dplyr::filter(Label == "HTSeq - Counts" | Label == "survival data")


# Define UI for application that draws a histogram
ui <- bootstrapPage(

    # Application title
    #titlePanel("Old Faithful Geyser Data"),

        navbarPage('UCSC Xena Analyzer', 
                   tabPanel('Data',
                            h3("Data selection"),
                            selectInput('selectedCohort', label = 'Selected Data: ', choices = unique(MyXenaData$XenaCohorts),
                                        width = "30%"), br(),
                            dataTableOutput("XenaData")
                   ),
                   tabPanel('Download and Process',
                            uiOutput('confirmCohort'),
                            actionButton('downloadCohort', 'Go!'), br(), hr(),
                            verbatimTextOutput('barCodeTypeSummary'),
                            
                            fluidRow(
                                column(6,
                                    dataTableOutput("SurvivalData")
                                ), column(6,
                                )
                            )
                   ),
                   
                   tabPanel('Analyze', 
                            div(
                                div(style = 'display:inline-block; width: 25%',
                                    selectizeInput('selectGene', 'Select Gene', choices = 1:2)
                                ), div(style = 'display:inline-block; width: 25%',
                                    selectizeInput('selectProbe', 'Select Probe', choices = 1:2)
                                )
                            )
                    )
        
        )
)

# Define server logic required to draw a histogram
server <- function(session, input, output) {
    print("running server...")
    
    get_survival_data <- function() {
        rds_file = paste0('cache/', input$selectedCohort, '_survival.rds')
        if (file.exists(rds_file)) {
            return(readRDS(rds_file))    
        }
        
        # limit to desired cohort
        blca <- MyXenaData %>% filter(XenaCohorts == input$selectedCohort)
        
        # Get the phenotype / clinical data
        cli_query = blca %>%
            filter(Label == "survival data") %>%  # select survival data
            XenaGenerate() %>%  # generate a XenaHub object
            XenaQuery() %>%     # generate the query
            XenaDownload()      # download the data
        
        # prepare (load) the data into R
        blca_survival <- XenaPrepare(cli_query)
        saveRDS(blca_survival, rds_file)
        blca_survival
    }
    
    get_count_data <- function() {
        rds_file = paste0('cache/', input$selectedCohort, '_counts.rds')
        if (file.exists(rds_file)) {
            return(readRDS(rds_file))    
        }
        
        # limit to desired cohort
        blca <- MyXenaData %>% filter(XenaCohorts == input$selectedCohort)
        
        # Get the RNA-seq data, including the "probe map"
        cli_query <- blca %>% filter(Label == 'HTSeq - Counts') %>%
            XenaGenerate() %>%  # generate a XenaHub object
            XenaQuery() %>%
            XenaDownload(download_probeMap = TRUE)
        
        # prepare (load) the data into R
        blca_counts <- XenaPrepare(cli_query)

        g <- grep('counts', names(blca_counts))
        
        X <- data.frame(blca_counts[[g]])
        rownames(X) <- X$Ensembl_ID
        X <- X[,-1] 
        blca_counts[[g]] <- X
        names(blca_counts)[g] <- 'counts'
        names(blca_counts)[setdiff(1:2,g)] <- 'probeMap'
        
        saveRDS(blca_counts, rds_file)
        blca_counts
    }
    
    process_expr <- function(X, S, keep_samples) {
        
        colnames(X) <- gsub('\\.', '-', colnames(X))
        qry <- paste0(keep_samples, '$', collapse = '|')
        
        g <- grep(qry, colnames(X))
        X <- X[,g]
        
        # We still need to match the expression data with the clinical data
        # Let's do that by first finding the samples that are common
        # between the expression and clinical data. We can use 
        # intersect(a,b) to return a vector containing the elements common
        # to vectors 'a' and 'b'
        
        common_samples <- intersect(colnames(X), S$sample)
        
        # we then use match(x, t) to get a vector of indices. The value
        # x[i] is the index of 't' containing the i^th value of 'x'
        
        mx <- match(common_samples, colnames(X))
        my <- match(common_samples, S$sample)
        
        X <- X[,mx]
        S <- S[my,]
        
        # Make sure that the samples match -- if they don't, this will produce an error
        stopifnot(all(colnames(X) == S$sample))
        
        #######################################################
        # Process the expression data
        #######################################################
        
        ###################################################
        # Data is log2(counts + 1), so we need to get
        # back to the scale of counts
        # Note: This was left off of the original script
        #   and added after class on 11/10
        ####################################################
        X <- round(2**X - 1)
        
        # first create a Digital Gene Expression (DGE) list object,
        # which contains counts and library size information
        dge <- DGEList(counts=X)
        
        # remove genes with low counts, since these should not
        # be considered in our downstream analysis. The default
        # min.count = 10 will be used and we require this min
        # count in at least min.prop = 10% of samples
        keep <- filterByExpr(dge,min.prop = .10 )
        dge <- dge[keep,,keep.lib.sizes=FALSE]
        
        # apply TMM normalization, which computes the normalization 
        # factors. The actual normalization is done in a later step
        dge <- calcNormFactors(dge, method = "TMM")
        
        
        # Calculate the log CPM values, using the normalization factors;
        # 3 counts are added to each observation to prevent log 0 values
        logCPM <- cpm(dge, log = TRUE, prior.count = 3)
        
        list(logCPM = logCPM, survival = S)
        
    }
    
    output$XenaData <- DT::renderDT(MyXenaData %>% dplyr::filter(XenaCohorts == input$selectedCohort))
    
    output$confirmCohort <- renderUI({
        HTML('<p>Click button to download:',
             '<span style = "color:red">', input$selectedCohort, '</span>', 
             '(', SAMPLES_TO_KEEP, ')</p>')
        })
    
    probeMap <- reactiveVal()
    
    observeEvent(input$downloadCohort, {

        showNotification(id = 'progress', 'getting survival data...', type = 'message')
        blca_survival <- get_survival_data()

        output$SurvivalData <- DT::renderDT(blca_survival)
        
        showNotification(id = 'progress', 'getting count data...', type = 'message')        
        blca_counts <- get_count_data()
        
        X <- blca_counts$counts
        
        probeMap(blca_counts$probeMap)
        
        # output bar code types
        { 
            # 'change '.' to '-' so sample ID format is consistent
            colnames(X) <- gsub('\\.', '-', colnames(X))
            
            # Note that the sample ID is a barcode that has a special meaning:
            #   https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
            # In particular, the 4th section describes the 'Sample' which is 
            #   either tumor (01 - 09) or normal (10-19). For details see:
            #   https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
            
            # Summarize these samples, using a handy R trick that takes 
            # advantage of the fact that `[]` is a function
            s <- strsplit(colnames(X), '-')   # split each string in a vector by a '-'
            
            t <- sapply(s, `[`, 4) %>% table()
            
            output$barCodeTypeSummary <- renderPrint(t)
            
            
        }
        
        showNotification(id = 'progress', 'processing data...', type = 'message')        
        
        rds_file = paste0('cache/', input$selectedCohort, 
                          '_', paste0(SAMPLES_TO_KEEP, collapse = '_'),
                          '.rds')
        
        if (file.exists(rds_file)) {
            L <- readRDS(rds_file)
        } else {
            L <- process_expr(X, blca_survival, SAMPLES_TO_KEEP)
            saveRDS(L, file = rds_file)
        }
        
        output$SurvivalData <- DT::renderDT(L$survival)
        #output$CountData <- DT::renderDT(L$logCPM)
           
        pm <- probeMap()
        probeMap(pm[pm$id %in% rownames(L$logCPM),])
        showNotification(id = 'progress', 'Getting available genes', type = 'message')

        genes <-  probeMap()$gene %>% unique() %>% sort()
        
        updateSelectizeInput(session, 'selectGene', 'Select Gene', 
                             selected = DEFAULT_GENE,
                             choices = genes, server = TRUE)
        
        showNotification(id = 'progress', 'Complete -- select gene on next page', type = 'message') 
   
    })

    observeEvent(input$selectGene, {
        pm <- probeMap()
        
        if (is.null(pm)) {
            return()
        }
        
        print(head(pm))
        
        keep <- pm$gene %in% input$selectGene
        probes <- pm$id[keep]
        
        updateSelectizeInput(session, 'selectProbe', 'Select Probe', 
                             choices = probes)
        
        
    }, ignoreNULL = TRUE, ignoreInit = TRUE)
    
}


# Run the application 
shinyApp(ui = ui, server = server)
