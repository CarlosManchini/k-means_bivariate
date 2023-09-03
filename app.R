library(shiny)
library(shinythemes)
library(ggplot2)
library(MASS)

ui <- pageWithSidebar( 
    headerPanel(HTML(paste(h3('K-means clustering'),
                           p("Trabalho - Análise Multivariada")))),
    sidebarPanel(theme = shinytheme("cerulean"),
                 fluidRow(align="center", 
                          div(tags$b("Distribuição normal bivariada"), style="text-indent:10px;font-size:120%;"),
                          withMathJax(helpText("$$f(x)= \\frac{1}{(2\\pi)^{p/2} |\\Sigma|^{1/2}}\\exp\\left\\lbrace-\\frac{1}{2} (x-\\mu)^\\top \\Sigma^{-1} (x-\\mu) \\right\\rbrace$$", style="text-indent:20px;font-size:70%;")),
                          br(),
                          numericInput('nzera', "Tamanho amostral", value=50,step=50, width = "50%"),
                          h5(tags$b("Vetores de médias")),
                          column(width = 6, align="center",
                                 style='border-right: 1px solid silver',#lightgray; margin-top:-.3em',
                                 numericInput('mu11',  withMathJax("$$\\mu_1$$"), value=0),
                                 numericInput('mu12', NULL, value=0)
                          ),
                          column(width=6, align="center",
                                 style='padding-bottom:9px;',
                                 numericInput('mu21', withMathJax("$$\\mu_2$$"), value=2),
                                 numericInput('mu22', NULL, value=2)
                          ),
                          HTML("<hr>"), br(),
                          h5(tags$b("Matriz de covariância")),
                          column(width = 6, align="center",
                                 numericInput('mx11', NULL, value=1),
                                 numericInput('mx12', NULL, value=0)
                          ),
                          column(width=6, align="center",
                                 numericInput('mx21', NULL, value=0),
                                 numericInput('mx22', NULL, value=1)
                          ),
                          textOutput("mpd"),
                          tags$head(tags$style("#mpd{color: green;font-size: 13px}")),
                          br(),
                          numericInput('seed', "Semente", value=333, width = "30%")#,
                          # hr(),
                          # actionButton("run", "Atualize")
                 )
    ),
    mainPanel(
        tabsetPanel(type="tabs",
                    tabPanel(tags$b("K-means"),br(), plotOutput('plot1')),
                    tabPanel(tags$b("Verdadeiro"),br(), plotOutput('trueplot'))
                    
        )
        
    )
)

server <- function(input, output, session) {
    
    output$trueplot <- renderPlot({
        n<-input$nzera
        mu1 <- c(input$mu11, input$mu12)
        mu2 <- c(input$mu21, input$mu22)
        Sigma <- matrix(c(input$mx11,input$mx12,input$mx21,input$mx22),2)
        
        set.seed(input$seed)
        xnorm1 <- mvrnorm(n,mu1,Sigma)
        xnorm2 <- mvrnorm(n,mu2,Sigma)
        
        amostra <- rbind(xnorm1, xnorm2)
        grupo <- as.data.frame(cbind(amostra, c(rep(1,n),rep(2,n))))
        
        
        ggplot(grupo, aes(V1,V2,color=factor(grupo[,3]))) + geom_point(size=2.5) +
            labs(title = " ")+
            scale_color_manual(values=c("skyblue3", "tomato1"))+
            theme(legend.position="none")
    })
    
    output$mpd <- renderText({
        if(matrixcalc::is.positive.definite(matrix(c(input$mx11,input$mx12,input$mx21,input$mx22),2))==TRUE) paste("Matriz positiva definida!")
    })
    
    output$plot1 <- renderPlot({
        n<-input$nzera
        mu1 <- c(input$mu11, input$mu12)
        mu2 <- c(input$mu21, input$mu22)
        Sigma <- matrix(c(input$mx11,input$mx12,input$mx21,input$mx22),2)
        
        set.seed(input$seed)
        xnorm1 <- mvrnorm(n,mu1,Sigma)
        xnorm2 <- mvrnorm(n,mu2,Sigma)
        
        amostra <- rbind(xnorm1, xnorm2)
        grupo <- cbind(amostra, c(rep(1,n),rep(2,n)))
        
        k<-sample(1:(n*2), n*2)
        amostra3 <- grupo[k,]
        
        # particao inicial
        centros <- sample(1:(n*2),2)
        centro1 <- amostra[centros[1],]
        centro2 <- amostra[centros[2],]
        
        tidy <- data.frame(amostra3,NA)
        colnames(tidy) <- c("x1","x2","true","estimado")
        
        # inicializacao
        pare <- FALSE
        while(pare == FALSE){
            
            # distancia euclidiana
            for(i in 1:(n*2)){
                euclidist1 <- sqrt((tidy[i,1]-centro1[1])^2 + (tidy[i,2]-centro1[2])^2)
                euclidist2 <- sqrt((tidy[i,1]-centro2[1])^2 + (tidy[i,2]-centro2[2])^2)
                
                # atribuicao aos clusters
                ifelse(euclidist1 < euclidist2, tidy[i,4]<-1, tidy[i,4]<-2)
            }
            
            c1old <- centro1
            c2old <- centro2
            
            # novos centroids
            G1 <- tidy[tidy$estimado==1,1:2]
            centro1 <- colMeans(G1)
            
            G2 <- tidy[tidy$estimado==2,1:2]
            centro2 <- colMeans(G2)
            
            # condicao de parada
            if(identical(c1old, centro1) & identical(c2old, centro2)) pare <- TRUE
            # Sys.sleep(1)
        }
        
        anal <- mean(tidy$true==tidy$estimado)
        if(anal < 0.3){
            TIDY <- ifelse(tidy$estimado==1,2,1)
            anal <- mean(tidy$true==TIDY)
        }
        print(ggplot(tidy, aes(x1,x2,color=factor(estimado))) + geom_point(size=2.5) + scale_color_manual(values=c("tomato1","skyblue3"))+
                  annotate("text",x=centro1[1],y=centro1[2],colour="red2",label="C1",size=8) +
                  annotate("text",x=centro2[1],y=centro2[2],colour="blue2",label="C2",size=8) +
                  labs(title = as.character(paste("Acurácia:", anal)),size=3) +
                  theme(legend.position="none"))
    })
    
}

shinyApp(ui = ui, server = server)
