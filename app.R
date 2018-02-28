#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Circular Dichroism simulator"),
   
   #Plot output
   fluidRow(
     tabsetPanel(type = "tabs",
                 tabPanel("Data", plotOutput("CDtrace"),
                          verbatimTextOutput("text")),
                 tabPanel("More information",
                          includeHTML("CD.html")))
       
     )
   ,
   fluidRow(
     sidebarPanel(
        selectInput("mode",
                   "Which mode do you want to use?",
                   choices = list("Prediction" = "predict",
                                  "Display" = "display")
        )
        ,
                   submitButton("Submit")
     )
     
     ,
     mainPanel(
         #conditional panel 
         conditionalPanel(
           condition = "input.mode == 'predict'",
           checkboxInput("guides", 
                         "Show secondary structure guides?", 
                         value = FALSE, 
                         width = NULL)
           ,
           numericInput("helix",
                        "% alpha helix",
                        min = 0,
                        max = 100,
                        step = 5,
                        value = 50)
           ,
           numericInput("sheet",
                        "% beta sheet",
                        min = 0,
                        max = 100,
                        step = 5,
                        value = 50)
           ,
           numericInput("coil",
                        "% random coil",
                        min = 0,
                        max = 100,
                        step = 5,
                        value = 0)
           )
         ,
         conditionalPanel(
           condition = "input.mode == 'display'",
           selectInput("protein",
                       "Pick a protein to show predicted data for",
                       choices = list("Lysozyme (1LYD)" = "lyso",
                                      "Ubiquitin (1UBQ)" = "ub",
                                      "BST2 (3MQB)" = "bst",
                                      "Hemoglobin (2HHB)" = "hemo",
                                      "Antibody (1IGT)" = "ab")
           )
           ,
           "Simulated numbers generated in YASARA using the PDB IDs indicated"
           )
         
         
         )
       ),
   br(),
   br(),
   "Shiny app created by C.E. Berndsen, 2018",
   br(),
   "Simulator based on work by Abriata, J. Chem. Ed. (2011)", tags$a(href="https://pubs.acs.org/doi/full/10.1021/ed200060t", "REF"), "using data from Greenfield and Fasman, Biochemistry (1969)", tags$a(href="https://pubs.acs.org/doi/10.1021/bi00838a031", "REF"),
   br()
   
   )


# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$CDtrace <- renderPlot({
       if(input$mode == "predict") {
       #check to make sure fractions add to 100%          
       validate(
         need(input$helix/100 + input$sheet/100 + input$coil/100 == 1, 
              "% of helix, sheet, and coil must = 100%"
              )
       )
     
     
       #generate basis set
       CDdat <- data.frame(lambda = seq(190, 250, by = 0.2))
       CDdat <- CDdat %>% mutate(helix = 1*10^8 * (2230060.04151075*lambda^0 +
                                                   -100548.516559741*lambda^1 +
                                                   2037.18080475746*lambda^2 + 
                                                   -24.4244919907991*lambda^3 +
                                                   0.19190243015954*lambda^4 +
                                                   -0.00103245782924168*lambda^5 +
                                                   0.00000385211889091252*lambda^6 +
                                                   -9.84175959744622E-09*lambda^7 +
                                                   1.64786777298595E-11*lambda^8 +
                                                   -1.63282751503442E-14*lambda^9 +
                                                   7.27089674019501E-18*lambda^10)) %>% 
       mutate(beta = 1*10^8 * (-677807.330017282*lambda^0 +
                                 30975.2707887604*lambda^1 +
                                 -636.143263740698*lambda^2 + 
                                 7.73164864362657*lambda^3 +
                                 -6.15861633716145E-02*lambda^4 +
                                 3.35943314432255E-04*lambda^5 +
                                 -1.27092416556044E-06*lambda^6 +
                                 3.29272089581372E-09*lambda^7 +
                                 -5.59118151750062E-12*lambda^8 +
                                 5.61899389095424E-15*lambda^9 +
                                 -2.53794637234403E-18*lambda^10)) %>%
       mutate(coil = 1*10^8 * (-580939.072386969*lambda^0 +
                                 25845.2673351998*lambda^1 +
                                 -516.713088253122*lambda^2 + 
                                 6.1134023680003*lambda^3 +
                                 -4.74021175198809E-02*lambda^4 +
                                 2.51692531821056E-04*lambda^5 +
                                 -9.26824208397782E-07*lambda^6 +
                                 2.33714935193268E-09*lambda^7 +
                                 -3.86247107852678E-12*lambda^8 +
                                 3.77764956561175E-15*lambda^9 +
                                 -1.6603998403172E-18*lambda^10))
     
       #Predict spectrum based on user input  
       CDdat <- CDdat %>% mutate(prediction = (input$helix)/100*helix + (input$sheet/100)*beta + (input$coil/100)*coil)
      
      # draw the plot
       if(input$guides == FALSE) {
         ggplot(CDdat, aes(x = lambda, y = prediction/1000), color = "red") +
           geom_jitter(alpha = 0.4) +
           scale_x_continuous(breaks = seq(190, 250, by = 5)) +
           labs(x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum") +
           geom_hline(yintercept = 0) +
           ylim(-50, 75) +
           theme_classic() +
           theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 16, face = "bold"))
       }
       else {
         ggplot() +
           geom_jitter(data = CDdat, aes(x = lambda, y = prediction/1000), fill = "red", alpha = 0.4) +
           geom_line(data = CDdat, aes(x = lambda, y = helix/1000), color = "red") +
           geom_line(data = CDdat, aes(x = lambda, y = beta/1000), color = "green") +
           geom_line(data = CDdat, aes(x = lambda, y = coil/1000), color = "purple") +
           geom_hline(yintercept = 0) +
           scale_x_continuous(breaks = seq(190, 250, by = 5)) +
           labs(x = "wavelength (nm)", y = "Ellipticity", title = "Predicted CD spectrum") +
           annotate("text", x = 235, y = 25, label = "Helix", color = "red", size = 5) +
           annotate("text", x = 235, y = 20, label = "Beta Sheet", color = "green", size = 5) +
           annotate("text", x = 235, y = 15, label = "Random Coil", color = "purple", size = 5) +
           ylim(-50, 75) +
           theme_classic() +
           theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 16, face = "bold"))
        
       }
       }
     else {
       #Generate the wavelength values
       CDdat <- data.frame(lambda = seq(190, 250, by = 0.2))
       
       #Generate the basis set from Abriata, L., J. Chem. Educ., 2011, 88 (9), pp 1268–1273 and Davidson, B. and Fasman, G. D., Biochemistry 1967 6 (6) 1616-1629
       CDdat <- CDdat %>% mutate(helix = 1*10^8 * (2230060.04151075*lambda^0 +
                                                     -100548.516559741*lambda^1 +
                                                     2037.18080475746*lambda^2 + 
                                                     -24.4244919907991*lambda^3 +
                                                     0.19190243015954*lambda^4 +
                                                     -0.00103245782924168*lambda^5 +
                                                     0.00000385211889091252*lambda^6 +
                                                     -9.84175959744622E-09*lambda^7 +
                                                     1.64786777298595E-11*lambda^8 +
                                                     -1.63282751503442E-14*lambda^9 +
                                                     7.27089674019501E-18*lambda^10)) %>% 
         mutate(beta = 1*10^8 * (-677807.330017282*lambda^0 +
                                   30975.2707887604*lambda^1 +
                                   -636.143263740698*lambda^2 + 
                                   7.73164864362657*lambda^3 +
                                   -6.15861633716145E-02*lambda^4 +
                                   3.35943314432255E-04*lambda^5 +
                                   -1.27092416556044E-06*lambda^6 +
                                   3.29272089581372E-09*lambda^7 +
                                   -5.59118151750062E-12*lambda^8 +
                                   5.61899389095424E-15*lambda^9 +
                                   -2.53794637234403E-18*lambda^10)) %>%
         mutate(coil = 1*10^8 * (-580939.072386969*lambda^0 +
                                   25845.2673351998*lambda^1 +
                                   -516.713088253122*lambda^2 + 
                                   6.1134023680003*lambda^3 +
                                   -4.74021175198809E-02*lambda^4 +
                                   2.51692531821056E-04*lambda^5 +
                                   -9.26824208397782E-07*lambda^6 +
                                   2.33714935193268E-09*lambda^7 +
                                   -3.86247107852678E-12*lambda^8 +
                                   3.77764956561175E-15*lambda^9 +
                                   -1.6603998403172E-18*lambda^10))
       
         #Predict spectrum based on user input  
         CDdat <- CDdat %>% 
         mutate(lyso = 0.7*helix + 0.07*beta + 0.17*coil) %>%
         mutate(ub = 0.17*helix + 0.40*beta + 0.43*coil) %>%
         mutate(bst = 0.92*helix + 0*beta + 0.08*coil) %>%
         mutate(hemo = 0.74*helix + 0*beta + 0.17*coil) %>%
         mutate(ab = 0.05*helix + 0.46*beta + 0.36*coil)
         
         CDdat <- CDdat %>%
           select(1:4, input$protein)
         
         colnames(CDdat) <- c("lambda", "helix", "beta", "coil", "display")  
         
         ggplot() +
         geom_jitter(data = CDdat, aes(x = lambda, y = display/1000), fill = "red", alpha = 0.4) +
         geom_line(data = CDdat, aes(x = lambda, y = helix/1000), color = "red") +
         geom_line(data = CDdat, aes(x = lambda, y = beta/1000), color = "green") +
         geom_line(data = CDdat, aes(x = lambda, y = coil/1000), color = "purple") +
         geom_hline(yintercept = 0) +
         scale_x_continuous(breaks = seq(190, 250, by = 5)) +
         labs(x = "wavelength (nm)", y = "Ellipticity") +
         annotate("text", x = 235, y = 25, label = "Helix", color = "red", size = 7) +
         annotate("text", x = 235, y = 20, label = "Beta Sheet", color = "green", size = 7) +
         annotate("text", x = 235, y = 15, label = "Random Coil", color = "purple", size = 7) +
         ylim(-50, 75) +
         theme_classic() +
         theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 16, face = "bold"))
         
         
         
     }
})
   
   
   output$text <- renderText({
     if(input$mode == "predict") {
       paste("Spectrum showing", input$helix, "% alpha helix", input$sheet, "% beta sheet", input$coil, "% random coil")
       
     }
     else{
     label <- ifelse(input$protein == 'lyso', "70% helix, 7% beta strand, 17% random coil",
                     ifelse(input$protein == 'ub', "17% helix, 40% beta strand, 43% random coil",
                            ifelse(input$protein == 'bst', "92% helix, 0% beta strand, 8% random coil",
                                   ifelse(input$protein == 'ab', "5% helix, 46% beta strand, 36% random coil",
                                          ifelse(input$protein == 'hemo', "74% helix, 0% beta strand, 17% random coil", "N/A")))))
     
     paste("Structure contains", label)
     }
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

