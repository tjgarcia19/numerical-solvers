#***********************************************************************************************************************
#   CMSC 150: Numerical and Symbolic Computation Project
#   Created by Trestan Janos G. Garcia
#
#   This application will solve Polynomial Regression, Quadratic Spline Interpolation and Simplex Method.
#
#***********************************************************************************************************************
source("solvers.R")

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   tags$html(tags$body(style="background-color:#22272e")),#background color of the entire app
   
   # Application title
   titlePanel(tags$h1("CSMC 150 Solvers", style="color:#c9cbd0")),
   tabsetPanel(
     tabPanel(tags$h5("Polynomial Regression Solver", style="color:#226096"),
              sidebarLayout(
                sidebarPanel(
                  tags$style(".well{background-color:#484d56}"),
                  fileInput(inputId = "file1", label = "Enter CSV File"),
                  textInput(inputId="order", label = "Enter the order of the polynomial"),
                  textInput(inputId="number1", label = "Enter a number"),
                  actionButton(inputId = "solve1", label="Solve"),
                  tableOutput("points1"),
                  width = 3
                  
                ),
                  
                mainPanel(
                  plotOutput("plotPts1"),
                  plotOutput("fxn1"),
                  verbatimTextOutput("fxnstr"),
                  tags$style("#fxnstr{background-color:#8f939d}"),
                  verbatimTextOutput("estimate1"),
                  tags$style("#estimate1{background-color:#8f939d}"),
                  width = 9
                )
              )
            ),
     tabPanel(tags$h5("Quadratic Spline Interpolation Solver", style="color:#226096"), 
              sidebarLayout(
                sidebarPanel(
                  fileInput(inputId = "file2", label = "Enter CSV File"),
                  textInput(inputId="number2", label = "Enter a number"),
                  actionButton(inputId = "solve2", label="Solve"),
                  tableOutput("points2"),
                  width = 3
                  ),
                
                mainPanel(
                  plotOutput("plotPts2"),
                  plotOutput("fxn2"),
                  verbatimTextOutput("fxns"),
                  tags$style("#fxns{background-color:#8f939d}"),
                  verbatimTextOutput("estimate2"),
                  tags$style("#estimate2{background-color:#8f939d}"),
                  width = 9
                )
              )
            ),
     tabPanel(tags$h5("Simplex Solver", style="color:#226096"),
              sidebarLayout(
                sidebarPanel(
                  textInput(inputId = "supply1", label="Supply in Denver"),
                  textInput(inputId = "supply2", label="Supply in Phoenix"),
                  textInput(inputId = "supply3", label="Supply in Dallas"),
                  textInput(inputId = "demand1", label="Demand in Sacramento"),
                  textInput(inputId = "demand2", label="Demand in Salt Lake"),
                  textInput(inputId = "demand3", label="Demand in Albuquerque"),
                  textInput(inputId = "demand4", label="Demand in Chicago"),
                  textInput(inputId = "demand5", label="Demand in New York"),
                  textInput(inputId = "cost1", label="Denver to Sacramento"),
                  textInput(inputId = "cost2", label="Denver to Salt Lake"),
                  textInput(inputId = "cost3", label="Denver to Albuquerque"),
                  textInput(inputId = "cost4", label="Denver to Chicago"),
                  textInput(inputId = "cost5", label="Denver to New York"),
                  textInput(inputId = "cost6", label="Phoenix to Sacramento"),
                  textInput(inputId = "cost7", label="Phoenix to Salt Lake"),
                  textInput(inputId = "cost8", label="Phoenix to Albuquerque"),
                  textInput(inputId = "cost9", label="Phoenix to Chicago"),
                  textInput(inputId = "cost10", label="Phoenix to New York"),
                  textInput(inputId = "cost11", label="Dallas to Sacramento"),
                  textInput(inputId = "cost12", label="Dallas to Salt Lake"),
                  textInput(inputId = "cost13", label="Dallas to Albuquerque"),
                  textInput(inputId = "cost14", label="Dallas to Chicago"),
                  textInput(inputId = "cost15", label="Dallas to New York"),
                  actionButton(inputId = "solve3", label="Solve"),
                  checkboxInput("initTab", "Show Initial Tableau"),
                  checkboxInput("perItab", "Show Tableau Per Iteration"),
                  checkboxInput("perIsoln", "Show Basic Solution Per Iteration"),
                  width = 3
                ),
                
                mainPanel(
                  verbatimTextOutput("warning"),
                  tags$style("#warning{background-color:#cf6679}"),
                  tags$div(tableOutput("init_tab"),style="background-color:#8f939d"),
                  verbatimTextOutput("tabPerIteration"),
                  tags$style("#tabPerIteration{background-color:#8f939d}"),
                  verbatimTextOutput("solnPerIteration"),
                  tags$style("#solnPerIteration{background-color:#8f939d}"),
                  tags$div(tableOutput("solution"),style="background-color:#8f939d"),
                  width = 9
            )
         )
      )
   )
)

server <- function(input, output) {
  
  #REGRESSION*******************************************************************************************************************************
  fileLoad = reactive({loadCSV(input$file1[1,4])})
  regression = reactive({PolynomialRegression(as.integer(input$order), fileLoad()$points)})
  
  observeEvent(input$solve1,{
    output$points1<-renderTable(fileLoad()$xyTab, digits = 0, width="100%", align = "c")
    output$plotPts1 <- renderPlot({
      par(bg="#28416e")
      plot(fileLoad()$points[[1]], fileLoad()$points[[2]], xlab = "x", ylab ="y",pch=20, col="red", main = "Plot of x and y");
    })
    output$fxn1 <- renderPlot({
      isolate({
        par(bg="#28416e")
        curve(regression()$polynomial_function(x),from=fileLoad()$points[[1]][1], to=fileLoad()$points[[1]][length(fileLoad()$points[[1]])],
              xlab="x", ylab="y",col="blue", main="Plot of the Function")
      })
    })
    output$fxnstr <- renderPrint({
      isolate({
        cat("Regression Function:\n")
        cat(regression()$polynomial_string, "\n")
      })
    })
    output$estimate1 <- renderPrint({
      isolate({
        cat("Estimated Value:", toString(regression()$polynomial_function(as.integer(input$number1))))
      })
    })
  })
  
  #SPLINE***********************************************************************************************************************************
  fileLoad2 = reactive({loadCSV(input$file2[1,4])})
  quadSpline = reactive({QuadraticSpline(as.integer(input$number2), fileLoad2()$points)})
  
  observeEvent(input$solve2,{
    output$points2<-renderTable(fileLoad2()$xyTab, digits = 0, width="100%", align = "c")
    output$plotPts2 <- renderPlot({
      par(bg="#28416e")
      plot(fileLoad2()$points[[1]], fileLoad2()$points[[2]], xlab = "x", ylab ="y",pch=20, col="red", main = "Plot of x and y");
    })
    output$fxn2 <- renderPlot({
      isolate({
        par(bg="#28416e")
        curve(quadSpline()$estimating_fxn(x),from=fileLoad2()$points[[1]][1], to=fileLoad2()$points[[1]][length(fileLoad2()$points[[1]])],
              xlab="x", ylab="y",col="blue", main="Plot of the Estimating Function")
      })
    })
    #output$fxns <- renderPrint(unlist(quadSpline()$interval_fxns))
    output$fxns <- renderPrint({isolate({
      i=1
      while(i<=quadSpline()$num_interval){
        cat("Interval",i,":",quadSpline()$interval_fxns[i],"\n")
        i=i+1
      }
    })
    })
    output$estimate2 <- renderPrint({ isolate({ cat("Estimated Value: ",toString(quadSpline()$estimate)) }) })
  })
  
  
  #SIMPLEX*************************************************************************************************************************************
  supplies <- reactive(list(as.integer(input$supply1), as.integer(input$supply2), as.integer(input$supply3)))
  demands <- reactive(list(as.integer(input$demand1), as.integer(input$demand2), as.integer(input$demand3), as.integer(input$demand4), as.integer(input$demand5)))
  costs <- reactive(list(as.integer(input$cost1),as.integer(input$cost2),as.integer(input$cost3),as.integer(input$cost4),as.integer(input$cost5),
                         as.integer(input$cost6),as.integer(input$cost7),as.integer(input$cost8),as.integer(input$cost9),as.integer(input$cost10),
                         as.integer(input$cost11),as.integer(input$cost12),as.integer(input$cost13),as.integer(input$cost14),as.integer(input$cost15)))
  simplex<- reactive(Simplex(unlist(costs()), unlist(supplies()), unlist(demands())))
  observe({
    if(input$initTab){output$init_tab <- renderTable(simplex()$init_tableau, spacing = "s", digits = 0)
    }else output$init_tab <- renderTable(NA, colnames = FALSE);
    
    if(input$perItab){
      output$tabPerIteration <- renderPrint({
        i=1
        while(i<=simplex()$num_iteration){
          cat("Iteration",i,":\n")
          print(simplex()$tableau_iteration[[i]])
          i=i+1
        }
      });
    }else output$tabPerIteration <- renderPrint({cat("Tableau per Iteration is hidden.")});
    
    if(input$perIsoln){
      output$solnPerIteration <- renderPrint({
      i=1
      while(i<=simplex()$num_iteration){
        cat("Iteration",i,":\n")
        print(simplex()$soln_iteration[[i]])
        i=i+1
      }
    })}else output$solnPerIteration <- renderPrint({cat("Basic Solution per Iteration is hidden.")});
  })
  
    observeEvent(input$solve3,{
      if(!simplex()$optimize){ output$warning<-renderText("ERROR: No Feasible Solution")
      }else output$warning<-renderText("The cost has been optimized!");
      output$solution <- renderTable(isolate({simplex()$solution}), rownames = TRUE, digits = 0)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
