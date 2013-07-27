#
#      Quick interface to explore ISIIS data
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

# This code is run the first time the page is loaded in a browser
# For subsequent requests, the code hasn't changed since the last page load, it is not re-run (it is probably cached somehow by Shiny)
transects <- list.files("../transects")


# Define UI for the application (menus, checkboxes, etc.)
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("VISUFRONT"),

  # Sidebar with selectors
  sidebarPanel(

    div(class="row-fluid",
      div(
        class="span6",
        checkboxGroupInput(inputId="transect",
                    label="Transect",
                    choices=transects,
                    selected=transects[1]
                    )
      ),
      div(
        class="span6",
        checkboxGroupInput(inputId="vars",
                           label="Variables",
                           choices=c("Temp.C.", "Salinity.PPT.", "Density", "Fluoro.volts.", "Oxygen.ml.l.", "Irrandiance.UE.cm."),
                           selected="Salinity.PPT."
                           )
      )
    ),
    selectInput(inputId="dist",
                label="Distance measure",
                choices=c("distanceFromStart", "distanceFromVlfr"),
                selected="distanceFromStart"
                ),
 
  
    sliderInput(inputId="xstep",
                label="Distance step (m)",
                min=100, max=2000,
                value=1000,
                step=100
                ),
    sliderInput(inputId="ystep",
                label="Depth step (m)",
                min=0.5, max=5,
                value=2.5,
                step=0.5
                ),

    div(style="height: 20px"),

    # Button to force manual reload of the plot
    submitButton("Redraw plot"),
    
    # make sure we finish the panel correctly
    div(style="height: 1px")
  ),

  # Show a plot of the selected data
  mainPanel(
    # verbatimTextOutput("sessionInfo"),
    plotOutput("dataPlot", height="auto")
    # textOutput("testText")
  )
))
