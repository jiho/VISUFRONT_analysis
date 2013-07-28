#
#
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

# This part of the code is run once, when the server is started

library("stringr")
library("plyr")
library("ggplot2")
library("reshape2")
library("grid")
library("gridExtra")  # for grid.arrange
source("../lib_plot.R")

# This part of the code is run on every page load
shinyServer(function(input, output) {

  # Extract data based on the criterions specified by the used in the UI
  get.data <- reactive({
    # read ISIIS hydrological data for this transect
    d <- adply(input$transect, 1, function(transect) {
      read.csv(str_c("../transects/", transect, "/isiis.csv"))
    })
    names(d)[1] <- "transect"
    d$transect <- input$transect[d$transect]

    dm <- melt(d, id.vars=c("transect", "Depth.m.", input$dist, "down.up"), measure.vars=input$vars)

    # interpolate every variable
    di <- ddply(dm, ~transect+variable, function(x) {
      x <- na.omit(x[which(x$down.up=="up"),])
      xi <- interp.dist(x=x[,input$dist], y=x$Depth.m., z=x$value, duplicate="mean", x.step=input$xstep, y.step=input$ystep)
    })
    di <- rename(di, c("x"=input$dist, "y"="Depth.m."))

    # return raw and interpolated data
    list(dm=dm, di=di)
  })

  # Dynamically set plot height
  plotHeight <- function(){
    length(input$vars) * length(input$transect) * 200 + 50
  }

  # Generate the plot of the extracted data
  output$dataPlot <- renderPlot({

    d <- get.data()

    if ( nrow(d$dm)==0 | length(input$vars)==0 | all(is.na(d$dm$value)) ) {
      stop("No data, try again")

    } else {
      plots <- dlply(d$di, ~variable, function(x) {
        ggplot(x, aes_string(x=input$dist, y="-Depth.m.")) +
          # geom_point(aes(fill=value), shape=21, colour=NA, na.rm=T) +
          geom_tile(aes(fill=value), na.rm=T) +
          geom_contour(aes(z=value), colour="white", alpha=0.9, bins=5, na.rm=T) +
          scale_fill_gradientn(name=x$variable[1], colours=spectral()) +
          scale_x_continuous(expand=c(0,0)) +
          scale_y_continuous(expand=c(0,0)) +
          facet_grid(transect~.)
      })

      do.call(grid.arrange, c(plots,list(ncol=1)))
    }
  }, height=plotHeight)

})
