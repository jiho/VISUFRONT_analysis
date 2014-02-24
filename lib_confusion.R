#
#      Functions to compute and investigate confusion in automatic
#      classification
#
#  (c) Copyright 2014 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------


confusion_matrix <- function(predicted, true) {
  #
  # Compute a confusion matrix (i.e. contingency table)
  # pred    predicted category
  # true    true category

  t <- table(predicted=predicted, true=true)
  return(t)
}

fortify.table <- function(x) {
  #
  # Convert a contingency table (i.e. confusion matrix) into a data.frame
  #
  # x   contingency table, returned by function table()

  return(as.data.frame(x))
}

library("ggplot2")
autoplot.table <- function(x, norm="none", trans=NULL) {
  #
  # Plot a contingency table (i.e. confusion matrix) as a heatmap
  #
  # x     contingency table, returned by function table()
  # norm  normalisation method: none, by row ("rows") or by column ("columns")
  # trans function used to transform the counts in the contingency table

  # normalise data
  norm <- match.arg(norm, c("none", "rows", "columns"))
  variableName <- "Freq"
  if ( norm == "rows") {
    n <- rowSums(x)
    x <- x / n
    variableName <- "Freq\nby row"
  }
  if ( norm == "columns") {
    n <- colSums(x)
    x <- t(t(x) / n)
    variableName <- "Freq\nby column"
  }

  # make table into a data.frame
  x <- fortify.table(x)

  # transform frequencies
  if ( ! is.null(trans) ) {
    if ( is.function(trans) ) {
      x$Freq <- trans(x$Freq)
      fun <- deparse(substitute(trans))
      library("stringr")
      variableName <- str_c(fun, "(", variableName, ")")
    } else {
      stop("Cannot find function ", fun)
    }
  }

  # make the plot
  library("ggplot2")
  p <- ggplot(x) +
        geom_tile(aes_string(x=names(x)[2], y=names(x)[1], fill="Freq")) +
        coord_fixed(1) + labs(fill=variableName) +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0))
  return(p)
}

confusion_stats <- function(x, sort.by=NULL) {
  #
  # Confusion statistics (recall, precision, etc.)
  #
  # x    confusion matrix (table class with prediction as line and observations as columns)

  # reduce to common categories
  rowCats <- rownames(x)
  colCats <- colnames(x)
  if ( any( ! c(rowCats %in% colCats, colCats %in% rowCats) ) ) {
    warning("Confusion statistics can only be computed for categories present in both lines and columns of the confusion matrix. Reducing data to common categories")
    commonCats <- intersect(rowCats, colCats)
    x <- x[rowCats %in% commonCats, colCats %in% commonCats]
  }

  # check if matrix is square
  if ( nrow(x) != ncol(x)) {
    stop("The confusion matrix needs to be square")
  }

  # compute base stats
  (tp <- diag(x))              # true positive
  (fp <- rowSums(x) - tp)      # false positive
  (fn <- colSums(x) - tp)      # false negative
  (tn <- sum(x) - tp - fp -fn) # true negative

  # store it
  stats <- data.frame(tp, fp, fn)

  # define a formater for percentages
  format_percent <- function(x, precision=1) {
    round(x * 100, precision)
  }

  # precision = quantify how "pure" the identified signals are
  stats$precision <- format_percent(tp / (tp + fp))

  # recall = capacity to get signals of a given origin
  stats$recall <- format_percent(tp / (tp + fn))

  # F1 score = combination of precision and recall
  # http://en.wikipedia.org/wiki/F1_score
  # the higher the better
  stats$F1 <- with(stats, (2 * precision * recall) / (precision + recall))

  if ( ! is.null(sort.by) ) {
    sort.by <- match.arg(sort.by, names(stats))
    stats <- stats[order(stats[,sort.by]),]
  }

  return(stats)
}

confusion_plot <- function(p, true, pred, x=c("probability", "order"), y=c("count", "proportion")) {
  #
  # Plot confusion in each predicted category
  # The y axis is the proportion/number of images incorrectly classified
  # The x axis has all images of the category ordered according to their probability of being in that category; either the probability itself is represented or the order within the images in the category (normalised between 0 and 1)
  #
  # p     matrix or data.frame containing the probabilities of belonging to each category, for all images
  # true  vector of the true category of each image
  # pred  vector of the predicted category of each image (usually the one with the highest proba in p)
  # x     nature of the x axis (Cf above)
  # y     nature of the y axis (Cf above)
  #

  # check arguments
  p <- as.data.frame(p)

  x <- match.arg(x)
  y <- match.arg(y)

  # check number of records
  if ( nrow(p) != length(true) | nrow(p) != length(pred) ) {
    stop("probs, true, and pred must have the same number of elements")
  }

  # check compatibilities of categories/levels
  predCategories <- sort(unique(na.omit(pred)))
  trueCategories <- sort(unique(na.omit(true)))
  probsCategories <- sort(names(p))
  if ( (! all(predCategories == trueCategories)) ) {
    stop("Levels in true and pred should be the same")
  }
  if ( (! all(predCategories == probsCategories)) ) {
    stop("Levels in true and pred should be the same as columns in probs")
  }

  # for each record, get the proba of being in the predicted class
  p$recordNb <- 1:nrow(p)
  pm <- melt(p, id.vars="recordNb", value.name="probability", variable.name="pred")
  d <- data.frame(recordNb=1:length(true), true=true, pred=pred)
  d <- join(na.omit(d), pm, by=c("recordNb", "pred"))

  # detect (in)correct predictions
  d$incorrect <- ( is.na(d$true) | ( d$true != d$pred ) )
  # NB: records with NA as the true category are not predictable by the model but sould still be counted as incorrect because they indeed are

  # within each predicted category, sort records by their probability of being in that category
  d <- ddply(d, ~pred, function(X) {
    # sort data by probability of being the class
    X <- arrange(X, probability, decreasing=TRUE)
    # compute the number of mistakes remaining for each level of probability
    # i.e. the cumulated number of mistakes when starting from the highest probability
    X$count <- cumsum(X$incorrect)
    X$proportion <- X$count / nrow(X)
    # compute the order of data from lowest to highest probability
    X$order <- seq(from=100, to=0, length.out=nrow(X))

    return(X)
  })

  # TODO add counts per category in facet labels

  plot <- ggplot(d) + geom_path(aes_string(x=x, y=y)) + facet_wrap(~pred)
  plot <- plot + labs(y=str_c(y, " of incorrect predictions"), x=str_c("prediction ", x))

  return(plot)
}


#--------------------------------------------------------------------------
# Test data

#
# set.seed(123)
# n <- 5
# cat <- letters[1:n]
# pred <- cat[ceiling(runif(50)*n)]
# n <- 5
# cat <- letters[1:n]
# true <- cat[ceiling(runif(50)*n)]
# x <- confusion_matrix(pred, true)
#
# norm <- rowSums(x)
# rowSums(x/norm)
#
# norm <- colSums(x)
# colSums(t(t(x)/norm))
#
# autoplot(x)
# autoplot(x)
# autoplot(x, norm="row")
# autoplot(x, norm="col")
# autoplot(x, norm="col", trans=log)
# autoplot(x, norm="col", trans=sqrt)
# autoplot(x, trans=sqrt)
#
