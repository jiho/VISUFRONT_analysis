#
#      Deal with files and variables used by zooprocess and PlanktonIdentifier
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

read.pid <- function(file, metadata=FALSE) {
  #
  # Read .pid (or _dat1.txt) files
  #
  # file      file path
  # metadata  boolean, wether to extract the metadata from the header
  #           and store it in the attributes of the returned data.frame

  library("stringr")

  # read every line as text
  d <- scan(file, what="character", skip=1, sep="\n", quiet=T, fileEncoding="ISO-8859-1", encoding="UTF-8")

  # get line number where the data table starts
  dataIdx <- which(str_detect(d, fixed("[Data]")))

  # read data table
  # NB: we skipped the first line in the scan() call above, so the line numbers in the file are those computed here + 1, hence the dataIdx+1
  dt <- read.table(file, skip=dataIdx+1, sep=";", header=T, as.is=T)
  names(dt)[1] <- "Item"

  # make validation column name uniform
  names <- names(dt)
  names <- str_replace_all(names, "pred_valid_Id_", "Valid")
  # force latest validation to have the name "Valid"
  validationColumns <- which(str_detect(names, "Valid"))
  names[max(validationColumns)] <- "Valid"
  names(dt) <- names

  if (metadata) {
    # remove pseudo blank lines
    d <- d[d!=" "]
    # get line number where the data table starts
    dataIdxWithoutBlanks <- which(str_detect(d, fixed("[Data]")))

    # read meta data in the header
    m <- list()
    count <- 0
    # walk every line and test for patterns
    for (i in 1:(dataIdxWithoutBlanks-1)) {
      if (str_detect(d[i], "^\\[")) {
        # item is group title with format '[title]'
        # create a new sub-list
        count <- count+1
        m[[count]] <- list()
        names(m)[count] <- str_replace_all(d[i], "\\[|\\]", "")
      } else {
        # item is meta data with format 'name = value'
        # extract name and value
        line <- str_split(d[i], "=")[[1]]
        metaName <- line[1]
        metaData <- str_trim(line[2])
        # detect numbers with regular expressions
        if (!is.na(metaData)) {
          if (str_detect(metaData, "^-*[0-9]+\\.?[0-9]*$")) {
            metaData <- as.numeric(metaData)
          }
        }
        # store it
        m[[count]][metaName] <- metaData
      }
    }

    # attach meta data as attribute
    attr(dt, "meta") <- m
  }

  return(dt)
}


read.pids.in.project <- function(project, ...) {
  #
  # Walk the hierachy of a zooprocess project and read the pid + dat1.txt files for all "profiles" in the project merging the information in both
  # pid = image characteristics
  # dat1 = identification
  #
  # project   path to the project directory

  library("stringr")
  library("plyr")

  # find pid and dat1 files in project
  # NB: one pid/dat1 per profile
  pidFiles <- list.files(str_c(project, "/PID_process/Pid_results"), pattern="pid$", full=TRUE)

  dat1Files <- str_replace(pidFiles, "pid$", "txt")
  dat1Files <- str_replace(dat1Files, fixed("/Pid_results/"), "/Pid_results/Dat1_validated/")

  # read all files
  D <- ldply(1:length(pidFiles), function(i, pid, dat, ...) {
    # read pid
    p <- read.pid(pid[i], ...)

    # look for dat1 and read identifications from it when present
    if (file.exists(dat[i])) {
      d <- read.pid(dat[i], ...)
      if (nrow(d) != nrow(p)) {
        stop("dat1.txt and dat1.pid do not have the same number of lines. Something is wrong")
      }
      p$Pred <- d$Pred
      p$Valid <- d$Valid
    } else {
      warning("No dat1.txt file for pid: ", basename(pid[i]), "\nIdentifications (prediction and validation) will be absent", call.=FALSE)
    }

    return(p)
  }, pid=pidFiles, dat=dat1Files, ...)

  return(D)
}


read.learning.set.in.project <- function(project) {
  #
  # Walk the hierachy of a zooprocess project and read the pid file corresponding to the learning set
  #
  # project   path to the project directory

  library("stringr")

  learningSetFiles <- list.files(str_c(project, "/PID_process/Learning_set/"), pattern="pid$")
  if (length(learningSetFiles) > 1) {
    choice <- menu(learningSetFiles, graphics=FALSE)
    learningSetFile <- learningSetFiles[choice]
  } else {
    learningSetFile <- learningSetFiles
  }
  learn <- read.csv(str_c(project, "/PID_process/Learning_set/", learningSetFile), sep=";", skip=1, stringsAsFactors=FALSE)
  names(learn)[1] <- "Item"

  return(learn)
}


get.prediction.variables <- function(x, var.removed=c("X", "Y", "XM", "YM", "BX", "XMg5", "YMg5", "Width", "Height", "Angle", "XStart", "YStart", "Compentropy", "Compmean", "Compslope", "CompM1", "CompM2", "CompM2")) {
  #
  # Select all variables useful for prediction and compute derived variables
  #
  # x             data.frame from which to select/compute the variables, usually read with read.pid
  # var.removed   names of variables to remove

  library("stringr")
  library("plyr")

  # select only variables useful for prediction
  # = remove variables provided in the arguments (usually variables denoting position or so, which are not meaningful for prediction) + identification variables
  x <- x[, ! names(x) %in% c(var.removed, "X.Item", "Item", "Tag", "Ident", "Status", "Pred", "Valid", "Label")]
  x <- x[, ! str_detect(names(x), "Valid")]

  # check that all are numeric
  classes <- laply(x, class)
  nonNum <- which(! classes %in% c("numeric", "integer"))
  if (length(nonNum) != 0) {
    warning("Column(s) ", str_c(names(x)[nonNum], collapse=", "), " is/are not numeric and will cause trouble during the prediction")
  }

  # add derived variables (see PkID for all possibilities)
  # TODO not all possibles derived variables are computed here, code a better mechanism to specify them
  x$ESD <- 2 * sqrt(x$Area / pi)
  x$Elongation <- x$Major / x$Minor
  x$Range <- x$Max - x$Min
  x$MeanPos <- (x$Mean - x$Max) / (x$Max - x$Min)
  x$CV <- 100 * (x$StdDev / x$Mean)
  x$SR <- 100 * (x$StdDev / x$Range)
  x$PerimFeret <- x$Perim. / x$Feret
  x$PerimMaj <- x$Perim. / x$Major
  x$Circexc <- (4 * pi * x$Area_exc) / x$Perim.^2

  return(x)
}


presort.category <- function(project, category) {
  #
  # Allow pre-sorting images in a given category by their probability of being correctly sorted
  #
  # project   path to the root of the zooprocess project
  # category  name of the category to presort
  #
  # project="uvp5_sn999_visufront_2013_isiis_cc4_05"
  # category="noise"

  library("stringr")
  library("plyr")
  suppressPackageStartupMessages(library("randomForest", quietly=TRUE))

  # checks
  if ( ! file.exists(project) ) {
    stop("Cannot find project ", project)
  }

  # get prediction folder
  pred <- list.files(str_c(project, "/PID_process/Sorted_vignettes/"))
  categoryDir <- str_c(project, "/PID_process/Sorted_vignettes/", pred, "/", category)
  if ( ! file.exists(categoryDir) ) {
    stop("Cannot find predicted images for category ", category)
  }


  message("Learn identifications from learning set")
  # read learning set
  learn <- read.learning.set.in.project(project)
  learnIDs <- learn$Ident
  learnVars <- get.prediction.variables(learn)

  # fit model on learning set
  m <- randomForest(x=learnVars, y=factor(learnIDs), ntree=300)


  message("Read image characteristics for images in \"", category, "\"")
  # get profile id(s)
  metaFile <- list.files(str_c(project, "/meta/"), full=TRUE)
  meta <- read.csv(metaFile, sep=";", stringsAsFactors=FALSE)
  profile <- meta$profileid
  # NB: there could be several, so we should make a loop or something at this point
  #     but so far we have only used one and we will only code for this situation

  # get all images in category of interest
  imgs <- list.files(categoryDir, pattern="jpg$")

  # get corresponding numbers
  imgNbs <- imgs
  imgNbs <- str_replace(imgNbs, fixed(str_c(profile, "_")), "")
  imgNbs <- str_replace(imgNbs, fixed(".jpg"), "")
  imgNbs <- as.numeric(imgNbs)
  if (any(is.na(imgNbs))) {
    stop("Cannot find image number in some images")
  }

  # read pid
  pid <- read.pid(str_c(project, "/PID_process/Pid_results/", profile, "_dat1.pid"), metadata=FALSE)

  # get corresponding lines
  pidCat <- pid[imgNbs,]

  message("Predict relevance score for those images")
  # predict probability to be in each category from the model fitted on the learning set
  pidCatVars <- get.prediction.variables(pidCat)
  probs <- predict(m, newdata=pidCatVars, type="prob")
  # get the proba to be in the category of interest
  pCat <- probs[,category]


  message("Pre-pend images names with relevance score")
  # create names
  source <- str_c(pidCat$Label, "_", row.names(pidCat), ".jpg")
  n <- 5    # number of digits to consider
  score <- formatC(pCat*10^5, format="d", width=5, flag=0)
  dest <- str_c(score, "--", source)

  source <- str_c(categoryDir, source, sep="/")
  dest <- str_c(categoryDir, dest, sep="/")

  # file.rename(source, dest)
  # !! looooooong !!

  nbImgs <- length(source)
  l_ply(1:nbImgs, function(i, s, d) {
    system(str_c("mv \"", s[i], "\" \"", d[i], "\""))
  }, s=source, d=dest, .progress="text")
  # seemd faster and has progress bar, yay!

  return(invisible(categoryDir))
}


revert.presorting <- function(project) {
  #
  # Presorting images is possible by pre-pending the name with the probability to be correctly sorted
  # For zooprocess to work correctly after manual sorting is done, the names must be the original ones
  # This removes the pre-pended probabilities for all images in the given project
  #
  # project   project path

  library("stringr")
  library("plyr")

  # checks
  if ( ! file.exists(project) ) {
    stop("Cannot find project ", project)
  }

  # list images to which the proba was added at the beginning of the name
  vignettesDir <- str_c(project, "/PID_process/Sorted_vignettes/")
  imgs <- list.files(vignettesDir, pattern="^[0-9]*--.*\\.jpg$", recursive=TRUE)
  nbImgs <- length(imgs)

  if (nbImgs == 0) {
    message("No images names to revert in project ", project)
  } else {
    message("Reverting ", nbImgs, " images to their original name")
    # remove the proba from the name
    imgsClean <- str_replace(imgs, "[0-9]*--", "")
    # make the changes
    l_ply(1:nbImgs, function(i) {
      system(str_c("mv \"", vignettesDir, "/", imgs[i], "\" \"", vignettesDir, "/", imgsClean[i], "\""))
    }, .progress="text")
  }

  return(invisible(nbImgs))
}

