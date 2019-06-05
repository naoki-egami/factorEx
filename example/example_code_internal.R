make.equation <- function(dvs, ivs, interactions = c()){
  return(paste(dvs, "~", paste(c(ivs, interactions), sep = "", collapse = " + ")))
}

prepDist <- function(refDist, cleanDist, include.facs = NULL,
                     name.ref = "Reference",
                     name.clean = NULL){
  if(!is.null(include.facs)){ # If user has supplied a list of factors to subset to, run this
    refDist <- refDist[refDist$factor %in% include.facs, ]
  }

  refFac <- paste(refDist$factor, refDist$levels, sep = "-")
  if(class(cleanDist) != "list"){useDist <- list(cleanDist)}
  if(class(cleanDist) == "list"){useDist <- cleanDist}
  finalDist <- list()
  for(i in 1:length(useDist)){
    cleanFac <- paste(useDist[[i]]$factor, useDist[[i]]$levels, sep = "-")
    finalDist[[i]] <- useDist[[i]][match(refFac, cleanFac), ]
  }
  if(is.null(name.clean)){print("ERROR: Need to give names of distributions to be cleaned")}
  if(!is.null(name.clean)){

    if(length(name.clean) != length(useDist)){print("ERROR: length of clean distribution name vector does not match number of distributions in list")}

    if(length(name.clean) == length(useDist)){
      names(finalDist) <- name.clean
      finalDist[[name.ref]] <- refDist
      return(finalDist)
    }
  }
}
