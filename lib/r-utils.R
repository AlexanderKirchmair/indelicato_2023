
## R helper functions

cutstr <- function (x, maxchar = 25, add = "...", add_incl = TRUE){
  ix <- nchar(x) > maxchar
  x[ix] <- substr(x[ix], 1, ifelse(add_incl, maxchar - nchar(add), 
                                   maxchar))
  x[ix] <- paste0(x[ix], add)
  x
}

baseext <- function(path, ...){
  path <- basename(path)
  ext <- gsub(x = path, pattern = ".*\\.", replacement = "")
  ext[!grepl(pattern = ".", x = path, fixed = TRUE)] <- ""
  ext
}

rownames2col <- function (data, col = id, keep = FALSE){
  col <- rlang::enquo(col)
  names <- rownames(data)
  if (is.null(names)) {
    warning("Warning: No rownames found.")
    return(data)
  }
  if (keep == FALSE) 
    rownames(data) <- NULL
  i <- rlang::as_name(col)
  if (i %in% colnames(data)) {
    data[, i] <- names
  }
  else {
    if (class(names) %in% unique(sapply(data, class))) {
      data <- cbind(names, data)
    }
    else {
      data <- data.frame(names, data)
    }
    colnames(data)[1] <- i
  }
  data
}


matScale <- function(data, rows = FALSE, cols = FALSE, FUN = scale, ...){
  data.org <- data
  ix <- sapply(as.data.frame(data), is.numeric)
  data <- data[, ix, drop = FALSE]
  names.org <- dimnames(data)
  if (rows == TRUE & cols == TRUE) 
    stop("Error: Do not scale rows and columns at once!")
  if (rows == TRUE) {
    data <- t(apply(data, 1, function(tmp) as.numeric(FUN(tmp, 
                                                          ...))))
  }
  if (cols == TRUE) {
    data <- apply(data, 2, function(tmp) as.numeric(FUN(tmp, 
                                                        ...)))
  }
  dimnames(data) <- names.org
  data <- cbind(data, data.org[, !ix, drop = FALSE])
  if (!is.null(colnames(data.org))) 
    data <- data[, colnames(data.org)]
  stopifnot(all.equal(dim(data.org), dim(data)))
  return(data)
}


roundup <- function (x, digits = 0, ...){
  ceiling(x * 10^digits, ...)/10^digits
}

rounddown <- function (x, digits = 0, ...){
  floor(x * 10^digits, ...)/10^digits
}


clusterData <- function (data, method = "hclust", rows = NULL, cols = NULL, 
          inf = NULL, na = NULL, ...){
  if (is.null(rows)) 
    rows <- nrow(data) < 1000
  if (is.null(cols)) 
    cols <- ncol(data) < 1000
  if (is.null(na)) 
    na <- ""
  if (is.null(inf)) 
    inf <- ""
  if (any(!is.finite(nat(data)))) {
    if (naf(inf == FALSE)) 
      data <- subInf(data)
    if (is.na(inf)) 
      data[!is.finite(data)] <- inf
  }
  if (any(is.na(data))) {
    if (is.numeric(na)) 
      data[is.na(data)] <- na
  }
  tmp <- list(rows = data, cols = t(data))
  res <- lapply(tmp[c(rows, cols)], function(data) {
    if (na == "omit") 
      data <- na.omit(data)
    clust <- NULL
    if (method == "hclust") 
      clust <- dendsort::dendsort(stats::hclust(stats::dist(data)))
    clust
  })
  res
}



writeTables <- function(data, file, rowNames = TRUE, adjwidths = TRUE, ...){
  
  if (!"list" %in% class(data)) {
    data <- setNames(list(data), gsub("\\..*$", "", basename(file)))
  }
  newnames <- cutstr(names(data), maxchar = 29)
  if (any(duplicated(newnames))) {
    newnames <- cutstr(names(data), maxchar = 26)
    newnames <- dedupl(newnames)
  }
  names(data) <- newnames
  wb <- openxlsx::createWorkbook()
  invisible(lapply(names(data), function(tmpname) {
    tmpdata <- data[[tmpname]]
    openxlsx::addWorksheet(wb, tmpname)
    openxlsx::writeData(wb, sheet = tmpname, x = tmpdata, 
                        rowNames = rowNames, ...)
    if (adjwidths == TRUE) {
      openxlsx::setColWidths(wb, sheet = tmpname, cols = 1:(ncol(tmpdata) + 
                                                              1), widths = "auto")
    }
  }))
  
  openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
}


dedup <- function (x, sep = "", ...){
  
  xorg <- x
  ix <- duplicated(x, ...)
  
  ix <- x %in% x[ix]
  
  i <- 1
  while (any(ix)) {
    x[ix] <- paste0(xorg[ix], sep, letters[i])
    ix <- duplicated(x, ...)
    i <- i + 1
  }
  
  stopifnot(!any(duplicated(x, ...)))
  x
}


nat <- function(x){
  x[is.na(x)] <- TRUE
  x
}

naf <- function(x){
  x[is.na(x)] <- FALSE
  x
}

