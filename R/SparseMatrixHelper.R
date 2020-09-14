#https://slowkow.com/notes/sparse-matrix/

#' @param x A sparse matrix from the Matrix package.
#' @param file A filename that ends in ".gz".
writeMMgz <- function(batch, fileName) {
  data.table::fwrite(
    x = data.frame(i = as.data.frame(batch %>% select(rowId))$rowId,
                   j = as.data.frame(batch %>% select(covariateId))$covariateId,
                   x = as.data.frame(batch %>% select(covariateValue))$covariateValue),
    file = fileName,
    append = TRUE,
    sep = " ",
    row.names = FALSE,
    col.names = FALSE
  )
}

startWritingMMgz <- function(fileName, maxX, maxY, rowNum, mType = "real", removeIfFileExist = TRUE){
  if(removeIfFileExist){
    if (file.exists(fileName)) 
      #Delete file if it exists
      file.remove(fileName)
  }
  
  writeLines(
    c(
      sprintf("%%%%MatrixMarket matrix coordinate %s general", mType),
      sprintf("%s %s %s", maxX, maxY, rowNum)
    ),
    gzfile(fileName)
  )
}

readBigMM<-function(file){
  if (is.character(file)) 
    file <- if (file == "") 
      stdin()
  else file(file)
  if (!inherits(file, "connection")) 
    stop("'file' must be a character string or connection")
  if (!isOpen(file)) {
    open(file)
    on.exit(close(file))
  }
  scan1 <- function(what, ...) scan(file, nmax = 1, what = what, 
                                    quiet = TRUE, ...)
  if (scan1(character()) != "%%MatrixMarket") 
    stop("file is not a MatrixMarket file")
  if (!(typ <- tolower(scan1(character()))) %in% "matrix") 
    stop(gettextf("type '%s' not recognized", typ), domain = NA)
  if (!(repr <- tolower(scan1(character()))) %in% c("coordinate", 
                                                    "array")) 
    stop(gettextf("representation '%s' not recognized", repr), 
         domain = NA)
  elt <- tolower(scan1(character()))
  if (!elt %in% c("real", "complex", "integer", "pattern")) 
    stop(gettextf("element type '%s' not recognized", elt), 
         domain = NA)
  sym <- tolower(scan1(character()))
  if (!sym %in% c("general", "symmetric", "skew-symmetric", 
                  "hermitian")) 
    stop(gettextf("symmetry form '%s' not recognized", sym), 
         domain = NA)
  nr <- scan1(integer(), comment.char = "%")
  nc <- scan1(integer())
  nz <- as.integer(scan1(numeric()))
  checkIJ <- function(els) {
    if (any(els$i < 1 | els$i > nr)) 
      stop("readMM(): row\t values 'i' are not in 1:nr", 
           call. = FALSE)
    if (any(els$j < 1 | els$j > nc)) 
      stop("readMM(): column values 'j' are not in 1:nc", 
           call. = FALSE)
  }
  if (repr == "coordinate") {
    switch(elt, real = , integer = {
      els <- scan(file, nmax = nz, quiet = TRUE, what = list(i = numeric(), 
                                                             j = numeric(), x = numeric()))
      els$i <- as.integer(els$i)
      els$j <- as.integer(els$j)
      checkIJ(els)
      switch(sym, general = {
        new("dgTMatrix", Dim = c(nr, nc), i = els$i - 
              1L, j = els$j - 1L, x = els$x)
      }, symmetric = {
        new("dsTMatrix", uplo = "L", Dim = c(nr, nc), 
            i = els$i - 1L, j = els$j - 1L, x = els$x)
      }, `skew-symmetric` = {
        stop("symmetry form 'skew-symmetric' not yet implemented for reading")
        new("dgTMatrix", uplo = "L", Dim = c(nr, nc), 
            i = els$i - 1L, j = els$j - 1L, x = els$x)
      }, hermitian = {
        stop("symmetry form 'hermitian' not yet implemented for reading")
      }, stop(gettextf("symmetry form '%s' is not yet implemented", 
                       sym), domain = NA))
    }, pattern = {
      els <- scan(file, nmax = nz, quiet = TRUE, what = list(i = integer(), 
                                                             j = integer()))
      checkIJ(els)
      switch(sym, general = {
        new("ngTMatrix", Dim = c(nr, nc), i = els$i - 
              1L, j = els$j - 1L)
      }, symmetric = {
        new("nsTMatrix", uplo = "L", Dim = c(nr, nc), 
            i = els$i - 1L, j = els$j - 1L)
      }, `skew-symmetric` = {
        stop("symmetry form 'skew-symmetric' not yet implemented for reading")
        new("ngTMatrix", uplo = "L", Dim = c(nr, nc), 
            i = els$i - 1L, j = els$j - 1L)
      }, hermitian = {
        stop("symmetry form 'hermitian' not yet implemented for reading")
      }, stop(gettextf("symmetry form '%s' is not yet implemented", 
                       sym), domain = NA))
    }, complex = {
      stop("element type 'complex' not yet implemented")
    }, stop(gettextf("'%s()' is not yet implemented for element type '%s'", 
                     "readMM", elt), domain = NA))
  }
  else stop(gettextf("'%s()' is not yet implemented for  representation '%s'", 
                     "readMM", repr), domain = NA)
}

