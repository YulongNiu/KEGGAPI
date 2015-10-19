##' Get a R matrix object if the weblink returned as a matrix.
##'
##' If the web return a matrix, use this function to extract it as a R matrix object. An empty element is used to represent a "NA" in KEGG. So the empty web element is set to "NA" according the "ncol" parameter.
##' @title Get R matrix from weblink
##' @param url The weblink.
##' @param ncol The column number of the matrix.
##' @return A R matrix
##' @examples
##' pathUrl <- 'http://rest.kegg.jp/list/pathway'
##' keggPath <- webTable(pathUrl, ncol = 2)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl getURL
##' @export
##'
##' 
webTable <- function(url, ncol) {

  webPage <-getURL(url)

  ## transfer webpage into a matrix
  webMat <- unlist(strsplit(webPage, split = '\n', fixed = TRUE))
  webMat <- sapply(webMat, strsplit, split = '\t', fixed = TRUE)

  ## transfer empty web elements to NAs
  webMat <- sapply(webMat, function(x) {
    lenSub <- ncol - length(x)

    if (lenSub > 0) {
      x <- c(x, rep(NA, lenSub))
    } else {}

    return(x)
  })

  webMat <- matrix(unlist(webMat), ncol = ncol, byrow = TRUE)

  return(webMat)
}



##' Miscellaneous - Retreive grepexpr content
##'
##' Retreive the content from the output of grepexpr() function.
##' @title Retreive grepexpr content
##' @param s Input string.
##' @param g A list that output from grepexpr()
##' @return A string
##' @examples
##' testStr <- 'wwwtax.cgi?mode=Info&id=593907'
##' regList <- gregexpr('\\d+', testStr)
##' getcontent(testStr, regList[[1]])
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##'
##' 
getcontent <- function(s, g) {
  
  subText <- substring(s, g, g + attr(g, 'match.length') - 1)

  return(subText)
}
