##' KEGG Database API - Get species KEGG/NCBI ID.
##'
##' Get the phylogenetic information of given species.
##' It supports both batch input and regular expression search.
##'
##' @title Get species list from KEGG.
##' @name getKEGGPhylo
##' @param speList The species list that is a vector like "c('hsa', 'eco')". The input "speList" should be consistent with the parameter "speType".
##' @param speType It supports five types: "KEGG", "Tnum", "regexpr", and "phylo".
##' KEGG type is a three or four letters, for exmaple "hsa" is the KEGG ID for Homo sapiens,
##' while the corresponding T number is "T01001".
##' The "regexpr" is used for regulare expression search with the Latin name ("Escherichia coli"), sub-species name ("K-12 MG1655"), and common name ("human").
##' The "phylo" uses phylogentic orders for search, and it supports "Domain" (either "Eukaryotes" or "Prokaryotes"), "Kingdom" ("Animals"), "phylum" ("Vertebrates"), and "class" ("Mammals").
##' But it does not support mixed "KEGG".
##' Attention: Mutiple KEGG species ID may correspond to one taxonomy ID, for exmaple "lph" and "lpo" to "91891"
##' @param whole Whether or not get the whole KEGG species list,
##' and the default value is FALSE.
##' @return Matrix of species information.
##' @examples
##' ## search species list from KEGG ID
##' getKEGGPhylo(c('hsa', 'eco'))
##' ## search species whose names include 'Escherichia coli'
##' getKEGGPhylo('Escherichia coli', speType = 'regexpr')
##' ## search species whose class is 'Mammals'
##' getKEGGPhylo('Mammals', speType = 'phylo')
##' ## get whole KEGG species information table
##' getKEGGPhylo(whole = TRUE)
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.kegg.jp/kegg/rest/}
##' @export
##'
##' 
getKEGGPhylo <- function(speList, speType = 'KEGG', whole = FALSE){

  if (!(speType %in% c('KEGG', 'regexpr', 'phylo', 'Tnum'))) {
    stop('"speType" now only supports "NCBI", "KEGG", "Tnum", "regexpr", and "phylo".')
  } else {}

  ## get whole KEGG species information
  speAnno <- webTable('http://rest.kegg.jp/list/organism', ncol = 4)

  colnames(speAnno) <- c('TID', 'KEGGID', 'LatinName', 'Phylo')

  if (whole) {
    return(speAnno)
  } else {
    if (speType == 'Tnum') {
      selSpeAnno <- speAnno[speAnno[, 1] %in% speList, , drop = FALSE]
    }
    else if (speType == 'KEGG') {
      selSpeAnno <- speAnno[speAnno[, 2] %in% speList, , drop = FALSE]
    }
    else if (speType == 'regexpr') {
      selSpeAnno <- speAnno[grep(speList, speAnno[, 3]), , drop = FALSE]
    }
    else if (speType == 'phylo') {
      selSpeAnno <- speAnno[grep(speList, speAnno[, 4]), , drop = FALSE]
    }
  }

  return(selSpeAnno)

}


##' KEGG Database API - Get the KEGG orthology list.
##'
##' Get the KEGG orthology list by a given KEGG KO ID.
##' @title Get KEGG orthology.
##' @param KOID The KEGG orthology ID.
##' @return A character vector of KEGG gene IDs
##' @examples
##' KOMat <- getKEGGKO('K02110')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.kegg.jp/kegg/rest/}
##' @export
##'
##' 
getKEGGKO <- function(KOID){

  ## get KO webpage
  url <- paste('http://rest.kegg.jp/link/genes/', KOID, sep = '')

  ## get KO list
  KOWebMat <- webTable(url, ncol = 2)

  geneIDs <- KOWebMat[, 2]
  
  return(geneIDs)

}


##' KEGG Database API - Get the whole pathway ID from KEGG database.
##'
##' Get the pathway ID and annoation of a given KEGG species ID.
##' @title List pathway of a given species ID
##' @param specID KEGSS species org code or T number , for example "hsa" or "T01001".
##' @return A matrix of pathway ID and annotation.
##' @examples
##' hasPath <- getKEGGPathAnno('hsa')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.kegg.jp/kegg/rest/}
##' @export
##' 
##'
getKEGGPathAnno <- function(specID){

  ## get KEGG pathway annotation list
  url <- paste('http://rest.kegg.jp/list/pathway/', specID, sep = '')
  pathAnno <- webTable(url, ncol = 2)

  colnames(pathAnno) <- c('pathID', 'Annotation')

  return(pathAnno)
}


##' KEGG Database API - Get the pathway and genes.
##'
##' Get the pathway and genes according to KEGG species ID.
##' @title List pathways and genes of a given KEGG species ID
##' @inheritParams getKEGGPathAnno
##' @return A List named with KEGG pathway IDs, and each element of the list contains the KEGG gene IDs.
##' @examples
##' hasPathGenes <- getKEGGPathGenes('hsa')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.kegg.jp/kegg/rest/}
##' @export
##'
getKEGGPathGenes <- function(specID){

  ## get KEGG pathway list
  url <- paste('http://rest.kegg.jp/link/', specID, '/pathway',  sep = '')

  pathAnno <- webTable(url, ncol = 2)

  pathInfo <- aggregate(pathAnno[, 2], list(pathAnno[, 1]), '[')
  pathList <- list()
  # transfer into list
  for (i in 1:nrow(pathInfo)){
    pathList[[i]] <- as.character(pathInfo[i, 2][[1]])
  }
  names(pathList) <- pathInfo[, 1]

  return(pathList)
}


##' KEGG Database API - Get the whole KEGG IDs from one species.
##'
##' Get the KEGG protein ID list and annotation.
##' @title Get whole KEGG IDs and annotation
##' @inheritParams getKEGGPathAnno
##' @return A matrix of KEGG IDs and annotation
##' @examples
##' ## KEGG org cord
##' getProID('eco')
##'
##' \dontrun{
##' ## KEGG T number
##' getProID('T00007')}
##'
##' ## KEGG T number with empty elements
##' getProID('T10004')
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.kegg.jp/kegg/rest/}
##' @export
##'
##' 
getProID <- function(specID){
  
  ## get KEGG ID annotation list
  url <- paste('http://rest.kegg.jp/list/', specID, sep = '')

  ## transfer webpage into a matrix
  speIDAnno <- webTable(url, ncol = 2)

  return(speIDAnno)
}



##' KEGG Database API - Convert IDs between KEGG databases and outside databases
##'
##' Convert gene identifiers or chemical substance identifiers between KEGG databases and oursite outside databases. For gene identifiers, this API provides functions to convert IDs/databases between KEGG databases and ncbi-gi/ncbi-geneid/uniprot. For chemical substance identifiers, it converts IDs/databases between drug/compound/glycan and pubchem/chebi. This API doesn't convert IDs between outside database. For example, the convert between pubchem and chebi IDs is not allowed. Try to use mutiple CPUs when conver large number of IDs. 
##'
##' The IDs and database convert is controlled by the argument "convertType".
##' @title KEGG convert function
##' @param targetDB "targetDB" and "convertType" are set correspondingly.
##' "convertType" --> "database"
##' For gene convert from KEGG to outside databases.
##' "sourceEntry" --> KEGG organism code, or T number.
##' "targetDB" --> "ncbi-gi", "ncbi-geneid", or "uniprot".
##' For gene convert from outside databases to KEGG.
##' "sourceEntry" --> "ncbi-gi", "ncbi-geneid", or "uniprot".
##' "targetDB" --> KEGG organism code, or T number.
##' For chemical substance convert from KEGG to outside databases.
##' "sourceEntry" --> "drug", "compound", or "glycan".
##' "targetDB" --> "pubchem", or "chebi".
##' For chemical substance convert from outside databases to KEGG.
##' "sourceEntry" --> "pubchem", or "chebi".
##' "targetDB" --> "drug", "compound", or "glycan".
##' 
##' 
##' "convertType" --> "identity"
##' For gene convert
##' "sourceEntry" --> KEGG organism code, T number, "genes", "ncbi-gi", "ncbi-geneid", or "uniprot". "genes" is set to convert outside identities to KEGG when the organism code is not known.
##' "targetDB" --> A vector (length can be bigger than 1) of identities. 
##' For chemical substance convert
##' "sourceEntry" --> "drug", "compound", "glycan", "pubchem", or "chebi".
##' "targetDB" --> A vector (length can be bigger than 1) of identities.
##' 
##' @param sourceEntry see "targetDB"
##' @param convertType set to be "database" or "identity".
##' @inheritParams getKEGGGeneSeq
##' @return A matrix that the first column is "targetDB"
##' @examples
##' ## convert database from KEGG to outside databases.
##' convKEGG('ncbi-geneid', 'eco')
##' convKEGG('pubchem', 'drug')
##' 
##' \dontrun{
##' ## convert database from outside databases to KEGG.
##' convKEGG('smu', 'uniprot')
##' convKEGG('glycan', 'chebi')}
##'
##' ## convert identities from KEGG to outside database.
##' ## mutiple organism convert.
##' convKEGG('ncbi-gi', c('hsa:10458', 'ece:Z5100'), convertType = 'identity', n = 2)
##' convKEGG('pubchem', 'cpd:C00004', convertType = 'identity', n = 2)
##'
##' ## convert identities from outside databases to KEGG.
##' ## the organism code is unknown.
##' convKEGG('genes', 'ncbi-geneid:3113320', convertType = 'identity', n = 2)
##' convKEGG('genes', 'ncbi-gi:54293358', convertType = 'identity', n = 2)
##' @importFrom foreach foreach %dopar%
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom ParaMisc CutSeqEqu
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.kegg.jp/kegg/rest/keggapi.html}
##' @export
##'
##' 
convKEGG <- function(targetDB, sourceEntry, convertType = 'database', n = 1) {

  if (convertType == 'identity') {
    
    ## register multiple core
    registerDoParallel(cores = n)
    
    ## cut the every 10 sourceEntry
    cutMat <- CutSeqEqu(length(sourceEntry), 10)
    sourceSeq <- apply(cutMat, 2, function(x) {
      eachSeq <- paste(sourceEntry[x[1] : x[2]], collapse = '+')
      return(eachSeq)
    })

    urlSeq <- paste0('http://rest.kegg.jp/conv/', targetDB, '/', sourceSeq)
    
    convRes <- foreach(i = 1:length(urlSeq), .combine = rbind) %dopar% {
      convRes <- webTable(urlSeq[i], ncol = 2)
      ## deal with NA
      ## rbind(NULL, NULL) is NULL
      if (is.na(convRes[1, 1])) {
        convRes <- NULL
      } else {}
      
      return(convRes)
    }

    ## remove no result
    if (is.null(convRes)) {
      convRes <- matrix(rep(NA, 2), nrow = 1)
    } else {}

    ## stop multiple core
    stopImplicitCluster()
    
  } else {
    url <- paste0('http://rest.kegg.jp/conv/', targetDB, '/', sourceEntry)
    convRes <- webTable(url, ncol = 2)
  }

  return(convRes)
}

