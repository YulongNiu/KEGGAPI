##' KEGG Database Additional API - Get the NCBI taxonomy ID from a given KEGG ID
##'
##' NCBI taxonomy ID is used as unique ID accoss KEGG and BioCyc databases. This functions is used to get the corresponding NCBI Taxonomy ID from KEGG.
##' @title Get NCBI Taxonomy ID From KEGG ID
##' @param KEGGID The KEGG support multiple species ID, for example c('hsa', 'eco').
##' @param n The number of CPUs or processors, and the default value is 4.
##' @return The corresponding NCBI Taxonomy ID in character vector.
##' @examples
##' # get human and Ecoli NCBI taxonomy ID with 2 threads
##' transPhyloKEGG2NCBI(c('hsa', 'eco', 'ath', 'smu'), n = 2)
##' 
##' \dontrun{
##' # transfer all KEGG species ID to NCBI taxonomy ID
##' wKEGGSpe <- getKEGGPhylo(whole = TRUE)
##' wNCBISpe <- transPhyloKEGG2NCBI(wKEGGSpe[, 2])
##' }
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl getURL
##' @importFrom doMC registerDoMC
##' @importFrom foreach foreach %dopar%
##' @export
##'
transPhyloKEGG2NCBI <- function(KEGGID, n = 4){

  registerDoMC(n)

  getSingleTax <- function (KEGGspeID) {
    # USE: get KEGGSpeID webpage
    # INPUT: 'KEGGID' is the KEGG species ID.
    # OUTPUT: The NCBI taxonomy ID.
    KEGGLink <- paste('http://www.genome.jp/kegg-bin/show_organism?org=', KEGGspeID, sep = '')
    KEGGWeb <- getURL(KEGGLink)

    # get Taxonomy ID. The taxonomy ID is in the web-link like 'http://www.ncbi.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=593907'
    taxIDLink <- gregexpr('wwwtax\\.cgi\\?mode=Info&id=\\d+', KEGGWeb)
    taxIDLink <- getcontent(KEGGWeb, taxIDLink[[1]])
    taxID <- gregexpr('\\d+', taxIDLink)
    taxID <- getcontent(taxIDLink, taxID[[1]])

    return(taxID)
  }

  NCBITax <- foreach(i = 1:length(KEGGID), .combine = c) %dopar% {
    print(paste0('It is running ', i, ' in a total number of ', length(KEGGID), '.'))
    taxID <- getSingleTax(KEGGID[i])
    names(taxID) <- KEGGID[i]
    return(taxID)
  }

  return(NCBITax)

}


##' KEGG Database Additional API - Get nucleotide acid and amino acid sequences according to the T numbers
##'
##' Get protein and gene sequences from KEGG T number in fasta format. As there is no direct API for retrieving the sequence from T number, for example "T10017:100009". The fasta sequence is extract from a webpage like "http://www.genome.jp/dbget-bin/www_bget?-f+-n+a+t10017:100009". The function singleTIDSeq() get a sequence one time, and the function getKEGGTIDGeneSeq() provides a parallel way to download sequences.
##' @title Get protein and gene sequences from T numbers
##' @param TID The T number ID for the protein or gene.
##' @param seqType  Choose nucleotide acid ('ntseq') or amino acid ('aaseq') seqences, and the default is amino acid sequences.
##' @return A BStringSet
##' @importFrom RCurl getURL
##' @importFrom Biostrings BStringSet
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @keywords internal
##'
##' 
singleTIDSeq <- function(TID, seqType = 'aaseq') {
  
  if (seqType == 'aaseq') {
    KEGGLink <- paste0('http://www.genome.jp/dbget-bin/www_bget?-f+-n+', 'a+', TID)
  }
  else if (seqType == 'ntseq') {
    KEGGLink <- paste0('http://www.genome.jp/dbget-bin/www_bget?-f+-n+', 'n+', TID)
  }
  KEGGWeb <- getURL(KEGGLink)

  splitPage <- unlist(strsplit(KEGGWeb, split = '\n', fixed = TRUE))

  ## get sequence name
  ## `nameInd` is also the start number (logic)
  nameInd <- grepl(TID, splitPage, fixed = TRUE)
  seqName <- splitPage[nameInd]
  seqNameStart <- gregexpr(TID, seqName)
  seqNameStart[[1]]
  seqName <- substring(seqName, seqNameStart)

  ## get sequence
  seqStart <- which(nameInd) + 1
  seqEnd <- which(grepl('</pre></div>', splitPage, fixed = TRUE)) - 1
  seq <- paste(splitPage[seqStart:seqEnd], collapse = '')
  seqBS <- BStringSet(seq)
  names(seqBS) <- seqName

  return(seqBS)
}



##' KEGG Database Additional API - Get mutiple nucleotide acid and amino acid sequences according to the T numbers
##'
##' Get mutiple protein and gene sequences from KEGG T number in fasta format. As there is no direct API for retrieving the sequence from T number, for example "T10017:100009". The fasta sequence is extract from a webpage like "http://www.genome.jp/dbget-bin/www_bget?-f+-n+a+t10017:100009".
##' @param TIDs A vector of T number IDs. for the protein or gene.
##' @inheritParams getKEGGGeneSeq
##' @return A BStringSet
##' @examples
##' tNumMultiSeqs <- getKEGGTIDGeneSeq(c('T10017:100009', 'T10017:100036', 'T10017:100044'), n = 2)
##' @importFrom foreach foreach %dopar%
##' @importFrom doMC registerDoMC
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.genome.jp/dbget-bin/www_bget?-f+-n+a+t10017:100009}
##' @references \url{http://www.genome.jp/dbget-bin/www_bget?-f+-n+n+t10017:100009}
##' @seealso getKEGGGeneSeq
##' @export
##'
##' 
getKEGGTIDGeneSeq <- function(TIDs, seqType = 'aaseq', n = 4) {

  # register mutiple cores
  registerDoMC(n)
  
  seqMulRes <- foreach(i = 1:length(TIDs), .combine = append) %dopar% {
    seqRes <- singleTIDSeq(TIDs[i], seqType = seqType)
    return(seqRes)
  }

  return(seqMulRes)
  
}


##' KEGG Database Additional API - Get the gene motif table
##'
##' Get the gene motif tables and additional information.
##' 
##' @title Get motif information from KEGG
##' @param geneID A single KEGG gene ID
##' @param hasAddInfo A logic element whether to show the additional information. The default value is "FALSE".
##' @return A list if "hasAddInfo" is TRUE, else a matrix. A special case retuns "NULL" instead of a matrix (other information is not effected), if there is no motif information.
##' @examples
##' getKEGGGeneMotif('brp:103873230')
##' 
##' # no motif information
##' getKEGGGeneMotif('hsa:4558', hasAddInfo = TRUE)
##' @importFrom XML readHTMLTable
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @export
##' 
getKEGGGeneMotif <- function(geneID, hasAddInfo = FALSE) {
  
  ## motif webpage URL
  url <- paste0('http://www.kegg.jp/ssdb-bin/ssdb_motif?kid=', geneID)

  ## read merged TableList
  readMotifList <- readHTMLTable(url)

  ## gene information
  Organism <- colnames(readMotifList[[1]])[2]

  ## GeneName
  GeneName <- geneID

  ## Definition
  Definition <- as.character(readMotifList[[1]][2, 2])

  ## motif Table
  motifTable <- readMotifList[[2]]

  if (!is.null(motifTable)) {
    motifTable <- apply(motifTable, 1:2, as.character)
  } else {}

  if (hasAddInfo) {
    motifList <- list(Organism = Organism,
                      GeneName = GeneName,
                      Definition = Definition,
                      motifTable = motifTable)
  } else {
    motifList <- motifTable
  }

  return(motifList)
}


##' KEGG Database Additional API - Get genes processing certain motif
##'
##' Get all the gene having certain motif.
##' 
##' @title Get motif list from KEGG.
##' @param motifName A single KEGG motif ID
##' @rdname KEGGMotifList
##' @return A matrix of KEGG genes and description
##' @examples
##' \dontrun{
##' getKEGGMotifList('pf:DUF3675')}
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl getURL
##' @export
##' 
getKEGGMotifList <- function(motifName) {

  ## motif list url
  url <- paste0('www.genome.jp/dbget-bin/get_linkdb?-t+genes+', motifName)

  ## process webpage
  webPage <- getURL(url)

  ## get the webpage contains gene information
  getGeneReg <- gregexpr('<a href=\"/dbget-bin/www_bget?.*$', webPage)
  webPage <- getcontent(webPage, getGeneReg[[1]])

  ## split webPage
  motifList <- unlist(strsplit(webPage, split = '\n', fixed = TRUE))
  motifList <- motifList[1:(length(motifList) - 3)]

  ## get gene names and description
  motifList <- lapply(motifList, function(x) {
    geneReg <- gregexpr('>.*</a>', x)
    geneVal <- getcontent(x, geneReg[[1]])
    geneValNchar <- nchar(geneVal)
    geneVal <- substr(geneVal, start = 2, stop = geneValNchar - 4)

    desReg <- gregexpr('</a>.*$', x)
    desVal <- getcontent(x, desReg[[1]])
    desValNchar <- nchar(desVal)
    desVal <- substring(desVal, 5)
    ## remove space
    desBlankReg <- gregexpr('[^ ].*[^ ]', desVal)
    desVal <- getcontent(desVal, desBlankReg[[1]])

    return(c(geneVal, desVal))
  })

  motifMat <- do.call(rbind, motifList)
  colnames(motifMat) <- c('GeneName', 'Description')
  
  return(motifMat)
}



##' @inheritParams getKEGGGeneMotif
##' @rdname KEGGMotifList
##' @return A vector of protein names including UniProt and SWISS-PROT. Not all these "protID" have KEGG IDs. KEGG uses UniProt and SWISS-PROT for gene UniProt annotation.
##' @examples
##' pfUniprot <- getKEGGMotifList2('pf:DUF3675')
##' \dontrun{
##' # convert uniprot ID to KEGG IDs
##' pfUniprot <- paste0('uniprot:', pfUniprot)
##' pfGeneID <- convKEGG('genes', pfUniprot, convertType = 'identity')
##' }
##' @seealso convKEGG
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @importFrom RCurl getURL
##' @export
##' 
getKEGGMotifList2 <- function(motifName) {

  ## motif list url
  url <- paste0('http://www.genome.jp/dbget-bin/get_linkdb?-t+9+', motifName)

  ## process webpage
  webPage <- getURL(url)

  ## get the webpage contains uniprot information
  getUnipReg <- gregexpr('<a href="/dbget-bin/www_bget?.*?</a>', webPage)
  webPage <- getcontent(webPage, getUnipReg[[1]])

  ## get each uniprot
  unipVec <- sapply(webPage, function(x) {
    getEachReg <- gregexpr('>.*</a>', x)
    getEach <- getcontent(x, getEachReg[[1]])
    getEachNchar <- nchar(getEach)

    getEach <- substr(getEach, start = 2, stop = getEachNchar - 4)
    return(getEach)
  })

  names(unipVec) <- NULL

  return(unipVec)
}


