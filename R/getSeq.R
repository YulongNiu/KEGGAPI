##' KEGG Database API - Get the nucleotide acid and amino acid sequences 
##'
##' Get protein and gene sequences in fasta format. This function supports mutiple querys.
##' @title Get protein and gene sequences
##' @param KEGGGeneIDs A vector of KEGG IDs. Seqences from different species could be combined together.
##' @param seqType Choose nucleotide acid (ntseq) or amino acid (aaseq) seqences, and the default is amino acid sequences.
##' @param n The number of CPUs or processors, and the default value is 1.
##' @inheritParams webTable
##' @return A BStringSet.
##' @examples
##' ## two amino acid seqences from different sepecies with 2 threads.
##' twoAASeqs <- getKEGGGeneSeq(c('mja:MJ_0011', 'hsa:10458'), n = 2)
##' \dontrun{
##' ## export fasta format files
##' require('Biostrings')
##' writeXStringSet(twoAASeqs, 'twoAASeqs.fasta')}
##'
##' \dontrun{
##' ## more examples
##' twoNTSeqs <- getKEGGGeneSeq(c('shy:SHJG_7159', 'shy:SHJG_7160'), 'ntseq')
##' mutilAASeqs <- getKEGGGeneSeq(c('eco:b0202', 'eco:b0203', 'eco:b0204',
##' 'eco:b0205', 'eco:b0206', 'eco:b0216', 'eco:b0244',
##' 'eco:b4626', 'eco:b3796', 'eco:b3797', 'eco:b3296',
##' 'eco:b3297'))}
##' 
##' \dontrun{
##' ## get the whole E.coli genome protein seqences
##' ecoProIDs <- getProID('eco')
##' ecoGenomePro <- getKEGGGeneSeq(ecoProIDs[, 1])}
##' @importFrom RCurl getURL
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%

##' @importFrom ParaMisc CutSeqEqu
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @references \url{http://www.kegg.jp/kegg/rest/keggapi.html}
##' @rdname getSeq
##' @seealso getKEGGTIDGeneSeq
##' @export
##'
##' 
getKEGGGeneSeq <- function(KEGGGeneIDs, seqType = 'aaseq', enforceURL = FALSE, n = 1){

  ## register multiple core
  registerDoParallel(cores = n)

  ## base url
  urlBase <- 'http://rest.kegg.jp/get/'

  
  ## cut the input 'KEGGGeneIDs' into 10
  cutMat <- CutSeqEqu(length(KEGGGeneIDs), 10)

  ## deal with ten sequences each time
  print(paste0('Input ', length(KEGGGeneIDs), ' gene IDs.'))
  proSeq <- foreach(i = cutMat[1, ], j = cutMat[2, ], .combine = append) %dopar% {
    print(paste0('It is getting ', j , ' in a total of ', length(KEGGGeneIDs), '.'))
    if (i == j){
      ## only one input KEGG ID
      seqUrl <- paste0(urlBase, KEGGGeneIDs[i], '/', seqType)
    } else {
      mergeID <- paste(KEGGGeneIDs[i:j], collapse = "+")
      seqUrl <- paste0(urlBase, mergeID, '/', seqType)
    }

    if (enforceURL) {
      webStr <- EnforceGetURL(seqUrl)
    } else {
      webStr <- getURL(seqUrl) 
    }

    geneSeq <- doTen(webStr, seqType)
    return(geneSeq)
  }

  ## stop multiple core
  stopImplicitCluster()
  
  return(proSeq)

}


##' @param tenWebSeq return results from function getURL()
##' @inheritParams getKEGGGeneSeq
##' @return  BStingSet or NULL (for the case one gene has no cording protein sequence)
##' @importFrom Biostrings BStringSet
##' @importFrom foreach foreach %do%
##' @rdname getSeq
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @keywords internal
##'
##' 
doTen <- function(tenWebSeq, seqType){
  ## USE: a temporary function to deal with the return web "10 seqence", and it also for less ten or without coding sequence. The basic idea to split sequences is by the marker of '(A)' for amino acid sequences and '(N)' for nucleotide sequences.
  if (nchar(tenWebSeq) != 0){
    splitTen <- unlist(strsplit(tenWebSeq, split = '\n', fixed = TRUE))
    if (seqType == 'aaseq') {
      namePoint <- which(grepl(' \\(A\\)', splitTen))
    }
      else if (seqType == 'ntseq') {
        namePoint <- which(grepl(' \\(N\\)', splitTen))
      }
    ## sequence start point
    startPoint <- namePoint + 1
    endPoint <- c(namePoint[-1] - 1, length(splitTen))
  } else {
    ## gene without coding sequence
    return(NULL)
  }

  tenSeqBS <- foreach(i = 1:length(namePoint), .combine = append) %do% {
    proSeqName <- splitTen[namePoint[i]]
    proSeqName <- substring(proSeqName, 2)
    proSeq <- paste(splitTen[startPoint[i] : endPoint[i]], collapse = '')
    ## to BStringSet
    proSeqBS <- BStringSet(proSeq)
    names(proSeqBS) <- proSeqName
    return(proSeqBS)
  }

  return(tenSeqBS)
}

