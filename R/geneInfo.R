##' KEGG Database API - Get gene information
##'
##' getKEGGGeneInfo(): Get gene information. This function supports multiple querys.
##'
##' ExtractProLocus(): Extract prokaryotic gene location.
##' 
##' @title Get gene informatino
##' @inheritParams getKEGGGeneSeq
##' @param enforceURL whether to enfoce get get url until the results returned. The default value is "FALSE".
##' @return A string vectors.
##' @examples
##' genes <- c('eco:b4600', 'ece:Z5100', 'eco:b3160', 'dra:DR_0001', 'dra:DR_A0001', 'dra:DR_B0001')
##' genesInfo <- getKEGGGeneInfo(genes, n = 2)
##' @importFrom RCurl getURL
##' @importFrom doParallel registerDoParallel stopImplicitCluster
##' @importFrom foreach foreach %dopar%
##' @importFrom ParaMisc CutSeqEqu
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname geneInfo
##' @export
##'
##' 
getKEGGGeneInfo <- function(KEGGGeneIDs, n = 1, enforceURL = FALSE) {

  ## register multiple core
  registerDoParallel(cores = n)

  ## base url
  urlBase <- 'http://rest.kegg.jp/get/'

  ## cut the input 'KEGGGeneIDs' into 10
  cutMat <- CutSeqEqu(length(KEGGGeneIDs), 10)

  ## deal with ten sequences each time
  print(paste0('Input ', length(KEGGGeneIDs), ' gene IDs.'))
  geneInfo <- foreach(i = cutMat[1, ], j = cutMat[2, ], .combine = c) %dopar% {
    print(paste0('It is getting ', j , ' in a total of ', length(KEGGGeneIDs), '.'))
    if (i == j) {
      ## only one input KEGG ID 
      infoUrl <- paste0(urlBase, KEGGGeneIDs[i])
    } else {
      mergeID <- mergeID <- paste(KEGGGeneIDs[i:j], collapse = "+")
      infoUrl <- paste0(urlBase, mergeID)
    }

    if (enforceURL) {
      webStr <- EnforceGetURL(infoUrl)
    } else {
      webStr <- getURL(infoUrl) 
    }

    geneAnno <- doTenInfo(webStr)
    
    return(geneAnno)
  }

  ## stop multiple core
  stopImplicitCluster()
  
  return(geneInfo)

}


##' @param tenWebInfo return results from function getURL().
##' @return a list each contains the annotation items (for the case without annotation returns NULL).
##' @importFrom stringr str_trim str_extract
##' @author Yulong Niu \email{niuylscu@@gmail.com}
##' @rdname geneInfo
##' @keywords internal
##'
##' 
doTenInfo <- function(tenWebInfo) {
  ## USE: a temporary function to deal with the return web "10 genes", and it also for less ten or without gene information.
  ## oneAnno <- c('ENTRY       DR_0001           CDS       T00025\nDEFINITION  (RefSeq) DNA polymerase III subunit beta\nORTHOLOGY   K02338  DNA polymerase III subunit beta [EC:2.7.7.7]\nORGANISM    dra  Deinococcus radiodurans\nPATHWAY     dra00230  Purine metabolism\n            dra00240  Pyrimidine metabolism\n            dra01100  Metabolic pathways\n            dra03030  DNA replication\n            dra03430  Mismatch repair\n            dra03440  Homologous recombination\nBRITE       KEGG Orthology (KO) [BR:dra00001]\n             Metabolism\n              Nucleotide metabolism\n               00230 Purine metabolism\n                DR_0001\n               00240 Pyrimidine metabolism\n                DR_0001\n             Genetic Information Processing\n              Replication and repair\n               03030 DNA replication\n                DR_0001\n               03430 Mismatch repair\n                DR_0001\n               03440 Homologous recombination\n                DR_0001\n            Enzymes [BR:dra01000]\n             2. Transferases\n              2.7  Transferring phosphorus-containing groups\n               2.7.7  Nucleotidyltransferases\n                2.7.7.7  DNA-directed DNA polymerase\n                 DR_0001\n            DNA replication proteins [BR:dra03032]\n             Prokaryotic Type\n              DNA Replication Elongation Factors\n               Elongation factors (bacterial)\n                DNA polymerase III holoenzyme\n                 DR_0001\n            DNA repair and recombination proteins [BR:dra03400]\n             Prokaryotic Type\n              SSBR (single strand breaks repair)\n               MMR (mismatch exicision repair)\n                DNA polymerase III holoenzyme\n                 DR_0001\nPOSITION    1:complement(1..1182)\nMOTIF       Pfam: DNA_pol3_beta DNA_pol3_beta_2 DNA_pol3_beta_3\nDBLINKS     NCBI-ProteinID: NP_293727\n            NCBI-GI: 15805043\n            NCBI-GeneID: 1799546\n            UniProt: Q9RYE8\nAASEQ       393\n            MMKANVTKKTLNEGLGLLERVIPSRSSNPLLTALKVETSEGGLTLSGTNLEIDLSCFVPA\n            EVQQPENFVVPAHLFAQIVRNLGGELVELELSGQELSVRSGGSDFKLQTGDIEAYPPLSF\n            PAQADVSLDGGELSRAFSSVRYAASNEAFQAVFRGIKLEHHGESARVVASDGYRVAIRDF\n            PASGDGKNLIIPARSVDELIRVLKDGEARFTYGDGMLTVTTDRVKMNLKLLDGDFPDYER\n            VIPKDIKLQVTLPATALKEAVNRVAVLADKNANNRVEFLVSEGTLRLAAEGDYGRAQDTL\n            SVTQGGTEQAMSLAFNARHVLDALGPIDGDAELLFSGSTSPAIFRARRWGRRVYGGHGHA\n            ARLRGLLRPLRGMSALAHHPESSPPLEPRPEFA\nNTSEQ       1182\n            gtgatgaaagccaatgtcaccaaaaagaccctgaacgagggcctgggcctgctcgaacgt\n            gtgattccgagccgttcgagcaatccgctgctgacggcgctgaaggtcgaaacgtcggaa\n            ggtggcctgacgctgagcggcaccaacctggaaatcgacctgtcgtgcttcgtgcctgcc\n            gaggtgcagcagcccgaaaacttcgtggtgccggcgcacctgttcgcgcaaatcgttcgc\n            aacctcggcggtgagctcgtcgaactcgaactgagcggccaggaactctcggtgcgctcg\n            ggcggctcagatttcaagctccagaccggtgacatcgaagcgtacccgccactctctttc\n            cccgcacaggccgatgtgagcctggacggcggcgaactgtcccgcgccttttccagcgtg\n            cgctacgcggcaagcaacgaggcgtttcaggcggtgtttcgcggcattaagcttgagcac\n            cacggcgagagcgcccgcgtggtggcgtccgacggttaccgggtggctatccgcgacttt\n            ccggcgagcggcgacggcaaaaacctgattattcccgcccgcagcgtggacgaactgatt\n            cgcgtgctcaaggacggcgaggcgcggttcacctacggcgacggcatgctcaccgtgacc\n            accgaccgcgtgaagatgaacctcaagctgctcgacggtgattttcccgactacgagcgg\n            gtcattcccaaggacatcaaacttcaggtgacactgcccgccaccgccctcaaggaagcg\n            gtcaaccgtgtggccgtgctggccgacaaaaacgccaacaaccgcgtcgagtttctggtg\n            tccgaaggcactctgcgcctcgctgcggagggcgactatggccgcgctcaggacacgctc\n            agcgtcacccagggcggcaccgagcaggcgatgagcctcgccttcaacgctcgccatgtg\n            ctcgatgcgctgggcccgattgacggagacgccgagctgctgttctccgggtccaccagc\n            cccgccattttccgcgcccgtaggtgggggaggcgggtatatggcggtcatggtcacgct\n            gcgcgtttaaggggccttctgaggccgttacggggcatgtctgccctggcccatcacccg\n            gaaagttcaccgccgcttgaaccgaggccagagttcgcgtga')
  ## doTenInfo(oneAnno)

  if (nchar(tenWebInfo) != 0){
    splitTen <- unlist(strsplit(tenWebInfo, split = '\n///\n', fixed = TRUE))
    geneAnno <- lapply(splitTen, function(x){
      ## The basic idea the item name is all capticalized (or with "_").
      eachInfo <- unlist(strsplit(x, split = '\n', fixed = TRUE))
      
      ## detect item names, and trim blanks after names were identified
      namePattern <- '^[A-Z_]+? +?'
      firstIdx <- which(grepl(namePattern, eachInfo))
      eleNum <- length(firstIdx)
      eachInfo <- str_trim(eachInfo)
      firstEle <- eachInfo[firstIdx]
      
      ## first word is item name, and the rest are contents
      ## initial names and contents in the name
      namesEle <- character(eleNum)
      namesCont <- character(eleNum)
      namesReg <- gregexpr(namePattern, firstEle)
      for (i in 1:eleNum) {
        namesEle[i] <- substr(firstEle[i],
                              start = namesReg[[i]],
                              stop = attr(namesReg[[i]], 'match.length'))
        namesCont[i] <- substring(firstEle[i],
                                  first = attr(namesReg[[i]], 'match.length') + 1)}
      namesEle <- str_trim(namesEle)
      namesCont <- str_trim(namesCont)
      
      ## firstIdx is the startP
      endP <- c(firstIdx[-1] - 1, length(eachInfo))
      eachList <- vector('list', eleNum)
      for (i in 1:eleNum) {
        eachList[[i]] <- eachInfo[firstIdx[i] : endP[i]]
        eachList[[i]][1] <- namesCont[i]
      }
      names(eachList) <- namesEle

      return(eachList)
    })
  } else {
    geneAnno <- NULL
  }

  if (!is.null(geneAnno)) {
    ## names of each anno
    geneAnnoLocus <- str_extract(sapply(geneAnno, function(x) x$ENTRY), '^[^ ]+')
    geneAnnoSpe <- str_extract(sapply(geneAnno, function(x) x$ORGANISM), '^[^ ]+')
    names(geneAnno) <- paste(geneAnnoSpe, geneAnnoLocus, sep = ':')
  } else {}

  return(geneAnno)
}





##' @param locusStr named location strings whose names are gene names.
##' @return a 5-column matrix, 1st gene names, 2ed is genome,  3rd is start, 4th is end, and the 5th is the strand.
##' @examples
##' locusVec <- c('complement(join(1631002..1632285,1652755..1652838))', 'complement(4658240..4658986)',
##' '3303448..3304455', '1:complement(1..1182)', '2:653..1435', 'MP1:520..1296')
##' names(locusVec) <- c('eco:b4600', 'ece:Z5100',
##' 'eco:b3160', 'dra:DR_0001', 'dra:DR_A0001', 'dra:DR_B0001')
##' ExtractProLocus(locusVec)
##' @importFrom stringr str_extract_all
##' @rdname geneInfo
##' @export
##'
##' 
ExtractProLocus <- function(locusStr) {
  ## The basic idea is:
  ## step 1: determine the genome or plasmid by ":". If the returned length of element is 1, then only one genome; else if only numeric element return, for example, "1" and "2", it means multiple genomes; the rest should be plasmid.
  ## step 2: extract locations from a pattern "11111..222222", one gene may have multiple locations.
  ## step 3: whether from the complement string.

  locusStr <- strsplit(locusStr, split = ':', fixed = TRUE)

  ## step 1
  locGenome <- rep('genome1', length(locusStr))
  locLen <- sapply(locusStr, length)
  mulLogic <- locLen > 1
  locGenome[mulLogic] <- sapply(locusStr[which(mulLogic)], '[[', 1)
  mulGenomeLogic <- grepl('^\\d+$', locGenome)
  locGenome[mulGenomeLogic] <- paste0('genome', locGenome[mulGenomeLogic])

  ## step 2
  locFT <- str_extract_all(locusStr, '\\d+..\\d+')
  repNum <- sapply(locFT, length)
  locFT <- str_extract_all(unlist(locFT), '\\d+', simplify = TRUE)

  ## step3
  locStrand <- rep('+', length(locusStr))
  strandLogic <- grepl('complement', locusStr)
  locStrand[strandLogic] <- '-'

  ## combine
  locMat <- cbind(rep(names(locusStr), times = repNum),
                  rep(locGenome, times = repNum),
                  locFT,
                  rep(locStrand, times = repNum),
                  deparse.level = 0)
  colnames(locMat) <- c('GeneName', 'Genome', 'Start', 'End', 'Strand')
  
  return(locMat)
}


