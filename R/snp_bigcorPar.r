#' Calculate pairwise distance among columns
#'
#' @param x The matrix for processing. Notably, the pairwise distance is calculated among columns.
#' @param size Number of columns in each block. Default value is 1000.
#' @param verbose
#' @param ncore Number of CPU to use
#' @param ... Other parameters for calculate distance, pass to cor()
#'
#' @return A distance matrix among columns
#' @export
#'
#' @examples
snp_bigcorPar <-
  function(x, size = 1000, verbose = TRUE, ncore = "all", ...) {
  library(ff, quietly = TRUE)
  require(doMC)
  if (ncore == "all") {
    ncore <- multicore:::detectCores()
    registerDoMC(cores = ncore)
  } else {
    registerDoMC(cores = ncore)
  }

  NCOL <- ncol(x)

  REST <- NCOL %% size
  LARGE <- NCOL - REST
  NBLOCKS <- NCOL %/% size

  corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))

  GROUP <- rep(1:NBLOCKS, each = size)
  if (REST > 0) {
    GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
  }
  SPLIT <- split(1:NCOL, GROUP)


  #
  # ## test if ncol(x) %% nblocks gives remainder 0
  # if (NCOL %% nblocks != 0){stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")}
  #
  # ## preallocate square matrix of dimension
  # ## ncol(x) in 'ff' single format
  # corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
  #
  # ## split column numbers into 'nblocks' groups
  # SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
  #
  ## create all unique combinations of blocks
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)

  ## iterate through each block combination, calculate correlation matrix
  ## between blocks and store them in the preallocated matrix on both
  ## symmetric sides of the diagonal
  results <- foreach(i = 1:nrow(COMBS)) %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
    flush.console()
    COR <- cor(x[, G1], x[, G2], ...)
    corMAT[G1, G2] <- COR
    corMAT[G2, G1] <- t(COR)
    COR <- NULL
  }

  gc()
  return(corMAT)
}
