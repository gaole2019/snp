

#' Title
#'
#' @param x
#' @param size
#' @param verbose
#' @param ncore
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
snp_bigcorPar2 <-
  function(x, size = 1000, verbose = TRUE, ncore = "all", ...) {

  library(ff)
  library(foreach)

  cl <- setup_parallel(ncore)

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

  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)

  # ✅ 并行只计算，不写 corMAT
  results <- foreach(i = 1:nrow(COMBS), .packages = "stats") %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]

    if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")

    COR <- cor(x[, G1], x[, G2], ...)

    list(G1 = G1, G2 = G2, COR = COR)
  }

  # ✅ 主进程写入 ff
  for (res in results) {
    corMAT[res$G1, res$G2] <- res$COR
    corMAT[res$G2, res$G1] <- t(res$COR)
  }

  if (!is.null(cl)) {
    parallel::stopCluster(cl)
  }

  gc()
  return(corMAT)
}

setup_parallel <- function(ncore) {
  if (ncore == "all") {
    ncore <- parallel::detectCores()
  }

  if (.Platform$OS.type == "windows") {
    library(doParallel)
    cl <- parallel::makeCluster(ncore)
    doParallel::registerDoParallel(cl)
    return(cl)
  } else {
    library(doMC)
    doMC::registerDoMC(cores = ncore)
    return(NULL)
  }
}


