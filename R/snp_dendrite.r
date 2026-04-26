#' Calculate morphological features of dendrites
#'
#' @param neuronlist list of dendrites
#'
#' @returns feature matrix
#' @export
#'
#' @examples
#'
snp_calculate_dendrite_feature <-
  function(neuronlist) {
    library(nat)
    library(geometry)
    library(igraph)

    l <- lapply(neuronlist, calculate)
    d <- do.call(rbind, l)
  }


calculate <-
  function(neuron) {
    data.frame(
      NumberOfNeurites = cal_NumberOfNeurites(neuron),
      NumOfSegments = neuron$NumSegs,
      NumOfBranchpoints = length(neuron$BranchPoints),
      NumOfEndpoints = length(neuron$EndPoints),
      Length = cal_total_cable(neuron),
      MeanSegmentLength = cal_total_cable(neuron) / neuron$NumSegs,
      NodeToRoot_EucDist_max = cal_root_EucDist(neuron, "max"),
      NodeToRoot_EucDist_mean = cal_root_EucDist(neuron, "mean"),
      NodeToRoot_PathDist_max = cal_root_PathDist(neuron, "max"),
      NodeToRoot_PathDist_mean = cal_root_PathDist(neuron, "mean"),
      Volume = cal_volume(neuron),
      Density = cal_total_cable(neuron) / cal_volume(neuron)
    )
  }

cal_NumberOfNeurites <-
  function(neuron) {
    tmpD <- neuron$d
    tmpID_root <- tmpD$PointNo[which(tmpD$Parent == -1)]
    return(length(which(tmpD$Parent == tmpID_root)))
  }

cal_total_cable <-
  function(neuron) {
    diffs <- neuron$d[, c("X","Y","Z")] -
      neuron$d[match(neuron$d$Parent, neuron$d$PointNo), c("X","Y","Z")]
    return(sum(sqrt(rowSums(diffs*diffs)), na.rm = TRUE))
  }

cal_root_EucDist <-
  function(neuron, option) {
    EucDist <- sqrt((neuron$d$X - neuron$d$X[1])^2 +
                    (neuron$d$Y - neuron$d$Y[1])^2 +
                    (neuron$d$Z - neuron$d$Z[1])^2)
    if(option == "max") return(max(EucDist))
    if(option == "mean") return(mean(EucDist))
  }

cal_root_PathDist <-
  function(neuron, option) {
    PathDist <- get_distance_to_soma(neuron)
    if(option == "max") return(max(PathDist))
    if(option == "mean") return(mean(PathDist))
  }

cal_volume <-
  function(neuron) {
    return(geometry::convhulln(neuron$d[,3:5], options="FA")$vol)
  }

# copied from nat
get_distance_to_soma <- function(n) {
  gw <- nat::as.ngraph(n, weights=TRUE)
  dst <- igraph::distances(gw, v = rootpoints(n))
  as.numeric(dst)
}
