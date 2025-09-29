#' Calculate the axon lengths in annotated brain regions
#'
#' @details Annotated brain regions can be any 3D array with
#' an unique index denoting a brain region or domain
#'
#' @param PathOfNeuron Full path of SWC file for each neuron
#' @param Hemisphere Which hemisphere does the neuron reside in. Note that neurons in right hemisphere were flipped to the left side.
#' @param Annotation 3D array with each element assigned a positive value
#' @param Resolution Voxel size of each annotation, e.g., 100, 50, 25, 10 um
#' @param NumOfCores Number of cores for parallel computation
#' @param OptionMetric Default value is "LengthOfAxon". ("LengthOfAxon", "NumberOfEndNodes", "NumberOfBranchNodes")
#' @param StructureIDOfTargets Structure ID of brain regions or domains in the annotation.
#' Column name of returned matrix is the same as StuctureIDOfTargets.
#'
#' @return projection matrix
#' @export
#'
#' @examples
snp_calculateProjections <-
  function(PathOfNeuron,
           Hemisphere,
           StructureIDOfTargets,
           Annotation,
           Resolution,
           NumOfCores,
           OptionMetric = "LengthOfAxon") {

    library(nat)
    library(doMC)
    library(foreach)

    NumNeuron <- length(PathOfNeuron)
    BasicStructureID <- StructureIDOfTargets
    NumOfTargets <- length(StructureIDOfTargets)

    print(paste0("Number of neurons: ", NumNeuron))
    print(paste0("Number of targets: ", NumOfTargets))

    # order of axis is important
    NumRow_Anno <- dim(Annotation)[1] # X
    NumColumn_Anno <- dim(Annotation)[2] # Y
    NumSlice_Anno <- dim(Annotation)[3] # Z

    doMC::registerDoMC(cores = NumOfCores)

    result <-
      foreach(iNeuron = 1:NumNeuron, .combine = rbind, .inorder = TRUE) %dopar% {
        vector_proj <- matrix(data = 0, nrow = 1, ncol = NumOfTargets)

        neuron <- nat::read.neuron(PathOfNeuron[iNeuron])
        tmpSWC <- neuron$d

        # flip neurons in the right hemisphere to the left side
        if (Hemisphere[iNeuron] == "Right")
          tmpSWC$Z <- 11400 - tmpSWC$Z

        # LengthOfAxon ----
        if ("LengthOfAxon" == OptionMetric) {
          anno_nodes <-
            apply(tmpSWC[,c("X","Y","Z")], 1,
                  function(x) {
                    if (x[1] > 1 & x[1] < NumRow_Anno*Resolution &
                        x[2] > 1 & x[2] < NumColumn_Anno*Resolution &
                        x[3] > 1 & x[3] < NumSlice_Anno*Resolution) {
                      Annotation[floor(x[1]/Resolution) + 1,
                                 floor(x[2]/Resolution) + 1,
                                 floor(x[3]/Resolution) + 1]
                    }
                    else {
                      -1
                    }
                  })


          numBasicStructureID <- length(BasicStructureID)
          for (iBasicStructureID in 1:numBasicStructureID) {
            index_node_current <-
              which(anno_nodes == BasicStructureID[iBasicStructureID])
            if (length(index_node_current) == 0) {
              vector_proj[1, iBasicStructureID] <- 0
            } else {
              vector_proj[1, iBasicStructureID] <-
                sum(sqrt(rowSums((tmpSWC[index_node_current,3:5] -
                                    tmpSWC[match(tmpSWC$Parent[index_node_current],
                                                 tmpSWC$PointNo),3:5])^2)), na.rm = T)
            }
          }
        } # end of "LengthOfAxon" metric


        # NumberOfEndNodes ----
        if ("NumberOfEndNodes" == OptionMetric) {
          index_end_points <- neuron$EndPoints

          for (iEnd in 1:length(index_end_points)) {
            coord_node_current <- c(
              neuron$d$X[index_end_points[iEnd]],
              neuron$d$Y[index_end_points[iEnd]],
              neuron$d$Z[index_end_points[iEnd]]
            )

            ## what if the node is out of the annotation.
            tmpRow_current <- floor(coord_node_current[1] / Resolution) + 1
            tmpColumn_current <- floor(coord_node_current[2] / Resolution) + 1
            tmpSlice_current <- floor(coord_node_current[3] / Resolution) + 1

            if (tmpRow_current < 1 || tmpRow_current > NumRow_Anno ||
                tmpColumn_current < 1 || tmpColumn_current > NumColumn_Anno ||
                tmpSlice_current < 1 || tmpSlice_current > NumSlice_Anno) {
              next
            }

            structureID_node_current <-
              Annotation[tmpRow_current, tmpColumn_current, tmpSlice_current]

            vector_proj[1, which(BasicStructureID == structureID_node_current)] <-
              vector_proj[1, which(BasicStructureID == structureID_node_current)] + 1
          }
        } # end of "NumberOfEndNodes" metric


        # NumberOfBranchNodes ----
        if ("NumberOfBranchNodes" == OptionMetric) {
          #neuron <- nat::read.neuron(PathOfNeuron[iNeuron])
          index_branch_points <- neuron$BranchPoints

          for (iEnd in 1:length(index_branch_points)) {
            coord_node_current <- c(
              neuron$d$X[index_branch_points[iEnd]],
              neuron$d$Y[index_branch_points[iEnd]],
              neuron$d$Z[index_branch_points[iEnd]]
            )

            ## what if the node is out of the annotation.
            tmpRow_current <- floor(coord_node_current[1] / Resolution) + 1
            tmpColumn_current <- floor(coord_node_current[2] / Resolution) + 1
            tmpSlice_current <- floor(coord_node_current[3] / Resolution) + 1

            if (tmpRow_current < 1 || tmpRow_current > NumRow_Anno ||
                tmpColumn_current < 1 || tmpColumn_current > NumColumn_Anno ||
                tmpSlice_current < 1 || tmpSlice_current > NumSlice_Anno) {
              next
            }

            structureID_node_current <-
              Annotation[tmpRow_current, tmpColumn_current, tmpSlice_current]

            vector_proj[1, which(BasicStructureID == structureID_node_current)] <-
              vector_proj[1, which(BasicStructureID == structureID_node_current)] + 1
          }
        } # end of "NumberOfBranchNodes" metric

        vector_proj
      }


    colnames(result) <- BasicStructureID
    return(result)
  }
