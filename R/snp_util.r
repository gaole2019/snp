#' convert Acronym to StructureID
#'
#' @description
#'
#' This function returns the structureID of specified regions (acronym)
#'
#' @param acronymList
#'
#' @return
#' @export
#'
#' @examples
#' snp::getStructureFromAcronym(c("l_Isocortex", "r_CP"))
snp_getStructureIDFromAcronym <-
  function(acronymList, flag_debug = FALSE) {
    NumAcronym <- length(acronymList)
    vec_structure <- matrix(nrow = NumAcronym, ncol = 1)
    for (iAcronym in 1:NumAcronym) {
      tmpAcronym <- acronymList[iAcronym]
      tmpHemisphere <- substr(tmpAcronym, 1, 2)
      tmpAcronymName <- substr(tmpAcronym, 3, nchar(tmpAcronym))

      if (length(match(tmpAcronymName, snp::shared_allen_anno$Acronym)) == 0) {
        vec_structure[iAcronym, 1] <- NA
      } else {
        vec_structure[iAcronym, 1] <-
          paste0(tmpHemisphere,
                 as.character(
                   snp::shared_allen_anno$CurrentID[
                     match(tmpAcronymName,
                           snp::shared_allen_anno$Acronym)]))
      }

      if (flag_debug) {
        print(vec_structure[iAcronym, 1])
      }
    }
    return(vec_structure)
  }


#' Convert StructureID to Acronym
#'
#' @param structureList
#'
#' @return
#' @export
#'
#' @examples
snp_getAcronymFromStructureID <-
  function(structureList) {
    NumStructure <- length(structureList)
    vec_acronym <- matrix(nrow = NumStructure, ncol = 1)
    for (iStructure in 1:NumStructure) {
      tmpStructure <- structureList[iStructure]
      tmpHemisphere <- substr(tmpStructure, 1, 2)
      tmpStructureID <- substr(tmpStructure, 3, nchar(tmpStructure))

      vec_acronym[iStructure, 1] <-
        paste0(tmpHemisphere,
               as.character(snp::shared_allen_anno$Acronym[match(
                 tmpStructureID,
                 snp::shared_allen_anno$CurrentID)]))
    }
    return(vec_acronym)
  }
