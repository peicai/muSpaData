#' Wei22_full
#' 
#' This function loads the Wei22_full dataset from ExperimentHub.
#' 
#' @param metadata A \code{logical} value indicating whether to return only
#' the ExperimentHub metadata or the entire dataset. Defaults to \code{FALSE}.
#' @return A spatial experiment object associated with Wei22_full.
#' @export
#' @examples
#' Wei22_full() 
Wei22_full <- function(metadata = FALSE) {
    # Calls internal function to load dataset
    .load_dataset("EH9612", metadata)
}

#' Wei22_example
#' 
#' This function loads the Wei22_example dataset from ExperimentHub.
#' 
#' @param metadata A \code{logical} value indicating whether to return only
#' the ExperimentHub metadata or the entire dataset. Defaults to \code{FALSE}.
#' @return A spatial experiment object associated with Wei22_example.
#' @export
#' @examples
#' Wei22_example() 
Wei22_example <- function(metadata = FALSE) {
    # Calls internal function to load dataset
    .load_dataset("EH9613", metadata)
}

#' Load a dataset from ExperimentHub
#' 
#' This function loads a specified dataset from ExperimentHub.
#' 
#' @param ehid A character string representing the ExperimentHub dataset ID.
#' @param metadata A \code{logical} value indicating whether to return only
#' the ExperimentHub metadata, which describes the overall dataset,
#' or to load the entire dataset. Defaults to \code{FALSE}.
#' @return A spatial experiment object associated with the specified dataset.
#' @internal
#' @examples
#' load_dataset("EH9612") # for Wei22_full
#' load_dataset("EH9613") # for Wei22_example
.load_dataset <- function(ehid, metadata = FALSE) {
    eh <- ExperimentHub()
    if (metadata) {
        return(eh[ehid])
    } else {
        return(eh[[ehid]])
    }
}
