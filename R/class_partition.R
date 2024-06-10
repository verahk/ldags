

#' Class: partition
#' @name partition
#' @description
#' An object of class `partition` is a list that defines for a partition of the 
#' outcome space of a set of categorical variables. 
NULL

#' @rdname partition
#' @param x (list)
#'  a list with 
#' @param nlev integer vector:
#'  the cardinality of each variable
#' @return an object of class `partition`
#' @examples
#' # generate a new partition
#' nlev <- c(2, 2)
#' P <- new_partition(list(c(1, 3), 2, 0), nlev)
#' P
new_partition <- function(x, nlev) {
  structure(x, 
            nlev = nlev, 
            class = "partition")
}

validate_partition <- function(x) {
  n <- prod(attr(x, "nlev"))
  outcomes <- sort(unlist(x))
  if (anyDups(outcomes) > 0) {
    stop("The partition contains duplicate outcomes")
  } else if (length(outcomes) != n) {
    stop("The partition contains to many or to few outcomes")
  } else if (! all(outcomes == seq_len(n)-1)) {
    stop("The outcomes should be encoded as 0, ..., prod(nlev)-1")
  }
}

#' @rdname partition
#' @param P an object of class `partition`
#' @return 
#' - `get_parts:` a vector with that assigns each outcome to a part
#' @export
#'
#' @examples
#' # list which part each outcome belongs to
#' get_parts(P)
get_parts <- function(P, lens = lengths(P)) {
  rep.int(seq_along(P), lens)[order(unlist(P))]
}

