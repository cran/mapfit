#' Padding zero to omitted diagonal elements
#' 
#' This function provides a sparse matrix whose all diagonal elements are not
#' omitted even if they have zero values.
#' 
#' @param A A sparse matrix
#' @return A sparse matrix whose all diagonal elements are not omitted even if 
#' they have zero values.
#' 

diag.padding.zero <- function(A) {
  x <- diag(A)
  diag(A)[x == 0] <- NA
  A@x[is.na(A@x)] <- 0
  A
}

