#' @title Model of Larval Production
#'
#' @description This function models the number of larvae present in an ant
#'  colony as function of time.
#'
#' @param t A vector of the times at which to compute the number of larvae.
#'
#' @param A The relative amplitude of the larval production cycle.
#'
#' @param M The average number of larvae in a colony at any one time.
#'
#' @param P The period of the larval production cycle.
#'
#' @param e The smoothness of the larval production cycle.
#'
#' @return A vector of the same length as \code{t}.
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @export
#'
lt <- function(t, A, M, P, e) {
  M * (1 + ((1 + sqrt(10 ^ (2 * e))) * sin((2 * pi * t) / P) /
              (1 + sqrt(10 ^ (2 * e) * sin((2 * pi * t) / P) ^ 2)))) +
    (1 - A) * (M - M * (1 + ((1 + sqrt(10 ^ (2 * e))) * sin((2 * pi * t) / P)) /
                          (1 + sqrt(10 ^ (2 * e) * sin((2 * pi * t) / P) ^ 2))))
}


#' @title Foraging Cost Function
#'
#' @description This function models the foraging cost associated with a given
#'  number of larvae present in an ant colony.
#'
#' @param l The number of larvae in the colony.
#'
#' @param k The maximum number of larvae that a colony can have at any given
#'  time.
#'
#' @param n A parameter that determines how the cost of foraging scales with the
#'  number of larvae to be fed.
#'
#' @return A vector of the same length as \code{l}.
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @export
#'
fc <- function(l, k, n) {
  (1 / 2) * (k ^ (1 - n)) * (l ^ n) * (1 + n)
}


#' @title Total Foraging Cost
#'
#' @description This function models the total foraging cost for an ant colony
#'  as a function of time. It corresponds to the composite of \code{\link{fc}}
#'  and \code{\link{lt}}.
#'
#' @param t A vector of the times at which to compute the number of larvae.
#'
#' @param A The relative amplitude of the larval production cycle.
#'
#' @param M The average number of larvae in a colony at any one time.
#'
#' @param P The period of the larval production cycle.
#'
#' @param e The smoothness of the larval production cycle.
#'
#' @return A vector of the same length as \code{t}.
#'
#' @param l The number of larvae in the colony.
#'
#' @param k The maximum number of larvae that a colony can have at any given
#'  time.
#'
#' @param n A parameter that determines how the cost of foraging scales with the
#'  number of larvae to be fed.
#'
#' @return A vector of the same length as \code{t}.
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @seealso \code{\link{fc}}, \code{\link{lt}}
#'
#' @export
#'
tc <- function(t, A, M, P, e, k, n) {
  fc(lt(t, A, M, P, e), k, n)
}


#' @title Approximate Local Gradient
#'
#' @description This function computes an approximation of the local 2D gradient
#'  of a matrix using the Sobel operator.
#'
#' @param m A matrix.
#'
#' @return A list of two matrices \code{r} and \code{c} representing the local
#'  gradient along the rows and columns of the matrix respectively.
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @export
#'
approxGradient <- function(m) {
  sobel_r <- matrix(c(-1, -2, -1, 0, 0, 0, 1, 2, 1), 3, 3)
  sobel_c <- t(sobel_r)

  r <- mmand::morph(m, sobel_r, "*")
  c <- mmand::morph(m, sobel_c, "*")

  list(r = r, c = c)
}
