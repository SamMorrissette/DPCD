#' Add two numbers together
#'
#' This function takes two numeric inputs and returns their sum.
#'
#' @param x A numeric value
#' @param y A numeric value
#' @return The sum of x and y
#' @export
#' @examples
#' add_numbers(2, 3)
#' add_numbers(10, -5)
add_numbers <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("Both inputs must be numeric")
  }
  x + y
}
