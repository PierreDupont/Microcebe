#' WildMap color palette
#'
#' R utility function to createa "WildMap" color palette. Currently supports up to 11 color levels.
#'
#' The \code{wildMap} function is used in advance of model building.
#'
#' @param n an integer denoting the number of colors levels to be returned. 
#'
#' @return It returns a vector of R colors in  “hexadecimal triplets” format.
#'
#' @author Pierre Dupont
#'
#' @examples
#' \donttest{
#' for( i in 3:12){
#' plot(1:i, 1:i, pch = 19, cex = 15, col = wildMap(n = i))
#' }#i
#' }
#'
#' @export
wildMap <- function(n){
  if (n < 3) {
    warning("minimal value for n is 3, returning WildMap palette with 3 different levels")
    return(WildMap(3))
  }
  
  if (n > 11) {
    warning(paste(" maximum alue for n is 11, returning WildMap palette with 11 different levels"))
    return(WildMap(11))
  }
  
  switch(n - 2,
         rgb(c(166, 220, 1  ),
             c(97 , 190, 102),
             c(26 , 130, 94), maxColorValue = 255), # n = 3
         rgb(c(166, 223, 128, 1  ),
             c(97 , 194, 205, 102),
             c(26 , 125, 193, 94), maxColorValue = 255), # n = 4
         rgb(c(140, 216, 223, 128, 1  ),
             c(81 , 179, 212, 205, 102), 
             c(10 , 101, 185, 193, 94), maxColorValue = 255), # n = 5
         rgb(c(140, 216, 226, 189, 90 , 1  ),
             c(81 , 179, 212, 214, 180, 102),
             c(10 , 101, 175, 219, 172, 94 ), maxColorValue = 255), # n = 6
         rgb(c(140, 216, 226, 213, 189, 90 , 1),
             c(81 , 179, 212, 221, 224, 180, 102),
             c(10 , 101, 175, 202, 219, 172, 94), maxColorValue = 255), # n = 7
         rgb(c(140, 191, 223, 230, 189, 128, 53 , 1),
             c(81 , 129, 194, 207, 214, 205, 151, 102),
             c(10 , 45 , 125, 170, 209, 193, 143, 94), maxColorValue = 255), # n = 8
         rgb(c(84 , 140, 191, 223, 230, 169, 128, 53 , 1),
             c(48 , 81 , 129, 194, 207, 194, 205, 151, 102), 
             c(5  , 10 , 45 , 125, 170, 190, 193, 143, 94), maxColorValue = 255), # n = 9
         rgb(c(84 , 140, 191, 223, 220, 169, 128, 53 , 1  , 0),
             c(48 , 81 , 129, 194, 207, 194, 205, 151, 102, 60),
             c(5  , 10 , 45 , 125, 170, 190, 193, 143, 94 , 48), maxColorValue = 255), # n = 10
         rgb(c(84 , 140, 191, 223, 220, 155, 127, 98, 47 , 1  , 0),
             c(48 , 81 , 129, 194, 207, 180, 175, 175, 141, 98, 60),
             c(5  , 10 , 45 , 125, 170, 175, 172, 163, 133, 90 , 48), maxColorValue = 255)) # n = 11
}
