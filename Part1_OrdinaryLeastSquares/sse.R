#Sum of squared errors function

sse <- function(par, x, y) {
  B1 = par[1]
  B0 = par[2]
  # fit our data with these parameters
  fitted = B1*x + B0
  # Find our SSE
  sse = sum((fitted - y)^2)
  return(sse)
}
