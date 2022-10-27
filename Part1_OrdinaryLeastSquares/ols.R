# Analytical solution to least squares regression

linear.coef <- function(x,y) {
  n = length(x) # sample size
  # calculate SSxy and SSxx
  SSxy = sum(x*y) - (sum(x)*sum(y))/n
  SSxx = sum(x^2) - sum(x)^2/n
  # Find the coefficients 
  B1 = SSxy / SSxx
  B0 = mean(y) - B1*mean(x)
  
  coeffs = data.frame(Intecept = B0, slope = B1)
  return(coeffs)
}