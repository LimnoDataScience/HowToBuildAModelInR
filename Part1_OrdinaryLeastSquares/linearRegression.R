simple.linear.coef <- function(x, y) {
  # Find length of x to get sample size. Assuming x and y have the same sample size.
  n <- length(x)
  # Calculate the error statistics Sxx, Syy, and Sxy
  sxx <- sum(x^2) - sum(x)^2 / n
  syy <- sum(y^2) - sum(y)^2 / n
  sxy <- sum(x * y) - (sum(x) * sum(y)) / n
  # Coefficients beta0 and beta1
  b1 <- sxy / sxx
  b0 <- mean(y) - b1 * mean(x)
  # Sum of standard error 
  sse <- syy - sxy^2 / sxx

  # Coefficient of determination R-squared
  r2 <- (syy - sse) /syy
  # R-squared adjusted
  r2.adj <- r2 - (1 - r2) * ((2 - 1) / (length(y) - 2))
  
  rsquare <- paste('Multiple R-squared: ', round(r2, 4), ', Adjusted R-squared: ', round(r2.adj, 4))
  
  coeffs <- data.frame(cbind(c(b0, b1)))
  colnames(coeffs) <- c('Estimate')
  rownames(coeffs) <- c('Intercept', 'x1')
  
  # Fit the line to the data with beta0 and beta1 found above
  fitted <- x * b1 + b0
  
  # The F-Statistic
  msr <- sum((fitted - mean(y))^2) / 1
  mse2 <- sum((y - fitted)^2) / (length(y) - 2)
  f <- msr / mse2
  # p-value
  p <- pf(f, 1, length(y) - 2, lower.tail = FALSE)
  p.val <- paste('p-value: ', format(p, digits = 3, scientific = TRUE))
  
  ### Return list of coefficients, r2 and p value
  regres <- list('Coefficients'=coeffs, rsquare, p.val)
  
  return(regres)
}

# Run function
simple.linear.coef(cars$speed, cars$dist)

# Compare to lm
cars.lm <- lm(dist ~ speed, data = cars)
summary(cars.lm)
plot(dist ~ speed, data = cars)
