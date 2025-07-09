# implement the formula of rcs components
get_rcs <- function(x, knots) {
  k <- length(knots)
  res <- data.frame(rcs.1 = x)
  for (j in 1:(k-2)) {
    kd <- (knots[k] - knots[1])^(2/3) # (knot_k-knot_1)^(2/3)
    
    vec1 <- x - knots[j] # X-knot_j
    
    val2 <- knots[k-1] - knots[j] # knot_{k-1}-knot_j
    vec2 <- x - knots[k] # X-knot_k
    
    val3 <- knots[k] - knots[j] # knot_k-knot_j
    vec3 <- x - knots[k-1] # X-knot_{k-1}
    val4 <- knots[k] - knots[k-1] # knot_k-knot_{k-1}
    
    vec_j <- pmax(vec1/kd, 0)^3 + (val2*pmax(vec2/kd, 0)^3 - val3*pmax(vec3/kd, 0)^3)/val4
    
    # View(Hmisc::rcspline.eval): line 111-114
    # vec_j <- pmax((x - knots[j])/kd, 0)^3 + 
    #   ((knots[k-1] - knots[j]) * pmax((x - knots[k])/kd, 0)^3 - 
    #      (knots[k] - knots[j]) * (pmax((x - knots[k-1])/kd, 0)^power))/(knots[k] - knots[k-1])
    
    res$new <- vec_j
    names(res)[ncol(res)] <- paste0("rcs.", j+1)
  }
  res <- as.matrix(res)
  return(res)
}

# --------------------------
# calculate the offset value
offset_rcs <- function(y_mean, x_mean, knots, coefs) {
  res <- get_rcs(x_mean, knots)
  sum <- sum(res*coefs[-1])
  offset <- coefs[1] + sum - y_mean
  return(offset)
}