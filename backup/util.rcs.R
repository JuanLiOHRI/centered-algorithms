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
# define the function for `rootSolve::uniroot.all`
fun_rcs <- function(x, object, y_mean, varname, knots = NULL) {
  if (is.null(knots)) { # single predictor
    df <- data.frame(var = x)
    names(df) <- varname
    result <- predict(object, newdata = df)-y_mean
    return(result)
  } else {
    coefs <- coef(object)
    coefs <- coefs[c(1,which(stringr::str_detect(names(coefs), varname)))]
    rcs <- c(1, as.vector(get_rcs(x, knots))) # add intercept
    result <- rcs %*% coefs - y_mean
    return(result)
  }
}

# --------------------------
# A blunt-force solution for searching for 2 roots
root.search <- function(model, rng1, rng2, varnames, y_mean, length.out = 101, err = 1e-4) {
  df <- data.frame(
    var1 = seq(rng1[1],rng1[2],length.out=length.out),
    var2 = seq(rng2[1],rng2[2],length.out=length.out))
  names(df) <- varnames
  df$y <- predict(model, newdata = df)
  df$diff <- df$y - y_mean
  start <- df %>% filter(diff >= 0)
  start <- start %>% filter(diff == min(diff))
  end <- df %>% filter(diff < 0)
  end <- end %>% filter(diff == max(diff))

  if (abs(start$diff) <= err) {
    return(start)
  } else if (abs(end$diff) <= err) {
    return(end)
  } else {
    root.search(model, rng1 = sort(c(start[,varnames[1]], end[,varnames[1]])), 
      rng2 = sort(c(start[,varnames[2]], end[,varnames[2]])), varnames, y_mean, length.out, err)
  }
}