pgrad <- function(x, t, sigj){

  tmp <- (t-x) / sigj
  n <- exp(1.702 * tmp)
  d <- n / (1+n)

  return(d)

}
