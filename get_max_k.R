get_max_k <- function(delta, p0, pk) {
  ceiling_val <- ceiling((pk - p0) / delta)
  k_max <- 2 * ceiling_val
  return(k_max)
}


k_max <- get_max_k(0.2, 0.0753, 0.6347)

print(k_max)