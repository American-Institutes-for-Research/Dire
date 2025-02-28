getOneD_E_sd <- function(xb, rr1, si, nodes) {
  Xb0 <- t(matrix(nodes, nrow=length(nodes), ncol=length(xb)))
  nodes.minus.XBi <- t(Xb0 - xb)
  pr <- rr1 * dnorm(nodes.minus.XBi, mean=0, sd=si)
  E_sum <- rep(0, length(xb))
  for(i in 1:length(nodes)) {
    E_sum <- E_sum + nodes[i] * pr[i,]
  }
  pr_sum <- apply(pr, 2, sum)
  E <- E_sum / pr_sum
  V_sum <- rep(0, length(xb))
  for(i in 1:length(nodes)) {
    V_sum <- V_sum + (nodes[i] - E)^2 * pr[i,]
  }
  V <- V_sum / pr_sum
  return(list(E=E, V=V))
}

