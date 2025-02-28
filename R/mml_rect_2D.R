# updates E_i, E_j, vari_vec, varj_vec, cov_vec
# uses Sigma, xbi, xbj, rr1i, rr1j, nodes, n
e_step_rect_grid <- function(values, i, j, e_only=FALSE, addGr=FALSE) {
  # temporary vectors to calculate integrals. One element per person in the sample
  pr_sumi <- pr_sumj <- E_i_sum <- E_j_sum <- var_i_sum <- var_j_sum <- cov_sum <- pr_sum <- rep(0, values$n)
  if(addGr) {
    gr_mat <- matrix(0, nrow=values$n, ncol=2)
  }
  Sigma <- values$Sigma[c(i,j),c(i,j)]
  nodes <- values$nodes
  Q <- length(nodes)
  Sigma_inv <- solve(Sigma)
  det_sigma_12 <- 1/(2*pi) * 1/sqrt(Sigma[1,1]*Sigma[2,2] - Sigma[1,2] * Sigma[2,1])
  for(ii in 1:Q) {
    resid_ii <- nodes[ii] - values$Xb[,i]
    for(jj in 1:Q) { # this is 5.2 generalized to covariance
      resid_jj <- nodes[jj] - values$Xb[,j]
      # drop normalizing terms because they will be the same in the numerator and denominator
      # BUT note that this does not have the det(Sigma) so can't be used to maxamize likelihood
      mvn0 <- det_sigma_12 * exp(dmvn_log_exp(resid_ii, resid_jj, Sigma_inv[1,1], Sigma_inv[2,2], Sigma_inv[1,2])) 
      # this is pr(theta|responses) * phi term
      pr_part <- values$rr1_list[[i]][ii,] * values$rr1_list[[j]][jj,] * mvn0
      # for denominator
      pr_sum <- pr_sum + pr_part
      if(addGr) {
        gr_mat <- gr_mat + pr_part * cbind(resid_ii, resid_jj) %*% Sigma_inv
      }
      E_i_sum <- E_i_sum + nodes[ii] * pr_part
      E_j_sum <- E_j_sum + nodes[jj] * pr_part
    }
  }
  if(addGr) {
    values$gr[,i] <- gr_mat[,1]/pr_sum
    values$gr[,j] <- gr_mat[,2]/pr_sum
  }
  values$E[,i] <- E_i_sum/pr_sum
  values$E[,j] <- E_j_sum/pr_sum
  if(e_only) {
    return(values)
  }
  for(ii in 1:Q) {
    resid_ii <- nodes[ii] - values$Xb[,i]
    resid_ii_post <- nodes[ii] - values$E[,i]
    for(jj in 1:Q) { # this is 5.2 generalized to covariance
      resid_jj <- nodes[jj] - values$Xb[,j]
      resid_jj_post <- nodes[jj] - values$E[,j]
      # drop normalizing terms because they will be the same in the numerator and denominator
      # BUT note that this does not have the det(Sigma) so can't be used to maxamize likelihood
      mvn0 <- det_sigma_12 * exp(dmvn_log_exp(resid_ii, resid_jj, Sigma_inv[1,1], Sigma_inv[2,2], Sigma_inv[1,2])) 
      # this is pr(theta|responses) * phi term
      pr_part <- values$rr1_list[[i]][ii,] * values$rr1_list[[j]][jj,] * mvn0
      # for denominator
      # don't sum this twice pr_sum <- pr_sum + pr_part
      # this is the numerator that calculates the generalization of 5.2 of Thomas 1993, Journal of Computational and Graphical Statistics (2)3 209-322
      cov_sum <- cov_sum + resid_ii_post * resid_jj_post * pr_part
      var_i_sum <- var_i_sum + resid_ii_post^2 * pr_part
      var_j_sum <- var_j_sum + resid_jj_post^2 * pr_part
    }
  }
  mm <- matmap_factory(values$k)
  # calculate the covariance term vectors
  values$cov[,i] <- var_i_sum/pr_sum
  values$cov[,j] <- var_j_sum/pr_sum
  values$cov[,mm(i,j)] <- cov_sum/pr_sum
  return(values)
}

# makes a crosswalk from the upper(/lower) triangular part (excluding the diagonal) to an index.
# this allows the upper/lower triangular part to be stored in an array
# test code:
# k <- 7
# dmat <- matrix(0, nrow=k, ncol=k)
# matmap <- matmap_factory(k)
# for(i in 1:k){
#   for(j in 1:k) {
#     dmat[i,j] <- matmap(i,j)
#   }
# }
# dmat
matmap_factory <- function(k) {
  presum <- -k + cumsum(0:k)
  function(i,j) {
    if(i == j) {
      return(i)
    }
    ii <- min(i,j)
    jj <- max(i,j)
    return(k*(ii-1) - presum[ii] + jj-ii)
  }
}

# uses: n, Sigma, Xbmat, rrlist, nodes,
# updates: E_i, E_j, xbi, xbj, vari_vec, varj_vec, cov_vec, rr1i, rrij, posterior_dist (E, V, VC, sigma_list)
e_step_nD_by2D <- function(values, addGr=FALSE) {
  n <- values$n
  k <- values$k
  if(addGr) {
    if(k > 2) {
      warning("Using addGr=TRUE is not advised when there are more than two dimensions.")
    }
    values$gr <- matrix(0, nrow=n, ncol=k)
  }
  for(i in 1:(k-1)) {
    for(j in (i+1):k) {
      values <- e_step_rect_grid(values, i, j, addGr=addGr)
    }
  }
  return(values)
}

#'@importFrom stats cov2cor
add_sigma_list <- function(values) {
  k <- values$k
  n <- values$n
  pre_sigma <- matrix(0, nrow=n, ncol=k^2)
  mm <- matmap_factory(values$k)
  psi <- 1
  for(i in 1:k) {
    for(j in 1:k) {
      if(i==j) {
        pre_sigma[, psi] <- values$cov[,i]
      } else {
        pre_sigma[ , psi] <- values$cov[,mm(i,j)]
      }
      psi <- psi + 1
    }
  }
  sigma_list <- lapply(split(pre_sigma,1:nrow(pre_sigma)), function(z) {
    m <- matrix(z, nrow=k)
    c2c <- cov2cor(m)
    for(i in 1:(k-1)) {
      for(j in (i+1):k) {
        if(c2c[i,j] > 0.95) {
          # limit correlation to 0.95
          m[i,j] <- m[j,i] <- 0.95 * sqrt(m[i,i]*m[j,j])
        }
      }
    }
    eig <- eigen(m)
    if(min(eig$values) < 0 | min(eig$values)/max(eig$values) < (.Machine$double.eps)^0.25) {
      m <- as(Matrix::nearPD(m, posd.tol=1e-4)$mat, "matrix")
    }
    chol_robust(m)
  } )
  names(sigma_list) <- values$id
  values$sigma_list <- sigma_list
  return(values)
}
