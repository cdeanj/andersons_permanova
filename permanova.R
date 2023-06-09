data(dune)
data(dune.env)

bray_curtis <- function(mat) {
  nr <- nrow(mat)
  bc <- matrix(nrow = nr, ncol = nr)
  
  for (i in 1:nr) {
    for (j in 1:nr) {
      vec_i <- mat[i,]
      vec_j <- mat[j,]
      
      C_ij <- sum(pmin(vec_i,vec_j))
      
      s_i <- sum(vec_i)
      s_j <- sum(vec_j)
      
      d <- 1-(2*C_ij)/(s_i+s_j)
      
      bc[i,j] <- d
    }
  }
  return(bc)
}

sst <- function(mat) {
  n <- nrow(mat)
  d_sqr <- 0
  for(i in 1:n) {
    end <- i
    for(j in 1:end) {
      d_sqr <- d_sqr + mat[i,j]^2
    }
  }
  ss_t <- d_sqr/n
  return(ss_t)
}

ssw <- function(mat) {
  grps <- sort(unique(colnames(mat)))
  n_grps <- length(grps)
  
  l <- as.list(rep(0, n_grps))
  names(l) <- grps
  
  c <- as.list(rep(0, n_grps))
  names(c) <- grps
  
  n <- nrow(mat)
  grp_i <- ''
  grp_j <- ''
  for(i in 1:n) {
    end <- i
    grp_i <- colnames(mat)[i]
    c[[grp_i]] <- c[[grp_i]] + 1
    for(j in 1:end) {
      grp_i <- colnames(mat)[i]
      grp_j <- colnames(mat)[j]
      if(grp_i == grp_j) {
        l[[grp_i]] <- l[[grp_i]] + mat[i,j]^2
      }
    }
  }
  ss_w <- 0
  for(k in 1:n_grps) {
    ss_w <- ss_w + l[[k]]/c[[k]]
  }
  return(ss_w[[1]])
}

ssa <- function(ss_t, ss_w) {
  return(ss_t - ss_w)
}

psuedo_f_ratio <- function(ss_a, ss_w, a, N) {
  return((ss_a/(a-1))/(ss_w/(N-a)))
}

r_squared <- function(ss_w, ss_t) {
  return(1-(ss_w/ss_t))
}

permutation_test <- function(mat, perms) {
  f_pi <- c()
  col_labs <- colnames(mat)
  for(i in 1:perms) {
    cols <- sample(col_labs)
    permuted_matrix <- mat
    colnames(permuted_matrix) <- cols
    pseudo_ss_t <- sst(permuted_matrix)
    pseudo_ss_w <- ssw(permuted_matrix)
    pseudo_ss_a <- ssa(pseudo_ss_t, pseudo_ss_w)
    pseudo_f <- psuedo_f_ratio(pseudo_ss_a, pseudo_ss_w, 4, 20)
    f_pi <- c(f_pi, pseudo_f)
  }
  return(f_pi)
}

p_value <- function(f_pi_list, f_ratio, perms) {
  return(length(which(f_pi_list > f_ratio))/perms)
}

bc_mat <- bray_curtis(dune)

ss_total <- sst(bc_mat)
ss_total

grp <- dune.env$Management
grp

colnames(bc_mat) <- grp

ss_within <- ssw(bc_mat)
ss_within

ss_between <- ssa(ss_total, ss_within)
ss_between

f_ratio <- psuedo_f_ratio(ss_between, ss_within, 4, 20)
f_ratio

r2 <- r_squared(ss_within, ss_total)
r2

f_pii <- permutation_test(bc_mat, 1000)
p_value(f_pii, f_ratio, 1000)
