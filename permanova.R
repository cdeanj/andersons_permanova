library(vegan)

data(dune)
data(dune.env)

bc_dist <- vegan::vegdist(dune, method = "bray")

pnova <- adonis2(bc_dist ~ Management, data = dune.env)
pnova

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