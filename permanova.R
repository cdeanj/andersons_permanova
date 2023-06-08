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
  n <- nrow(mat)
  d_sqr <- 0
  for(i in 1:n) {
    end <- i
    for(j in 1:end) {
      d_sqr <- d_sqr + mat[i,j]^2
    }
  }
  ss_w <- d_sqr/n
  return(ss_w)
}