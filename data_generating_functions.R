library(MASS)

genB_rand_200 <- function(p = 200, a = 0.3, b = 0.7, m = NULL, z = 0.05, s = 0.5){
  set.seed(12345)
  D = rep(1,p)
  plower = p*(p-1)/2
## off-diagonals
  T = diag(p)
  T[upper.tri(T)] = 0
  T[lower.tri(T)] = (ifelse(runif(plower)<s, -1, 1) *
		   ifelse(runif(plower)<z,  1, 0) *
		   runif(plower, a, b))
  L = diag(1.0/sqrt(D)) %*% T   # cholesky factor
  if(is.null(m)||(m==1)){
      Border = sample(1:p)
  } else {
      M = floor(p/m)
      Border = c()
      pos = 1
      while((pos+M+1)<p){
          Border = append(Border,sample(pos:(pos+M-1),M))
          pos = pos+M
      }
      Border = append(Border,sample(pos:p,p-pos+1))
  }
  B = L[Border,Border]
  adj = (B!=0)
  return(list(adj = adj,B = B))
}

genB_rand_100 <- function(p = 100, a = 0.3, b = 0.7, m = NULL, z = 0.05, s = 0.5){
  set.seed(12345)
  D = rep(1,p)
  plower = p*(p-1)/2
  ## off-diagonals
  T = diag(p)
  T[upper.tri(T)] = 0
  T[lower.tri(T)] = (ifelse(runif(plower)<s, -1, 1) *
                       ifelse(runif(plower)<z,  1, 0) *
                       runif(plower, a, b))
  L = diag(1.0/sqrt(D)) %*% T   # cholesky factor
  if(is.null(m)||(m==1)){
    Border = sample(1:p)
  } else {
    M = floor(p/m)
    Border = c()
    pos = 1
    while((pos+M+1)<p){
      Border = append(Border,sample(pos:(pos+M-1),M))
      pos = pos+M
    }
    Border = append(Border,sample(pos:p,p-pos+1))
  }
  B = L[Border,Border]
  adj = (B!=0)
  return(list(adj = adj,B = B))
}

genB_rand_140 <- function(p = 40, a = 0.3, b = 0.7, m = NULL, z = 0.05, s = 0.5){
  set.seed(12345)
  D = rep(1,p)
  plower = p*(p-1)/2
  ## off-diagonals
  T = diag(p)
  T[upper.tri(T)] = 0
  T[lower.tri(T)] = (ifelse(runif(plower)<s, -1, 1) *
                       ifelse(runif(plower)<z,  1, 0) *
                       runif(plower, a, b))
  L = diag(1.0/sqrt(D)) %*% T   # cholesky factor
  if(is.null(m)||(m==1)){
    Border = sample(1:p)
  } else {
    M = floor(p/m)
    Border = c()
    pos = 1
    while((pos+M+1)<p){
      Border = append(Border,sample(pos:(pos+M-1),M))
      pos = pos+M
    }
    Border = append(Border,sample(pos:p,p-pos+1))
  }
  B = L[Border,Border]
  adj = (B!=0)
  return(list(adj = adj,B = B))
}

genB_mult <- function(Yeast1, a = 0.3, b = 0.7, m = NULL){
  p = dim(Yeast1)[1]
  set.seed(12345)
  adj = topSort(Yeast1)
  if(is.null(m)||(m==1)){
      Border = 1:p
      } else {
      M = floor(p/m)
      cat('M = ', M, '\n')
      Border = c()
      pos = 1
      while((pos+M+1)<p){
          cat('pos = ', pos,'\n')
          cat('pos + M -1 = ', (pos + M - 1),'\n')
          Border = append(Border,sample(pos:(pos+M-1),M))
          pos = pos+M
      }
      Border = append(Border,sample(pos:p,p-pos+1))
  }
  adj = adj[Border,Border]
  adj = t(adj)
  Bmat = matrix(runif(p^2,a,b),nrow = p, ncol = p)
  signmat = matrix(0,nrow = p, ncol = p)
  signmat = ifelse(runif(p^2)<0.5,-1,1)
  B = adj*Bmat*signmat
  diag(B) = rep(1,p)
  return(list(adj = adj,B = B))
}

genB_mult_Yeast1 <- function(a = 0.3, b = 0.7, m = NULL){
  load('amat.Rdata')
  p = dim(Yeast1)[1]
  set.seed(12345)
  adj = topSort(Yeast1)
  if(is.null(m)||(m==1)){
    Border = 1:p
  } else {
    M = floor(p/m)
    cat('M = ', M, '\n')
    Border = c()
    pos = 1
    while((pos+M+1)<p){
      cat('pos = ', pos,'\n')
      cat('pos + M -1 = ', (pos + M - 1),'\n')
      Border = append(Border,sample(pos:(pos+M-1),M))
      pos = pos+M
    }
    Border = append(Border,sample(pos:p,p-pos+1))
  }
  adj = adj[Border,Border]
  adj = t(adj)
  Bmat = matrix(runif(p^2,a,b),nrow = p, ncol = p)
  signmat = matrix(0,nrow = p, ncol = p)
  signmat = ifelse(runif(p^2)<0.5,-1,1)
  B = adj*Bmat*signmat
  diag(B) = rep(1,p)
  return(list(adj = adj,B = B))
}

genB_mult_Yeast2 <- function(a = 0.3, b = 0.7, m = NULL){
  load('amat.Rdata')
  p = dim(Yeast2)[1]
  set.seed(12345)
  adj = topSort(Yeast2)
  if(is.null(m)||(m==1)){
    Border = 1:p
  } else {
    M = floor(p/m)
    cat('M = ', M, '\n')
    Border = c()
    pos = 1
    while((pos+M+1)<p){
      cat('pos = ', pos,'\n')
      cat('pos + M -1 = ', (pos + M - 1),'\n')
      Border = append(Border,sample(pos:(pos+M-1),M))
      pos = pos+M
    }
    Border = append(Border,sample(pos:p,p-pos+1))
  }
  adj = adj[Border,Border]
  adj = t(adj)
  Bmat = matrix(runif(p^2,a,b),nrow = p, ncol = p)
  signmat = matrix(0,nrow = p, ncol = p)
  signmat = ifelse(runif(p^2)<0.5,-1,1)
  B = adj*Bmat*signmat
  diag(B) = rep(1,p)
  return(list(adj = adj,B = B))
}

genB_mult_Yeast3 <- function(a = 0.3, b = 0.7, m = NULL){
  load('amat.Rdata')
  p = dim(Yeast3)[1]
  set.seed(12345)
  adj = topSort(Yeast3)
  if(is.null(m)||(m==1)){
    Border = 1:p
  } else {
    M = floor(p/m)
    cat('M = ', M, '\n')
    Border = c()
    pos = 1
    while((pos+M+1)<p){
      cat('pos = ', pos,'\n')
      cat('pos + M -1 = ', (pos + M - 1),'\n')
      Border = append(Border,sample(pos:(pos+M-1),M))
      pos = pos+M
    }
    Border = append(Border,sample(pos:p,p-pos+1))
  }
  adj = adj[Border,Border]
  adj = t(adj)
  Bmat = matrix(runif(p^2,a,b),nrow = p, ncol = p)
  signmat = matrix(0,nrow = p, ncol = p)
  signmat = ifelse(runif(p^2)<0.5,-1,1)
  B = adj*Bmat*signmat
  diag(B) = rep(1,p)
  return(list(adj = adj,B = B))
}

genB_mult_Ecoli1 <- function(a = 0.3, b = 0.7, m = NULL){
  load('amat.Rdata')
  p = dim(Ecoli1)[1]
  set.seed(12345)
  adj = topSort(Ecoli1)
  if(is.null(m)||(m==1)){
    Border = 1:p
  } else {
    M = floor(p/m)
    cat('M = ', M, '\n')
    Border = c()
    pos = 1
    while((pos+M+1)<p){
      cat('pos = ', pos,'\n')
      cat('pos + M -1 = ', (pos + M - 1),'\n')
      Border = append(Border,sample(pos:(pos+M-1),M))
      pos = pos+M
    }
    Border = append(Border,sample(pos:p,p-pos+1))
  }
  adj = adj[Border,Border]
  adj = t(adj)
  Bmat = matrix(runif(p^2,a,b),nrow = p, ncol = p)
  signmat = matrix(0,nrow = p, ncol = p)
  signmat = ifelse(runif(p^2)<0.5,-1,1)
  B = adj*Bmat*signmat
  diag(B) = rep(1,p)
  return(list(adj = adj,B = B))
}

genB_mult_Ecoli2 <- function(a = 0.3, b = 0.7, m = NULL){
  load('amat.Rdata')
  p = dim(Ecoli2)[1]
  set.seed(12345)
  adj = topSort(Ecoli2)
  if(is.null(m)||(m==1)){
    Border = 1:p
  } else {
    M = floor(p/m)
    cat('M = ', M, '\n')
    Border = c()
    pos = 1
    while((pos+M+1)<p){
      cat('pos = ', pos,'\n')
      cat('pos + M -1 = ', (pos + M - 1),'\n')
      Border = append(Border,sample(pos:(pos+M-1),M))
      pos = pos+M
    }
    Border = append(Border,sample(pos:p,p-pos+1))
  }
  adj = adj[Border,Border]
  adj = t(adj)
  Bmat = matrix(runif(p^2,a,b),nrow = p, ncol = p)
  signmat = matrix(0,nrow = p, ncol = p)
  signmat = ifelse(runif(p^2)<0.5,-1,1)
  B = adj*Bmat*signmat
  diag(B) = rep(1,p)
  return(list(adj = adj,B = B))
}

genB <- function(Yeast1, a = 0.3, b = 0.7, m = NULL){
  set.seed(12345)
    p = dim(Yeast1)[1]
    if(is.null(m)){
        m = floor(p/2)
    }
    set.seed(12345)
    adj = topSort(Yeast1)
    Border = c(sample(1:m,m),sample((m+1):p,p-m))
    adj = adj[Border,Border]
    adj = t(adj)
    Bmat = matrix(runif(p^2,a,b),nrow = p, ncol = p)
    signmat = matrix(0,nrow = p, ncol = p)
    signmat = ifelse(runif(p^2)<0.5,-1,1)
    B = adj*Bmat*signmat
    #diag(B) = runif(p,2,5)
    diag(B) = rep(1,p)
    return(list(adj = adj,B = B))
}

genB_Yeast1_informative <- function(a = 0.3, b = 0.7, m = 5){
  load('amat.Rdata')
  set.seed(12345)
  p = dim(Yeast1)[1]
  if(is.null(m)){
    m = floor(p/2)
  }
  adj = Yeast1
  Border1 = sample(c(11,9,4,16,24),5)
  Border2 = sample(setdiff(1:50,c(11,9,4,16,24)),45)
  Border = c(Border1,Border2)
  adj = adj[Border,Border]
  adj = t(adj)
  Bmat = matrix(runif(p^2,a,b),nrow = p, ncol = p)
  signmat = matrix(0,nrow = p, ncol = p)
  signmat = ifelse(runif(p^2)<0.5,-1,1)
  B = adj*Bmat*signmat
  diag(B) = rep(1,p)
  return(list(adj = adj,B = B))
}

genB_Yeast1_informative_sec <- function(Yeast1, a = 0.3, b = 0.7, m = 5){
  set.seed(12345)
  p = dim(Yeast1)[1]
  if(is.null(m)){
    m = floor(p/2)
  }
  adj = Yeast1
  Bset1 = c(11,9,4,13,5,7,12,2,16,24,33,26,37,38,21)
  Bset2 = setdiff(1:50,Bset1)
  Border1 = sample(Bset1,length(Bset1))
  Border2 = sample(Bset2,length(Bset2))
  Border = c(Border1,Border2)
  adj = adj[Border,Border]
  adj = t(adj)
  Bmat = matrix(runif(p^2,a,b),nrow = p, ncol = p)
  signmat = matrix(0,nrow = p, ncol = p)
  signmat = ifelse(runif(p^2)<0.5,-1,1)
  B = adj*Bmat*signmat
  diag(B) = runif(p,2,5)
  diag(B) = rep(1,p)
  return(list(adj = adj,B = B))
}

genB_Yeast1_noninformative <- function(a = 0.3, b = 0.7, m = 5){
  load('amat.Rdata')
  set.seed(12345)
  p = dim(Yeast1)[1]
  if(is.null(m)){
    m = floor(p/2)
  }
  adj = Yeast1
  Border1 = sample(c(11,9,10,16,13),5)
  Border2 = sample(setdiff(1:50,c(11,9,10,16,13)),45)
  Border = c(Border1,Border2)
  adj = adj[Border,Border]
  adj = t(adj)
  Bmat = matrix(runif(p^2,a,b),nrow = p, ncol = p)
  signmat = matrix(0,nrow = p, ncol = p)
  signmat = ifelse(runif(p^2)<0.5,-1,1)
  B = adj*Bmat*signmat
  diag(B) = rep(1,p)
  return(list(adj = adj,B = B))
}

genB_Yeast1_noninformative_sec <- function(Yeast1, a = 0.3, b = 0.7, m = 5){
  set.seed(12345)
  p = dim(Yeast1)[1]
  if(is.null(m)){
    m = floor(p/2)
  }
  adj = Yeast1
  Bset1 = 1:15
  Bset2 = setdiff(1:50,Bset1)
  Border1 = sample(Bset1,length(Bset1))
  Border2 = sample(Bset2,length(Bset2))
  Border = c(Border1,Border2)
  Border = c(Border1,Border2)
  adj = adj[Border,Border]
  adj = t(adj)
  Bmat = matrix(runif(p^2,a,b),nrow = p, ncol = p)
  signmat = matrix(0,nrow = p, ncol = p)
  signmat = ifelse(runif(p^2)<0.5,-1,1)
  B = adj*Bmat*signmat
  diag(B) = runif(p,2,5)
  diag(B) = rep(1,p)
  return(list(adj = adj,B = B))
}



