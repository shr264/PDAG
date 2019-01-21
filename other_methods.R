source("helper_functions.R")
library(sparsebn)
library(pcalg)

pcalg_custom <- function(X,l=NULL, a = NULL, m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL,eps = 10^(-4),maxitr = 100, init=NULL){
  (n = dim(X)[1])
  (p = dim(X)[2])
  V = sapply(1:p,toString)
  time = proc.time()[3]
  pc.fit <- pc(suffStat = list(C = cor(X), n = n),
               indepTest = gaussCItest, ## indep.test: partial correlations
               alpha=a, labels = V, verbose = FALSE)
  B = (as.matrix(as(pc.fit, "amat")))
  diag(B) = 1
  time = time - proc.time()[3]
  return(list(B=as.matrix(B),itr = 1, time = time))
}

ccdr_paper_t <- function(X,l,a=NULL,m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL,eps = 10^(-4),maxitr = 100, init=NULL){
  (n = dim(X)[1])
  (p = dim(X)[2])
  (S = (1/n)*t(X)%*%X)
  if(is.null(init)){
    B = Bold = diag(p)
  } else {B = Bold = init}
  for(i in 1:p){
    B[i,i] = Bii(S,B,i)
  }
  time = proc.time()[3]
  data = sparsebnData(X, type = "continuous")
  dags <- estimate.dag(data, lambdas = l)
  dags.fit = estimate.parameters(dags, data = data)
  Btemp = diag(p)-t(as.matrix(dags.fit[[1]]$coefs))
  diag(Btemp) = diag(B)
  B = Btemp
  time = time - proc.time()[3]
  return(list(B=as.matrix(B),itr = 1, time = time))
}