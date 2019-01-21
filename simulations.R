library(MASS)
library(Matrix)
library(ggm)
library(sparsebn)
library(ROC)
library(lassoshooting)
library(pcalg)
library(parallel)
library(igraph)
library(caret)
library(partitionDAG)
source('metrics_functions.R')
source('data_generating_functions.R')
source('other_methods.R')
load('amat.Rdata')
set.seed(12345)

get_metrics_by_method = function(method, nlambda, Btype, m, n, seed, m1=NULL,m2=NULL,m3=NULL,m4=NULL,m5=NULL,m6=NULL,m7=NULL,m8=NULL,m9=NULL){
  lambda = seq(0.001,3, length= nlambda)^3
  alpha = exp(-seq(0.1,16, length=nlambda))
  B = get(Btype)(m = m)$B
  p = dim(B)[1]
  omega = t(B) %*% B            # omega
  sigma = solve(omega)
  sigma[abs(sigma)<1e-10] = 0   # set numerical error to zero
  set.seed(23456 + seed) #seed for generating data
  X = mvrnorm(n, mu=rep(0, p), Sigma=sigma) # observations
  X = scale(X,center = TRUE, scale=TRUE)

  Bhatlist = vector("list", nlambda)
  time_vec = rep(0,nlambda)
  for(i in 1:length(lambda)){
    B_ = get(method)(X=X,l=lambda[i], a=alpha[i], m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,m6=m6,m7=m7,m8=m8,m9=m9)
    Bhatlist[[i]] = B_$B
    time_vec[i] = B_$time
  }

  return(get_avg_metrics2(B, Bhatlist, n, p, method, Btype, mean(time_vec), seed, debug=FALSE))
}

generate_tables = function(Methods,Btypes,Ns,Seeds, m = NULL, m1=NULL, m2=NULL, m3=NULL){
  for (AA in Methods){
    for (BB in Btypes){
      for (CC in Ns){
        for (DD in Seeds){
          if(exists('table1')){
            table1 = rbind(table1,get_metrics_by_method(method = AA,
                                                        nlambda = 30,
                                                        Btype = BB,
                                                        m = m,
                                                        n = CC,
                                                        seed = DD,
                                                        m1=m1,
                                                        m2=m2,
                                                        m3=m3))
          } else {
            table1 = get_metrics_by_method(method = AA,
                                           nlambda = 30,
                                           Btype = BB,
                                           m = m,
                                           n = CC,
                                           seed = DD,
                                           m1=m1,
                                           m2=m2,
                                           m3=m3)
          }
        }
      }
    }
  }
  return(table1)
}

# Table 1
generate_tables(Methods = c('pcalg_custom','ccdr_paper_t','partial2'),
                Btypes = c('genB_mult_Yeast1','genB_mult_Yeast2','genB_mult_Yeast3',
                           'genB_mult_Ecoli1','genB_mult_Ecoli2'),
                Ns = c(40,50,100,200),
                Seeds = 1:20,
                m = 2,
                m1 = 25)

# Table 2
generate_tables(Methods = c('pcalg_custom','ccdr_paper_t',
                            'partial2','partial3','partial4'),
                Btypes = c('genB_mult_Yeast1','genB_mult_Yeast2','genB_mult_Yeast3',
                           'genB_mult_Ecoli1','genB_mult_Ecoli2'),
                Ns = c(40,50,100,200),
                Seeds = 1:20,
                m = 4,
                m1 = 12,
                m2 = 24,
                m3 = 36)

# Table 3
generate_tables(Methods = c('partial2'),
                Btypes = c('genB_Yeast1_informative','genB_Yeast1_noninformative'),
                Ns = c(40,50,100,200),
                Seeds = 1:20,
                m1 = 5)

# Table 4
generate_tables(Methods = c('pcalg_custom','ccdr_paper_t',
                            'partial2','partial3','partial4'),
                Btypes = c('genB_rand_100'),
                Ns = c(50,75,100,200),
                Seeds = 1:20,
                m = 4,
                m1 = 25,
                m2 = 50,
                m3 = 75)

# Table 5
generate_tables(Methods = c('pcalg_custom','ccdr_paper_t',
                            'partial2','partial3','partial4'),
                Btypes = c('genB_rand_200'),
                Ns = c(100,150,200,400),
                Seeds = 1:20,
                m = 4,
                m1 = 50,
                m2 = 100,
                m3 = 150)

