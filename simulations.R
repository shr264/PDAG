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
library(foreach)
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

ncores <- 20

# Test on p = 40

values = expand.grid(list(
  Methods = c('pcalg_custom',
              'partial2','partial3','partial4',
              'pcalg_addBG2','pcalg_addBG3','pcalg_addBG4',
              'lingam_custom'),
  nlambda = c(30),
  Btypes = c('genB_rand_40'),
  Ns = c(50,75,100,200),
  Seeds = 1:20,
  m = 4,
  m1 = 10,
  m2 = 20,
  m3 = 30), stringsAsFactors = FALSE)

#mapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m=values$m, m1=values$m1, m2=values$m2, m3=values$m3)

table0 = mcmapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m=values$m, m1=values$m1, m2=values$m2, m3=values$m3,
                  mc.cores=ncores)

print(table0)

# Table 1
# generate_tables(Methods = c('pcalg_custom','ccdr_paper_t','partial2','pcalg_addBG2', 'lingam_custom'),
#                 Btypes = c('genB_mult_Yeast1','genB_mult_Yeast2','genB_mult_Yeast3',
#                            'genB_mult_Ecoli1','genB_mult_Ecoli2'),
#                 Ns = c(40,50,100,200),
#                 Seeds = 1:20,
#                 m = 2,
#                 m1 = 25)

print('table1')

values = expand.grid(
  list(
    Methods = c('pcalg_custom','ccdr_paper_t','partial2','pcalg_addBG2', 'lingam_custom'),
    nlambda = c(30),
    Btypes = c('genB_mult_Yeast1','genB_mult_Yeast2','genB_mult_Yeast3',
              'genB_mult_Ecoli1','genB_mult_Ecoli2'),
    Ns = c(40,50,100,200),
    Seeds = 1:20,
    m = c(2),
    m1 = c(25)), stringsAsFactors = FALSE)


#mapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m=values$m, m1=values$m1)

table1 = mcmapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m=values$m, m1=values$m1,
                  mc.cores=ncores)

save(table1, file = "table1.RData")

# Table 2
# generate_tables(Methods = c('pcalg_custom','ccdr_paper_t',
#                             'partial2','partial3','partial4',
#                             'pcalg_addBG2','pcalg_addBG3','pcalg_addBG4',
#                             'lingam_custom'),
#                 Btypes = c('genB_mult_Yeast1','genB_mult_Yeast2','genB_mult_Yeast3',
#                            'genB_mult_Ecoli1','genB_mult_Ecoli2'),
#                 Ns = c(40,50,100,200),
#                 Seeds = 1:20,
#                 m = 4,
#                 m1 = 12,
#                 m2 = 24,
#                 m3 = 36)

values = expand.grid(list(
Methods = c('pcalg_custom','ccdr_paper_t',
             'partial2','partial3','partial4',
             'pcalg_addBG2','pcalg_addBG3','pcalg_addBG4',
             'lingam_custom'), 
nlambda = c(30), 
Btypes = c('genB_mult_Yeast1','genB_mult_Yeast2','genB_mult_Yeast3',
             'genB_mult_Ecoli1','genB_mult_Ecoli2'),
  Ns = c(40,50,100,200),
  Seeds = 1:20,
  m = 4,
  m1 = 12,
  m2 = 24,
  m3 = 36), stringsAsFactors = FALSE)

#mapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m=values$m, m1=values$m1, m2=values$m2, m3=values$m3)

print('table2')

table2 = mcmapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m=values$m, m1=values$m1, m2=values$m2, m3=values$m3,
                  mc.cores=ncores)

save(table2, file = "table2.RData")

# Table 3
# generate_tables(Methods = c('partial2'),
#                 Btypes = c('genB_Yeast1_informative','genB_Yeast1_noninformative'),
#                 Ns = c(40,50,100,200),
#                 Seeds = 1:20,
#                 m1 = 5)

print('table3')

values = expand.grid(list(
  Methods = c('partial2'),
   nlambda = c(30),
   Btypes = c('genB_Yeast1_informative','genB_Yeast1_noninformative'),
   Ns = c(40,50,100,200),
   Seeds = 1:20,
   m1 = 5), stringsAsFactors = FALSE)

#mapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m1=values$m1)

table3 = mcmapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m1=values$m1,
                  mc.cores=ncores)

save(table3, file = "table3.RData")

# Table 4
# generate_tables(Methods = c('pcalg_custom',
#                             'partial2','partial3','partial4',
#                             'pcalg_addBG2','pcalg_addBG3','pcalg_addBG4',
#                             'lingam_custom'),
#                 Btypes = c('genB_rand_100'),
#                 Ns = c(50,75,100,200),
#                 Seeds = 1:20,
#                 m = 4,
#                 m1 = 25,
#                 m2 = 50,
#                 m3 = 75)

print('table4')

values = expand.grid(list(
  Methods = c('pcalg_custom',
              'partial2','partial3','partial4',
              'pcalg_addBG2','pcalg_addBG3','pcalg_addBG4',
              'lingam_custom'),
  nlambda = c(30),
  Btypes = c('genB_rand_100'),
  Ns = c(50,75,100,200),
  Seeds = 1:20,
  m = 4,
  m1 = 25,
  m2 = 50,
  m3 = 75), stringsAsFactors = FALSE)

#mapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m=values$m, m1=values$m1, m2=values$m2, m3=values$m3)

table4 = mcmapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m=values$m, m1=values$m1, m2=values$m2, m3=values$m3,
                  mc.cores=ncores)

save(table4, file = "table4.RData")

# # Table 5
# generate_tables(Methods = c('pcalg_custom',
#                             'partial2','partial3','partial4',
#                             'pcalg_addBG2','pcalg_addBG3','pcalg_addBG4',
#                             'lingam_custom'),
#                 Btypes = c('genB_rand_200'),
#                 Ns = c(100,150,200,400),
#                 Seeds = 1:20,
#                 m = 4,
#                 m1 = 50,
#                 m2 = 100,
#                 m3 = 150)

print('table5')

values = expand.grid(list(
  Methods = c('pcalg_custom',
              'partial2','partial3','partial4',
              'pcalg_addBG2','pcalg_addBG3','pcalg_addBG4',
              'lingam_custom'),
  Btypes = c('genB_rand_200'),
  Ns = c(100,150,200,400),
  Seeds = 1:20,
  m = 4,
  m1 = 50,
  m2 = 100,
  m3 = 150), stringsAsFactors = FALSE)

#mapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m=values$m, m1=values$m1, m2=values$m2, m3=values$m3)

table5 = mcmapply(get_metrics_by_method, method=values$Methods, nlambda=values$nlambda, Btype=values$Btypes, n=values$Ns, seed=values$Seeds, m=values$m, m1=values$m1, m2=values$m2, m3=values$m3,
                  mc.cores=ncores)

save(table5, file = "table5.RData")


