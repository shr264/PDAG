addBgKnowledge_helper = function(B, m1=NULL, m2=NULL, m3=NULL){
  for(i in (m1+1):m2){
    for(j in (m2+1):m3){
      if(B[i,j] == 1 & B[j,i] == 1){
        B <- addBgKnowledge(gInput = B, x = j, y = i)
      }
    }
  }
  return(B) 
}

addBgKnowledge_helper_t = function(B, m1=NULL, m2=NULL, m3=NULL){
  for(i in (m1+1):m2){
    for(j in (m2+1):m3){
      if(B[i,j] == 1 & B[j,i] == 1){
        B <- addBgKnowledge(gInput = B, x = i, y = j)
      }
    }
  }
  return(B) 
}


dfs <- function(node,p,graph){
  nodes_to_visit = setdiff(1:p,node)
  paths = which(graph[node,]>0)
  if(length(paths)==0){
    return(paths)
  } else {
    while(length(nodes_to_visit)>0){
      currentnode = nodes_to_visit[1]
      nodes_to_visit = setdiff(nodes_to_visit,currentnode)
      paths = unique(append(paths,which(graph[,]>0)))
    }
    return(paths)
  }
}

softhresh <- function(x,l){
  sign(x)*max(abs(x)-l,0)
}

Bii <- function(S,B,i){
  temp = sum(S[i,-i]*B[i,-i])
  out = (-2*temp + sqrt(4*temp^2+8*S[i,i]))/(4*S[i,i])
  out
}

Bki <- function(S,B,k,i,l){
  out = softhresh(-2*sum(S[i,-i]*B[k,-i])/(4*S[i,i]),l/(4*S[i,i]))
  out
}

qBki <- function(S,Bki,B,k,i,l){
  out1 = S[i,i]*Bki^2 + 2*Bki*sum(S[i,-i]*B[k,-i]) + l*abs(Bki)
  out2 = 0
  if(out1<out2){return(list(q = out1, b = Bki))
  } else {return(list(q = out2, b = 0))}
}

penloglik <- function(B,S,l){
  Omega = t(B)%*%B
  out = sum(diag(Omega%*%S)) - log(det(Omega)) + l*sum(abs(B))
  out
}

qBki_ll <- function(S,Bki,B,k,i,l){
  B[k,i] = 0
  old_ll = penloglik(B,S,l)
  B[k,i] = Bki
  B[i,k] = 0
  new_ll = penloglik(B,S,l)
  if(new_ll<old_ll){return(list(q = new_ll, b = Bki))
  } else {return(list(q = old_ll, b = 0))}
}

Bik <- function(S,B,i,k,l){
  out = softhresh(-2*sum(S[k,-k]*B[i,-k])/(4*S[k,k]),l/(4*S[k,k]))
  out
}

qBik <- function(S,Bik,B,i,k,l){
  out1 = S[k,k]*Bik^2 + 2*Bik*sum(S[k,-k]*B[i,-k]) + l*abs(Bik)
  out2 = 0
  if(out1<out2){return(list(q = out1, b = Bik))
  } else {return(list(q = out2, b = 0))}
}

qBik_ll <- function(S,Bik,B,i,k,l){
  B[i,k] = 0
  old_ll = penloglik(B,S,l)
  B[i,k] = Bik
  B[k,i] = 0
  new_ll = penloglik(B,S,l)
  if(new_ll<old_ll){return(list(q = new_ll, b = Bik))
  } else {return(list(q = old_ll, b = 0))}
}

Lii <- function(S,L,i){
  temp = sum(S[i,-i]*L[i,-i])
  out = (-2*temp + sqrt(4*temp^2+8*S[i,i]))/(4*S[i,i])
  out
}

Lij <- function(S,L,i,j){
  out = -sum(S[j,-j]*L[i,-j])/S[j,j]
  out
}

CSCS2 <- function( Y, lambda, L=NULL, maxitr=100, tol=1e-4, warmstart=FALSE) {
  
  #### inputs
  ## Y: n by p matrix of data
  ## lambda: l1-penalty parameter
  ## maxitr: maximum number of iterations allowed (diagonal/off-diagonal optimization)
  ## tol: if maximum difference of iterates is less than tol, consider converged
  ## warmstart: warmstarting actually made the runs slower (overhead may be too expensive)
  
  #### outputs
  ## L: lower triangular matrix of cholesky factor computed with CSCS algorithm
  ## itr_log: (p-1)-vector of number of iterations
  ## eps_log: (p-1)-vector of number maximum difference for considering convergence
  
  
  n = nrow(Y)
  p = ncol(Y)
  
  if (is.null(L)) L = diag(p)
  
  S = (t(Y)%*%Y)/n
  
  itr_log = eps_log = NULL
  
  L[1, 1] = 1/sqrt(S[1,1])
  for (k in 2:p){ ## Loop in Algorithm 2
    
    ## nu_k vector (equation 2.5)
    nuk_old = nuk_new = c(rep(0, k-1), 1)
    r = 0
    
    repeat {      ## Loop in Algorithm 1
      
      r = r + 1
      
      km1_ = 1:(k-1)    ## 1, ..., k-1 is off-diagonal elements indices
      
      ## Update off-diagonal terms
      hk = lassoshooting(XtX    =  S[km1_, km1_, drop=FALSE],
                         Xty    = -nuk_old[k] * S[km1_, k],
                         lambda =  0.5*lambda)
      nuk_new[km1_] = hk$coefficients
      
      ## Update diagonal term
      sumterm = sum(nuk_new[km1_] * S[k, km1_])
      nuk_new[k] = (-sumterm + sqrt(sumterm^2 + 4*S[k, k]))/(2*S[k, k])
      
      ## Check convergence
      maxdiff = max(abs(nuk_new - nuk_old))
      if (maxdiff < tol || r >= maxitr){
        L[k, 1:k] = nuk_new
        eps_log = c(eps_log, maxdiff)
        itr_log = c(itr_log, r)
        break
      } else {
        nuk_old = nuk_new
      }
    }
  }
  
  list(L=L, itr=itr_log, eps=eps_log)
  
}

