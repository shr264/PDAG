require(MASS)

get_metric = function(L, Lhat, i, debug = FALSE){
  if(debug){
    cat('L = ', L, '\n')
    cat('Lhat = ', Lhat, '\n')}
  L = (L == i)
  Lhat = (Lhat == i)
  tp = sum((L == 1)*(Lhat==1))
  tn = sum((L == 0)*(Lhat==0))
  fn = sum((L == 1)*(Lhat==0))
  fp = sum((L == 0)*(Lhat==1))
  if(debug){
    cat('tp = ', tp, '\t')
    cat('tn = ', tn, '\t')
    cat('fn = ', fn, '\t')
    cat('fp = ', fp, '\n')}
  list(tp = tp, tn = tn, fn = fn, fp = fp, fpr = fp/(fp+tn), tpr = tp/(tp+fn))
}

cpDagtoDag = function(Bhat, B){
  if(is.null(B)){
    return(Bhat)
  }
  if(sum((Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B != 0))>0){
    Bhat[t(Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B != 0)] = 0
  }
  
  if(sum(Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B == 0)>0){
    unif = runif(1)
    if(unif < 0.5){x = c(1,0)
    } else{x = c(0,1)}
    
    Bhat[Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B == 0] = x
  }
  
  return(Bhat)
}

cpDagtoDagW = function(Bhat, B){
  if(is.null(B)){
    return(Bhat)
  }
  if(sum((Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B != 0))>0){
    Bhat[(Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B != 0)] = 0
  }
  
  if(sum(Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B == 0)>0){
    unif = runif(1)
    if(unif < 0.5){x = c(1,0)
    } else{x = c(0,1)}
    
    Bhat[Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B == 0] = x
  }
  
  return(Bhat)
}

cpDagtoDagBest = function(Bhat, B){
  if(is.null(B)){
    return(Bhat)
  }
  if(sum((Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B != 0))>0){
    Bhat[(Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B == 0)] = 0
  }
  
  return(Bhat)
}

cpDagtoDagWorst = function(Bhat, B){
  if(is.null(B)){
    return(Bhat)
  }
  if(sum((Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B != 0))>0){
    Bhat[t(Bhat == t(Bhat) & row(Bhat) != col(Bhat) & Bhat  != 0 & B == 0)] = 0
  }
  
  return(Bhat)
}


getL = function(B, B_truth=NULL){
  if(!is.null(B_truth)){
    B = cpDagtoDagWorst(B,B_truth) 
  }
  B = (abs(B) != 0)*1
  L = lower.tri(B)*B - t(B*upper.tri(B))
  L = (L[lower.tri(L)])
  return(L)
}

getLhatlist = function(Bhatlist, B=NULL, debug = FALSE){
  Lhatlist = vector("list",length(Bhatlist))
  if(debug){
    message(Lhatlist)
  }
  for (i in 1:length(Bhatlist)){
    Lhatlist[[i]] = getL(Bhatlist[[i]], B_truth=B)
  }
  return(Lhatlist)
}

get_avg_metrics2 = function(B,
                            Bhatlist,
                            n,
                            p,
                            method,
                            Btype,
                            mean_time,
                            seed,
                            debug=FALSE){
  
  # check whether and NAs in return Bhatlist
  
  i = 1
  while(i <= length(Bhatlist)){
    if( anyNA(Bhatlist[[i]]) | is.null(Bhatlist[[i]]) ){
      Bhatlist = Bhatlist[-i]
    }
    else{
      i = i + 1
    }
  }
  
  if(length(Bhatlist)==0){
    auc = rep(0,17)
    names_auc = rep('auc',17)
    for (j in c(-1,0,1)){
      auc[j+2] = NA
      names_auc[j+2] = paste('norm_auc',toString(j),sep='')
      auc[3+j+2] = NA
      names_auc[3+j+2] = paste('auc',toString(j),sep='')
      auc[6+j+2] = NA
      names_auc[6+j+2] = paste('max_auc',toString(j),sep='')}
    auc[10] = NA
    names_auc[10] = 'avg_norm_auc'
    auc[11] = NA
    names_auc[11] = 'avg_auc'
    auc[12] = NA
    names_auc[12] = 'avg_max_auc'
    auc[13] = n
    names_auc[13] = 'n'
    auc[14] = p
    names_auc[14] = 'p'
    auc[15] = method
    names_auc[15] = 'method'
    auc[16] = Btype
    names_auc[16] = 'Btype'
    auc[17] = mean_time
    names_auc[17] = 'mean_time'
    auc[18] = seed
    names_auc[18] = 'seed'
    names(auc) = names_auc
    
    return(auc)
  }
  L = getL(B)
  Lhatlist = getLhatlist(Bhatlist, B)
  auc = rep(0,17)
  names_auc = rep('auc',17)
  for (j in c(-1,0,1)){
    fpr_ = rep(0,length(Bhatlist))
    tpr_ = rep(0,length(Bhatlist))
    for (i in 1:length(Bhatlist)){
      fpr_tpr = get_metric(L, Lhatlist[[i]], j)
      fpr_[i] = fpr_tpr$fpr
      tpr_[i] = fpr_tpr$tpr
    }
    if(debug){
      cat('j = ', j, '\t')
      cat('fpr: ', fpr_, '\t')
      cat('tpr: ', tpr_, '\n')
    }
    ord.fp = order(fpr_)
    auc[j+2] = trapezint(fpr_[ord.fp],tpr_[ord.fp],min(fpr_),max(fpr_))/(max(fpr_)-min(fpr_))
    names_auc[j+2] = paste('norm_auc',toString(j),sep='')
    auc[3+j+2] = trapezint(fpr_[ord.fp],tpr_[ord.fp],min(fpr_),max(fpr_))
    names_auc[3+j+2] = paste('auc',toString(j),sep='')
    auc[6+j+2] = (max(fpr_)-min(fpr_))
    names_auc[6+j+2] = paste('max_auc',toString(j),sep='')
  }
  auc[10] = mean(auc[1:3])
  names_auc[10] = 'avg_norm_auc'
  auc[11] = mean(auc[4:6])
  names_auc[11] = 'avg_auc'
  auc[12] = mean(auc[7:9])
  names_auc[12] = 'avg_max_auc'
  auc[13] = n
  names_auc[13] = 'n'
  auc[14] = p
  names_auc[14] = 'p'
  auc[15] = method
  names_auc[15] = 'method'
  auc[16] = Btype
  names_auc[16] = 'Btype'
  auc[17] = mean_time
  names_auc[17] = 'mean_time'
  auc[18] = seed
  names_auc[18] = 'seed'
  names(auc) = names_auc
  
  return(auc)
}
