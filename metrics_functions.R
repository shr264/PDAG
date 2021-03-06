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

getL = function(B){
  B = (abs(B) != 0)*1
  L = lower.tri(B)*B - t(B*upper.tri(B))
  L = (L[lower.tri(L)])
  return(L)
}

getLhatlist = function(Bhatlist, debug = FALSE){
  Lhatlist = vector("list",length(Bhatlist))
  if(debug){
    message(Lhatlist)
  }
  for (i in 1:length(Bhatlist)){
    Lhatlist[[i]] = getL(Bhatlist[[i]])
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
  L = getL(B)
  Lhatlist = getLhatlist(Bhatlist)
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
