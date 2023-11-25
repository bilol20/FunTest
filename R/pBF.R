pBF.test = function(X,Y,argval, R = 200, kernel = c("L2","log","exp")){
  n = nrow(X)
  m = nrow(Y)
  V = rbind(X,Y)
  N = nrow(V)
  S = inp_matrix(V, grid)

  if(kernel == "L2"){
    L = Gram_calc_phi1(S)
  }
  if(kernel == "log"){
    L = Gram_calc_log(S)
  }
  if(kernel == "exp"){
    L = Gram_calc_exp(S)
  }
  T = sapply(1:N, function(t){
    kmmd_cpp(L[[t]],n,m)
  })
  T = 0.5*mean(T[1:n])+0.5*mean(T[n+1:m])
  W = replicate(R, sample(N))-1

  #T1 = T_perm(W,L,n,m)

  T1 = T_perm(W,L,n,m)

  pval = (sum(T1>T)+1)/(R+1)
  return(list(Stat = T, p.value = pval))
}


