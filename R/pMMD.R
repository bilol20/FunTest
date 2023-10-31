inp = function(x,y,argval){
  return(sum(x[-1]*y[-1]*diff(argval)))
}


l2norm = function(x){sum(x^2)}

ker_phi1 = function(x,y) return(-sqrt(l2norm(x-y))/2)
ker_log = function(x,y) return(-log(1+l2norm(x-y)))
ker_exp = function(x,y) return(-1+exp(-l2norm(x-y)/2))



kmmd = function(D,n,m){
  T = sum(D[1:n,1:n])/n^2+sum(D[n+1:m,n+1:m])/m^2-2*sum(D[1:n,n+1:m])/n/m
  return(T)
}

rowcolsam = function(A,l){
  n1 = nrow(A)
  m1 = ncol(A)
  A1 = matrix(0, nrow = n1, ncol = m1)
  for(i in 1:n1){
    for(j in 1:m1){
      A1[i,j] = A[l[i],l[j]]
    }
  }
  return(A1)
}


Gram.matrix = function(x,ker){
  N = length(x)
  A = matrix(0,N,N)
  for(i in 1:N){
    for(j in 1:i){
      A[j,i] = A[i,j] = ker(x[i]-x[j])
    }
  }
  return(A)
}


pMMD.test = function(X,Y,argval, R = 200, kernel = c("L2","log","exp")){
  if(kernel == "L2"){
    ker = ker_phi1
  }
  if(kernel == "log"){
    ker = ker_log
  }
  if(kernel == "exp"){
    ker = ker_exp
  }
  n = nrow(X)
  m = nrow(Y)
  V = rbind(X,Y)
  N = nrow(V)
  S = sapply(1:N, function(x){
    sapply(1:N, function(y){
      inp(V[x,],V[y,],argval)
    })
  })
  #D = unlist(lapply(1:N, function(j){
  #  x1 = S[j,1:n]
 #   y1 = S[j,n+1:m]
 #   return(dist(c(x1,y1))^2)
 # }))
 # sigma = median(D)
  #D = unlist(lapply(1:N, function(j){
  #  x1 = S[j,1:n]
  #   y1 = S[j,n+1:m]
  #   return(dist(c(x1,y1))^2)
  # }))
  # sigma = median(D)
  #L = lapply(1:N, function(k){
  #  x1 = S[k,1:n]
  # y1 = S[k,n+1:m]
  #  v = c(x1,y1)
  #  D = sapply(1:N, function(i){
  #    sapply(1:N, function(j){
  #      ker(v[i]-v[j])
  #    })
  #  })
  #  return(D)
  #})
  L = lapply(1:N, function(k){
    v = S[k,]
    D = Gram.matrix(v,ker)
    return(D)
  })
  T = sapply(1:N, function(t){
    kmmd(L[[t]],n,m)
  })
  T = 0.5*mean(T[1:n])+0.5*mean(T[n+1:m])
  T1 = unlist(lapply(1:R, function(i){
    l = sample(N)
    F = sapply(1:N, function(t){
      D1 = L[[l[t]]]
      D1 = rowcolsam(D1,l)
      return(kmmd(D1,n,m))
    })
    return(0.5*mean(F[1:n])+0.5*mean(F[n+1:m]))
  }
  ))
  pval = (sum(T1>T)+1)/(R+1)
  return(list(Stat = T, p.value = pval))
}

