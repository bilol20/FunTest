inp = function(x,y,argval){
  return(sum(x[-1]*y[-1]*diff(argval)))
}
norm = function(x){ sum(abs(x))}

ker = function(x,y) return(0.5*(norm(x)+norm(y)-norm(x-y)))

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

pMMD.test = function(X,Y,argval, R = 200){
  n = nrow(X)
  m = nrow(Y)
  V = rbind(X,Y)
  N = nrow(V)
  S = sapply(1:N, function(x){
    sapply(1:N, function(y){
      inp(V[x,],V[y,],argval)
    })
  })
  D = unlist(lapply(1:N, function(j){
    x1 = S[j,1:n]
    y1 = S[j,n+1:m]
    return(dist(c(x1,y1))^2)
  }))
 # sigma = median(D)
  L = lapply(1:N, function(k){
    x1 = S[k,1:n]
    y1 = S[k,n+1:m]
    v = c(x1,y1)
    D = sapply(1:N, function(i){
      sapply(1:N, function(j){
        ker(v[i],v[j])
      })
    })
    return(D)
  })
  T = mean(sapply(1:N, function(t){
    kmmd(L[[t]],n,m)
  }))
  T1 = unlist(lapply(1:R, function(i){
    l = sample(N)
    mean(sapply(1:N, function(t){
      D1 = L[[l[t]]]
      D1 = rowcolsam(D1,l)
      return(kmmd(D1,n,m))
    }))
  }
  ))
  pval = (sum(T1>T)+1)/(R+1)
  return(list(Stat = T, p.value = pval))
}

