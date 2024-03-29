
band.selection = function(data,time,q=0.5)
{
  s = unlist(median_huristic(data,time))
  s = s[-which(s==0)]
  sqrt(quantile(s,q)/2)
}

phsic.test = function(data, time, alpha = 0.05, R = 100)
{
  l = as.numeric(lapply(data, nrow))
  p = length(data)
  x = rep(nrow(data[[1]]),p)
  if(all(l == x) & alpha>0 & R>1){
    b = band.selection(data,time, 0.5)
    D = lapply(data, proj_cpp, argval = time)
    T_0 = stat(D, b)
    T = permutation(D,R,b)
    q = quantile(T,1-alpha)
    e = list(
      Estimate = T_0,
      p.value = (sum(T>T_0)+1)/(R+1),
      permutation.repetition = R,
      sample.size = l[1])
    return(e)
  }else{
    if(q == 0) print("ERROR: Please provide a valid order of the quantile.")
    if(alpha == 0) print("ERROR: Significance level can not be 0.")
    if(R<2) print("ERROR: Please provide a positive value for the numer of permutations.")
    if(all(l!=x)) print("ERROR:The size of the two samples are different.")
  }
}

phsic.MultiTest = function(data, time, q = c(0.5,0.4,0.3,0.2,0.1,0.05), alpha = 0.05,
                           R = 100, plot = FALSE)
{
  l = as.numeric(lapply(data, nrow))
  p = length(data)
  x = rep(nrow(data[[1]]),p)
  if(all(l == x) & alpha>0 & R>1){
    n = length(q)
    b = numeric(n)
    T_0 = numeric(n)
    S = numeric(n)
    p = numeric(n)
    for(i in 1:n){b[i] = band.selection(data,time, q[i])}
    D = lapply(data, proj_cpp, argval = time)
    for(i in 1:n) T_0[i] = stat(D, b[i])
    A = multi_permutation(D,R,b)
    T = apply(A,1,sum)
    q1 = (sum(T>sum(T_0))+1)/(R+1)
    T = apply(A,1,max)
    q2 = (sum(T>max(T_0))+1)/(R+1)
    for(i in 1:n){p[i] = (sum(A[,i]>T_0[i])+1)/(R+1)}
    e = list(
      Estimates = T_0,
      Sum.Estimate = sum(T_0),
      Max.Estimate = max(T_0),
      p.values = p,
      p.value.sum = q1,
      p.value.max = q2,
      permutation.repetition = R,
      sample.size = l[1])
    return(e)
  }else{
    if(q == 0) print("ERROR: Please provide a valid order of the quantile.")
    if(alpha == 0) print("ERROR: Significance level can not be 0.")
    if(R<2) print("ERROR: Please provide a positive value for the numer of permutations.")
    if(all(l!=x)) print("ERROR:The size of the two samples are different.")
  }
}


sigma_multi = function(data,time,q=0.5)
{
  l = median_huristic_2(data,time)
  d = length(l)
  l1 = list(d)
  for(i in 1:d){
    s = unlist(l[[i]])
    s = s[-which(s==0)]
    l1[[i]] = sqrt(quantile(s,q)/2)
  }
  return(l1)
}

Test_multi = function(data, time, alpha = 0.05, R = 100)
{
  l = as.numeric(lapply(data, nrow))
  p = length(data)
  x = rep(nrow(data[[1]]),p)
  if(all(l == x) & alpha>0 & R>1){
    b = sigma_multi(data,time, 0.5)
    D = lapply(data, proj_cpp, argval = time)
    T_0 = stat_2(D, b)
    T = permutation_2(D,R,b)
    q = quantile(T,1-alpha)
    e = list(
      Estimate = T_0,
      p.value = (sum(T>T_0)+1)/(R+1),
      permutation.repetition = R,
      sample.size = l[1])
    return(e)
  }else{
    if(q == 0) print("ERROR: Please provide a valid order of the quantile.")
    if(alpha == 0) print("ERROR: Significance level can not be 0.")
    if(R<2) print("ERROR: Please provide a positive value for the numer of permutations.")
    if(all(l!=x)) print("ERROR:The size of the two samples are different.")
  }
}
