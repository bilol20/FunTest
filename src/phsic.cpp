#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double inp(NumericVector x, NumericVector y, NumericVector time){
  int n = x.length();
  double s = 0;
  for(int i = 0;i<n-1;i++){
    s += x(i)*y(i)*(time(i+1)-time(i));
  }
  return(s);
}

//[[Rcpp::export]]
List median_huristic(List list, NumericVector time){
  int d = list.length();
  NumericMatrix A = list[0];
  int n = A.nrow();
  int i,j,k,l;
  List s(n-2);
  for(i = 0; i<n-2;i++){
    NumericMatrix B(n);
    for(j = i+1; j<n-1; j++){
      for(k = j+1; k<n; k++){
        for(l = 0; l<d;l++){
          NumericMatrix A = list[l];
          B(j,k) = B(j,k) + inp(A(i,_),A(j,_)-A(k,_),time)*inp(A(i,_),A(j,_)-A(k,_),time);
        }
      }
    }
    s[i] = B;
  }
  return(s);
}


//[[Rcpp::export]]
NumericMatrix proj_cpp(NumericMatrix x, NumericVector argval){
  int n = x.nrow();
  NumericMatrix A(n);
  for(int i =0; i<n; i++){
    for(int j =0; j<n;j++){
      A(i,j) = inp(x(i,_),x(j,_),argval);
    }
  }
  return(A);
}


//[[Rcpp::export]]
NumericVector perm(int n) {
  Function f("sample");
  return( f(Named("x") = n, Named("size") = n));
}

//[[Rcpp::export]]
double summat(NumericMatrix A)
{
  double s = 0;
  int p = A.ncol();
  int n = A.nrow();
  for(int i=0;i<n;i++){
    for(int j=0; j<p;j++){
      s+=A(i,j);
    }
  }
  return(s);
}

//[[Rcpp::export]]
double sumvec(NumericVector A)
{
  double s = 0;
  int p = A.length();
  for(int j=0; j<p;j++){
    s+=A(j);
  }
  return(s);
}


//[[Rcpp::export]]
NumericVector colSum_cpp(NumericMatrix A)
{
  int p = A.ncol();
  int n = A.nrow();
  NumericVector s(p);
  for(int j=0;j<p;j++){
    for(int i=0; i<n;i++){
      s[j]=s[j]+A(i,j);
    }
  }
  return(s);
}

//[[Rcpp::export]]
NumericMatrix operation1(NumericMatrix l, NumericMatrix A)
{
  int p = A.ncol();
  int n = A.nrow();
  NumericMatrix s(n,p);
  for(int i=0; i<n;i++){
    for(int j=0;j<p;j++){
      s(i,j)=l(i,j)*A(i,j);
    }
  }
  return(s);
}

//[[Rcpp::export]]
NumericVector operation2(NumericVector l, NumericVector A)
{
  int p = A.length();
  NumericVector s(p);
  for(int j=0;j<p;j++){
    s[j]=l(j)*A(j);
  }
  return(s);
}


//[[Rcpp::export]]
double g_cpp(List K)
{
  int d = K.length();
  NumericMatrix A = K[0];
  int len = A.nrow();
  int p = A.ncol();
  NumericMatrix s1(len,p);
  for(int i=0; i<len; i++){
    for(int j=0; j<p; j++){
      s1(i,j) = 1;
    }
  }
  double term2 = 1;
  NumericVector s(p,2.0/len);
  for (int j=0; j<d;j++) {
    s1 = operation1(s1, K[j]);
    term2 =  (term2 * summat(K[j]))/(len*len);
    s = operation2(s , colSum_cpp(K[j]))/len;
  }
  double term1 = summat(s1)*1/(len*len);
  double term3 = sumvec(s);
  double T = term1+term2-term3;
  return(T);
}

//[[Rcpp::export]]
NumericMatrix gaussian_kernel(NumericVector x, double b)
{
  int n = x.length();
  NumericMatrix A(n);
  for(int  i= 0; i<n; i++)
  {
    for(int j = 0; j<n; j++)
    {
      double s = -(x[i]-x[j])*(x[i]-x[j])/(2*b*b);
      A(i,j) = std::exp(s);
    }
  }
  return(A);
}

//[[Rcpp::export]]
double stat(List D, double b)
{
  int d = D.length();
  NumericMatrix A = D[0];
  int N = A.nrow();
  NumericVector q(N);
  for(int i = 0; i<N; i++)
  {
    List K(d);
    for(int j=0; j<d; j++)
    {
      NumericMatrix B = D[j];
      K[j] = gaussian_kernel( B(_,i), b);
    }
    q(i) = g_cpp(K);
  }
  double g = sumvec(q)/N ;
  return(g);
}

//[[Rcpp::export]]
NumericMatrix rowcolsam(NumericMatrix A)
{
  int N = A.ncol();
  NumericVector l = perm(N);
  NumericMatrix B(N);
  for(int i = 0; i<N;i++){
    for(int j = 0; j<N; j++){
      B(i,j) = A(l[i]-1,l[j]-1);
    }
  }
  return(B);
}



//[[Rcpp::export]]
NumericVector permutation(List D, int R, double b)
{
  int d = D.length();
  NumericMatrix A = D[0];
  int N = A.nrow();
  NumericVector q(N);
  NumericVector T(R);
  List l(d);
  List K(d);
  List D1(d);
  D1[0] = D[0];
  for(int k = 0; k<R; k++)
  {
    for(int j=1; j<d; j++)
    {
      D1[j] = rowcolsam(D[j]);
    }
    T(k) = stat(D1,b);
  }
  return(T);
}


//[[Rcpp::export]]
NumericMatrix multi_permutation(List D, int R, NumericVector b)
{
  int d = D.length();
  NumericMatrix A = D[0];
  int N = A.nrow();
  NumericVector q(N);
  int n = b.length();
  NumericMatrix T(R,n);
  List l(d);
  List K(d);
  List D1(d);
  for(int k = 0; k<R; k++)
  {
    for(int j=0; j<d; j++)
    {
      D1[j] = rowcolsam(D[j]);
    }
    for(int i = 0; i<n; i++)
    {
      T(k,i) = stat(D1,b[i]);
    }
  }
  return(T);
}










