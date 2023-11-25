#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
float inp_cpp(const NumericVector& x, const NumericVector& y, const NumericVector& argval){

  int n = x.length();
  float s = 0;

  for(int i = 0; i<n-1; i++)
    s = s + x(i+1)*y(i+1)*(argval(i+1)-argval(i));

  return(s);
}

//[[Rcpp::export]]
NumericMatrix inp_matrix(const NumericMatrix& V, const NumericVector& argval){
  int N = V.nrow();
  NumericMatrix A(N,N);

  for(int i = 0; i<N; i++)
    for(int j = 0; j<N; j++)
      A(i,j) = inp_cpp(V(i,_),V(j,_), argval);

  return(A);
}

//[[Rcpp::export]]
float l2norm_cpp(const NumericVector& x){
  int n = x.length();
  float s = 0;
  for(int i = 0; i<n; i++)
    s = s + (x(i)*x(i));

  return(s);
}

//[[Rcpp::export]]
double ker_phi1_cpp(double x, double y){
  double s = x-y;
  float a = -(std::sqrt((s*s)))/2;
  return(a);
}

//[[Rcpp::export]]
double ker_log_cpp(double x, double y){
  double s = x-y;
  float a = -std::log(1+(s*s));
  return(a);
}

//[[Rcpp::export]]
double ker_exp_cpp(double x, double y){
  double s = x-y;
  float a = s*s;
  a = -1+std::exp(-a/2);
  return(a);
}

//[[Rcpp::export]]
NumericMatrix Gram_matrix_cpp_phi1(const NumericVector& x){
  int N = x.length();
  NumericMatrix A(N,N);

  for(int i=0; i<N; i++){
    for(int j=0; j<i; j++){
      A(i,j) = ker_phi1_cpp(x(i),x(j));
      A(j,i) = A(i,j);
    }
  }

  return(A);
}

//[[Rcpp::export]]
NumericMatrix Gram_matrix_cpp_exp(NumericVector x){
  int N = x.length();
  NumericMatrix A(N,N);

  for(int i=0; i<N; i++){
    for(int j=0; j<i; j++){
      A(i,j) = ker_exp_cpp(x(i),x(j));
      A(j,i) = A(i,j);
    }
  }

  return(A);
}

//[[Rcpp::export]]
NumericMatrix Gram_matrix_cpp_log(NumericVector x){
  int N = x.length();
  NumericMatrix A(N,N);

  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      A(i,j) = ker_log_cpp(x(i),x(j));
      A(j,i) = A(i,j);
    }
  }

  return(A);
}

//[[Rcpp::export]]
List Gram_calc_phi1(NumericMatrix S){
  int N = S.nrow();
  List L(N);

  for(int i = 0; i<N; i++)
    L[i] = Gram_matrix_cpp_phi1(S(i,_));

  return(L);
}

//[[Rcpp::export]]
List Gram_calc_exp(NumericMatrix S){
  int N = S.nrow();
  List L(N);

  for(int i = 0; i<N; i++)
    L[i] = Gram_matrix_cpp_exp(S(i,_));

  return(L);
}

//[[Rcpp::export]]
List Gram_calc_log(NumericMatrix S){
  int N = S.nrow();
  List L(N);

  for(int i = 0; i<N; i++)
    L[i] = Gram_matrix_cpp_log(S(i,_));

  return(L);
}


//[[Rcpp::export]]
double kmmd_cpp(const NumericMatrix& D, const int& n, const int& m){

  int N = n+m;
  double s = 0;

  for(int i=0; i<n; i++)
    for(int j=0; j<i; j++)
      s = s + 2*D(i,j)/(n*n);

  for(int i=n; i<N; i++)
    for(int j=n; j<i; j++)
      s = s + 2*D(i,j)/(m*m);

  for(int i=n; i<N; i++)
    for(int j=0; j<n; j++)
      s = s - 2*D(i,j)/(n*m);


  return(s);
}

//[[Rcpp::export]]
NumericMatrix test(const NumericMatrix& M, const IntegerVector& ind){
  int col = M.ncol();
  NumericMatrix M2(col,col);

  for(int j = 0; j<col; j++)
    for(int i = 0; i<col; i++)
      M2(i,j) = M(ind[i], ind[j]);

  return(M2);
}

//[[Rcpp::export]]
double op1(const List L, const IntegerVector& l, int& n, int& m){

  int k = l.length();
  NumericVector s(k);

  double T=0;

  for(int j = 0; j<k; j++)
    s(j) = kmmd_cpp(test(L[l(j)],l),n,m);

  for(int i=0;i<n;i++)
    T = T + 0.5*s(i)/n;

  for(int j=n;j<k;j++)
    T = T + 0.5*s(j)/m;

  return(T);
}

//[[Rcpp::export]]
NumericVector T_perm(const IntegerMatrix& W, List L, int n, int m){
  int R = W.ncol();

  NumericVector S(R);

  for(int i = 0; i<R; i++)
    S(i) = op1(L, W(_,i), n,m);

  return(S);
}

