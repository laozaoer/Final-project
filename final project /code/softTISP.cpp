#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

inline arma::vec softThresholding(arma::vec z, double lambda){
  int n=z.n_elem;
  arma::vec x(n);
  for(int i=0;i<n;i++){
    if(z(i) > lambda)
      x(i) = z(i)-lambda;
    else if(z(i) < (-lambda))
      x(i) = z(i)+lambda;
    else 
      x(i) = 0;
  }
  return x;
}



  // [[Rcpp::export]]
  arma::vec softTISP(arma::mat X0, arma::vec y0,double lambda){
    int n=X0.n_rows;
    int p=X0.n_cols;
   
   arma::mat S(n,p,arma::fill::ones);
   arma::mat Xcent=X0-S*diagmat(mean(X0));
   arma::mat Xstand;
   arma::mat colnorm(p,p,arma::fill::zeros);
   colnorm=diagmat(stddev(X0,1,0));
   colnorm=colnorm.i()/sqrt(n);
   //for(int i=0;i<=p-1;i++){
   //colnorm(i,i)=pow(sqrt(sum(Xcent.col(i)%Xcent.col(i))),-1);
   //}
   arma::mat X=Xcent*colnorm;
    arma::vec s;
    s.ones(n);
    arma::vec y = y0-mean(y0)*s;
    
    arma::mat Xt=X.t();
    arma::mat Sigma=Xt*X;
    arma::vec eigenvalue=eig_sym(Sigma);
    double m=max(eigenvalue);
    arma::mat I;
    I.eye(p,p);
    arma::mat A=I-Sigma/m;
    arma::vec c=Xt*y/m;
    lambda=lambda/m;
    
    arma::vec beta0;
    beta0.zeros(p);
    arma::vec t;
    arma::vec beta;
    int k=0;
    do{
      beta=softThresholding(A*beta0+c,lambda);
      t=beta0;
      beta0=beta;
      k++;
    } while (norm(beta-t,2)>1e-3);
    beta=colnorm*beta;
    return beta;
  }

// [[Rcpp::export]]
double softTISPcrossValidation(arma::mat X, arma::vec y, int K, arma::vec Lambda){
  int N = X.n_rows;
  int size=N/K;
  int testnum=N-size;
  int M=Lambda.n_elem;
  arma::mat Xcopy;
  arma::mat Xtrain;
  arma::mat Xtest;
  arma::vec ytrain;
  arma::vec ytest;
  arma::vec betahat;
  arma::vec yhat;
  arma::mat AveForecastError(M,K);
  double lambda;
  for(int j=0;j<=M-1;j++ ){
    for(int i=1;i<=K;i++){
      Xtest=X.rows((i-1)*size,i*size-1);
      ytest=y.subvec((i-1)*size,i*size-1);  //where i=K is included. because
      //without loss of generality,the N is the interger multiple of K.
      
      if(i==1){
        Xtrain=X.rows(size,N-1);
        ytrain=y.subvec(size,N-1);
      }
      else if(i==K){
        Xtrain=X.rows(0,(K-1)*size-1);
        ytrain=y.subvec(0,(K-1)*size-1);
      }
      else{
        Xcopy=X;
        Xcopy.shed_rows((i-1)*size,i*size-1);
        Xtrain=Xcopy;
        ytrain=join_vert(y.subvec(0,(i-1)*size-1),y.subvec(i*size,N-1));
      }
      betahat=softTISP(Xtrain, ytrain,Lambda(j));
      yhat=Xtest*betahat;
      AveForecastError(j,i-1)=sum((ytest-yhat)%(ytest-yhat))/testnum;
    }
  }
  arma::vec perform=mean(AveForecastError,1);
  double minimun=min(perform);
  for(int m=0;m<=M-1;m++){
    if(perform(m)==minimun){
      lambda=Lambda(m);
      break;}
    else
      continue;
  }
  return lambda;
}




