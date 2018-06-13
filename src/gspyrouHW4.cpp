

#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

//' Implementation of Viterbi Algorithm on a Symmetric Continuous-time Markov Chain HMM using Rcpp
//'
//' Function that implements a Viterbi algorithm, given likelihood of observed data,timestamp of each observation and a parameter theta.
//'
//' @param ts size m Numeric vector containing the timestamp of each observation
//' @param theta double type value which represents the parameter theta in the model of transition probability
//' @param obs n x m Numeric matrix containing the likelihood of the observed data
//'
//' @return the Viterbi path(path)
//'
//' @examples
//' obs <- matrix(c(0.88,0.10,0.88,0.10,0.02,0.30,0.02,0.30,0.10,0.60),2,5)
//' theta <- log(2)
//' ts <- c(1,2,3,4,5)
//' ctmcViterbi(ts,theta,obs)
//' @export
// [[Rcpp::export]]
IntegerVector ctmcViterbi(NumericVector ts, double theta, NumericMatrix obs)
{
  int m = (int)ts.size();
  int n = (int)obs.nrow();
  if ( obs.ncol() != m )
    stop("The input matrix does not conform to the other parameters");
  IntegerVector viterbiPath(m);
  // creating a delta matrix, as in the lecture notes and fill it with -1
  // also we create a phi matrix
  NumericMatrix delta(n, m);
  std::fill(delta.begin(), delta.end(), -1);
  IntegerMatrix phi(n, m);
  NumericMatrix final_matrix(n,m); // matrix for final values
  
  double same_s,not_same_s,temp,size_p = n,pi = 1/double(n); // these will hold values and help us calculate the path
  
  for(int i=0; i < n; ++i) {
    delta(i,0) = pi * obs(i, 0);
    if (delta(i,0) > final_matrix(0,0)) {
      final_matrix(0,0) = delta(i,0);
      final_matrix(1,0) = i;
    }
  }
  
  for(int t=1; t < m; ++t) {
    same_s =  1.0 - (((size_p-1)/size_p)*(1.0 - exp(-1*theta*(ts[t] - ts[t-1])))); //choice 1
    not_same_s = (1/size_p) * (1.0 - exp(-1*theta*(ts[t] - ts[t-1]))); // choise 2
    
    for(int i=0; i < n; ++i) {
      if(i==final_matrix(1,t-1))
      {temp = delta(i,t-1)*same_s*obs(i, t);}
      else{
        temp = max(delta(final_matrix(1,t-1),t-1)*not_same_s*obs(i, t),delta(i,t-1)*same_s*obs(i, t));}
      
      if (temp > delta(i,t)) {
        delta(i,t) = temp;
        if (i==final_matrix(1,t-1)) {
          phi(i,t) = i;
        }
        else if (delta(i,t-1)*same_s*obs(i, t) < delta(final_matrix(1,t-1),t-1)*not_same_s*obs(i, t)){
          phi(i,t) = final_matrix(1,t-1);
        }
        else {
          phi(i,t) = i;
        }
        if (delta(i,t) > final_matrix(0,t)) {
          final_matrix(0,t) = delta(i,t);
          final_matrix(1,t) = i;
        }
      }
      
    }
  }
  double ml = -1;
  
  for(int i=0; i < n; ++i) {
    if ( ml < delta(i,m-1) ) {
      ml = delta(i,m-1);
      viterbiPath[m-1] = i;
    }
  }  
  for(int i=m-1; i > 0; --i) {
    viterbiPath[i-1] = phi(viterbiPath[i],i);
  }
  
  return viterbiPath;
}

/* end of Viterbi Algorithm / it should return the Viterbi path*/


/* Now we start with the two loops
* forwardLoop,backwardLoop and then 
* we will implement them inside the whole
* forwardbackward algorithm
*/

void forwardLoop(NumericVector& ts, NumericMatrix& obs, 
                 double& pi, NumericMatrix& alpha,
                 double& theta,int& m, int& n) {
  
  double sum_t,c_value; //c_value is going to hold the "constant" value for each case
  // while sum_t will hold the total sum at this point
  
  //define some values for alpha, as initial ones
  for(int i=0; i < n; ++i) {
    alpha(i,0) = pi * obs(i,0);
  }
  
  for(int t=1; t < m; ++t) {
    c_value = (1.0 - exp(-1*theta*(ts[t] - ts[t-1]))),sum_t = 0;;
    // total sum for alphas for each t and all n / one a time
    for(int i=0; i < n; ++i) {
      sum_t += alpha(i,t-1);
    }
    for (int i=0; i<n; ++i) {
      alpha(i,t) = ((c_value/n)*sum_t + alpha(i,t-1)*(1-c_value))*obs(i, t);
    }
  }
}

void backwardLoop(NumericVector& ts, NumericMatrix& obs, 
                  NumericMatrix& beta, 
                  double& theta,int& m, int& n) {
  
  double sum_t,c_value; //c_value is going to hold the "constant" value for each case
  // while sum_t will hold the total sum at this point - same logic with forwardloop
  //at first we initialize
  for(int i=0; i < n; ++i) 
    beta(i,m-1) = 1;
  
  // now we are going to fill the beta's
  for(int t=m-2; t >=0; --t) {
    c_value = (1.0 - exp(-1*theta*(ts[t+1] - ts[t]))),sum_t = 0;;
    
    for(int i=0; i < n; ++i) {
      sum_t += (beta(i,t+1)*obs(i, t+1));
    }
    for (int i=0; i<n; ++i) {
      beta(i,t) = (c_value/n)*sum_t + beta(i,t+1)*(1-c_value)*obs(i, t+1);
    }
  }
}

//' Implementation of ForwardBackward Algorithm on a Symmetric Continuous-time Markov Chain HMM
//'
//' Implementing the ForwardBackward algorithm,given likelihood of observed data,timestamp of each observation and a parameter theta.
//'
//' @param ts size m Numeric vector containing the timestamp of each observation
//' @param theta double type value which represents the parameter theta in the model of transition probability
//' @param obs n x m Numeric matrix containing the likelihood of the observed data
//'
//' @return a n x m matrix containing the conditional probabilities of each observation
//'
//' @examples
//' obs <- matrix(c(0.88,0.10,0.88,0.10,0.02,0.30,0.02,0.30,0.10,0.60),2,5)
//' theta <- log(2)
//' ts <- c(1,2,3,4,5)
//' ctmcForwardBackward(ts,theta,obs)
//' @export
// [[Rcpp::export]]
NumericMatrix ctmcForwardBackward(NumericVector ts, double theta, NumericMatrix obs) {
  int m = (int)ts.size();  
  int n = (int)obs.nrow();
  if ( obs.ncol() != m )
    stop("The input matrix does not conform to the other parameters");
  NumericMatrix condProb(n, m);
  // initialize pi, it will be equal to 1/n for every case as we have uniform initial states
  double pi = 1/double(n);
  //alpha and beta will store the values from the two loops we defined previously
  NumericMatrix alpha(n, m),beta(n, m);
  
  //the loops that we defined before - we run them
  forwardLoop(ts, obs, pi, alpha, theta,m,n);
  backwardLoop(ts, obs, beta, theta,m,n);
  // now we have alpha,beta and we can calculate the probabilities
  for(int t=0; t < m; ++t) {
    double sum = 0;
    
    for(int i=0; i < n; ++i) 
      sum += ( alpha(i,t) * beta(i,t) );
    for(int i=0; i < n; ++i) 
      condProb(i,t) = alpha(i,t) * beta(i,t) / sum;
  }
  return condProb; // return the matrix of prob's
}

