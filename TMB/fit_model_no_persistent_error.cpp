#include <TMB.hpp>  // Links in the TMB libraries
# include <vector> // to try the .insert in sim section

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_INTEGER(nrec);
  DATA_VECTOR(L); // obsv length
  DATA_VECTOR(age);
  
  PARAMETER(log_sigma); // for residuals
  PARAMETER(log_Linf_hyper); // hyper population means
  PARAMETER(log_K_hyper);
  PARAMETER(log_L1_hyper);
  PARAMETER(log_sigma_PE); // process error variance
  PARAMETER_VECTOR(log_PE); // process error; random effects

  // back transform parameters
  Type sigma = exp(log_sigma);
  Type sigma_PE = exp(log_sigma_PE);
  vector<Type> PE = exp(log_PE); // process error

  Type Linf = exp(log_Linf_hyper);
  Type L1 = exp(log_L1_hyper);
  Type K = exp(log_K_hyper);

  vector<Type> L_hat(nrec); // store predicted growth at age
  vector<Type> L_hat_tmp(nrec); // for storing temporary L_hat for logsitic smoothing function
  vector<Type> resid(nrec);
  vector<Type> log_resid(nrec); 
  vector<Type> delta(nrec); // store increments

  Type w2 = 0; // weights for logistic smoothing function
  Type w1 = 1;
  Type alpha = 2; // slope for logistic smoothing function, fixed
        
  // grow the fish
  for(int i=0; i<nrec; i++){
    if(age(i) == 1){
      L_hat(i) = L1;
    }else{
      // calculate L_hat_tmp to then use in weighted average
      delta(i-1) = (Linf - L_hat(i-1)) * (1-exp(-K));
      L_hat_tmp(i) = L_hat(i-1) + (delta(i-1) * PE(i));

      // set the weights for logistic smoothing function, based on distance
      // from Linf (for that individual) and the calculated L_hat (for this time step)
      w2 = 1/(1+exp(-alpha*((0.99*Linf) - L_hat_tmp(i)))); // condition
      w1 = 1 - w2;
      // calculate L_hat (real) using weighted average to prevent value from
      // exceeding Linf
      L_hat(i) = (w2 * L_hat_tmp(i)) + (w1 * Linf);
    }
  }

  // Calculate residuals
  resid = L - L_hat;
  log_resid = log(L) - log(L_hat);
  
  // NLL for the observed and expected increments (true length)
  Type nll_length = -sum(dnorm(log(L), (log(L_hat)- pow(sigma,2)/2),
                              sigma, true));
  // NLL for process error
  // back-transformed expected value to be 0 -> avg length, not median length

  Type nll_pe = -sum(dnorm(log(PE),
                           -(pow(sigma_PE,2)/2), sigma_PE, true));

  // Get the full log likelihood
  Type f = nll_length + nll_pe;

  REPORT(PE);
  REPORT(L);
  REPORT(L_hat);
  REPORT(L_hat_tmp);
  REPORT(resid);
  REPORT(age);
  REPORT(f);
  REPORT(nll_pe);
  REPORT(sigma_PE);
  REPORT(sigma);
  REPORT(delta);
  REPORT(Linf);
  REPORT(K);
  REPORT(L1);
  REPORT(log_PE);
  
  return f;
  
}
