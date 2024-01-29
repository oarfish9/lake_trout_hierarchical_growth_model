#include <TMB.hpp>  // Links in the TMB libraries
# include <vector> // to try the .insert in sim section

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_IVECTOR(ind_idx);  // nrec long with values indicating fish id
  //DATA_IVECTOR(pop_idx);  // nrec long with values indicating pop id
  DATA_IVECTOR(pop_idx_ind);  // nind long with values indicating pop id, created on R side
  DATA_IVECTOR(ind_idx_sim); // for simulating the balanced dataset
  DATA_IVECTOR(pop_idx_ind_sim); // for simulating the balanced dataset
  DATA_INTEGER(npop);
  DATA_INTEGER(npop_sim);
  DATA_INTEGER(nrec);
  DATA_INTEGER(nind);
  DATA_VECTOR(L); // obsv length
  DATA_VECTOR(age);
  DATA_INTEGER(maxage); // for simulation, maximum age you want to simulate to
  DATA_VECTOR(age_sim);
  
  PARAMETER(log_sigma); // for residuals
  PARAMETER(log_sigma_log_Linf_devs); // for individual-level deviations
  PARAMETER(log_sigma_log_K_devs);
  PARAMETER(log_sigma_log_L1_devs);
  PARAMETER(log_sigma_log_Linf_hyper_devs); // for population-level deviations
  PARAMETER(log_sigma_log_K_hyper_devs);
  PARAMETER(log_sigma_log_L1_hyper_devs);
  PARAMETER(log_Linf_hyper); // hyper population means
  PARAMETER(log_K_hyper);
  PARAMETER(log_L1_hyper);
  PARAMETER(theta1); // for among-individuals corr matrix, Linf/K
  PARAMETER(theta2);
  PARAMETER(theta3);
  PARAMETER(theta4);
  PARAMETER(theta5);
  PARAMETER(theta6);
  PARAMETER(log_sigma_PE); // process error variance
  PARAMETER_VECTOR(log_Linf_devs); // nind, for individual-level deviations; random effects
  PARAMETER_VECTOR(log_K_devs);
  PARAMETER_VECTOR(log_L1_devs);
  PARAMETER_VECTOR(log_Linf_hyper_devs); // npop, for population-level deviations; random effects
  PARAMETER_VECTOR(log_K_hyper_devs);
  PARAMETER_VECTOR(log_L1_hyper_devs);
  PARAMETER_VECTOR(log_PE); // process error; random effects

  // back transform parameters
  Type sigma = exp(log_sigma);
  Type sigma_PE = exp(log_sigma_PE);
  vector<Type> PE = exp(log_PE); // process error
  
  vector<Type> Linf_pop(nind); // vector for storing estimated population Linfs from individual level deviations
  vector<Type> K_pop(nind); // vector for storing estimated population Ks from individual level deviations
  vector<Type> L1_pop(nind); // vector for storing estimated population L1s from individual level deviations

  vector<Type> L_hat(nrec); // store predicted growth at age
  vector<Type> L_hat_tmp(nrec); // for storing temporary L_hat for logsitic smoothing function
  vector<Type> resid(nrec);
  vector<Type> log_resid(nrec); 
  vector<Type> delta(nrec); // store increments

  
  // add in individual- and population-level deviations to the hyper-population mean for each
  // vb parameter
  for(int i=0; i<nind; i++){
    Linf_pop(i) = exp(log_Linf_hyper + log_Linf_hyper_devs(pop_idx_ind(i)) + log_Linf_devs(i));
    L1_pop(i) = exp(log_L1_hyper + log_L1_hyper_devs(pop_idx_ind(i)) + log_L1_devs(i));
    K_pop(i) = exp(log_K_hyper + log_K_hyper_devs(pop_idx_ind(i)) + log_K_devs(i));
  }

  Type w2 = 0; // weights for logistic smoothing function
  Type w1 = 1;
  Type alpha = 2; // slope for logistic smoothing function, fixed
        
  // grow the fish
  for(int i=0; i<nrec; i++){
    if(age(i) == 1){
      L_hat(i) = L1_pop(ind_idx(i));
    }else{
      // calculate L_hat_tmp to then use in weighted average
      delta(i-1) = (Linf_pop(ind_idx(i)) - L_hat(i-1)) * (1-exp(-K_pop(ind_idx(i))));
      L_hat_tmp(i) = L_hat(i-1) + (delta(i-1) * PE(i)); 

      // set the weights for logistic smoothing function, based on distance
      // from Linf (for that individual) and the calculated L_hat (for this time step)
      w2 = 1/(1+exp(-alpha*((0.99*Linf_pop(ind_idx(i))) - L_hat_tmp(i)))); // condition
      w1 = 1 - w2;
      // calculate L_hat (real) using weighted average to prevent value from
      // exceeding Linf
      L_hat(i) = (w2 * L_hat_tmp(i)) + (w1 * Linf_pop(ind_idx(i)));
    }
  }

  // Calculate residuals
  resid = L - L_hat;
  log_resid = log(L) - log(L_hat);
  
  // backtransform sd parameters for individual growth parameters
  Type sigma_log_Linf_devs = exp(log_sigma_log_Linf_devs);
  Type sigma_log_K_devs = exp(log_sigma_log_K_devs);
  Type sigma_log_L1_devs = exp(log_sigma_log_L1_devs);

  // backtransform sd parameters for pop-level distribution
  Type sigma_log_Linf_hyper_devs = exp(log_sigma_log_Linf_hyper_devs);
  Type sigma_log_K_hyper_devs = exp(log_sigma_log_K_hyper_devs);
  Type sigma_log_L1_hyper_devs = exp(log_sigma_log_L1_hyper_devs);


  /////////////////////////////////////////////////////////////////
  //// Multivariate Normal Calculations //////
  
  // CALCULATE MULTIVARIATE CORRELATION OF K/Linf/L1 at individual level//
  
  using namespace density;
  
  int n=3; //three parameters
  
  // combine deviations into a matrix for later calculations
  matrix<Type> all_devs(n,nind);
  all_devs.row(0)=log_Linf_devs;
  all_devs.row(1)=log_K_devs;
  all_devs.row(2)=log_L1_devs;
  
  // construct unstructured density object for trivariate normal distribution
  vector<Type> thetas(n); // vector of correlation parameters
  thetas(0) = theta1; // Linf/K
  thetas(1) = theta2; // Linf/L1
  thetas(2) = theta3; // K/L1
  vector<Type> sds(n); // vector of variance parameters
  sds(0) = sigma_log_Linf_devs;
  sds(1) = sigma_log_K_devs;
  sds(2) = sigma_log_L1_devs;

  // Use UNSTRUCTURED_CORR_t to calculate MVN density object
  // as a function of thetas to use in likelihood calculation
  // for L1, K, and Linf
  UNSTRUCTURED_CORR_t<Type> my_mvn(thetas);
  Type nll_mvn = 0;
  
  // use the MVN object to calculate the probability density of each trio of deviations from the mean
  for(int i=0; i<nind; i++){
    nll_mvn += VECSCALE(my_mvn, sds)(all_devs.col(i));
  }

  // BVN calculation for population parameters

  // combine deviations into a matrix for later calculations
  
  matrix<Type> all_pops(n,npop);
  all_pops.row(0) = log_Linf_hyper_devs;
  all_pops.row(1) = log_K_hyper_devs;
  all_pops.row(2) = log_L1_hyper_devs;

  vector<Type> thetas_hyper(n*(n-1)/2);
  thetas_hyper(0) = theta4;
  thetas_hyper(1) = theta5;
  thetas_hyper(2) = theta6;
  vector<Type> sds_hyper(n);
  sds_hyper(0) = sigma_log_Linf_hyper_devs;
  sds_hyper(1) = sigma_log_K_hyper_devs;
  sds_hyper(2) = sigma_log_L1_hyper_devs;

  UNSTRUCTURED_CORR_t<Type> my_mvn_hyper(thetas_hyper);
  Type nll_mvn_hyper = 0;

  for(int i=0; i<npop; i++){
    nll_mvn_hyper += VECSCALE(my_mvn_hyper, sds_hyper)(all_pops.col(i));
  }
  
  // NLL for the observed and expected increments (true length)
  Type nll_length = -sum(dnorm(log(L), (log(L_hat)- pow(sigma,2)/2),
                              sigma, true));
  // NLL for process error
  // back-transformed expected value to be 0 -> avg length, not median length

  Type nll_pe = -sum(dnorm(log(PE),
                           -(pow(sigma_PE,2)/2), sigma_PE, true));

  // Get the full log likelihood
  Type f = nll_length + nll_pe + nll_mvn + nll_mvn_hyper;


    // Simulation block 11-9-2021
  SIMULATE
    {
 

      // set up the simulated data
      int nind_sim = maxage*npop_sim; // 258 individuals
      int nrec_ind_sim = maxage * (maxage + 1)/2;
      int nrec_sim = nrec_ind_sim*npop_sim; // length for 1 individual of each age (1:43) 946 * npop (for each population)

      
      // initialize simulated deviations vectors with differet nind values
      vector<Type> sim_log_Linf_devs(nind_sim);
      vector<Type> sim_log_K_devs(nind_sim);
      vector<Type> sim_log_L1_devs(nind_sim);
      vector<Type> sim_log_Linf_hyper_devs(npop_sim);
      vector<Type> sim_log_K_hyper_devs(npop_sim);
      vector<Type> sim_log_L1_hyper_devs(npop_sim);
      // now simulate random effects
      
 
      // create variance covariance  matrices
      // first, variances for individuals (bivariate)
      matrix<Type> sddiag(n,n); // (3x3)
      sddiag.fill(0);
      sddiag(0,0) = sds(0);
      sddiag(1,1) = sds(1);
      sddiag(2,2) = sds(2);

      //now multiply by correlation matrices extracted from the estimates
      matrix<Type> ind_vcv(n,n);
      ind_vcv = sddiag * my_mvn.cov() * sddiag;

      // generate random deviations for individual vb parameters
      // the MVN will generate 3 observations: the Linf, K, and L1 for an individual
      // store them for each individual and put them in the deviation vectors

      matrix<Type> simdevs_ind(n, nind_sim); // to store the devs for K and Linf
      for(int i=0; i<nind_sim; i++){
        simdevs_ind.col(i) = MVNORM(ind_vcv).simulate();
      }
      
      sim_log_Linf_devs = simdevs_ind.row(0); // log_Linf_dev for the individual
      sim_log_K_devs = simdevs_ind.row(1); // log_K_dev for the individual
      sim_log_L1_devs = simdevs_ind.row(2); // log_L1_dev for the individual
 
      // generate MVN estimates at the population level

      matrix<Type> sddiag_hyper(n,n); // (3x3)
      sddiag_hyper.fill(0);
      sddiag_hyper(0,0) = sds_hyper(0);
      sddiag_hyper(1,1) = sds_hyper(1);
      sddiag_hyper(2,2) = sds_hyper(2);

      matrix<Type> pop_vcv(n,n);
      pop_vcv = sddiag_hyper * my_mvn_hyper.cov() * sddiag_hyper;

      // bivariate for K/Linf at population level
      matrix<Type> hyperdevs(n, npop_sim);
      for(int i=0; i<npop_sim; i++){
        hyperdevs.col(i) = MVNORM(pop_vcv).simulate();
      }

      sim_log_Linf_hyper_devs = hyperdevs.row(0);
      sim_log_K_hyper_devs = hyperdevs.row(1);
      sim_log_L1_hyper_devs = hyperdevs.row(2);

      // simulate process error (last edit 4-20)
      vector<Type> sim_log_PE_mu(nrec_sim);
      sim_log_PE_mu = -pow(sigma_PE,2)/2;
      vector<Type> sim_log_PE(nrec_sim);
      sim_log_PE = rnorm(sim_log_PE_mu, sigma_PE);
      vector<Type> simPE = exp(sim_log_PE);

      // now simulate growth of the fish
      //first, add in deviations
      vector<Type> Linf_pop_sim(nind_sim);
      vector<Type> K_pop_sim(nind_sim);
      vector<Type> L1_pop_sim(nind_sim);


      for(int i=0; i<nind_sim; i++){
        Linf_pop_sim(i) = exp(log_Linf_hyper + sim_log_Linf_hyper_devs(pop_idx_ind_sim(i)) + sim_log_Linf_devs(i));
        K_pop_sim(i) = exp(log_K_hyper + sim_log_K_hyper_devs(pop_idx_ind_sim(i)) + sim_log_K_devs(i));
        L1_pop_sim(i) = exp(log_L1_hyper + sim_log_L1_hyper_devs(pop_idx_ind_sim(i)) + sim_log_L1_devs(i));
      }

      // initialize simulated vectors to store L_hat, delta, and L_hat_tmp for sim
      vector<Type> sim_L_hat(nrec_sim);
      vector<Type> sim_L_hat_tmp(nrec_sim);
      vector<Type> sim_delta(nrec_sim);

      w2 = 0;
      w1 = 1;
      alpha = 2;
      

      for(int i=0; i<nrec_sim; i++){
        if(age_sim(i) == 1){
          sim_L_hat(i) = L1_pop_sim(ind_idx_sim(i));
        }else{
          sim_delta(i-1) = (Linf_pop_sim(ind_idx_sim(i)) - sim_L_hat(i-1)) * (1-exp(-K_pop_sim(ind_idx_sim(i))));
          sim_L_hat_tmp(i) = sim_L_hat(i-1) + (sim_delta(i-1) * simPE(i));
          w2 = 1/(1+exp(-alpha*((0.99*Linf_pop_sim(ind_idx_sim(i))) - sim_L_hat_tmp(i))));
          w1 = 1-w2;
          sim_L_hat(i) = (w2*sim_L_hat_tmp(i)) + (w1*Linf_pop_sim(ind_idx_sim(i)));
        }
      }

      // add in observation error
      vector<Type> log_L_mu(nrec_sim); // 4-20
      log_L_mu = log(sim_L_hat) - (pow(sigma,2)/2); // 4-20
      //vector<Type> log_L_mu = log(sim_L_hat) - pow(sigma,2)/2;
      vector<Type> L_sim(nrec_sim);
      L_sim = exp(rnorm(log_L_mu, sigma));
      //L_sim = rnorm(log_L_mu, sigma);

      REPORT(L_sim);
      REPORT(age_sim);
      REPORT(sim_L_hat);
    }

  REPORT(sds);
  REPORT(sds_hyper);
  REPORT(all_devs);
  REPORT(sigma_log_L1_devs);
  REPORT(sigma_log_L1_hyper_devs);
  REPORT(sigma_log_K_hyper_devs);
  REPORT(sigma_log_Linf_hyper_devs);
  REPORT(my_mvn.cov());
  REPORT(my_mvn_hyper.cov());
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
  REPORT(nll_mvn);
  REPORT(delta);
  REPORT(Linf_pop);
  REPORT(K_pop);
  REPORT(L1_pop);
  REPORT(log_L1_hyper);
  REPORT(log_Linf_hyper);
  REPORT(log_K_hyper);
  REPORT(log_Linf_devs);
  REPORT(log_K_devs);
  REPORT(log_L1_devs);
  REPORT(log_Linf_hyper_devs);
  REPORT(log_K_hyper_devs);
  REPORT(log_L1_hyper_devs);
  REPORT(log_PE);


  
  return f;
  
}
