#include <TMB.hpp>  // Links in the TMB libraries
# include <vector> // added 8-4-2022

template<class Type>
Type objective_function<Type>::operator() ()
{
 
  
  DATA_IVECTOR(ind_idx);
  DATA_IVECTOR(pop_idx_ind);
  DATA_IVECTOR(ind_idx_sim);
  DATA_IVECTOR(pop_idx_ind_sim);
  DATA_INTEGER(npop);
  DATA_INTEGER(npop_sim);
  DATA_INTEGER(nrec);
  DATA_INTEGER(nind);
  DATA_VECTOR(L);
  DATA_VECTOR(age);
  DATA_INTEGER(maxage);
  DATA_VECTOR(age_sim);
  
  PARAMETER(log_sigma);
  PARAMETER(log_sigma_log_g_devs);
  PARAMETER(log_sigma_log_h_devs); 
  PARAMETER(log_sigma_log_L1_devs);
  PARAMETER(log_g_hyper); // hyper-population mean reproductive investment
  PARAMETER(log_h_hyper); // hyper-population mean rate of energy aquisition
  PARAMETER(log_L1_hyper); // hyper-population mean length for age == 1
  PARAMETER_VECTOR(log_g_devs); // for MVN among-individuals
  PARAMETER_VECTOR(log_h_devs); // for MVN among-individuals
  PARAMETER_VECTOR(log_L1_devs); // for MVN
  PARAMETER(alpha); // intercept for calculation of age at maturity (T)
  PARAMETER(log_beta); // slope for calculation of age at maturity (T)
  PARAMETER(theta1); // for corr matrix
  PARAMETER(theta2);
  PARAMETER(theta3);
  PARAMETER_VECTOR(log_PE); // process error
  PARAMETER(log_sigma_PE);

  // back transform parameters
  Type beta = exp(log_beta);
  Type sigma = exp(log_sigma);
  Type sigma_PE = exp(log_sigma_PE);
  vector<Type> PE = exp(log_PE); // process error

  vector<Type> g(nind); // different value for g for each individual
  vector<Type> h(nind); // different value for h for each individual
  vector<Type> L1(nind); // different value for L1 for each individual
  vector<Type> Linf(nind); // different value for Linf for each individual
  vector<Type> T(nind); // different value for T for each individual, calculated later
  vector<Type> K(nind);// different value for K for each individual, calculated later
  vector<Type> L_hat(nrec); // store predicted growth at age
  vector<Type> resid(nrec);
  vector<Type> log_resid(nrec);
  vector<Type> lin(nrec); // linear increment calculation
  vector<Type> vB(nrec); // vB increment calculation


  // calculate parameters

  // Add in individual- and population-level deviations to the hyper-population means
  for(int i=0; i<nind; i++){
    g(i) = exp(log_g_hyper + log_g_devs(i));
    h(i) = exp(log_h_hyper +  log_h_devs(i));
    L1(i) = exp(log_L1_hyper + log_L1_devs(i));
  }
  // calculate Linf, K, and T from parameters
  Linf = (3*h)/g;
  K = log(1+(g/3));
  T = alpha + beta*h; // calculating without error term to begin with
  //T = beta*h;

  Type w2 = 0;
  Type w1 = 1;

  // grow the fish
  for(int i=0; i<nrec; i++){
    if(age(i) == 1){
      L_hat(i) = L1(ind_idx(i)); // set length at age==1 for all fish to L1
      //delta(i) = -99; // no increment for age 0 to 1
    }else{
      // weights
      w2 = 1/(1+exp(-3*(age(i) - T(ind_idx(i)))));
      w1 = 1 - w2;
                   
      // while age(i) <<  T, calculate linear increment
      //lin(i-1) = h(ind_idx(i))*(age(i) - t1); // linear growth
      lin(i-1) = h(ind_idx(i));
      // while age(i) >> T, calculate vB increment
      vB(i-1) = (Linf(ind_idx(i)) - L_hat(i-1)) * (1-exp(-K(ind_idx(i)))); // vB growth
      // add increment to previous length, add process error, and store as predicted length
      L_hat(i) = L_hat(i-1) + (((lin(i-1) * w1) + (vB(i-1) * w2)) * PE(i)); // for every last age of a fish, no delta and should leave as 0
    }
  }

  // Calculate residuals
  resid = L - L_hat;
  log_resid = log(L) - log(L_hat);
  
  // backtransform sd parameters for individual growth parameters
  Type sigma_log_g_devs = exp(log_sigma_log_g_devs);
  Type sigma_log_h_devs = exp(log_sigma_log_h_devs);
  Type sigma_log_L1_devs = exp(log_sigma_log_L1_devs);

  /////////////////////////////////////////////////////////////////
  // CALCULATE CORRELATION OF g,h,L1 //
  
  using namespace density;
  
  int n=3; //three parameters
  
  // combine deviations into a matrix for later calculations
  matrix<Type> all_devs(n,nind);
  all_devs.row(0)=log_g_devs;
  all_devs.row(1)=log_h_devs;
  all_devs.row(2)=log_L1_devs;
  
  // construct unstructured density object for trivariate normal distribution
  vector<Type> thetas(n); // vector of parameters
  thetas(0) = theta1; // g/h
  thetas(1) = theta2; // g/L1
  thetas(2) = theta3; // h/L1
  vector<Type> sds(n);
  sds(0) = sigma_log_g_devs;
  sds(1) = sigma_log_h_devs;
  sds(2) = sigma_log_L1_devs;

  // Use UNSTRUCTURED_CORR_t to calculate MVN density object
  // as a function of thetas to use in likelihood calculation
  // for L1, g, and h
  UNSTRUCTURED_CORR_t<Type> my_mvn(thetas);
  // use the MVN object to calculate the probability density of each trio of deviations from the mean
  Type nll_mvn = 0;
  
  // use the MVN object to calculate the probability density of each trio of deviations from the mean
  for(int i=0; i<nind; i++){
    nll_mvn += VECSCALE(my_mvn, sds)(all_devs.col(i));
  }
  
  // NLL for the observed and expected increments
  Type nll_length = -sum(dnorm(log(L), (log(L_hat)-(pow(sigma,2)/2)),
                         sigma, true));
  // NLL FOR PROCESS ERROR
  Type nll_pe = -sum(dnorm(log(PE),
                           -(pow(sigma_PE,2)/2), sigma_PE, true));

  // Get the full log likelihood
  Type f = nll_length + nll_mvn + nll_pe;

  SIMULATE
    {
      int nind_sim = maxage*npop_sim;
      int nrec_ind_sim = maxage * (maxage + 1)/2;
      int nrec_sim = nrec_ind_sim * npop_sim;

      // initialize simulated deviation vectors
      vector<Type> sim_log_g_devs(nind_sim);
      vector<Type> sim_log_h_devs(nind_sim);
      vector<Type> sim_log_L1_devs(nind_sim);

      // now simulate random effects
      matrix<Type> sddiag(n,n); //(3x3)
      sddiag.fill(0);
      sddiag(0,0) = sds(0);
      sddiag(1,1) = sds(1);
      sddiag(2,2) = sds(2);

      // extract variance/covariance matrix
      matrix<Type> ind_vcv(n,n);
      ind_vcv = sddiag * my_mvn.cov() * sddiag;

      // now multiply by vcov matrix
      matrix<Type> simdevs_ind(n, nind_sim);
      for(int i=0; i<nind_sim; i++){
        simdevs_ind.col(i) = MVNORM(ind_vcv).simulate();
      }

      // fill in deviations
      sim_log_g_devs = simdevs_ind.row(0);
      sim_log_h_devs = simdevs_ind.row(1);
      sim_log_L1_devs = simdevs_ind.row(2);

      // simulate process error
      vector<Type> sim_log_PE_mu(nrec_sim); // 4-20-2022
      sim_log_PE_mu = -pow(sigma_PE,2)/2; // 4-20-2022
      vector<Type> sim_log_PE(nrec_sim); // earlier
      sim_log_PE = rnorm(sim_log_PE_mu, sigma_PE);
      //sim_log_PE = rnorm(-pow(sigma_PE,2)/2, sigma_PE); // earlier
      vector<Type> simPE = exp(sim_log_PE); // earlier

      // simulate growth of fish
      vector<Type> g_sim(nind_sim);
      vector<Type> h_sim(nind_sim);
      vector<Type> L1_sim(nind_sim);

      for(int i=0; i<nind_sim; i++){
        g_sim(i) = exp(log_g_hyper + sim_log_g_devs(i));
        h_sim(i) = exp(log_h_hyper + sim_log_h_devs(i));
        L1_sim(i) = exp(log_L1_hyper + sim_log_L1_devs(i));
      }

      // initialize simulated vectors to store L_hat, delta, and L_hat_tmp
      vector<Type> sim_L_hat(nrec_sim);
      vector<Type> sim_lin(nrec_sim);
      vector<Type> sim_vB(nrec_sim);

      // grow the fish
        // calculate Linf, K, and T from parameters
      vector<Type> Linf_sim = (3*h_sim)/g_sim;
      vector<Type> K_sim = log(1+(g_sim/3));
      vector<Type> T_sim = alpha + beta*h_sim;

      w2 = 0;
      w1 = 1;

      // grow the fish
      for(int i=0; i<nrec_sim; i++){
        if(age_sim(i) == 1){
          sim_L_hat(i) = L1_sim(ind_idx_sim(i)); // set length at age==1 for all fish to L1
        }else{
          // weights
          w2 = 1/(1+exp(-3*(age_sim(i) - T_sim(ind_idx_sim(i)))));
          w1 = 1 - w2;    
          // while age(i) <<  T, calculate linear increment
          sim_lin(i-1) = h_sim(ind_idx_sim(i));
          // while age(i) >> T, calculate vB increment
          sim_vB(i-1) = (Linf_sim(ind_idx_sim(i)) - sim_L_hat(i-1)) * (1-exp(-K_sim(ind_idx_sim(i)))); // vB growth
          // add increment to previous length, add process error, and store as predicted length
          sim_L_hat(i) = sim_L_hat(i-1) + (((sim_lin(i-1) * w1) + (sim_vB(i-1) * w2)) * simPE(i)); // for every last age of a fish, no delta and should leave as 0
        }
      }

      // add in obsv error
      vector<Type> log_L_mu(nrec_sim); // 4-20
      log_L_mu = log(sim_L_hat) - (pow(sigma,2)/2); // 4-20
      //vector<Type> log_L_mu = log(sim_L_hat) - pow(sigma,2)/2;
      vector<Type> L_sim(nrec_sim);
      L_sim = exp(rnorm(log_L_mu, sigma)); // 4-20
      //L_sim = rnorm(log_L_mu, sigma);

      REPORT(simPE);
      REPORT(sim_log_g_devs);
      REPORT(sim_log_h_devs);
      REPORT(sim_log_L1_devs);
      REPORT(L_sim);
      REPORT(age_sim);
      REPORT(sim_L_hat);
    }

  Type g_hyper = exp(log_g_hyper);
  Type h_hyper = exp(log_h_hyper);
  Type L1_hyper = exp(log_L1_hyper);
 
  REPORT(thetas);
  REPORT(f);
  REPORT(sds);
  REPORT(all_devs);
  REPORT(my_mvn.cov());
  REPORT(PE);
  REPORT(resid);
  REPORT(age);
  REPORT(g_hyper);
  REPORT(h_hyper);
  REPORT(L1_hyper);
  REPORT(Linf);
  REPORT(vB);
  REPORT(lin);
  REPORT(K);
  REPORT(T);
  REPORT(L);
  REPORT(L_hat);
  REPORT(g);
  REPORT(h);
  REPORT(L1);
  REPORT(alpha);
  REPORT(beta);
  //REPORT(w1);
  //REPORT(w2);
  REPORT(sigma);
  REPORT(sigma_PE);

  return f;
  
}
