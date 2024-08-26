//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//


data{

int<lower = 0> n; // number of observations

vector[n] Ei; // vector of escapement

vector[n] Ri; // vector of recruits

real max_r;  // max observed recruitment
}


parameters {
  
  real<lower=0> sigma_y; // standard deviation of recruitment 
  
  real<lower=1, upper=10> r; // growth rate
  
  real<lower = 0.1 * max_r, upper =10 * max_r> k; // carrying capacity
}

model {
  
  sigma_y ~ cauchy(0, 1);
  
  r ~ uniform(1,10);
  
  k ~ uniform(0.1 * max_r,10 * max_r);
  
  for (i in 1:n){
    Ri[i] ~ lognormal(log(Ei[i]) + log(r) - log(1 + (r - 1) * Ei[i] / k), sigma_y); 
  }
}


generated quantities {
  vector[n] log_lik;
  for (i in 1:n){
    log_lik[i] = lognormal_lpdf(Ri[i] | log(Ei[i]) + log(r) - log(1 + (r - 1) * Ei[i] / k), sigma_y);
  }
}





