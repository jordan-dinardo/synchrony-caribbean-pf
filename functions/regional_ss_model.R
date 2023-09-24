SS_trend_ord_multi<-"functions {
  real induced_dirichlet_lpdf(vector c, vector alpha, real phi) {
    int K = num_elements(c) + 1;
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);
    
    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];
    
    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;
    
    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }
    
    return   dirichlet_lpdf(p | alpha)
           + log_determinant(J);
  }
}


data{
  int<lower=1> N;//number of observations (REEF surveys)
  int y[N]; //abundance category for each survey
  int<lower=0> N_site; //number of sites
  int<lower=1,upper=N_site> site[N]; // vector of site identities
  int<lower=0> N_dv; //number of divers
  int<lower=1,upper=N_dv> diver[N]; // vector of diver identities
  int<lower=0> N_dmy; //number of site day clusters
  int<lower=1,upper=N_dmy> dmy[N]; // vector of site day cluster identities
  int N_R; //number of regions
  int<lower=1,upper=N_R> R[N]; //Region id
  int<lower=1> NC_1;  // number of correlations
  int Z; // columns in the covariate matrix
  matrix[N,Z] X; // design matrix X
  int K; //ordinal levels
  int TT; // timespan
  int<lower=0> N_yr; //number of years
  int yr_index[N_yr]; //index of years
  int<lower=1,upper=N_yr> year_id[N]; // vector of year
  int<lower=1,upper=N_R*N_yr> year_id_R[N]; // vector of year-region-id
}
parameters {
  ordered[K-1] c; //cutpoints
  real x0[N_R]; //initial popn size

  //deviations from intercept
  vector[Z] beta; //effort coefficients - RVC
  vector[N_site] a_site; //deviation between sites
  vector[N_dv] a_dv; //deviation between divers
  vector[N_dmy] a_dmy; //deviation between site day clusters
 
  //variance on the deviance components
  real<lower = 0> sd_site;
  real<lower = 0> sd_dv;
  real<lower = 0> sd_dmy;
  vector<lower = 0>[N_R] sd_r;
  vector<lower = 0>[N_R] sd_q;
  cholesky_factor_corr[N_R] Lcorr;
  matrix[TT-1,N_R] z_q;  // standardized process deviations
  matrix[TT,N_R] z_r;// measurement errors
}

transformed parameters{
  matrix[TT,N_R] a_yr; //year
  matrix[TT,N_R] x;  // states by region
  matrix[TT-1,N_R] q_mat;  //process deviations by region
 
  // compute actual group-level effects 
  
  for(t in 1:TT-1){
     q_mat[t,] = (diag_pre_multiply(sd_q, Lcorr) * to_vector(z_q[t,]))';
  }
  
 for(r in 1:N_R){
    x[1,r] = x0[r];
	}
 for(t in 2:TT){
    x[t,] = x[t-1,] + q_mat[t-1,];
      }
 
 for(r in 1:N_R){   
  for(q in 1:N_yr){
    a_yr[q,r] = x[yr_index[q],r] + z_r[q,r]*sd_r[r]; 
  }
  }
}  

model{
  //priors
  c ~ induced_dirichlet(rep_vector(1, K), 0); //prior on ordered cut-points
  beta ~ normal(0,2); //covariates
  
  //correlation matrix
  Lcorr ~ lkj_corr_cholesky(2.0); // prior for cholesky factor; 2.0 instead of 1.0 (OG) might reduce the likelihood of strong correlations

  //variance terms
  sd_q ~ inv_gamma(3,0.5);
  sd_r ~ inv_gamma(3,0.5);
  sd_site ~ inv_gamma(4, 2);
  sd_dv ~ inv_gamma(4, 2);
  sd_dmy ~ inv_gamma(4, 2);
  
  //varying intercepts
  a_site ~ std_normal();
  a_dv ~ std_normal();
  a_dmy ~ std_normal();
  
  for(r in 1:N_R){
      x0[r] ~ normal(0,5);
  }

   to_vector(z_q) ~ std_normal();
    to_vector(z_r) ~ std_normal();
  
  y ~ ordered_logistic(to_vector(a_yr)[year_id_R]+a_site[site]*sd_site+a_dmy[dmy]*sd_dmy+a_dv[diver]*sd_dv+X*beta,c);
}



generated quantities {
  vector<lower=-1,upper=1>[NC_1] cor_1; 
  corr_matrix[N_R] Cor_1 = multiply_lower_tri_self_transpose(Lcorr);
   
   for (k in 1:N_R){
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
}
"
