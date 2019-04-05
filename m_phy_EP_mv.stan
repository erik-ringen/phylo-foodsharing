data{
// The first part of the data block are integers that denote the length of the data and the number of societies. Other objects in the datablock will be appended with these integers to denote that they are vectors of length N. These integers will also be useful for executing loops later in the Stan program.
int N;
int N_society;
int precip_pred_num_missing;
int temp_pred_num_missing;
int npp_pred_num_missing;
int comm_size_num_missing;

int y[N]; // outcome variable; daily food sharing
int food_share[N]; // same as above, but for the imputation model
int society[N];  // index for society

// Predictors
real strat[N];
real trade_food[N];
real hunt_s[N];
real hus_s[N];
real precip_pred_s[N];
real temp_pred_s[N];
real npp_pred_s[N];
real comm_size_s[N];
int food_store[N];
int labor_share[N];

// Indices for missing values
int food_store_missing[N];
int labor_share_missing[N];
int precip_pred_missing[N];
int temp_pred_missing[N];
int npp_pred_missing[N];
int comm_size_missing[N];

// Distance matrices
matrix[N_society,N_society] phy_dist;   // pairwise patristic distance between any two societies
matrix[N_society,N_society] time_dist;   // pairwise temporal distance between any two societies
} // end of data block

parameters{
vector[11] b; // parameter for every predictor +1 for intercept
vector[N_society] phy_z;  // matrix of varying effects for phylogeny, unscaled
vector[N_society] ep_z;    // matrix of varying effects for EP, unscaled
real<lower=0> eta_phy;    // GP parameter for the maximum covariance between any two societies
real<lower=0> eta_ep;      // GP parameter for the maximum covariance between any two societies
real<lower=0> rho_phy;    // GP parameter for the rate at which covariance declines with distance
real<lower=0> rho_ep;      // GP parameter for the maximum covariance between any two societies

vector[11] a; // intercepts for imputation model
matrix[11,N] res_z; // observation-level varying effects, unscaled and uncorrelated. We will use these to predict missing values
cholesky_factor_corr[11] L_Rho; // cholesky decomposition of the correlation matrix for missing data model
vector<lower=0>[11] sigma_res; // scale parameters for observation level varying effects

// Imputation parameters
// Continous variables
vector[precip_pred_num_missing] precip_pred_impute; // vector of missing values to be imputed
vector[temp_pred_num_missing] temp_pred_impute; 
vector[npp_pred_num_missing] npp_pred_impute;
vector[comm_size_num_missing] comm_size_impute;

} // end parameters block

transformed parameters{
// The parameters are in this block are 'transformed' in that they are not directly estimated, but rather are functions of the data and the parameters declared in the previous block.  

// Vector(s) of mixtures of observed and imputed observations for continous misssing
vector[N] precip_pred_merged;
vector[N] temp_pred_merged;
vector[N] npp_pred_merged;
vector[N] comm_size_merged;

// Transformed parameters for GP functions
matrix[N_society, N_society] L_phy;   // cholesky decomposition of the phy covariance matrix
matrix[N_society, N_society] L_time; // cholesky decomposition of the EP covariance matrix
vector[N_society] phy_v; // varying effects for phylogeny, scaled by the covariance matrix
vector[N_society] ep_v;  // varing effects for EP, scaled by the covariance matrix

// Observation-level, correlated varying effects
matrix[N,11] res_v;

res_v = (diag_pre_multiply(sigma_res, L_Rho) * res_z)';
////////////////////// Estimating Covariance Matrices ////////////////////////
// Phylo Cov
{
  matrix[N_society,N_society] sigma_phy_dist; // the covariance matrix
  for ( i in 1:(N_society-1) )
  for ( j in (i+1):N_society ) {
  sigma_phy_dist[i,j]  = eta_phy*exp(-(rho_phy*phy_dist[i,j])); // Cov function
  sigma_phy_dist[j,i] = sigma_phy_dist[i,j];  // makes the lower tri of the matix = the upper tri, symmetry
  }
  for ( q in 1:N_society )
  sigma_phy_dist[q,q] = eta_phy + 0.0001;     // small numerical constant to diag for computational reasons
  L_phy = cholesky_decompose(sigma_phy_dist); // cholesky decomposition of cov matrix for efficiency
}
phy_v = L_phy * phy_z;  // Matrix multiplying a vector of uncorrelated, unscaled varying effects by the cholesky decomposition of the covariance matrix returns a vector of correlated, scaled varying effects

// Time Cov, same code as above
{
  matrix[N_society,N_society] sigma_time;
  for ( i in 1:(N_society-1) )
  for ( j in (i+1):N_society ) {
  sigma_time[i,j]  = eta_ep*exp(-(rho_ep*time_dist[i,j]));
  sigma_time[j,i] = sigma_time[i,j];
  }
  for ( q in 1:N_society )
  sigma_time[q,q] = eta_ep + 0.0001;
  L_time = cholesky_decompose(sigma_time);
}
ep_v = L_time * ep_z;
///////////////////////////////////////////////////////////////////////
// Merge missing and observed for continous missing
for ( i in 1:N ) {
precip_pred_merged[i] = precip_pred_s[i];
if ( precip_pred_missing[i] > 0 ) precip_pred_merged[i] = precip_pred_impute[precip_pred_missing[i]];
}

for ( i in 1:N ) {
temp_pred_merged[i] = temp_pred_s[i];
if ( temp_pred_missing[i] > 0 ) temp_pred_merged[i] = temp_pred_impute[temp_pred_missing[i]];
}

for ( i in 1:N ) {
npp_pred_merged[i] = npp_pred_s[i];
if ( npp_pred_missing[i] > 0 ) npp_pred_merged[i] = npp_pred_impute[npp_pred_missing[i]];
}

for ( i in 1:N ) {
comm_size_merged[i] = comm_size_s[i];
if ( comm_size_missing[i] > 0 ) comm_size_merged[i] = comm_size_impute[comm_size_missing[i]];
}
} // end transformed parameters block

model {
vector[N] p; // containter for linear model, prob of daily food sharing

/////////////////////////// Priors /////////////////////////
b ~ normal(0,2);

// Unscaled, uncorrelated random effects
phy_z ~ normal(0,1);
ep_z ~ normal(0,1);

// Priors for GP functions
eta_phy ~ exponential(0.5);
eta_ep ~ exponential(0.5);
rho_phy ~ exponential(0.5);
rho_ep ~ exponential(0.5);

// Parameters missing data imputation
a ~ normal(0,2);
sigma_res ~ exponential(0.5);
L_Rho ~ lkj_corr_cholesky(2);
to_vector(res_z) ~ normal(0,1);

////////////////// Linear Model ////////////////////
// Due to the multiple discrete missing variables, this section of the program is long, and perhaps confusing. For a more stripped down example of discrete data imputation in Stan see: https://gist.github.com/rmcelreath/9406643583a8c99304e459e644762f82. The code gets ugly because, while we can treat missing continous values as data in the model block, we cannot do the same for discrete missingness. Instead we have to marginalzie over the uncertainty, and this gets convoluted when there is more than one discrete variable with missing values.

for ( i in 1:N ) {
// Food storage missing, labor share not missing
if ( food_store_missing[i]==1 && labor_share_missing[i]==0 ) {
 
 vector[2] lp; // log probabilities to marginalize over
 real pr_store = inv_logit(a[3] + res_v[i,3]); // probability food storage present
 
 // Model where storage is present, weighted by the probability of food storage
 lp[1] = log(pr_store) + bernoulli_logit_lpmf(y[i] | b[1] + phy_v[society[i]] + ep_v[society[i]] + b[2]*hunt_s[i] + b[3]*1 + b[4]*trade_food[i] + b[5]*strat[i] + b[6]*hus_s[i] + b[7]*precip_pred_merged[i] + b[8]*temp_pred_merged[i] + b[9]*npp_pred_merged[i] + b[10]*labor_share[i] + b[11]*comm_size_merged[i]);
 
  // Model where storage is absent, weightted by the probability of no food storage
 lp[2] = log(1-pr_store) + bernoulli_logit_lpmf(y[i] | b[1] + phy_v[society[i]] + ep_v[society[i]] + b[2]*hunt_s[i] + b[3]*0 + b[4]*trade_food[i] + b[5]*strat[i] + b[6]*hus_s[i] + b[7]*precip_pred_merged[i] + b[8]*temp_pred_merged[i] + b[9]*npp_pred_merged[i] + b[10]*labor_share[i] + b[11]*comm_size_merged[i]);
 
 // Marginalize over each possibility
 target += log_sum_exp(lp);
}

// Labor missing, food storage observed
else if ( food_store_missing[i]==0 && labor_share_missing[i]==1 ) {
 
 vector[2] lp; // log probabilities to marginalize over
 real pr_labor = inv_logit(a[10] + res_v[i,10]); // probability labor sharing present
 
 // Model where labor is present
 lp[1] = log(pr_labor) + bernoulli_logit_lpmf(y[i] | b[1] + phy_v[society[i]] + ep_v[society[i]] + b[2]*hunt_s[i] + b[3]*food_store[i] + b[4]*trade_food[i] + b[5]*strat[i] + b[6]*hus_s[i] + b[7]*precip_pred_merged[i] + b[8]*temp_pred_merged[i] + b[9]*npp_pred_merged[i] + b[10]*1 + b[11]*comm_size_merged[i]);
 
  // Model where labor is absent
 lp[2] = log(1-pr_labor) + bernoulli_logit_lpmf(y[i] | b[1] + phy_v[society[i]] + ep_v[society[i]] + b[2]*hunt_s[i] + b[3]*food_store[i] + b[4]*trade_food[i] + b[5]*strat[i] + b[6]*hus_s[i] + b[7]*precip_pred_merged[i] + b[8]*temp_pred_merged[i] + b[9]*npp_pred_merged[i] + b[10]*0 + b[11]*comm_size_merged[i]);
 
 // Marginalize over each possibility
 target += log_sum_exp(lp);
}

// Both missing
    else if ( food_store_missing[i]==1 && labor_share_missing[i]==1 ) {
      vector[4] lp; // log-lik containers to marginalize over
      
        // Pr (Labor Present, Storage Present) = Pr (Labor Present) * Pr (Storage Present)
        // Pr (Labor Present, Storage Absent) = Pr (Labor Present) * (1 - Pr(Storage Present))
        // Pr (Labor Absent, Storagee Present) = (1 - Pr(Labor Present)) * Pr(Storage Present)
        // Pr (Labor Absent, Storage Absent) = (1 - Pr(Labor Present)) * (1 - Pr(Storage Present))
          
          lp[1] = log( inv_logit(a[10] + res_v[i,10]) * inv_logit(a[3] + res_v[i,3]) ) + bernoulli_logit_lpmf(y[i] | b[1] + phy_v[society[i]] + ep_v[society[i]] + b[2]*hunt_s[i] + b[3]*1 + b[4]*trade_food[i] + b[5]*strat[i] + b[6]*hus_s[i] + b[7]*precip_pred_merged[i] + b[8]*temp_pred_merged[i] + b[9]*npp_pred_merged[i] + b[10]*1 + b[11]*comm_size_merged[i]);
          
          lp[2] = log( inv_logit(a[10] + res_v[i,10]) * (1-inv_logit(a[3] + res_v[i,3])) ) + bernoulli_logit_lpmf(y[i] | b[1] + phy_v[society[i]] + ep_v[society[i]] + b[2]*hunt_s[i] + b[3]*0 + b[4]*trade_food[i] + b[5]*strat[i] + b[6]*hus_s[i] + b[7]*precip_pred_merged[i] + b[8]*temp_pred_merged[i] + b[9]*npp_pred_merged[i] + b[10]*1 + b[11]*comm_size_merged[i]);
          
          lp[3] = log( (1-inv_logit(a[10] + res_v[i,10])) * inv_logit(a[3] + res_v[i,3]) ) + bernoulli_logit_lpmf(y[i] | b[1] + phy_v[society[i]] + ep_v[society[i]] + b[2]*hunt_s[i] + b[3]*1 + b[4]*trade_food[i] + b[5]*strat[i] + b[6]*hus_s[i] + b[7]*precip_pred_merged[i] + b[8]*temp_pred_merged[i] + b[9]*npp_pred_merged[i] + b[10]*0 + b[11]*comm_size_merged[i]);
          
          lp[4] = log( (1-inv_logit(a[10] + res_v[i,10])) * (1-inv_logit(a[3] + res_v[i,3]))) + bernoulli_logit_lpmf(y[i] | b[1] + phy_v[society[i]] + ep_v[society[i]] + b[2]*hunt_s[i] + b[3]*0 + b[4]*trade_food[i] + b[5]*strat[i] + b[6]*hus_s[i] + b[7]*precip_pred_merged[i] + b[8]*temp_pred_merged[i] + b[9]*npp_pred_merged[i] + b[10]*0 + b[11]*comm_size_merged[i]);
          
          target += log_sum_exp(lp);
    }
          
          // Neither are missing
          else if (food_store_missing[i]==0 && labor_share_missing[i]==0) {
            
p[i] = b[1] + phy_v[society[i]] + ep_v[society[i]] + b[2]*hunt_s[i] + b[3]*food_store[i] + b[4]*trade_food[i] + b[5]*strat[i] + b[6]*hus_s[i] + b[7]*precip_pred_merged[i] + b[8]*temp_pred_merged[i] + b[9]*npp_pred_merged[i] + b[10]*labor_share[i] + b[11]*comm_size_merged[i];

y[i] ~ bernoulli_logit(p[i]);
}
}

// Modelling all other variables
for (i in 1:N) {
// Because observation level random effects are highly correlation with the variance term of Gaussian distributions, we fix the resiudal sd to 0.5
food_share[i] ~ bernoulli_logit(a[1] + res_v[i,1]); // the outcome again, repeated here just for the imputation model
hunt_s[i] ~ normal(a[2] + res_v[i,2], 0.5);
trade_food[i] ~ normal(a[4] + res_v[i,4], 0.5);
strat[i] ~ normal(a[5] + res_v[i,5], 0.5);
hus_s[i] ~ normal(a[6] + res_v[i,6], 0.5);
precip_pred_merged[i] ~ normal(a[7] + res_v[i,7], 0.5);
temp_pred_merged[i] ~ normal(a[8] + res_v[i,8], 0.5);
npp_pred_merged[i] ~ normal(a[9] + res_v[i,9], 0.5);
comm_size_merged[i] ~ normal(a[11] + res_v[i,11], 0.5);

// Food Storage
if (food_store_missing[i]==1) {
  vector[2] lp_store;
  
 lp_store[1] = log(inv_logit(a[3] + res_v[i,3])) + bernoulli_logit_lpmf( 1 | a[3] + res_v[i,3] );
 lp_store[2] = log(1 - inv_logit(a[3] + res_v[i,3])) + bernoulli_logit_lpmf( 0 | a[3] + res_v[i,3] );
 target += log_sum_exp(lp_store);
}

else{
  food_store[i] ~ bernoulli_logit(a[3] + res_v[i,3]);
}

// Labor Sharing
if (labor_share_missing[i]==1) {
  vector[2] lp_labor;
  
 lp_labor[1] = log(inv_logit(a[10] + res_v[i,10])) + bernoulli_logit_lpmf( 1 | a[10] + res_v[i,10] );
 lp_labor[2] = log(1 - inv_logit(a[10] + res_v[i,10])) + bernoulli_logit_lpmf( 0 | a[10] + res_v[i,10] );
 target += log_sum_exp(lp_labor);
}

else{
  labor_share[i] ~ bernoulli_logit(a[10] + res_v[i,10]);
}

}

} // end model block

          generated quantities {
          // Here we will calculate the log likelihood. Discrete missing values are also problematic in this case, but to get around this problem we will use Stan' built-in random number generator. Each iteration of the sampler, discrete values will be generated based on the probability of food storage/labor sharing being present. This procedure propogates uncertainty in estimates and allows us to calculate log lik for model comparison.
          
          vector[N] store_gen;
          vector[N] labor_gen;
          vector[N] p;
          vector[N] log_lik;
          vector[N] y_hat;      // expected values for each observation; fixed effects only
          
          
          for (i in 1:N) {
            labor_gen[i] = labor_share[i];
            store_gen[i] = food_store[i];
           
          if (labor_share_missing[i]==1)
          labor_gen[i] = bernoulli_rng( inv_logit(a[10] + res_v[i,10]) );
          
          if (food_store_missing[i]==1)
          store_gen[i] = bernoulli_rng( inv_logit(a[3] + res_v[i,3]) );
          
            p[i] = b[1] + phy_v[society[i]] + ep_v[society[i]] + b[2]*hunt_s[i] + b[3]*store_gen[i] + b[4]*trade_food[i] + b[5]*strat[i] + b[6]*hus_s[i] + b[7]*precip_pred_merged[i] + b[8]*temp_pred_merged[i] + b[9]*npp_pred_merged[i] + b[10]*labor_gen[i] + b[11]*comm_size_merged[i];
            
            log_lik[i] = bernoulli_lpmf( y[i] | inv_logit(p[i]) );
            
            // Fixed effects only
            y_hat[i] = b[1] + b[2]*hunt_s[i] + b[3]*store_gen[i] + b[4]*trade_food[i] + b[5]*strat[i] + b[6]*hus_s[i] + b[7]*precip_pred_merged[i] + b[8]*temp_pred_merged[i] + b[9]*npp_pred_merged[i] + b[10]*labor_gen[i] + b[11]*comm_size_merged[i];
          }
          
          } // end generated quantities
