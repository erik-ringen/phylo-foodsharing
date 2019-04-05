data{
// The first part of the data block are integers that denote the length of the data and the number of societies. Other objects in the datablock will be appended with these integers to denote that they are vectors of length N. These integers will also be useful for executing loops later in the Stan program.
int N;
int N_society;
int precip_pred_num_missing;
int temp_pred_num_missing;
int npp_pred_num_missing;
int comm_size_num_missing;

int food_share[N]; // same as above, but for the imputation model
int society[N];  // index for society

// Predictors
int strat[N];
int trade_food[N];
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
matrix[N_society,11] phy_z;  // matrix of varying effects for phylogeny, unscaled
vector<lower=0>[11] eta_phy;    // GP parameter for the maximum covariance between any two societies
vector<lower=0>[11] rho_phy;    // GP parameter for the rate at which covariance declines with distance

vector[11] a; // intercepts for imputation model

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
matrix[N_society,11] phy_v; // varying effects for phylogeny, scaled by the covariance matrix

////////////////////// Estimating Covariance Matrices ////////////////////////
for (k in 1:11) {
matrix[N_society, N_society] L_phy;   // cholesky decomposition of the phy covariance matrix
// Phylo Cov
{
  matrix[N_society,N_society] sigma_phy_dist; // the covariance matrix
  for ( i in 1:(N_society-1) )
  for ( j in (i+1):N_society ) {
  sigma_phy_dist[i,j]  = eta_phy[k]*exp(-(rho_phy[k]*phy_dist[i,j])); // Cov function
  sigma_phy_dist[j,i] = sigma_phy_dist[i,j];  // makes the lower tri of the matix = the upper tri, symmetry
  }
  for ( q in 1:N_society )
  sigma_phy_dist[q,q] = eta_phy[k] + 0.0001;     // small numerical constant to diag for computational reasons
  L_phy = cholesky_decompose(sigma_phy_dist); // cholesky decomposition of cov matrix for efficiency
}
phy_v[,k] = L_phy * phy_z[,k];  // Matrix multiplying a vector of uncorrelated, unscaled varying effects by the cholesky decomposition of the covariance matrix returns a vector of correlated, scaled varying effects
}
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
/////////////////////////// Priors /////////////////////////
// Unscaled, uncorrelated random effects
to_vector(phy_z) ~ normal(0,1);

// Priors for GP functions
eta_phy ~ exponential(0.5);
rho_phy ~ exponential(0.5);

a ~ normal(0,2);

////////////////// Linear Model ////////////////////
// Due to the multiple discrete missing variables, this section of the program is long, and perhaps confusing. For a more stripped down example of discrete data imputation in Stan see: https://gist.github.com/rmcelreath/9406643583a8c99304e459e644762f82. The code gets ugly because, while we can treat missing continous values as data in the model block, we cannot do the same for discrete missingness. Instead we have to marginalzie over the uncertainty, and this gets convoluted when there is more than one discrete variable with missing values.

// Modelling all other variables
for (i in 1:N) {
// Because observation level random effects are highly correlation with the variance term of Gaussian distributions, we fix the resiudal sd to 0.5
food_share[i] ~ bernoulli_logit(a[1] + phy_v[i,1]); // the outcome again, repeated here just for the imputation model
hunt_s[i] ~ normal(a[2] + phy_v[i,2], 0.5);
trade_food[i] ~ bernoulli_logit(a[4] + phy_v[i,4]);
strat[i] ~ bernoulli_logit(a[5] + phy_v[i,5]);
hus_s[i] ~ normal(a[6] + phy_v[i,6], 0.5);
precip_pred_merged[i] ~ normal(a[7] + phy_v[i,7], 0.5);
temp_pred_merged[i] ~ normal(a[8] + phy_v[i,8], 0.5);
npp_pred_merged[i] ~ normal(a[9] + phy_v[i,9], 0.5);
comm_size_merged[i] ~ normal(a[11] +  phy_v[i,11], 0.5);

// Food Storage
if (food_store_missing[i]==1) {
  vector[2] lp_store;
  
 lp_store[1] = log(inv_logit(a[3] + phy_v[i,3])) + bernoulli_logit_lpmf( 1 | a[3] + phy_v[i,3]);
 lp_store[2] = log(1 - inv_logit(a[3] + phy_v[i,3])) + bernoulli_logit_lpmf( 0 | a[3] + phy_v[i,3]);
 target += log_sum_exp(lp_store);
}

else{
  food_store[i] ~ bernoulli_logit(a[3] + phy_v[i,3]);
}

// Labor Sharing
if (labor_share_missing[i]==1) {
  vector[2] lp_labor;
  
 lp_labor[1] = log(inv_logit(a[10] + phy_v[i,10])) + bernoulli_logit_lpmf( 1 | a[10] + phy_v[i,10]);
 lp_labor[2] = log(1 - inv_logit(a[10] + phy_v[i,10])) + bernoulli_logit_lpmf( 0 | a[10] + phy_v[i,10]);
 target += log_sum_exp(lp_labor);
}

else{
  labor_share[i] ~ bernoulli_logit(a[10] + phy_v[i,10]);
}

}

} // end model block

