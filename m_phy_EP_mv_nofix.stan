data{
// The first part of the data block are integers that denote the length of the data and the number of societies. Other objects in the datablock will be appended with these integers to denote that they are vectors of length N. These integers will also be useful for executing loops later in the Stan program.
int N;
int N_society;

int y[N]; // outcome variable; daily food sharing
int food_share[N]; // same as above, but for the imputation model
int society[N];  // index for society

// Distance matrices
matrix[N_society,N_society] phy_dist;   // pairwise patristic distance between any two societies
matrix[N_society,N_society] time_dist;   // pairwise temporal distance between any two societies
} // end of data block

parameters{
vector[1] b; // intercept
vector[N_society] phy_z;  // matrix of varying effects for phylogeny, unscaled
vector[N_society] ep_z;    // matrix of varying effects for EP, unscaled
real<lower=0> eta_phy;    // GP parameter for the maximum covariance between any two societies
real<lower=0> eta_ep;      // GP parameter for the maximum covariance between any two societies
real<lower=0> rho_phy;    // GP parameter for the rate at which covariance declines with distance
real<lower=0> rho_ep;      // GP parameter for the maximum covariance between any two societies

} // end parameters block

transformed parameters{
// The parameters are in this block are 'transformed' in that they are not directly estimated, but rather are functions of the data and the parameters declared in the previous block.  

// Transformed parameters for GP functions
matrix[N_society, N_society] L_phy;   // cholesky decomposition of the phy covariance matrix
matrix[N_society, N_society] L_time; // cholesky decomposition of the EP covariance matrix
vector[N_society] phy_v; // varying effects for phylogeny, scaled by the covariance matrix
vector[N_society] ep_v;  // varing effects for EP, scaled by the covariance matrix

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

for (i in 1:N) {
p[i] = b[1] + phy_v[society[i]] + ep_v[society[i]];
y[i] ~ bernoulli_logit(p[i]);
}

} // end model block

generated quantities {
  vector[N] p; 
// Here we will calculate the log likelihood. Discrete missing values are also problematic in this case, but to get around this problem we will use Stan' built-in random number generator. Each iteration of the sampler, discrete values will be generated based on the probability of food storage/labor sharing being present. This procedure propogates uncertainty in estimates and allows us to calculate log lik for model comparison.
          
          vector[N] log_lik;
          
          for (i in 1:N) {
            p[i] = b[1] + phy_v[society[i]] + ep_v[society[i]];
            log_lik[i] = bernoulli_lpmf( y[i] | inv_logit(p[i]) );
          }
          
} // end generated quantities
