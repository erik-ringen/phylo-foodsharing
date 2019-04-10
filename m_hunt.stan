data{
// The first part of the data block are integers that denote the length of the data and the number of societies. Other objects in the datablock will be appended with these integers to denote that they are vectors of length N. These integers will also be useful for executing loops later in the Stan program.
int N;
int N_society;

int y[N]; // outcome variable; daily food sharing
int society[N];  // index for society

// Predictors
real hunt_s[N];

// Distance matrices
matrix[N_society,N_society] phy_dist;   // pairwise patristic distance between any two societies
matrix[N_society,N_society] time_dist;   // pairwise temporal distance between any two societies
} // end of data block

parameters{
vector[2] b; // parameter for every predictor +1 for intercept
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
}

model {
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

////////////////// Linear Model ////////////////////
for (i in 1:N) {
y[i] ~ bernoulli( inv_logit( b[1] + phy_v[society[i]] + ep_v[society[i]] + b[2]*hunt_s[i] ) );
}

} // end model block
