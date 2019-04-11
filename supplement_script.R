# Supplemental Script for: The evolution of daily food sharing: A Bayesian Phylogenetic Analysis
# Erik J. Ringen, Pavel Duda, and Adrian V. Jaeggi

# This script relies on the following packages, all of which can be downloaded directly from cran except for 'rethinking'
# 1) rethinking (see instructions for installation at http://xcelab.net/rm/software/)
# 2) rstan (a dependency of rethinking)
# 2) phytools (used for importing the SCCS phylogeny)
# 3) tidyr (used for data wrangling)
# 4) geosphere (used for converting lat/lon into pairwise geographic distances)

# In order to replicate the figures from the main text, the following packages are also needed:
# 5) ggridges (stacked density plots, e.g., Fig 2A)
# 6) wesanderson (nice colors based on Wes Anderson films, just cosmetic)
# 7) mapdata (for reproducing the map in Fig. 1)
################################################################################################
################################################################################################
# Once all packages (and their dependencies) are installed, begin:
library(rethinking)
library(rstan)
library(phytools)
library(tidyr)
library(geosphere)
library(ggridges)
library(wesanderson)
library(mapdata)

setwd("")
d <- read.csv("supp_script_data.csv") # load the data
str(d) # structure of dataframe
 
# The dataframe contains the following variables (see Methods in the main text for citations and descriptions of variable transformations, where applicable):

# id = Standard Cross-Cultural Sample (SCCS) ID number; this is standardized across SCCS studies to easily combine variables from different studies
# socname = The name of each society; the spelling of some names may vary slightly across different datasets
# daily = the presence (1) or absence (0) of daily food sharing norm
# inter.trade = intercommunity trade, original scale
# food_storage = food storage, original scale
# soc_strat = social stratification, original scale
# DepGath_v203 = dependence on gathering (0-100%), binned into deciles
# DepGath_v204 = dependence on hunting (0-100%), binned into deciles
# DepGath_v205 = dependence on fishing (0-100%), binned into deciles
# DepGath_v206 = dependence on animal husbandry (0-100%), binned into deciles
# DepGath_v207 = dependence on agriculture (0-100%), binned into deciles
# Revised.latitude = latitude of the estimated center point of each societies' range, downloaded from DPLACE
# Revised.longitude = longitude of the estimated center point of each societies' range, downloaded from DPLACE
# focyear = the 'ethnographic present' or time foci for each society (codes are based on ethnographies written -15 years to +10 years around the ethnographic present)
# temp_pred = predictability of temperature
# npp_pred = predictability of net primary productivity (NPP)
# precip_pred = predictability of precipitation
# comm_size = mean size of local community
# labor_share = the presence (1) or absence (0) of daily labor sharing norm

############################################################################################
######################### Map plot (Figure 1) ##############################################
map_col <- ifelse(d$daily == 1, "#00000099", "#FF000099")
svg(filename="map.svg", 
    width=5, 
    height=8,
    pointsize=8)
map(database='worldHires', fill=TRUE, col="white", bg="lightblue", border=NA) # plots worldmap
points( d$Revised.longitude , d$Revised.latitude , xlab="longitude" , ylab="latitude" ,
        xlim=c(-180,180), ylim=c(-90,90), cex=0.8, col=map_col, pch=16)
legend(-180, 25, legend=c("Present", "Absent"),
       col=c("#00000099", "#FF000099"), pch=16, cex=0.8, bg="lightblue")

dev.off()

############################################################################################
############################################################################################
############### Importing SCCS Phylogeny ###################################################
sccs_tree <- read.nexus("Time-calibrated SCCS supertree.nex")

d$socname <- as.character(d$socname)
setdiff(d$socname, sccs_tree$tip.label) # checking for discrepancies between phylo tree names and dataframe names

# matching names to SCCS phylogeny, "MISSING" before socname indicates that genetic evidence was unavailable for this population, so phylogenetic position is based on lingustic data and genetic data from a closey-related group
d[c(1,11,13,27,30,43,55,61),"socname"] <- c("Kung_Bushmen", "Pastoral_Fulani", "Otoro_Nuba", "Punjabi_West", "Khalka_Mongols", "MISSING_Kimam","Copper_Eskimo", "MISSING_Huron")

# all names in the dataframe now match phylogeny if this returns "character(0)"
setdiff(d$socname, sccs_tree$tip.label)

# dropping tips of phylogeny that are not included in the sample
sccs_tree2 <- drop.tip(sccs_tree, subset(sccs_tree$tip.label, !(sccs_tree$tip.label %in% d$socname)))

# Plotting phylogeny
plot(sccs_tree2, type="fan", cex=0.8) 

d <- d[match(sccs_tree2$tip.label,d$socname),] # putting dataframe in same order as the phylogeny, THIS IS VERY IMPORTANT!

# Calculating phyogenetic distance; this returns a matrix of the patristic distances between each society in the phyloge
# Patristic distance is the sum of the branch lengths between each pair, and is thus 2x the time since the populations split; the units in this tree are in thousands of years (e.g., 280 = a patristic distance of 280,000 years and thus 140,000 years since divergence).
phy_dist <- cophenetic.phylo(sccs_tree2) 

# We will now scale these distances to make them unitless. This isn't strictly necessary, but given that we are working with 'distances' in ethngraphic present (EP) as well, standardizing the time makes it easier to set priors that are equally regularizing for both phylogeny and EP
phy_dist <- phy_dist / max(phy_dist)
max(phy_dist) # now the maximum distance between each society is scaled to equal 1

##########################################################################################
##########################################################################################
##################  Working with ethnographic presents ###################################
simplehist(d$focyear, xlab="Ethnographic Present") # distribution of dates

# Creating a distance matrix for ethnographic present
time_dist <- as.matrix( dist(d$focyear, method='euclidean'), nrow=length(d$id), ncol=length(d$id) )

# Scaling time marix
time_dist <- time_dist / max(time_dist)

##########################################################################################
##########################################################################################
##### Creating a distance matrix for geographic location (not used in main text analyses)#
long_lat <- c("Revised.longitude", "Revised.latitude")
geo_dist <- distm(d[,long_lat]) # pairwise geographic distance matrix; great circle distance
geo_dist <- geo_dist / max(geo_dist)
colnames(geo_dist) <- colnames(phy_dist)
rownames(geo_dist) <- colnames(geo_dist)

###########################################################################################
###########################################################################################
######## Identifying and indexing missing values ##########################################
# Multiple columns in our data frame contain NAs
sapply(d, function(y) sum(length(which(is.na(y))))) # count of missing values in each column

# Rather than drop all rows with a missing value (which would throw away a lot of other information in those rows), we will perform Bayesian imputation, which replaces missing values with a parameter and propogates uncertainty while fitting the model. We will make the conventional assumption that values are 'Missing Completely at Random' (MCAR), which means that we believe that observations are missing at random with respect to their own values (e.g., we do not belive that larger community sizes are more likely to be missing).

# For each column with discrete missing values we need a binary variable indicating whether its is missing. For variables that we are treating as continous, we need an index from 1:N_missing in order to assign a unique parameter to each observation. We do not do this for discrete missing values because Hamiltonian MCMC cannot sample discrete parameters, and thus we use a different strategy for imputing discrete variables.

# Indexing missing values for food storage
food_store_missing <- ifelse( is.na(d$food_storage), 1, 0 )

# Indexing missing values for temperature predictability
temp_pred_missing <- ifelse( is.na(d$temp_pred), 1, 0 )
# Creating an index of missing value positions 
temp_pred_missing <- sapply( 1:length(temp_pred_missing) ,
                             function(n) temp_pred_missing[n]*sum(temp_pred_missing[1:n]) )

# Indexing missing values for NPP predictability
npp_pred_missing <- ifelse( is.na(d$npp_pred), 1, 0 )
# Creating an index of missing value positions 
npp_pred_missing <- sapply( 1:length(npp_pred_missing) , 
                            function(n) npp_pred_missing[n]*sum(npp_pred_missing[1:n]) )

# Indexing missing values for precipitation predictability
precip_pred_missing <- ifelse( is.na(d$precip_pred), 1, 0 )
# Creating an index of missing value positions 
precip_pred_missing <- sapply( 1:length(precip_pred_missing) , 
                               function(n) precip_pred_missing[n]*sum(precip_pred_missing[1:n]) )

# Indexing missing values for community size
comm_size_missing <- ifelse( is.na(d$comm_size), 1, 0 )
# Creating an index of missing value positions 
comm_size_missing <- sapply( 1:length(comm_size_missing) , 
                             function(n) comm_size_missing[n]*sum(comm_size_missing[1:n]) )

# Indexing missing values for labor sharing
labor_share_missing <- ifelse( is.na(d$labor_share), 1, 0 )

##########################################################################################
##########################################################################################
################# Prepping variables for analysis #######################################

# For continous predictors, we standardize the observations by subtracting the mean from each value and dividing by 2 standard deviations. The former helps in model fitting, setting priors, and interpreting regression coefficients by making the mean = 0. The latter makes each variable unitless, enhancing interpretability of a one unit increase/decrease. We divide by two standard deviations rather than the conventional one standard deviation in order to make the effect sizes of discrete predictors and cont(inous predictors more directly comparable. See Gelman (2008) 'Scaling regression inputs by dividing by two standard deviations' for the logic underyling this standardization. Keep in mind that this scaling means that one unit in our regression coefficients = 2 SD.

# To visualize the difference in this transformation, let's simulate!
y_sim <- rnorm( n=1000 , mean = 13 , sd = 3.5 ) # an arbitrary gaussian variable
y_sim_s1 <- (y_sim - mean(y_sim)) / sd(y_sim) # conventional standardization
y_sim_s2 <- (y_sim - mean(y_sim)) / (sd(y_sim)*2) # 2 SD standardization
dens(y_sim_s1, col="cornflowerblue", ylim=c(0,1.3))
dens(y_sim_s2, col="coral", add=TRUE)
# Note that in the 2 SD standardization (orange density), a one unit increase implies twice as much variance as a one SD standardization.

# Because Stan will not accept NA values, we will replace NAs with a placeholder. I prefer to use absurd values like '-9999' so that if something goes wrong, it's easier to catch.

# Hunting
hunt_s <- (d$DepHunt_v204 - mean(d$DepHunt_v204, na.rm=TRUE)) / (sd(d$DepHunt_v204, na.rm=TRUE)*2)

# Animal Hus
hus_s <- (d$DepAnHusb_v206 - mean(d$DepAnHusb_v206, na.rm=TRUE)) / (sd(d$DepAnHusb_v206, na.rm = TRUE)*2)

# Temp Predictability
temp_pred_s <- (d$temp_pred - mean(d$temp_pred, na.rm=TRUE)) / (sd(d$temp_pred, na.rm=TRUE)*2)
temp_pred_s <- ifelse(temp_pred_missing > 0, -9999, temp_pred_s) # replacing NAs with placeholder

# NPP Predictability
npp_pred_s <- (d$npp_pred - mean(d$npp_pred, na.rm=TRUE)) / (sd(d$npp_pred, na.rm=TRUE)*2)
npp_pred_s <- ifelse(npp_pred_missing > 0, -9999, npp_pred_s)

# Precip Predictability
precip_pred_s <- (d$precip_pred - mean(d$precip_pred, na.rm=TRUE)) / (sd(d$precip_pred, na.rm=TRUE)*2)
precip_pred_s <- ifelse(precip_pred_missing > 0, -9999, precip_pred_s)

# Comm Size
comm_size_s <- (d$comm_size - mean(d$comm_size, na.rm = TRUE)) / (sd(d$comm_size,na.rm=TRUE)*2)
comm_size_s <- ifelse(comm_size_missing > 0, -9999, comm_size_s)

# Food Storage
food_store <- ifelse(d$food_storage > 1, 1, 0)
food_store <- ifelse(food_store_missing > 0, -9999, food_store)

# External Trade
trade_food <- ifelse(d$inter.trade > 4, 1, 0)

# Stratification
strat <- ifelse(d$soc_strat > 1, 1, 0)

# Labor Sharing
labor_share <- ifelse(labor_share_missing > 0, -9999, d$labor_share)

# Housekeeping for Stan
N <- nrow(d) # total number of observations
N_society <- length( unique(d$id) ) # same as N in this case
# For indexing random effects, Stan requires a vector from 1:N
society <- match(d$id, unique(d$id))

# Putting data into a list
data_list <- list(
  N = N,
  N_society = N_society,
  society = society,
  y = d$daily, # outcome variable, daily food sharing
  food_share = d$daily, # outcome again, for imputation model this time
  phy_dist = phy_dist,
  time_dist = time_dist,
  geo_dist = geo_dist,
  hunt_s = hunt_s,
  hus_s = hus_s,
  food_store = food_store,
  food_store_missing = food_store_missing,
  trade_food = trade_food,
  strat = strat,
  labor_share = labor_share,
  labor_share_missing = labor_share_missing,
  temp_pred_s = temp_pred_s,
  temp_pred_missing = temp_pred_missing,
  temp_pred_num_missing = length(temp_pred_s[temp_pred_s == -9999]),
  npp_pred_s = npp_pred_s,
  npp_pred_missing = npp_pred_missing,
  npp_pred_num_missing = length(npp_pred_s[npp_pred_s == -9999]),
  precip_pred_s = precip_pred_s,
  precip_pred_missing = precip_pred_missing,
  precip_pred_num_missing = length(precip_pred_s[precip_pred_s == -9999]),
  comm_size_s = comm_size_s,
  comm_size_missing = comm_size_missing,
  comm_size_num_missing = length(comm_size_s[comm_size_s == -9999])
)

#### Are our predictors multicolinear? Check using vif ##################################
library(car)
vif_df <- data.frame(
  y = data_list$y,
  hunt_s = hunt_s,
  food_store = food_store,
  trade_food = trade_food,
  strat = strat,
  hus_s = hus_s,
  precip_pred_s = precip_pred_s,
  temp_pred_s = temp_pred_s,
  npp_pred_s = npp_pred_s,
  comm_size_s = comm_size_s,
  labor_share = labor_share
)

# Counting whether any values are missing
missing_cases <- precip_pred_missing + food_store_missing + labor_share_missing + comm_size_missing
vif_df <- vif_df[missing_cases == 0, ]

# Checking the variance inflation factors, all below the commonly-used threshold of 10
vif( glm( y ~ hunt_s + food_store + trade_food + strat + hus_s + precip_pred_s + temp_pred_s + npp_pred_s + comm_size_s + labor_share, data=vif_df, family="binomial") )

##### Model fitting ######################################
n_chains_i <- 3

#### Fitting fixed effects + phylogney + EP
set.seed(250) # setting seed for reproducibility
fit_m_phy_EP <- stan(file="m_phy_EP_mv.stan", data=data_list, chains=n_chains_i, cores=n_chains_i, control=list(adapt_delta=0.99), iter=2000, init="0")

#### Fitting fixed effects + phylogney
set.seed(250)
fit_m_phy <- stan(file="m_phy_mv.stan", data=data_list, chains=n_chains_i, cores=n_chains_i, control=list(adapt_delta=0.99), iter=2000, init="0")

#### Fitting fixed effects + EP
set.seed(250)
fit_m_EP <- stan(file="m_EP_mv.stan", data=data_list, chains=n_chains_i, cores=n_chains_i, control=list(adapt_delta=0.99), iter=2000, init="0")

#### Fitting fixed effects only
set.seed(250)
fit_m_fixed <- stan(file="m_fixed_mv.stan", data=data_list, chains=n_chains_i, cores=n_chains_i, control=list(adapt_delta=0.99), iter=2000, init="0")


#### Fitting phy + EP, no fixed effects
set.seed(250)
fit_m_nofixed <- stan(file="m_phy_EP_mv_nofix.stan", data=data_list, chains=n_chains_i, cores=n_chains_i, control=list(adapt_delta=0.96), iter=2000, init="0")

#### Phylogenetic Mediation model for all variables, THIS ONE TAKES A WHILE TO RUN
set.seed(250)
fit_m_mediation <- stan(file="m_phy_EP_mediation.stan", data=data_list, chains=n_chains_i, cores=n_chains_i, control=list(adapt_delta=0.96), iter=2000, init="0")

############################################################################################
########################### Checking model convergence #####################################

# Should be a high number of effective samples, Rhat values very near 1
precis(fit_m_phy_EP)
precis(fit_m_phy)
precis(fit_m_EP)
precis(fit_m_fixed)

# Chains should have converged nicely. High values of eta and rho should not be an issue: they are only weakly identified, but their ratio is stable
tracerplot(fit_m_phy_EP, pars=c("b", "eta_phy", "rho_phy", "eta_ep", "rho_ep"))

# Summary of parameter estimates
precis(fit_m_phy_EP, depth=2, pars=c("b"))
post <- extract.samples(fit_m_phy_EP)
apply(post$b, 2, HPDI, prob=0.9)

precis(fit_m_phy_EP, pars=c("eta_phy", "eta_ep","rho_phy", "rho_ep"))
round(HPDI(post$eta_phy, prob=0.9),2)
round(HPDI(post$eta_ep, prob=0.9),2)
round(HPDI(post$rho_phy, prob=0.9),2)
round(HPDI(post$rho_ep, prob=0.9),2)

#############################################################################################
########################### Model Comparison using WAIC #####################################
rethinking::compare(fit_m_phy_EP, fit_m_phy, fit_m_EP, fit_m_fixed)

#############################################################################################
########################### Cohen's d (Figure 2A) ###########################################
post <- extract.samples(fit_m_phy_EP)
plot_cols <- c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest2"), "mediumorchid1")
plot_cols[7] <- "palegreen3"

cohens_d <- post$b[,c(2:11)] * (sqrt(3)/pi)
colnames(cohens_d) <- c("Hunting", "Food Storage", "External Trade", "Social Stratification", "Animal Husbandry", "Precipitation Pred.", "Temperature Pred.", "NPP Pred.", "Labor Sharing", "Community Size")

cohens_df <- gather(as.data.frame(cohens_d), key="var", value="est")
cohens_df <- gather(as.data.frame(cohens_d), key="var", value="est") 

# We want to sort by the absolute value of the median d
median_d <- abs( apply(cohens_d, 2, median) )

cohens_df$var <- factor(cohens_df$var, levels=names( sort(median_d) ))

# Now we want the posterior probability that each parameter is in our predicted direction

post_probs <- c(
  length(post$b[,2][post$b[,2] > 0]),
  length(post$b[,7][post$b[,7] < 0]),
  length(post$b[,5][post$b[,5] < 0]),
  length(post$b[,11][post$b[,11] < 0]),
  length(post$b[,4][post$b[,4] < 0]),
  length(post$b[,8][post$b[,8] < 0]),
  length(post$b[,3][post$b[,3] < 0]),
  length(post$b[,6][post$b[,6] < 0]),
  length(post$b[,10][post$b[,10] > 0]),
  length(post$b[,9][post$b[,9] < 0])
)

post_probs <- round(post_probs / length(post$b[,1]),2)*100
post_probs <- paste0(as.character(as.factor(post_probs)), "%")
post_probs[9] <- ">99%"

post_probs <- data.frame(
  labs = as.character(post_probs),
  x = c(2.3, rep(-2.3, 7), 2.3, -2.3),
  y = seq(2,11),
  jitter = c(rep(0.45, 9), 0.45)
)

svg(filename="Cohens_d.svg", 
    width=6, 
    height=8, 
    pointsize=12)

ggplot(cohens_df, aes(x=est, y=var)) + geom_hline(yintercept = c(1:10), alpha=0.6) + geom_density_ridges2(aes(fill=var, color=var), alpha=0.7, scale=0.8, rel_min_height=0.01) + ylab("") + xlab("Cohen's d") + theme_bw(base_size=20) + scale_fill_manual(values=plot_cols) + scale_color_manual(values=plot_cols) + scale_y_discrete(expand = c(0.00, 0)) + theme(legend.title = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", axis.text.y=element_text(vjust=0), axis.ticks.y = element_blank()) + geom_vline(xintercept = 0, linetype="dashed", lwd=1) + annotate("text", x=post_probs$x, y=post_probs$y-post_probs$jitter, label=as.character(post_probs$labs), color=plot_cols, size=5.5) + annotate("blank", x = 0, y=11) + scale_x_continuous(limits=c(-2.9,2.9))

dev.off()

########################### R^2 (Figure 3) ###############################################
#### R^2, conditional on Fixed Effects ####
post <- extract.samples(fit_m_phy_EP)
post2 <- extract.samples(fit_m_nofixed)
post3 <- extract.samples(fit_m_fixed)

fixed_var <- apply(post$y_hat, 1, var)
phy_var <- apply(post$phy_v, 1, var)
EP_var <- apply(post$ep_v, 1, var)
latent_var <- pi^2/3

phy_var2 <- apply(post2$phy_v, 1, var)
EP_var2 <- apply(post2$ep_v, 1, var)

fixed_var2 <- apply(post3$y_hat, 1, var)

r2_phy <- phy_var / (fixed_var + phy_var + EP_var + latent_var)
r2_EP <- EP_var / (fixed_var + phy_var + EP_var + latent_var)
r2_fixed <- fixed_var / (fixed_var + phy_var + EP_var + latent_var)

r2_phy2 <- phy_var2 / (phy_var2 + EP_var2 + latent_var)
r2_EP2 <- EP_var2 / (phy_var2 + EP_var2 + latent_var)

r2_fixed2 <- fixed_var2 / (fixed_var2 + latent_var)

r2_df <- data.frame(
  r2 = c(r2_phy, r2_phy2, r2_EP, r2_EP2, r2_fixed, r2_fixed2),
  var = rep(c("Phylogeny", "Ethnographic Present", "Fixed Effects"), times=c(length(r2_phy)*2, length(r2_phy)*2, length(r2_phy)*2)),
  type = c(rep("adjusted", length(r2_phy)), rep("unadjusted", length(r2_phy)), rep("adjusted", length(r2_phy)), rep("unadjusted", length(r2_phy)), rep("adjusted", length(r2_phy)), rep("unadjusted", length(r2_fixed2)))
)

r2_df$var <- factor(r2_df$var, levels=c("Ethnographic Present", "Phylogeny", "Fixed Effects"))
r2_df$type <- factor(r2_df$type, levels=c("unadjusted", "adjusted"))

svg(filename="r2.svg", 
    width=8, 
    height=6, 
    pointsize=12)

ggplot(r2_df, aes(x=r2, y=var)) + geom_hline(yintercept=seq(1:3), alpha=0.6) + geom_density_ridges2(aes(fill=var, color=var, linetype=type, alpha=type), lwd=0.8,scale=0.9, rel_min_height=0.01) + scale_x_continuous(limits=c(0,0.8), expand = c(0.00, 0.05)) + scale_y_discrete(expand = c(0.00, 0))  + theme_bw(base_size=20) + annotate("blank", x = 0, y=4) + xlab( expression(R^2)) + ylab("") + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_fill_manual(values=c("cornflowerblue", "orange", "coral")) + scale_color_manual(values=c("cornflowerblue", "orange", "coral")) + scale_linetype_manual(values=c("dashed", "solid")) + scale_alpha_manual(values=c(0,0.3))
 
dev.off()

#### Phylogenetic Mediation ####
post_med <- extract.samples(fit_m_mediation)

## Food Sharing
food_phy <- apply(post_med$phy_v[,,1], 1, var) / (apply(post_med$phy_v[,,1], 1, var) + pi^2/3 )

## Hunting
hunting_phy <- apply(post_med$phy_v[,,2], 1, var) / (apply(post_med$phy_v[,,2], 1, var) + 0.5^2 )

## Food Storage
store_phy <- apply(post_med$phy_v[,,3], 1, var) / (apply(post_med$phy_v[,,3], 1, var) + pi^2/3 )

## Trade
trade_phy <- apply(post_med$phy_v[,,4], 1, var) / (apply(post_med$phy_v[,,4], 1, var) + pi^2/3 )

## Strat
strat_phy <- apply(post_med$phy_v[,,5], 1, var) / (apply(post_med$phy_v[,,5], 1, var) + pi^2/3 )

## Animal Hus
hus_phy <- apply(post_med$phy_v[,,6], 1, var) / (apply(post_med$phy_v[,,6], 1, var) + 0.5^2 )

## Precip
precip_phy <- apply(post_med$phy_v[,,7], 1, var) / (apply(post_med$phy_v[,,7], 1, var) + 0.5^2 )

## Temp
temp_phy <- apply(post_med$phy_v[,,8], 1, var) / (apply(post_med$phy_v[,,8], 1, var) + 0.5^2 )

## NPP
npp_phy <- apply(post_med$phy_v[,,9], 1, var) / (apply(post_med$phy_v[,,9], 1, var) + 0.5^2 )

## labor
labor_phy <- apply(post_med$phy_v[,,10], 1, var) / (apply(post_med$phy_v[,,10], 1, var) + pi^2/3 )

## comm size
comm_phy <- apply(post_med$phy_v[,,11], 1, var) / (apply(post_med$phy_v[,,11], 1, var) + 0.5^2 )

phy_df <- data.frame(
  food_phy = food_phy,
  hunting_phy = hunting_phy,
  store_phy = store_phy,
  trade_phy = trade_phy,
  strat_phy = strat_phy,
  hus_phy = hus_phy,
  precip_phy = precip_phy,
  temp_phy = temp_phy,
  npp_phy = npp_phy,
  labor_phy = labor_phy,
  comm_phy = comm_phy
)

names(phy_df) <- c("Food Sharing", "Hunting", "Food Storage", "External Trade", "Social Stratification", "Animal Husbandry", "Precipitation Pred.", "Temperature Pred.", "NPP Pred.", "Labor Sharing", "Community Size")
phy_med <- apply(phy_df, 2, median)

phy_df <- gather(phy_df, key="var", value="est")
phy_df$var <- factor(phy_df$var, levels=names( sort(phy_med) ))

ggplot(phy_df, aes(x=est, y=var)) + geom_hline(yintercept=seq(1:11), alpha=0.6) + geom_density_ridges2(fill="orange", color="orange", lwd=0.8,scale=0.9, rel_min_height=0.01, alpha=0.5) + scale_x_continuous(limits=c(0,1), expand = c(0.00, 0.05)) + scale_y_discrete(expand = c(0.00, 0))  + theme_bw(base_size=18) + annotate("blank", x = 0, y=12) + xlab( expression(R^2)) + ylab("") + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("Phylogenetic Signal")

#########################################################################################
################## GP Covariance Functions ##############################################
post <- extract.samples(fit_m_phy_EP)
pred_seq <- seq(from=0, to=1, length.out = 25)

svg(filename="GP_cov.svg", 
    width=8, 
    height=6, 
    pointsize=8)

par(mfrow=c(1,2), mar=c(5.1, 5.1, 4.1, 2.1))

phy_GP <- matrix(NA, nrow=length(post$lp__), ncol=length(pred_seq))
for (i in 1:nrow(phy_GP)) {
  phy_GP[i,] = post$eta_phy[i] * exp(-(post$rho_phy[i]*pred_seq))
}

plot(NULL, xlim=c(0,1), ylim=c(0,5), xlab="Distance", ylab="Covariance", xaxt='n', cex.lab=1.8, cex.axis=1.8)
axis(1, at=c(0,1), labels=c("Min", "Max"), cex.lab=1.8, cex.axis=1.8)
shade(apply(phy_GP, 2, HPDI, prob=0.90), pred_seq, col=col.alpha("orange", 0.3))
shade(apply(phy_GP, 2, HPDI, prob=0.60), pred_seq, col=col.alpha("orange", 0.3))
shade(apply(phy_GP, 2, HPDI, prob=0.30), pred_seq, col=col.alpha("orange", 0.3))
mtext("Phylogenetic Covariance Function", cex=1.8)

ep_GP <- matrix(NA, nrow=length(post$lp__), ncol=length(pred_seq))
for (i in 1:nrow(ep_GP)) {
  ep_GP[i,] = post$eta_ep[i] * exp(-(post$rho_ep[i]*pred_seq))
}

plot(NULL, xlim=c(0,1), ylim=c(0,5), xlab="Distance", ylab="Covariance", xaxt='n', cex.lab=1.8, cex.axis=1.8)
axis(1, at=c(0,1), labels=c("Min", "Max"), cex.lab=1.8, cex.axis=1.8)
shade(apply(ep_GP, 2, HPDI, prob=0.90), pred_seq, col=col.alpha("cornflowerblue", 0.3))
shade(apply(ep_GP, 2, HPDI, prob=0.60), pred_seq, col=col.alpha("cornflowerblue", 0.3))
shade(apply(ep_GP, 2, HPDI, prob=0.30), pred_seq, col=col.alpha("cornflowerblue", 0.3))
mtext("EP Covariance Function", cex=1.8)

dev.off()
#########################################################################################
############# Additional tests/sanity checks ############################################
#########################################################################################

############################ Visuaizing our priors ######################################
svg(filename="priors.svg", 
    width=8, 
    height=6, 
    pointsize=8)

par(mfrow=c(1,2), mar=c(5.1, 5.1, 4.1, 2.1))

### Cohen's d for main effects
curve( dnorm( x*(sqrt(3)/pi), mean=0, sd=1), from=-6, to=6, lwd=3, xlab="Cohen's d", ylab="Density", cex.lab=1.8, cex.axis=1.8)
mtext("Main Effects Prior", cex=1.8)

### GP Function Prior
rho_prior <- rexp(n=200, 0.5)
rho_prior <- rho_prior[rho_prior > 0]
eta_prior <- rexp(n=200, 0.5)
eta_prior <- eta_prior[eta_prior > 0]
curve( ( median(eta_prior)*exp(-median(rho_prior)*x) ), from=0, to=1 ,lwd=0, ylab="Covariance", xlab="Distance", ylim=c(0,20), cex.lab=1.8, cex.axis=1.8)
for (i in 1:200) {
  curve( (eta_prior[i]*exp(-rho_prior[i]*x)), from=0, to=1, add=TRUE, col=col.alpha("black", 0.2) )
}
mtext("GP Covariance Prior", cex=1.8)

dev.off()

########################################################################
### How well can we predict location from phylogeny? ##################

xy_phy <- t(combn(colnames(phy_dist), 2))
phy_df <- data.frame(xy_phy, dist=phy_dist[xy_phy])

xy_geo <- t(combn(colnames(geo_dist), 2))
geo_df <- data.frame(xy_geo, dist=geo_dist[xy_geo])

society <- coerce_index(phy_df$X1)

beta_normalize <- function(x) {
  x_ <- ((x - min(x)) / (max(x) - min(x)))
  (x_ * (length(x_) - 1) + 0.5) / length(x_)
}

data_list_spatial <- list(
  N = N,
  N_dist = nrow(phy_df),
  society = society,
  geo_dist = beta_normalize(geo_df$dist),
  phy_dist = beta_normalize(phy_df$dist)
)

m_spatial_phy <- "
data{
int N;
int N_dist;
int society[N_dist];
vector[N_dist] geo_dist;
vector[N_dist] phy_dist;
}

parameters{
real b[2]; // intercept
real<lower=0> phi[2];
vector<lower=0>[2] sigma;
vector<lower=0>[2] sigma_soc;
matrix[2, N_dist] obs_z;
matrix[N,2] soc_z;
cholesky_factor_corr[2] L_Rho;
}

transformed parameters{
vector<lower=0,upper=1>[N_dist] mu_phy;
vector<lower=0,upper=1>[N_dist] mu_geo;
vector[N_dist] p_phy;
vector[N_dist] q_phy;
vector[N_dist] p_geo;
vector[N_dist] q_geo;
matrix[N_dist,2] obs_v;
matrix[N,2] soc_v;

obs_v = (diag_pre_multiply(sigma, L_Rho) * obs_z)';
soc_v[,1] = sigma_soc[1] * soc_z[,1];
soc_v[,2] = sigma_soc[2] * soc_z[,2];

for  (i in 1:N_dist) {
mu_phy[i] = inv_logit(b[1] + obs_v[i,1] + soc_v[society[i],1]);
mu_geo[i] = inv_logit(b[2] + obs_v[i,2] + soc_v[society[i],2]);

p_phy[i] = mu_phy[i] * phi[1];
q_phy[i] = (1 - mu_phy[i]) * phi[1];

p_geo[i] = mu_geo[i] * phi[2];
q_geo[i] = (1 - mu_geo[i]) * phi[2];
}
}

model{
// priors
b ~ normal(0,2);
phi ~ cauchy(0,2);
sigma ~ cauchy(0,2);
sigma_soc ~ cauchy(0,2);
L_Rho ~ lkj_corr_cholesky(2);
to_vector(obs_z) ~ normal(0,1);
to_vector(soc_z) ~ normal(0,1);

for (i in 1:N_dist) {
phy_dist[i] ~ beta( p_phy[i], q_phy[i] );
geo_dist[i] ~ beta( p_geo[i], q_geo[i] );
}

}

"

fit_m_spatial_phy <- stan(model_code=m_spatial_phy, data=data_list_spatial, iter=1500, chains=2, cores=2, control=list(adapt_delta=0.98))
post_spatial <- extract.samples(fit_m_spatial_phy)

cor_spatial_phy <- array(0, dim=c(length(post_spatial$lp__), 2, 2))
for (i in 1:nrow(cor_spatial_phy)) {
  cor_spatial_phy[i,,] = post_spatial$L_Rho[i,,] %*% t(post_spatial$L_Rho[i,,])
}

dens(cor_spatial_phy[,1,2])
median(cor_spatial_phy[,1,2])
HPDI(cor_spatial_phy[,1,2], prob=0.9)
#########################################################################
### Pairs plot of hunting and all other predictors ######################
hunt_df <- data.frame(
  food_store = food_store,
  trade_food = trade_food,
  strat = strat,
  hus_s = hus_s,
  precip_pred_s = precip_pred_s,
  temp_pred_s = temp_pred_s,
  npp_pred_s = npp_pred_s,
  labor_share = labor_share,
  comm_size_s = comm_size_s
)

# Dealing with missingness
hunt_df <- apply(hunt_df, 2, function(x) ifelse(x == -9999, NA, x))

plot_cols <- c("#A8DFA8", "#92D3E5", "#F6C957", "#EC9BFF", "#57C0B3", "#EDC0D7", "#E5C4BB", "#A1B9E3", "#FBAE57")

names <- c("Food Storage", "External Trade", "Social Stratification", "Animal Husbandry", "Precipitation Pred.", "Temperature Pred.", "NPP Pred.", "Labor Sharing", "Community Size")
binary_pred <- c(1,1,1,0,0,0,0,1,0)  # indicates whether need to make boxplot

par(cex=1.6, pty="s", mfrow=c(3,3), mar=c(2.7,1,2,1))

for (i in 1:9) {
spearman <- cor.test(d$DepHunt_v204, hunt_df[,i], method="spearman") # calculate spearmans' rho

if (binary_pred[i] == 0)
  plot(jitter(hunt_s) ~ jitter(hunt_df[,i]), ylab="Hunting", xlab=names[i], pch=16, col=plot_cols[i])

else
  boxplot(hunt_s ~ hunt_df[,i], ylab="Hunting", xlab=names[i], col=col.alpha(plot_cols[i], 0.4))

mtext(paste("rho =", round(as.numeric(spearman$estimate),2)), cex=0.7) # plot the correlation
}
#########################################################################
#### Is hunting correlated with sharing in a bivariate model? ###########
set.seed(79)
fit_m_hunt <- stan(file="m_hunt.stan", data=data_list, chains=3, cores=3, control=list(adapt_delta=0.96), iter=2000, init="0" )

post_hunt <- extract.samples(fit_m_hunt)
d_hunt <- post_hunt$b[,2] * (sqrt(3)/pi)
length(post_hunt$b[,2][post_hunt$b[,2] > 0]) / length(post_hunt$b[,2]) # posterior probability

pred_seq <- seq(from=-1,to=1, length.out=25)

hunt <- matrix(0, nrow=250, ncol=length(pred_seq))
for (i in 1:250) {
  hunt[i,] = logistic(post_hunt$b[i,1] + post_hunt$b[i,2]*pred_seq)
}

dev.off()
par(cex=1.4, mar=c(5,5,2,2))

plot(NULL, xlab="Hunting z-score", ylab="Pr (Food Sharing)", xlim=c(-1,1), ylim=c(0,1), xaxt='n', axes=F, ann=F)
axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), line=1)
axis(1, at=seq(-1,1,length.out = 5), labels=seq(-1,1,length.out = 5)*2, line=1)
shade(apply(hunt, 2, PI, 0.90), pred_seq, col=col.alpha(plot_cols[1], 0.18))
shade(apply(hunt, 2, PI, 0.72), pred_seq, col=col.alpha(plot_cols[1], 0.18))
shade(apply(hunt, 2, PI, 0.54), pred_seq, col=col.alpha(plot_cols[1], 0.18))
shade(apply(hunt, 2, PI, 0.36), pred_seq, col=col.alpha(plot_cols[1], 0.18))
shade(apply(hunt, 2, PI, 0.18), pred_seq, col=col.alpha(plot_cols[1], 0.18))
lines(x=pred_seq, y=apply(hunt,2,median), col=plot_cols[1], lwd=3)
mtext("Bivariate Model", cex=1.5, line=0.5)
mtext(side=1, "Hunting z-score", cex=1.5, line=3.5)
mtext(side=2, "Pr (Food Sharing)", cex=1.5, line=3.5)

#########################################################################, 
#### What if we substitude just hunting for a summary of hunting, gathering, and fishing?
hgf_s <- d$DepHunt_v204 + d$DepGath_v203 + d$DepFish_v205
hgf_s <- (hgf_s - mean(hgf_s)) / (2*sd(hgf_s))

data_list2 <- data_list
data_list2$hunt_s <- hgf_s # substituting the hunt_s variable for hgf_s so we don't have to re-write the model

set.seed(79)
fit_m_phy_EP_hgf <- stan(file="m_phy_EP_mv.stan", data=data_list2, chains=3, cores=3, control=list(adapt_delta=0.95), iter=2000, init="0")
post2 <- extract.samples(fit_m_phy_EP_hgf)
median(post2$b[,2] * (sqrt(3)/pi))

length(post2$b[,2][post2$b[,2] > 0]) / length(post2$b[,2]) # posterior probability

pred_seq <- seq(from=-1,to=1, length.out=25)

hunt <- matrix(0, nrow=250, ncol=length(pred_seq))
for (i in 1:250) {
  hunt[i,] = logistic(post2$b[i,1] + post2$b[i,2]*pred_seq)
}

par(cex=1.4)
plot(NULL, xlab="Hunting z-score", ylab="Pr (Food Sharing)", xlim=c(-1,1), ylim=c(0,1), xaxt='n', axes=F, ann=F)
axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), line=1)
axis(1, at=seq(-1,1,length.out = 5), labels=seq(-1,1,length.out = 5)*2, line=1)
shade(apply(hunt, 2, PI, 0.90), pred_seq, col=col.alpha(plot_cols[1], 0.18))
shade(apply(hunt, 2, PI, 0.72), pred_seq, col=col.alpha(plot_cols[1], 0.18))
shade(apply(hunt, 2, PI, 0.54), pred_seq, col=col.alpha(plot_cols[1], 0.18))
shade(apply(hunt, 2, PI, 0.36), pred_seq, col=col.alpha(plot_cols[1], 0.18))
shade(apply(hunt, 2, PI, 0.18), pred_seq, col=col.alpha(plot_cols[1], 0.18))
lines(x=pred_seq, y=apply(hunt,2,median), col=plot_cols[1], lwd=3)
mtext("Hunting, Gathering, and Fishing", cex=1.5, line=0.5)
mtext(side=1, "Hunting z-score", cex=1.5, line=3.5)
mtext(side=2, "Pr (Food Sharing)", cex=1.5, line=3.5)

#########################################################################
###### Plotting phylogenetic tree with maximum likelihood ASE ###########

daily <- d$daily
names(daily) <- sccs_tree2$tip.label
obj<- contMap(sccs_tree2, daily, plot=FALSE)

# change color palette --> see http://blog.phytools.org/2014/05/changing-color-ramp-in-contmap-or.html
n<- length(obj$cols)
obj$cols[1:n]<- colorRampPalette(c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"), space="Lab")(n)

plot(obj,legend=FALSE,lwd=6,outline=FALSE,fsize=c(0.9,1), sig=0, cex.lab=0.5)
add.color.bar(leg=29.8,cols=obj$cols,title="",prompt=FALSE,x=20,
              y=53.5,lwd=7,fsize=1.2,subtitle="Pr (Daily Food Share)", outline=FALSE)

#########################################################################################
#################### Posterior predictive plots (Figure 2B) #############################
post <- extract.samples(fit_m_phy_EP)
pred_seq <- seq(from=-1, to=1, length.out = 10)
plot_cols <- c(wes_palette("Darjeeling1"), wes_palette("GrandBudapest2"), "mediumorchid1")
plot_cols[7] <- "palegreen3"

### Posterior predictions
### Continous predictors will be projected on sequence from -2SD to +2SD ([-1,1] units given our standardization)
### Binary predictors will be projected onto a sequence c(0,1)
{
  labor <- matrix(0, nrow=250, ncol=2)
  for (i in 1:250) {
    labor[i,] = logistic(post$b[i,1] + post$b[i,10]*c(0,1))
  }
  food_store <- matrix(0, nrow=250, ncol=2)
  for (i in 1:250) {
    food_store[i,] = logistic(post$b[i,1] + post$b[i,3]*c(0,1))
  }
  npp_pred <- matrix(0, nrow=250, ncol=length(pred_seq))
  for (i in 1:250) {
    npp_pred[i,] = logistic(post$b[i,1] + post$b[i,9]*pred_seq)
  }
  trade_food <- matrix(0, nrow=250, ncol=2)
  for (i in 1:250) {
    trade_food[i,] = logistic(post$b[i,1] + post$b[i,4]*c(0,1))
  }
  strat <- matrix(0, nrow=250, ncol=2)
  for (i in 1:250) {
    strat[i,] = logistic(post$b[i,1] + post$b[i,5]*c(0,1))
  }
  hus <- matrix(0, nrow=250, ncol=length(pred_seq))
  for (i in 1:250) {
    hus[i,] = logistic(post$b[i,1] + post$b[i,6]*pred_seq)
  }
  temp_pred <- matrix(0, nrow=250, ncol=length(pred_seq))
  for (i in 1:250) {
    temp_pred[i,] = logistic(post$b[i,1] + post$b[i,8]*pred_seq)
  }
  comm_size <- matrix(0, nrow=250, ncol=length(pred_seq))
  for (i in 1:250) {
    comm_size[i,] = logistic(post$b[i,1] + post$b[i,11]*pred_seq)
  }
  precip_pred <- matrix(0, nrow=250, ncol=length(pred_seq))
  for (i in 1:250) {
    precip_pred[i,] = logistic(post$b[i,1] + post$b[i,7]*pred_seq)
  }
  hunt <- matrix(0, nrow=250, ncol=length(pred_seq))
  for (i in 1:250) {
    hunt[i,] = logistic(post$b[i,1] + post$b[i,2]*pred_seq)
  }
} # end predictions

svg(filename="pred_plot.svg", 
    width=6, 
    height=8, 
    pointsize=8)

par(mfrow=c(4,3), mar=c(4, 5, 0, 2), oma=c(1,3.5,1,0.5), cex.lab=2.4, pty='s', xaxs="r", las=1) 
{
  ##########################
  plot(NULL, xlab="", ylab="", xlim=c(-1,1), ylim=c(0,1), yaxt="n",  xaxt="n", axes=F, ann=F)
  axis(1, at=seq(-1,1,length.out = 5), labels=seq(-1,1,length.out = 5)*2, cex.axis=2.2, line=1)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=2.2, line=1, tck=-.02)
  shade(apply(npp_pred, 2, PI, 0.90), pred_seq, col=col.alpha(plot_cols[10], 0.18))
  shade(apply(npp_pred, 2, PI, 0.72), pred_seq, col=col.alpha(plot_cols[10], 0.18))
  shade(apply(npp_pred, 2, PI, 0.54), pred_seq, col=col.alpha(plot_cols[10], 0.18))
  shade(apply(npp_pred, 2, PI, 0.36), pred_seq, col=col.alpha(plot_cols[10], 0.18))
  shade(apply(npp_pred, 2, PI, 0.18), pred_seq, col=col.alpha(plot_cols[10], 0.18))
  lines(x=pred_seq, y=apply(npp_pred,2,median), col=plot_cols[10], lwd=3)
  mtext("NPP Pred. z-score", cex=1.5, line=0.5)
  ##############################
  plot(NULL, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxt='n',  yaxt='n' , axes=F, ann=F)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=2.2, line=1)
  axis(1, at=c(0.1,0.9), labels=c("Absent", "Present"), cex.axis=2.2, line=1, tck=-.02)
  shade(apply(labor, 2, PI, 0.90), c(0,1), col=col.alpha(plot_cols[9], 0.18))
  shade(apply(labor, 2, PI, 0.72), c(0,1), col=col.alpha(plot_cols[9], 0.18))
  shade(apply(labor, 2, PI, 0.54), c(0,1), col=col.alpha(plot_cols[9], 0.18))
  shade(apply(labor, 2, PI, 0.36), c(0,1), col=col.alpha(plot_cols[9], 0.18))
  shade(apply(labor, 2, PI, 0.18), c(0,1), col=col.alpha(plot_cols[9], 0.18))
  lines(x=c(0,1), y=apply(labor,2,median), col=plot_cols[9], lwd=3)
  mtext("Labor Sharing", cex=1.5, line=0.5)
  ##########################
  plot(NULL, xlab="", ylab="", xlim=c(-1,1), ylim=c(0,1),  xaxt="n", yaxt="n", axes=F, ann=F)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=2.2, line=1)
  axis(1, at=seq(-1,1,length.out = 5), labels=seq(-1,1,length.out = 5)*2, cex.axis=2.2, line=1, tck=-.02)
  shade(apply(hus, 2, PI, 0.90), pred_seq, col=col.alpha(plot_cols[8], 0.18))
  shade(apply(hus, 2, PI, 0.72), pred_seq, col=col.alpha(plot_cols[8], 0.18))
  shade(apply(hus, 2, PI, 0.54), pred_seq, col=col.alpha(plot_cols[8], 0.18))
  shade(apply(hus, 2, PI, 0.36), pred_seq, col=col.alpha(plot_cols[8], 0.18))
  shade(apply(hus, 2, PI, 0.18), pred_seq, col=col.alpha(plot_cols[8], 0.18))
  lines(x=pred_seq, y=apply(hus,2,median), col=plot_cols[8], lwd=3)
  mtext("Animal Hus. z-score", cex=1.5, line=0.5)
  ###########################
  plot(NULL, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxt='n',  yaxt="n", axes=F, ann=F)
  axis(1, at=c(0.1,0.9), labels=c("Absent", "Present"), cex.axis=2.2, line=1)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=2.2, line=1, tck=-.02)
  shade(apply(food_store, 2, PI, 0.90), c(0,1), col=col.alpha(plot_cols[7], 0.18))
  shade(apply(food_store, 2, PI, 0.72), c(0,1), col=col.alpha(plot_cols[7], 0.18))
  shade(apply(food_store, 2, PI, 0.54), c(0,1), col=col.alpha(plot_cols[7], 0.18))
  shade(apply(food_store, 2, PI, 0.36), c(0,1), col=col.alpha(plot_cols[7], 0.18))
  shade(apply(food_store, 2, PI, 0.18), c(0,1), col=col.alpha(plot_cols[7], 0.18))
  lines(x=c(0,1), y=apply(food_store,2,median), col=plot_cols[7], lwd=3)
  mtext("Food Storage", cex=1.5, line=0.5)
  ##############################
  plot(NULL, xlab="", ylab="", xlim=c(-1,1), ylim=c(0,1), xaxt='n', yaxt="n",  axes=F, ann=F)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=2.2, line=1, tck=-.02)
  axis(1, at=seq(-1,1,length.out = 5), labels=seq(-1,1,length.out = 5)*2, cex.axis=2.2, line=1)
  shade(apply(temp_pred, 2, PI, 0.90), pred_seq, col=col.alpha(plot_cols[6], 0.18))
  shade(apply(temp_pred, 2, PI, 0.72), pred_seq, col=col.alpha(plot_cols[6], 0.18))
  shade(apply(temp_pred, 2, PI, 0.54), pred_seq, col=col.alpha(plot_cols[6], 0.18))
  shade(apply(temp_pred, 2, PI, 0.36), pred_seq, col=col.alpha(plot_cols[6], 0.18))
  shade(apply(temp_pred, 2, PI, 0.18), pred_seq, col=col.alpha(plot_cols[6], 0.18))
  lines(x=pred_seq, y=apply(temp_pred,2,median), col=plot_cols[6], lwd=3)
  mtext("Temp Pred. z-score", cex=1.5, line=0.5)
  ##########################
  plot(NULL, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxt='n',  yaxt="n", axes=F, ann=F)
  axis(1, at=c(0.1,0.9), labels=c("Absent", "Present"), cex.axis=2.2, line=1)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=2.2, line=1, tck=-.02)
  shade(apply(trade_food, 2, PI, 0.90), c(0,1), col=col.alpha(plot_cols[5], 0.18))
  shade(apply(trade_food, 2, PI, 0.72), c(0,1), col=col.alpha(plot_cols[5], 0.18))
  shade(apply(trade_food, 2, PI, 0.54), c(0,1), col=col.alpha(plot_cols[5], 0.18))
  shade(apply(trade_food, 2, PI, 0.36), c(0,1), col=col.alpha(plot_cols[5], 0.18))
  shade(apply(trade_food, 2, PI, 0.18), c(0,1), col=col.alpha(plot_cols[5], 0.18))
  lines(x=c(0,1), y=apply(trade_food,2,median), col=plot_cols[5], lwd=3)
  mtext("External Trade", cex=1.5, line=0.5)
  ##############################
  plot(NULL, xlab="", ylab="", xlim=c(-1,1), ylim=c(0,1), xaxt='n', yaxt="n",  axes=F, ann=F)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=2.2, line=1)
  axis(1, at=seq(-1,1,length.out = 5), labels=seq(-1,1,length.out = 5)*2, cex.axis=2.2, line=1, tck=-.02)
  shade(apply(comm_size, 2, PI, 0.90), pred_seq, col=col.alpha(plot_cols[4], 0.18))
  shade(apply(comm_size, 2, PI, 0.72), pred_seq, col=col.alpha(plot_cols[4], 0.18))
  shade(apply(comm_size, 2, PI, 0.54), pred_seq, col=col.alpha(plot_cols[4], 0.18))
  shade(apply(comm_size, 2, PI, 0.36), pred_seq, col=col.alpha(plot_cols[4], 0.18))
  shade(apply(comm_size, 2, PI, 0.18), pred_seq, col=col.alpha(plot_cols[4], 0.18))
  lines(x=pred_seq, y=apply(comm_size,2,median), col=plot_cols[4], lwd=3)
  mtext("Comm. Size z-score", cex=1.5, line=0.5)
  #########################
  plot(NULL, xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxt='n',  yaxt='n', axes=F, ann=F)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=2.2, line=1)
  axis(1, at=c(0.1,0.9), labels=c("Egalitarian", "Stratified"), cex.axis=2.2, line=1, tck=-.02)
  shade(apply(strat, 2, PI, 0.90), c(0,1), col=col.alpha(plot_cols[3], 0.18))
  shade(apply(strat, 2, PI, 0.72), c(0,1), col=col.alpha(plot_cols[3], 0.18))
  shade(apply(strat, 2, PI, 0.54), c(0,1), col=col.alpha(plot_cols[3], 0.18))
  shade(apply(strat, 2, PI, 0.36), c(0,1), col=col.alpha(plot_cols[3], 0.18))
  shade(apply(strat, 2, PI, 0.18), c(0,1), col=col.alpha(plot_cols[3], 0.18))
  lines(x=c(0,1), y=apply(strat,2,median), col=plot_cols[3], lwd=3)
  mtext("Social Strat.", cex=1.5, line=0.5)
  ##############################
  plot(NULL, xlab="", ylab="", xlim=c(-1,1), ylim=c(0,1), xaxt="n", yaxt="n",  axes=F, ann=F)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=2.2, line=1)
  axis(1, at=seq(-1,1,length.out = 5), labels=seq(-1,1,length.out = 5)*2, cex.axis=2.2, line=1, tck=-.02)
  shade(apply(precip_pred, 2, PI, 0.90), pred_seq, col=col.alpha(plot_cols[2], 0.18))
  shade(apply(precip_pred, 2, PI, 0.72), pred_seq, col=col.alpha(plot_cols[2], 0.18))
  shade(apply(precip_pred, 2, PI, 0.54), pred_seq, col=col.alpha(plot_cols[2], 0.18))
  shade(apply(precip_pred, 2, PI, 0.36), pred_seq, col=col.alpha(plot_cols[2], 0.18))
  shade(apply(precip_pred, 2, PI, 0.18), pred_seq, col=col.alpha(plot_cols[2], 0.18))
  lines(x=pred_seq, y=apply(precip_pred,2,median), col=plot_cols[2], lwd=3)
  mtext("Precip Pred. z-score", cex=1.5, line=0.5)
  #############################
  plot(NULL, xlab="", ylab="", xlim=c(-1,1), ylim=c(0,1), xaxt='n', yaxt="n",  axes=F, ann=F)
  axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=2.2, line=1)
  axis(1, at=seq(-1,1,length.out = 5), labels=seq(-1,1,length.out = 5)*2, cex.axis=2.2, line=1, tck=-.02)
  shade(apply(hunt, 2, PI, 0.90), pred_seq, col=col.alpha(plot_cols[1], 0.18))
  shade(apply(hunt, 2, PI, 0.72), pred_seq, col=col.alpha(plot_cols[1], 0.18))
  shade(apply(hunt, 2, PI, 0.54), pred_seq, col=col.alpha(plot_cols[1], 0.18))
  shade(apply(hunt, 2, PI, 0.36), pred_seq, col=col.alpha(plot_cols[1], 0.18))
  shade(apply(hunt, 2, PI, 0.18), pred_seq, col=col.alpha(plot_cols[1], 0.18))
  lines(x=pred_seq, y=apply(hunt,2,median), col=plot_cols[1], lwd=3)
  mtext("Hunting z-score", cex=1.5, line=0.5)
}
dev.off()
########################## Posterior predictive checks ###################################
y_pred <- post$y_hat + post$phy_v + post$ep_v
y_pred <- logistic(y_pred)

y_sims <- matrix(NA, length(post$a[,1]), data_list$N)
for (i in 1:length(post$a[,1]))
  for (j in 1:data_list$N) {
  y_sims[i,j] = rbinom(1, size=1, prob=y_pred[i,j])
}

y_prop <- rowSums(y_sims)

hist(y_prop, breaks=10,  lty="blank", col=col.alpha("skyblue", 0.7), main="", ylab="", xlab="Predicted Frequency of Food Sharing", xlim=c(0,73), yaxt='n', yaxs='i', cex=1.5, cex.axis=1.5, cex.lab=1.5)
abline(v=sum(data_list$y), lwd=3, lty=2, col="orange" )
text(x=sum(data_list$y) + 30, y=500, expression("Observed Frequency =" ~ over(22,73)), col="orange", cex=1.4)

dev.off()
