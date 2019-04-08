## This analysis relies on a fit model and data from Koster & McElreath et al. (2019), available at https://osf.io/2kzb6/
#devtools::install_github("rmcelreath/cchunts")

setwd("")
library(rethinking)
library(cchunts)
library(geosphere)
library(jsonlite)
load("model_fix_17092018.Rdata") # this needs to be downloaded from https://osf.io/2kzb6/

post <- extract.samples(mfit)

# load map indices and labels
data(cch_mapkey)
key <- cch_mapkey

labels <- as.character( key$CODE )

parlist <- c("af")

pl2 <- paste( "af[1," , key$Site.number , "]" , sep="" )
zz <- precis( mfit , 3 , pars=parlist )

mean_return <- zz[pl2,"mean"]
sd_return <- zz[pl2,"sd"]

### Loading climate data ###
clim <- read.csv("eco_data.csv")
clim2 <- read.csv("eco_data2.csv")

EA_data <- jsonlite::fromJSON("EA.json", flatten = T)$features
EA_coord <- EA_data$geometry.coordinates

EA_matching <- data.frame(
  id = EA_data$id,
  lon = sapply(EA_coord, "[[", 1),
  lat = sapply(EA_coord, "[[", 2)
)

# Converting from 360 lon to -180-180 lon
EA_matching$lon <- ifelse(EA_matching$lon > 180, -360 + EA_matching$lon, EA_matching$lon)

## This loop will find the nearest EA population for each Koster population, minimizing geographic distance 
EA_match <- rep(NA, length(key$CODE))
for (i in 1:length(key$CODE)) {
  
    all_loc <- rbind(
      cbind(key[i,"Longitude"], key[i,"Latitude"]), cbind(EA_matching[,"lon"], EA_matching[,"lat"]) 
    )
    
    geo_dist <- distm(all_loc)
    rownames(geo_dist) <- c(key$CODE[i], as.character(EA_matching$id))
    diag(geo_dist) <- Inf
    EA_match[i] <- rownames(geo_dist)[which(geo_dist[1,] == min(geo_dist[1,]), arr.ind = TRUE)[1]]
  }

#########################################################################################################
precip_pred <- data.frame( precip = clim[clim$var_id == "PrecipitationPredictability", "code"],
                           id = clim[clim$var_id == "PrecipitationPredictability", "soc_id"] )

temp_pred <- data.frame( temp = clim[clim$var_id == "TemperaturePredictability", "code"],
                         id = clim[clim$var_id == "TemperaturePredictability", "soc_id"] )
                         
npp_pred <- data.frame( npp = clim2[clim2$var_id == "NetPrimaryProductionPredictability", "code"],
                        id = clim2[clim2$var_id == "NetPrimaryProductionPredictability", "soc_id"] )

d_clim <- merge(EA_matching, precip_pred, by.y="id")
d_clim <- merge(d_clim, temp_pred, by.y="id")
d_clim <- merge(d_clim, npp_pred, by.y="id")

d_clim <- d_clim[complete.cases(d_clim),]

### Giving each Koster pop the nearest EA pop value
precip_pred_s <- rep(NA, length(key$CODE))
temp_pred_s <- rep(NA, length(key$CODE))
npp_pred_s <- rep(NA, length(key$CODE))

for (i in 1:length(key$CODE)) {
  precip_pred_s[i] = precip_pred[precip_pred$id == EA_match[i],"precip"]
  temp_pred_s[i] = temp_pred[temp_pred$id == EA_match[i],"temp"]
  npp_pred_s[i] = npp_pred[npp_pred$id == EA_match[i],"npp"]
}

forager_data <- list(
  N = length(key$CODE),
  mean_return = mean_return,
  sd_return = sd_return,
  precip_pred_s = (precip_pred_s - mean(precip_pred_s)) / (sd(precip_pred_s)*2),
  temp_pred_s = (temp_pred_s - mean(temp_pred_s)) / (sd(temp_pred_s)*2),
  npp_pred_s = (npp_pred_s - mean(npp_pred_s)) / (sd(npp_pred_s)*2)
)

returns_model <- "
data{
int N; // Koster pops
real mean_return[N];
real sd_return[N];
real precip_pred_s[N];
real temp_pred_s[N];
real npp_pred_s[N];
}

parameters{
vector[N] returns; // true values, logit scale
vector[4] b;
real<lower=0> sigma_returns; // variance parameter
}

model{
// priors
b ~ normal(0,1);
sigma_returns ~ exponential(1);

for (i in 1:N) {
mean_return[i] ~ normal(returns[i], sd_return[i]); // incorporating measurement error
returns[i] ~ normal(b[1] + b[2]*precip_pred_s[i] + b[3]*temp_pred_s[i] + b[4]*npp_pred_s[i], sigma_returns);
}
}

"

fit_returns <- stan( model_code = returns_model , data=forager_data, chains=3, cores=3, iter=3000, init="0", control=list(adapt_delta=0.96) )
post_returns <- extract.samples(fit_returns)

cohens <- as.data.frame(post_returns$b[,2:4]) / post_returns$sigma_returns # cohen's d

length(post_returns$b[,2][post_returns$b[,2] > 0]) / length(post_returns$b[,2]) # posterior probability of pos effect of returns
length(post_returns$b[,3][post_returns$b[,3] > 0]) / length(post_returns$b[,3]) # posterior probability of pos effect on returns
length(post_returns$b[,4][post_returns$b[,4] > 0]) / length(post_returns$b[,4]) # posterior probability

pred_seq <- seq(from=-1,to=1, length.out=35)

precip <- matrix(0, nrow=250, ncol=length(pred_seq))
for (i in 1:250) {
  precip[i,] = logistic(post_returns$b[i,1] + post_returns$b[i,2]*pred_seq)
}

temp <- matrix(0, nrow=250, ncol=length(pred_seq))
for (i in 1:250) {
  
  temp[i,] = logistic(post_returns$b[i,1] + post_returns$b[i,3]*pred_seq)
}

npp <- matrix(0, nrow=250, ncol=length(pred_seq))
for (i in 1:250) {
  npp[i,] = logistic(post_returns$b[i,1] + post_returns$b[i,4]*pred_seq)
}

par(cex=1.6, cex.lab=1.6, mfrow=c(1,3), mar=c(3, 4.5, 2.5, 2), oma=c(1,1,1,0.5))
plot(NULL, xlab="", ylab="Pr (Hunting Success)", xlim=c(-1,1), ylim=c(0.5,1), xaxt='n', yaxt="n", xaxs="r")
axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.6)
axis(1, at=seq(-1,1,length.out = 5), labels=seq(-1,1,length.out = 5)*2, cex.axis=1.6)
mtext("Precip Pred. z-score", cex=1.4)
shade(apply(precip, 2, PI, 0.90), pred_seq, col=col.alpha("#57C0B3", 0.18))
shade(apply(precip, 2, PI, 0.72), pred_seq, col=col.alpha("#57C0B3", 0.18))
shade(apply(precip, 2, PI, 0.54), pred_seq, col=col.alpha("#57C0B3", 0.18))
shade(apply(precip, 2, PI, 0.36), pred_seq, col=col.alpha("#57C0B3", 0.18))
shade(apply(precip, 2, PI, 0.18), pred_seq, col=col.alpha("#57C0B3", 0.18))
lines(x=pred_seq, y=apply(precip,2,median), col=col.alpha("#57C0B3", 0.18), lwd=3)

plot(NULL, xlab="", ylab="Pr (Hunting Success)", xlim=c(-1,1), ylim=c(0.3,1), xaxt='n', yaxt="n", xaxs="r")
axis(2, at=c(0, 0.3, 1), labels=c("0", "0.5", "1"), cex.axis=1.6)
axis(1, at=seq(-1,1,length.out = 5), labels=seq(-1,1,length.out = 5)*2, cex.axis=1.6)
mtext("Temp Pred. z-score", cex=1.4)

shade(apply(temp, 2, PI, 0.90), pred_seq, col=col.alpha("#EDC0D7", 0.18))
shade(apply(temp, 2, PI, 0.72), pred_seq, col=col.alpha("#EDC0D7", 0.18))
shade(apply(temp, 2, PI, 0.54), pred_seq, col=col.alpha("#EDC0D7", 0.18))
shade(apply(temp, 2, PI, 0.36), pred_seq, col=col.alpha("#EDC0D7", 0.18))
shade(apply(temp, 2, PI, 0.18), pred_seq, col=col.alpha("#EDC0D7", 0.18))
lines(x=pred_seq, y=apply(temp,2,median), col=col.alpha("#EDC0D7", 0.18), lwd=3)

plot(NULL, xlab="", ylab="Pr (Hunting Success)", xlim=c(-1,1), ylim=c(0.5,1), xaxt='n', yaxt="n", xaxs="r")
axis(2, at=c(0, 0.5, 1), labels=c("0", "0.5", "1"), cex.axis=1.6)
axis(1, at=seq(-1,1,length.out = 5), labels=seq(-1,1,length.out = 5)*2, cex.axis=1.6)
mtext("NPP Pred. z-score", cex=1.4)
shade(apply(npp, 2, PI, 0.90), pred_seq, col=col.alpha("#E5C4BB", 0.18))
shade(apply(npp, 2, PI, 0.72), pred_seq, col=col.alpha("#E5C4BB", 0.18))
shade(apply(npp, 2, PI, 0.54), pred_seq, col=col.alpha("#E5C4BB", 0.18))
shade(apply(npp, 2, PI, 0.36), pred_seq, col=col.alpha("#E5C4BB", 0.18))
shade(apply(npp, 2, PI, 0.18), pred_seq, col=col.alpha("#E5C4BB", 0.18))
lines(x=pred_seq, y=apply(npp,2,median), col=col.alpha("#E5C4BB", 0.18), lwd=3)
