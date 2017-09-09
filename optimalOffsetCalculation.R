library("xlsx")

dat <- read.xlsx("~/Research/Mamre interaction for Bob.xlsx",sheetName="PRELIMINARY ANALYSIS")

# variance of snouts within mamre population comparable to average variance among populations,
# but not variance of mean traits among populations.
var(dat$TONGUE_LENGTH_EXTENDED)
snout_var
var(snout)

# same result for nectar tubes
var(dat$NECTAR_TUBE_LENGTH)
tubes_var
var(tubes)

# these metrics are highly correlated and volume produces higher resolution information
# so i'll use volume instead of height (height was reported in paper)
height <- dat$NECTAR_CONSUMED_.mm./dat$NECTAR_HEIGHT_BEFORE
volume <- dat$nectar_volume_consumed/dat$nectar_volume_before
cor(height,volume)

vrbl <- dat$NECTAR_TUBE_LENGTH-dat$TONGUE_LENGTH_EXTENDED
vrbl2 <- vrbl^2
pldp <- dat$POLLEN_GRAINS_DEPOSITED

opt.tongue.model <- lm(volume ~ vrbl + vrbl2)
coeffs <- -opt.tongue.model$coefficients
getOffset <- function(x){
  y <- as.numeric(coeffs[2]*x + coeffs[3]*x^2)
  return(y)
}
optimal_diff_for_pollinator <- optim(0,getOffset)$par

lm(pldp ~ vrbl + vrbl2)
# clearly there's an issue with the pollen deposited regression...
# so instead, i'll use the following as a stand-in for the optimal offset
max(vrbl) # but the resulting offset is associated with a small number of grains deposited
# instead it might make more sense to use the offset associated with the max number of grains deposited
vrbl[which(pldp==max(pldp))]
# note that these approximations are both very crude and lead to very different answers

# i'll try to fit a gaussian curve instead
gaussianfit <- function(c){
  sum((pldp - c[1]*exp(-(vrbl-c[2])^2/2*c[3]^2))^2)
}
sol <- optim( runif(3,0,1), gaussianfit )
for(i in 1:50){
  try <- optim( runif(3,0,1), gaussianfit )
  if(try$value < sol$value) sol <- try
}
optimal_diff_for_plant <- sol$par[2]

# what's interesting is looking at the differences in mean traits across populations
# typically flies have longer tongues than anceps have tubes
max(tubes-snout)