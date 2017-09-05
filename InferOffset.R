library(gdata)

# read in data
offset.data <- read.xls("PlantPollinatorOffset.xlsx")

# Lets begin with inferring the offset for the fly

# We start by pulling out the relevant data from offset.data
offset.differences <- offset.data$TONGUE.TUBE
offset.differences.squared <- offset.differences^2
percent.nectar.consumed <- offset.data$NECTAR_VOLUME_CONSUMED/offset.data$NECTAR_VOLUME_BEFORE

# Perform a quadratic regression of trait differences onto percent nectar consumed
opt.tongue.model <- lm(percent.nectar.consumed ~ offset.differences + offset.differences.squared)

plot

# Numerically solve for the apex of the resulting parabola
coeffs <- -opt.tongue.model$coefficients
getOffset <- function(x){
  y <- as.numeric(coeffs[2]*x + coeffs[3]*x^2)
  return(y)
}
optimal_diff_for_pollinator <- optim(0,getOffset)$par


# Using 

gaussianfit <- function(c){
  sum((off$P - c[1]*exp(-(off$D-c[2])^2/2*c[3]^2))^2)
}
sol <- optim( runif(3,0,1), gaussianfit )
for(i in 1:50){
  try <- optim( runif(3,0,1), gaussianfit )
  if(try$value < sol$value) sol <- try
}
optimal_diff_for_plant <- -sol$par[2]