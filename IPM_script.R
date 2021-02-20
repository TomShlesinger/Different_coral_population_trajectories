###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###----------------------- script for building Integral Projection Models (IPMs) ----------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


## This code accompanies the paper titled:
## Different population trajectories of two reef-building corals with similar life-history traits
## By Tom Shlesinger and Rob van Woesik
## Journal of Animal Ecology, 2021


## First, install and/or load the following packages:
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(MASS)) install.packages("MASS") 
if(!require(fields)) install.packages("fields") 
if(!require(IPMpack)) install.packages("IPMpack")
if(!require(glmmADMB)) install.packages("glmmADMB")



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###---------------------- Reading in the data and estimating reproductive values ----------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


setwd() #set your working directory

## Growth, survival, and recruitment data
data <- read.csv("Data.csv") 
## In the species column, "dip" == Dipsastraea favus, "plat" == Platygyra lamellina


## Reproduction data 
repro_mat <- read.csv("Repro_data.csv")
repro_mat

## Estimates of reproduction: Creating a simulated list of corals and their reproductive values 
## based on the size-grouped data in 'repro_mat'.

# For each size group in 'repro_mat', 'n' number of individuals will be created and their size will be 
# randomly sampled from a Gaussian distribution based on the group's mean diameter ('diameter_avg') and SD ('diameter_sd'). 
# According to the proportion of fertile colonies, 'colony_percent', within each size group a binary value is randomly assigned
# to 'reproductive' across a certain group's individuals to indicate whether they are fertile (='1') or infertile (='0').
# The percentage of fertile polyps will be also randomly sampled from a Gaussian distribution based on 'polyp_percent' and 'polyp_percent_sd'.

# Initializing an empty data frame to store the simulated individuals and their reproductive values  
repro_mat_sam <- data.frame(matrix(ncol = 4, nrow = sum(repro_mat$n))) 
colnames(repro_mat_sam) <- c("species", "size", "reproductive", "polyp_prop")

# Insert a list of species names according to the sample sizes in 'repro_mat':
repro_mat_sam$species <- c(rep("plat", sum(repro_mat$n[which(repro_mat$species=="plat")])), 
                          rep("dip", sum(repro_mat$n[which(repro_mat$species=="dip")])))

# Loop to create simulated population for the reproduction estimates: 
k <- 1 # Initializing a temporary variable to use in the second loop 
set.seed(1212)
for (i in 1:nrow(repro_mat)){  
  n <- repro_mat$n[i] # The sample size of the group in row 'i' 
  x <- repro_mat$n[i] * (repro_mat$colony_percent[i] / 100) # the number of fertile colonies in size group 'i'
  reproductive <- sample(c(rep(1, x), rep(0, n-x))) # A randomly shuffled list of length 'n', including '1' for reproductive colonies and '0' for non-reproductive colonies based on the number of fertile colonies 'x'.
  for (j in 1:n) {
    repro_mat_sam$size[k] <- round(rnorm(1, mean = repro_mat$diameter_avg[i], sd = repro_mat$diameter_sd[i]), 2) 
    repro_mat_sam$polyp_prop[k] <- round(rnorm(1, mean = repro_mat$polyp_percent[i], sd = repro_mat$polyp_percent_sd[i])/100, 2)
    repro_mat_sam$reproductive[k] <- reproductive[j] 
    k <- k + 1
  }   
}  

repro_mat_sam$log_size <- log(repro_mat_sam$size)

# Since the proportions cannot be larger than 1 or smaller than 0, correct the few proportions that fell outside this range:
repro_mat_sam$polyp_prop[which(repro_mat_sam$polyp_prop > 1)] <- 1 
repro_mat_sam$polyp_prop[which(repro_mat_sam$polyp_prop < 0)] <- 0 

## For estimates of the maximum potential of oocytes per coral: 
## Calculating the surface area of each coral based on its measured diameter and multiply it by the number of polyps per cm^2 and by the maximum number of oocytes per polyp 

# Calculate coral surface area assuming a symmetrical hemisphere morphology: 2 * pi * radius^2 
repro_mat_sam$surface_area <- 2 * pi * (repro_mat_sam$size * 0.5)^2 

# Multiplying coral surface area by the number of polyps per one centimeter.
# The number of polyps per cm will be randomly drawn from a Gaussian distribution based on the measurements (2.69 +/- 0.521 for Platygyra and 0.378 +/- 0.15 for Dipsastraea) 
# while limiting each randomly generated number of polyps to no more than two SD (mostly to avoid the possibility of Dipsastraea getting a few negative values).
# Max number of oocytes per polyp: 324 for Dipsastraea, 288 for Platygyra.

set.seed(1212)

repro_mat_sam_plat <- repro_mat_sam %>% filter(species == "plat") %>% mutate(polyps_per_cm = rnorm(n(), 2.69, 0.521))
repro_mat_sam_plat$polyps_per_cm[repro_mat_sam_plat$polyps_per_cm < 1.648] <- 1.648
repro_mat_sam_plat$polyps_per_cm[repro_mat_sam_plat$polyps_per_cm > 3.732] <- 3.732
repro_mat_sam_plat$oocytes_per_colony <- round(repro_mat_sam_plat$surface_area * repro_mat_sam_plat$polyps_per_cm) * 288

repro_mat_sam_dip <- repro_mat_sam %>% filter(species == "dip") %>% mutate(polyps_per_cm = rnorm(n(), 0.378, 0.15))
repro_mat_sam_dip$polyps_per_cm[repro_mat_sam_dip$polyps_per_cm < 0.078] <- 0.078
repro_mat_sam_dip$polyps_per_cm[repro_mat_sam_dip$polyps_per_cm > 0.678] <- 0.678
repro_mat_sam_dip$oocytes_per_colony <- round(repro_mat_sam_dip$surface_area * repro_mat_sam_dip$polyps_per_cm) * 324

repro_mat_sam <- rbind(repro_mat_sam_plat, repro_mat_sam_dip)



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###------------------------------------ Functions for vital rates -------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


# 1a. Survival function with size-independent mortality set to 99% (limiting survivorship to avoid "immortality")
surv.fun <- function(x, parameters) { 
  u <- exp(parameters$surv.int + parameters$surv.slope * x)
  return(0.99*u/(1+u))
}


# 2. Growth function - Calculating probability density function of size at t+3 (xp) from size at t (x)
growth.fun <- function(xp, x, parameters) { 			
  mu <- parameters$growth.int + parameters$growth.slope * x
  sig <- parameters$growth.sd
  p_den_growth <- dnorm(xp, mean = mu, sd = sig)
  return(p_den_growth)
}


# 3. Colony fertility  function
p.colony.x <- function(x, parameters) {
  u <- exp(parameters$colony.int + parameters$colony.slope * x)
  return(u/(1+u))
}


# 4. Polyp fertility function
p.polyp.x <- function(x, parameters) {
  u <- exp(parameters$polyp.int + parameters$polyp.slope * x)
  return(u/(1+u))
}


# 5. Reproduction function (including the probability of colonies being fertile, percentage of fertile polyps per colony,
#    potential number of oocytes, ratio of observed recruits at 2018 and potential oocytes since 2015, and recruits size)      
reproduction.fun <- function(xp, x, parameters) { 	
  repro <- p.colony.x(x, parameters) * p.polyp.x(x, parameters) * 
    exp(parameters$oocyte.int + parameters$oocyte.slope * x) *
    parameters$establishment.prob * 
    dnorm(xp, mean <- parameters$recruit.size.mean, sd <- parameters$recruit.size.sd)
  return(repro)
}


# Data frame for the model estimates used for the IPM 
parameters <- data.frame(
  colony.slope=NA,
  colony.int=NA,
  polyp.slope=NA,
  polyp.int=NA,
  oocyte.slope=NA,
  oocyte.int=NA,
  establishment.prob=NA,
  recruit.size.sd=NA,
  recruit.size.mean=NA,
  surv.slope=NA,
  surv.int=NA,
  growth.sd=NA,
  growth.slope=NA,
  growth.int=NA
  
)



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###----------------------------------------- Dipsastraea favus ----------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


# For Dipsastraea favus
growth_dip <- data %>% filter(species == "dip")
reproduction_dip <- repro_mat_sam %>% filter(species == "dip")
recruits_dip <- data %>% filter(species == "dip" & is.na(log_diam_2015) == T) 



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###----------------------------------- regressions for vital rates ------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


## 1. Growth regression 
growth_reg <- lm(log_diam_2018 ~ log_diam_2015, data = growth_dip)
summary(growth_reg)

# plot the regression 
mydata <- data.frame(log_diam_2015 = seq(min(growth_dip$log_diam_2015, na.rm = T), max(growth_dip$log_diam_2015, na.rm = T), length = 100))
P1 <- predict(growth_reg, newdata = mydata, type='response')
plot(x = growth_dip$log_diam_2015, y = growth_dip$log_diam_2018, xlab = "Size at 2015", ylab = "Size at 2018", main = "Growth")
lines(mydata$log_diam_2015, P1, lwd = 2, col = "red")

# Save the model estimates
parameters$growth.int <- coefficients(growth_reg)[1]
parameters$growth.slope <- coefficients(growth_reg)[2]
parameters$growth.sd <- summary(growth_reg)$sigma



## 2. survival regression
# Since in both species the intercept was not significantly different from zero, it might be better to force it to be 0.
survival_reg <- glm(survive ~ 0 + log_diam_2015, data=growth_dip, family=binomial())
summary(survival_reg)

# plot the regression 
par(mfrow=c(1,1))
P2 <- predict(survival_reg, newdata = mydata, type='response')
plot(x = growth_dip$log_diam_2015, y = growth_dip$survive, xlab = "Size at 2015", ylab = "Survival", main = "Survival")
lines(mydata$log_diam_2015, P2, lwd = 2, col = "red")

# Save the model estimates
parameters$surv.int <- 0
parameters$surv.slope <- coefficients(survival_reg)[1]



## 3. Colony reproductive percentage regression
colony_reg <- glm(reproductive ~ log_size, data=reproduction_dip, family=binomial())
summary(colony_reg)

# plot the regression 
mydata2 <- data.frame(log_size = seq(min(reproduction_dip$log_size, na.rm = T), max(reproduction_dip$log_size, na.rm = T), length = 100))
P3 <- predict(colony_reg, newdata = mydata2, type='response')
plot(x = reproduction_dip$log_size, y = reproduction_dip$reproductive, xlab = "Size at 2015", ylab = "Reproductive colonies", main = "Colony fertility")
lines(mydata2$log_size, P3, lwd = 2, col = "red")

# Save the model estimates
parameters$colony.int <- coefficients(colony_reg)[1]
parameters$colony.slope <- coefficients(colony_reg)[2]



## 4. Polyp reproductive percentage regression
polyp_reg <- glm(polyp_prop ~ log_size, data=reproduction_dip, family=binomial())
summary(polyp_reg)

# plot the regression 
P4 <- predict(polyp_reg, newdata = mydata2, type='response')
plot(x = reproduction_dip$log_size, y = reproduction_dip$polyp_prop, xlab = "Size at 2015", ylab = "Fertile polyps percentage", main = "Polyp fertility")
lines(mydata2$log_size, P4, lwd = 2, col = "red")

# Save the model estimates
parameters$polyp.int <- coefficients(polyp_reg)[1]
parameters$polyp.slope <- coefficients(polyp_reg)[2]



## 5a. oocytes regression
oocyte_reg <- glm.nb(oocytes_per_colony ~ log_size, data=reproduction_dip)
summary(oocyte_reg)

# plot the regression 
P5 <- predict(oocyte_reg, newdata = mydata2, type = "response")
plot(x = reproduction_dip$log_size, y = reproduction_dip$oocytes_per_colony, xlab = "Size at 2015", ylab = "Number of oocytes", main = "Oocytes per coral")
lines(mydata2$log_size, P5, lwd = 2, col = "red")

# Save the model estimates
parameters$oocyte.int <- coefficients(oocyte_reg)[1]
parameters$oocyte.slope <- coefficients(oocyte_reg)[2] 


## 5b. oocytes regression for estimates when oocytes are reduced to half
reproduction_dip$oocytes_per_colony2 <- reproduction_dip$oocytes_per_colony/2

oocyte_reg2 <- glm.nb(oocytes_per_colony2 ~ log_size, data=reproduction_dip)
summary(oocyte_reg2)

# plot the regression 
P5 <- predict(oocyte_reg2, newdata = mydata2, type = "response")
plot(x = reproduction_dip$log_size, y = reproduction_dip$oocytes_per_colony2, xlab = "Size at 2015", ylab = "Number of oocytes", main = "Oocytes per coral")
lines(mydata2$log_size, P5, lwd = 2, col = "red")

# Save the model estimates
oocyte.int_dip2 <- coefficients(oocyte_reg2)[1]
oocyte.slope_dip2 <- coefficients(oocyte_reg2)[2] 


# 6. size distribution of recruits
parameters$recruit.size.mean <- mean(recruits_dip$log_diam_2018)
parameters$recruit.size.sd <- sd(recruits_dip$log_diam_2018)


# 7. establishment probability, i.e., the ratio of recruits observed at 2018 and the total number of potential oocytes produced along 3 years 

# Multiplying coral surface area by a random number of polyps drawn from a Gaussian distribution based on the measurements 
# above while limiting each randomly generated number of polyps to no more than two SD (mostly to avoid the possibility of dipsastraea getting a few negative values).
# Max number of oocytes per polyp: 324 for Dipsastraea.

# Calculate coral surface area assuming a symmetrical hemisphere morphology: 2 * pi * radius^2 
growth_dip$surface_area <- 2 * pi * (exp(growth_dip$log_diam_2015) * 0.5)^2 

set.seed(1212)
total_ooc <- growth_dip %>% mutate(polyps_per_cm = rnorm(n(), 0.378, 0.15))
total_ooc$polyps_per_cm[total_ooc$polyps_per_cm < 0.078] <- 0.078
total_ooc$polyps_per_cm[total_ooc$polyps_per_cm > 0.678] <- 0.678
total_ooc$oocytes_per_colony <- round(total_ooc$surface_area * total_ooc$polyps_per_cm) * 324

parameters$establishment.prob <- length(recruits_dip$log_diam_2018) / sum(total_ooc$oocytes_per_colony * 3, na.rm=TRUE)

parameters_dip <- parameters



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###----------------------------------------------- IPMs -----------------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


# The integration limits (min.size and max.size) span the range of sizes observed in the data set, and then some.
# n - number of cells (mesh-points) in the IPM matrix
# h - width of the cells or "step size", i.e., the distance between mesh-points 
# y - mesh points; the centers of the cells defining the matrix and the points at which the matrix is evaluated for the midpoint rule of numerical integration

# max size for Dipsastrae favus = 50 cm (3.912023 on log scale), which is twice the size of the largest colony measured (24.51 cm)
min.size_dip <- 0.9 * min(c(growth_dip$log_diam_2015, growth_dip$log_diam_2018, recruits_dip$log_diam_2018), na.rm=T)
max.size_dip <- 3.912023 

n <- 300
h_dip <- (max.size_dip - min.size_dip) / n
y_dip <- min.size_dip + (1:n) * h_dip - h_dip/2

## IPM matrices and component kernels
# The function outer() evaluates the matrix at all pairwise combinations of the two vectors y and y and returns matrices representing 
# the kernel components for growth/survival and fecundity.
# For the numerical integration, we are using the midpoint rule to estimate the area under a curve. 
G_dip <- h_dip * outer(y_dip, y_dip, growth.fun, parameters = parameters_dip) 	# growth kernel
S_dip <- surv.fun(y_dip, parameters = parameters_dip) 							# survival kernel
P_dip <- G_dip * S_dip								# placeholder; redefine P on the next line
for(i in 1:n) {
  P_dip[,i] <- G_dip[,i] * S_dip[i]  		# growth/survival kernel
}
F_dip <- h_dip * outer(y_dip, y_dip, reproduction.fun, parameters = parameters_dip) 	# reproduction/recruitment kernel
K_dip <- P_dip + F_dip 															                                  # full IPM kernel



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###------------------------------------------ basic analyses ------------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


# Obtaining lambda value and sensitivity and elasticity matrices 
lam_dip <- Re(eigen(K_dip)$values[1]) # The dominant eigenvalue gives the asymptotic population growth rate (lamda). 
lam_dip

sensMatrix_dip<-sens(K_dip)
elasMatrix_dip<-elas(K_dip)


## Partitioned elasticities
elasPMatrix_dip <- P_dip * sensMatrix_dip / lam_dip
elasFMatrix_dip <- F_dip * sensMatrix_dip / lam_dip

# Relative contribution of the survival/growth sub-kernel vs. the reproduction/recruitment sub-kernel 
c(sum(elasMatrix_dip), sum(elasPMatrix_dip),sum(elasFMatrix_dip))


## Perturbation analysis:
# Assessing the sensitivity of each component model parameter individually, by adjusting each coefficient 
# by +/-1%, and see how it affect the population growth rate (lambda)
nparams <- length(parameters_dip)       # total number of component model params
lam.sn_dip <- rep(NA, nparams*2)
all.params.orig <- parameters_dip

for (i in 1:nparams) {         # tweak each parameter in turn
  for (j in 1:2) {             # by 1% down and then 1% up
    all.params <- all.params.orig
    all.params[i] <- all.params[i]*(1+2*(j-1.5)/100)
    G2 <- h_dip * outer(y_dip, y_dip, growth.fun, parameters = all.params) 	# growth kernel
    S2 <- surv.fun(y_dip, parameters = all.params) 							# survival 
    P2 <- G2 * S2								# placeholder; redefine P on the next line
    for(z in 1:n) {
      P2[,z] <- G2[,z]*S2[z]  		# growth/survival kernel
    }
    F2 <- h_dip * outer(y_dip, y_dip, reproduction.fun, parameters = all.params) 	# reproduction kernel
    K2 <- P2 + F2 	
    lam_new <- Re(eigen(K2)$values[1])
    lam.sn_dip[2*(i-1)+j] <- lam_new
  }
}

## Since the survival intercept is zero, the above loop does not work on it so perturbating manually for -/+ 0.01:
all.params <- all.params.orig
all.params[11] <-  -0.01
G3 <- h_dip * outer(y_dip, y_dip, growth.fun, parameters = all.params) 	# growth kernel
S3 <- surv.fun(y_dip, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
    P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_dip * outer(y_dip, y_dip, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new <- Re(eigen(K3)$values[1])
lam.sn_dip[21] <- lam_new

all.params <- all.params.orig
all.params[11] <-  0.01
G3 <- h_dip * outer(y_dip, y_dip, growth.fun, parameters = all.params) 	# growth kernel
S3 <- surv.fun(y_dip, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
  P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_dip * outer(y_dip, y_dip, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new <- Re(eigen(K3)$values[1])
lam.sn_dip[22] <- lam_new
lam.sn_dip



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###-------------------------------------- Population projections --------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


## Projections of 500 simulated populations based on the value of lambda from the IPM and the abundance at 2018

#Perturbating all IPM componenents by +/-1% to estimate variability in lambda
all.params <- parameters_dip * 1.01
G3 <- h_dip * outer(y_dip, y_dip, growth.fun, parameters = all.params) 	# growth kernel
S3 <- surv.fun(y_dip, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
  P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_dip * outer(y_dip, y_dip, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new[1] <- Re(eigen(K3)$values[1])

all.params <- parameters_dip * 0.99
G3 <- h_dip * outer(y_dip, y_dip, growth.fun, parameters = all.params) 	# growth kernel
S3 <- surv.fun(y_dip, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
  P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_dip * outer(y_dip, y_dip, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new[2] <- Re(eigen(K3)$values[1])

lamVar_dip <- sum(abs(lam_dip - lam_new)) 

## Loop to project lambda through time 
n.years <- 16                                              # Number of time steps (including the starting population)
popstart_dip <- sum(!is.na(growth_dip$log_diam_2018))          # Initial population size (at 2018)
mean.lambda <- lam_dip                                     # Population growth rate (i.e., lambda)
sigma.lambda <- lamVar_dip                                 # Variation of the growth rate
n2 <- 500                                                 # Number of simulated populations
sim_pop <- data.frame(matrix(ncol = n.years, nrow = n2))   # Initializing an empty data frame to store the simulated populations
colnames(sim_pop) <- c(seq(2018,2063,3))

set.seed(1212)
for (u in 1:n2){
  lambda_rand <- rnorm(n.years-1, mean.lambda, sigma.lambda)
  pop_dip <- numeric(n.years)
  pop_dip[1] <- popstart_dip
  for (t in 1:(n.years-1)){
    pop_dip[t+1] <- pop_dip[t] * lambda_rand[t]
  }
  sim_pop[u,] <- pop_dip   
}

sim_pop$simulation <- seq(1,n2)
long_data_dip1 <- sim_pop %>% gather(year, abundance, 1:16) 
long_data_dip1$simulation <- as.factor(long_data_dip1$simulation)
long_data_dip1$year <- as.integer(long_data_dip1$year)


## Projections of 500 simulated populations with halved fecundity
parameters_dip2 <- parameters_dip
parameters_dip2$oocyte.slope <- oocyte.slope_dip2
parameters_dip2$oocyte.int <- oocyte.int_dip2

G_dip2 <- h_dip * outer(y_dip, y_dip, growth.fun, parameters = parameters_dip2) 	# growth kernel
S_dip2 <- surv.fun(y_dip, parameters = parameters_dip2) 							# survival kernel
P_dip2 <- G_dip2 * S_dip2								# placeholder; redefine P on the next line
for(i in 1:n) {
  P_dip2[,i] <- G_dip2[,i] * S_dip2[i]  		# growth/survival kernel
}
F_dip2 <- h_dip * outer(y_dip, y_dip, reproduction.fun, parameters = parameters_dip2) 	# reproduction kernel
K_dip2 <- P_dip2 + F_dip2 		
lam_dip2 <- Re(eigen(K_dip2)$values[1])
lam_dip2

#Perturbating all IPM componenents by +/-1% to estimate variability in lambda
all.params <- parameters_dip2 * 1.01
G3 <- h_dip * outer(y_dip, y_dip, growth.fun, parameters = all.params) 	# growth kernel
S3 <- surv.fun(y_dip, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
  P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_dip * outer(y_dip, y_dip, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new[1] <- Re(eigen(K3)$values[1])

all.params <- parameters_dip2 * 0.99
G3 <- h_dip * outer(y_dip, y_dip, growth.fun, all.params) 	# growth kernel
S3 <- surv.fun(y_dip, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
  P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_dip * outer(y_dip, y_dip, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new[2] <- Re(eigen(K3)$values[1])

lamVar_dip2 <- sum(abs(lam_dip2 - lam_new)) 

## Loop to project lambda through time 
n.years <- 16                                              # Number of time steps (including the starting population)
popstart_dip <- sum(!is.na(growth_dip$log_diam_2018))      # Initial population size (at 2018)
mean.lambda <- lam_dip2                                    # Population growth rate (i.e., lambda)
sigma.lambda <- lamVar_dip2                                # Variation of the growth rate
n2 <- 500                                                  # Number of simulated populations
sim_pop <- data.frame(matrix(ncol = n.years, nrow = n2))   # Initializing an empty data frame to store the simulated populations
colnames(sim_pop) <- c(seq(2018,2063,3))

set.seed(1212)
for (u in 1:n2){
  lambda_rand <- rnorm(n.years-1, mean.lambda, sigma.lambda)
  pop_dip <- numeric(n.years)
  pop_dip[1] <- popstart_dip
  for (t in 1:(n.years-1)){
    pop_dip[t+1] <- pop_dip[t] * lambda_rand[t]
  }
  sim_pop[u,] <- pop_dip   
}

sim_pop$simulation <- seq(n2+1, n2*2)
long_data_dip2 <- sim_pop %>% gather(year, abundance, 1:16) 
long_data_dip2$simulation <- as.factor(long_data_dip2$simulation)
long_data_dip2$year <- as.integer(long_data_dip2$year)

long_data_dip1$fecundity <- "full"
long_data_dip2$fecundity <- "half"

long_data_dip <- rbind(long_data_dip1, long_data_dip2)
long_data_dip$abundance <- round(long_data_dip$abundance)

## Assessing differences between projections using a GLMM 
projections_glmm <- glmmadmb(abundance ~ fecundity, random = ~ 1|simulation, data = long_data_dip, family = "nbinom")
summary(projections_glmm)



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###---------------------------------------- Platygyra lamellina ---------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


# For Platygyra lamellina
growth_plat <- data %>% filter(species == "plat")
reproduction_plat <- repro_mat_sam %>% filter(species == "plat")
recruits_plat <- data %>% filter(species == "plat" & is.na(log_diam_2015) == T) 



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###----------------------------------- regressions for vital rates ------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


## 1. Growth regression 
growth_reg <- lm(log_diam_2018 ~ log_diam_2015, data = growth_plat)
summary(growth_reg)

# plot the regression 
mydata <- data.frame(log_diam_2015 = seq(min(growth_plat$log_diam_2015, na.rm = T), max(growth_plat$log_diam_2015, na.rm = T), length = 100))
P1 <- predict(growth_reg, newdata = mydata, type='response')
plot(x = growth_plat$log_diam_2015, y = growth_plat$log_diam_2018, xlab = "Size at 2015", ylab = "Size at 2018", main = "Growth")
lines(mydata$log_diam_2015, P1, lwd = 2, col = "red")

# Save the model estimates
parameters$growth.int <- coefficients(growth_reg)[1]
parameters$growth.slope <- coefficients(growth_reg)[2]
parameters$growth.sd <- summary(growth_reg)$sigma



## 2. survival regression
# Since in both species the intercept was not significantly different from zero, it might be better to force it to be 0.
survival_reg <- glm(survive ~ 0 + log_diam_2015, data=growth_plat, family=binomial())
summary(survival_reg)

# plot the regression 
par(mfrow=c(1,1))
P2 <- predict(survival_reg, newdata = mydata, type='response')
plot(x = growth_plat$log_diam_2015, y = growth_plat$survive, xlab = "Size at 2015", ylab = "Survival", main = "Survival")
lines(mydata$log_diam_2015, P2, lwd = 2, col = "red")

# Save the model estimates
parameters$surv.int <- 0
parameters$surv.slope <- coefficients(survival_reg)[1]



## 3. Colony reproductive percentage regression
colony_reg <- glm(reproductive ~ log_size, data=reproduction_plat, family=binomial())
summary(colony_reg)

# plot the regression 
mydata2 <- data.frame(log_size = seq(min(reproduction_plat$log_size, na.rm = T), max(reproduction_plat$log_size, na.rm = T), length = 100))
P3 <- predict(colony_reg, newdata = mydata2, type='response')
plot(x = reproduction_plat$log_size, y = reproduction_plat$reproductive, xlab = "Size at 2015", ylab = "Reproductive colonies", main = "Colony fertility")
lines(mydata2$log_size, P3, lwd = 2, col = "red")

# Save the model estimates
parameters$colony.int <- coefficients(colony_reg)[1]
parameters$colony.slope <- coefficients(colony_reg)[2]



## 4. Polyp reproductive percentage regression
polyp_reg <- glm(polyp_prop ~ log_size, data=reproduction_plat, family=binomial())
summary(polyp_reg)

# plot the regression 
P4 <- predict(polyp_reg, newdata = mydata2, type='response')
plot(x = reproduction_plat$log_size, y = reproduction_plat$polyp_prop, xlab = "Size at 2015", ylab = "Fertile polyps percentage", main = "Polyp fertility")
lines(mydata2$log_size, P4, lwd = 2, col = "red")

# Save the model estimates
parameters$polyp.int <- coefficients(polyp_reg)[1]
parameters$polyp.slope <- coefficients(polyp_reg)[2]



## 5a. oocytes regression
oocyte_reg <- glm.nb(oocytes_per_colony ~ log_size, data=reproduction_plat)
summary(oocyte_reg)

# plot the regression 
P5 <- predict(oocyte_reg, newdata = mydata2, type = "response")
plot(x = reproduction_plat$log_size, y = reproduction_plat$oocytes_per_colony, xlab = "Size at 2015", ylab = "Number of oocytes", main = "Oocytes per coral")
lines(mydata2$log_size, P5, lwd = 2, col = "red")

# Save the model estimates
parameters$oocyte.int <- coefficients(oocyte_reg)[1]
parameters$oocyte.slope <- coefficients(oocyte_reg)[2] 



## 5b. oocytes regression for estimates when oocytes are reduced to half
reproduction_plat$oocytes_per_colony2 <- reproduction_plat$oocytes_per_colony/2

oocyte_reg2 <- glm.nb(oocytes_per_colony2 ~ log_size, data=reproduction_plat)
summary(oocyte_reg2)

# plot the regression 
P5 <- predict(oocyte_reg2, newdata = mydata2, type = "response")
plot(x = reproduction_plat$log_size, y = reproduction_plat$oocytes_per_colony2, xlab = "Size at 2015", ylab = "Number of oocytes", main = "Oocytes per coral")
lines(mydata2$log_size, P5, lwd = 2, col = "red")

# Save the model estimates
oocyte.int_plat2 <- coefficients(oocyte_reg2)[1]
oocyte.slope_plat2 <- coefficients(oocyte_reg2)[2] 



# 6. size distribution of recruits
parameters$recruit.size.mean <- mean(recruits_plat$log_diam_2018)
parameters$recruit.size.sd <- sd(recruits_plat$log_diam_2018)


# 7. establishment probability, i.e., the ratio of recruits observed at 2018 and the total number of potential oocytes produced along 3 years 

# Multiplying coral surface area by a random number of polyps drawn from a Gaussian distribution based on the measurements 
# above while limiting each randomly generated number of polyps to no more than two SD (mostly to avoid the possibility of dipsastraea getting a few negative values).
# Max number of oocytes per polyp: 288 for Platygyra. 

# Calculate coral surface area assuming a symmetrical hemisphere morphology: 2 * pi * radius^2 
growth_plat$surface_area <- 2 * pi * (exp(growth_plat$log_diam_2015) * 0.5)^2 

set.seed(1212)

total_ooc <- growth_plat %>% mutate(polyps_per_cm = rnorm(n(), 2.69, 0.521))
total_ooc$polyps_per_cm[total_ooc$polyps_per_cm < 1.648] <- 1.648
total_ooc$polyps_per_cm[total_ooc$polyps_per_cm > 3.732] <- 3.732
total_ooc$oocytes_per_colony <- round(total_ooc$surface_area * total_ooc$polyps_per_cm) * 288

parameters$establishment.prob <- length(recruits_plat$log_diam_2018) / sum(total_ooc$oocytes_per_colony * 3, na.rm=TRUE)

parameters_plat <- parameters



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###----------------------------------------------- IPMs -----------------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


# The integration limits (min.size and max.size) span the range of sizes observed in the data set, and then some.
# n - number of cells (mesh-points) in the IPM matrix
# h - width of the cells or "step size", i.e., the distance between mesh-points 
# y - mesh points; the centers of the cells defining the matrix and the points at which the matrix is evaluated for the midpoint rule of numerical integration

# max size for platygyra = 100 cm (4.60517 on log scale)
min.size_plat <- 0.9 * min(c(growth_plat$log_diam_2015, growth_plat$log_diam_2018, recruits_plat$log_diam_2018), na.rm=T)
max.size_plat <- 4.60517 

n <- 300
h_plat <- (max.size_plat - min.size_plat) / n
y_plat <- min.size_plat + (1:n) * h_plat - h_plat/2

## IPM matrices and component kernels
# The function outer() evaluates the matrix at all pairwise combinations of the two vectors y and y and returns matrices representing 
# the kernel components for growth/survival and fecundity.
# For the numerical integration, we are using the midpoint rule to estimate the area under a curve. 

G_plat <- h_plat * outer(y_plat, y_plat, growth.fun, parameters = parameters_plat) 	# growth kernel
S_plat <- surv.fun(y_plat, parameters = parameters_plat) 							# survival 
P_plat <- G_plat * S_plat								# placeholder; redefine P on the next line
for(i in 1:n) {
  P_plat[,i] <- G_plat[,i] * S_plat[i]  		# growth/survival kernel
}
F_plat <- h_plat * outer(y_plat, y_plat, reproduction.fun, parameters = parameters_plat) 	# reproduction/recruitment kernel
K_plat <- P_plat + F_plat 															                                  # full IPM kernel



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###------------------------------------------ basic analyses ------------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


# Obtaining lambda value and sensitivity and elasticity matrices 
lam_plat <- Re(eigen(K_plat)$values[1]) # The dominant eigenvalue gives the asymptotic population growth rate (lamda).
lam_plat

sensMatrix_plat <- sens(K_plat)
elasMatrix_plat <- elas(K_plat)

## Partitioned elasticities
elasPMatrix_plat <- P_plat * sensMatrix_plat / lam_plat
elasFMatrix_plat <- F_plat * sensMatrix_plat / lam_plat

# Relative contribution of the survival/growth sub-kernel vs. the reproduction/recruitment sub-kernel 
c(sum(elasMatrix_plat), sum(elasPMatrix_plat),sum(elasFMatrix_plat))


## Perturbation analysis:
# Assessing the sensitivity of each component model parameter individually, by adjusting each coefficient 
# by +/-1%, and see how it effect the population growth rate (lambda)
nparams <- length(parameters_plat)       # total number of component model params
lam.sn_plat <- rep(NA, nparams*2)
all.params.orig <- parameters_plat

for (i in 1:nparams) {         # tweak each param in turn
  for (j in 1:2) {             # by 1% down and then 1% up
    all.params <- all.params.orig
    all.params[i] <- all.params[i]*(1+2*(j-1.5)/100)
    G2 <- h_plat * outer(y_plat, y_plat, growth.fun, parameters = all.params) 	# growth kernel
    S2 <- surv.fun(y_plat, parameters = all.params) 							# survival 
    P2 <- G2 * S2								# placeholder; redefine P on the next line
    for(z in 1:n) {
      P2[,z] <- G2[,z]*S2[z]  		# growth/survival kernel
    }
    F2 <- h_plat * outer(y_plat, y_plat, reproduction.fun, parameters = all.params) 	# reproduction kernel
    K2 <- P2 + F2 	
    lam_new <- Re(eigen(K2)$values[1])
    lam.sn_plat[2*(i-1)+j] <- lam_new
  }
}

## Since the survival intercept is zero, the above loop does not work on it so perturbating manually for -/+ 0.01
all.params <- all.params.orig
all.params[11] <-  -0.01
G3 <- h_plat * outer(y_plat, y_plat, growth.fun, parameters = all.params) 	# growth kernel
S3 <- surv.fun(y_plat, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
  P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_plat * outer(y_plat, y_plat, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new <- Re(eigen(K3)$values[1])
lam.sn_plat[21] <- lam_new

all.params <- all.params.orig
all.params[11] <-  0.01
G3 <- h_plat * outer(y_plat, y_plat, growth.fun, parameters = all.params) 	# growth kernel
S3 <- surv.fun(y_plat, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
  P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_plat * outer(y_plat, y_plat, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new <- Re(eigen(K3)$values[1])
lam.sn_plat[22] <- lam_new
lam.sn_plat


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
###-------------------------------------- Population projections --------------------------------------
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###


## Projections of 500 simulated populations based on the value of lambda from the IPM and the abundance at 2018

#Perturbating all IPM componenents by +/-1% to estimate variability in lambda
all.params <- parameters_plat * 1.01
G3 <- h_plat * outer(y_plat, y_plat, growth.fun, parameters = all.params) 	# growth kernel
S3 <- surv.fun(y_plat, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
  P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_plat * outer(y_plat, y_plat, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new[1] <- Re(eigen(K3)$values[1])

all.params <- parameters * 0.99
G3 <- h_plat * outer(y_plat, y_plat, growth.fun, parameters = all.params) 	# growth kernel
S3 <- surv.fun(y_plat, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
  P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_plat * outer(y_plat, y_plat, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new[2] <- Re(eigen(K3)$values[1])

lamVar_plat <- sum(abs(lam_plat - lam_new)) 

## Loop to project lambda through time 
n.years <- 16                                              # Number of time steps (including the starting population)
popstart_plat <- sum(!is.na(growth_plat$log_diam_2018))    # Initial population size (at 2018)
mean.lambda <- lam_plat                                    # Population growth rate (i.e., lambda)
sigma.lambda <- lamVar_plat                                # Variation of the growth rate
n2 <- 500                                                  # Number of simulated populations
sim_pop <- data.frame(matrix(ncol = n.years, nrow = n2))   # Initializing an empty data frame to store the simulated populations
colnames(sim_pop) <- c(seq(2018,2063,3))

set.seed(1212)
for (u in 1:n2){
  lambda_rand <- rnorm(n.years-1, mean.lambda, sigma.lambda)
  pop_plat <- numeric(n.years)
  pop_plat[1] <- popstart_plat
  for (t in 1:(n.years-1)){
    pop_plat[t+1] <- pop_plat[t] * lambda_rand[t]
  }
  sim_pop[u,] <- pop_plat   
}

sim_pop$simulation <- seq(1,n2)
long_data_plat1 <- sim_pop %>% gather(year, abundance, 1:16) 
long_data_plat1$simulation <- as.factor(long_data_plat1$simulation)
long_data_plat1$year <- as.integer(long_data_plat1$year)



## Projections of 500 simulated populations with halved fecundity
parameters_plat2 <- parameters_plat
parameters_plat2$oocyte.slope <- oocyte.slope_plat2
parameters_plat2$oocyte.int <- oocyte.int_plat2

G_plat2 <- h_plat * outer(y_plat, y_plat, growth.fun, parameters = parameters_plat2) 	# growth kernel
S_plat2 <- surv.fun(y_plat, parameters = parameters_plat2) 							# survival kernel
P_plat2 <- G_plat2 * S_plat2								# placeholder; redefine P on the next line
for(i in 1:n) {
  P_plat2[,i] <- G_plat2[,i] * S_plat2[i]  		# growth/survival kernel
}
F_plat2 <- h_plat * outer(y_plat, y_plat, reproduction.fun, parameters = parameters_plat2) 	# reproduction kernel
K_plat2 <- P_plat2 + F_plat2 		
lam_plat2 <- Re(eigen(K_plat2)$values[1])
lam_plat2

#Perturbating all IPM componenents by +/-1% to estimate variability in lambda
all.params <- parameters_plat2 * 1.01
G3 <- h_plat * outer(y_plat, y_plat, growth.fun, parameters = all.params) 	# growth kernel
S3 <- surv.fun(y_plat, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
  P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_plat * outer(y_plat, y_plat, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new[1] <- Re(eigen(K3)$values[1])

all.params <- parameters_plat2 * 0.99
G3 <- h_plat * outer(y_plat, y_plat, growth.fun, parameters = all.params) 	# growth kernel
S3 <- surv.fun(y_plat, parameters = all.params) 							# survival 
P3 <- G3 * S3								# placeholder; redefine P on the next line
for(z in 1:n) {
  P3[,z] <- G3[,z]*S3[z]  		# growth/survival kernel
}
F3 <- h_plat * outer(y_plat, y_plat, reproduction.fun, parameters = all.params) 	# reproduction kernel
K3 <- P3+F3 	
lam_new[2] <- Re(eigen(K3)$values[1])

lamVar_plat2 <- sum(abs(lam_plat2 - lam_new)) 

## Loop to project lambda through time 
n.years <- 16                                              # Number of time steps (including the starting population)
popstart_plat <- sum(!is.na(growth_plat$log_diam_2018))    # Initial population size (at 2018)
mean.lambda <- lam_plat2                                   # Population growth rate (i.e., lambda)
sigma.lambda <- lamVar_plat2                               # Variation of the growth rate
n2 <- 500                                                  # Number of simulated populations
sim_pop <- data.frame(matrix(ncol = n.years, nrow = n2))   # Initializing an empty data frame to store the simulated populations
colnames(sim_pop) <- c(seq(2018,2063,3))

set.seed(1212)
for (u in 1:n2){
  lambda_rand <- rnorm(n.years-1, mean.lambda, sigma.lambda)
  pop_plat <- numeric(n.years)
  pop_plat[1] <- popstart_plat
  for (t in 1:(n.years-1)){
    pop_plat[t+1] <- pop_plat[t] * lambda_rand[t]
  }
  sim_pop[u,] <- pop_plat   
}

sim_pop$simulation <- seq(n2+1, n2*2)
long_data_plat2 <- sim_pop %>% gather(year, abundance, 1:16) 
long_data_plat2$simulation <- as.factor(long_data_plat2$simulation)
long_data_plat2$year <- as.integer(long_data_plat2$year)

long_data_plat1$fecundity <- "full"
long_data_plat2$fecundity <- "half"

long_data_plat <- rbind(long_data_plat1, long_data_plat2)
long_data_plat$abundance <- round(long_data_plat$abundance)

## Assessing differences between projections using a GLMM 
projections_glmm <- glmmadmb(abundance ~ fecundity, random = ~ 1|simulation, data = long_data_plat, family = "nbinom")
summary(projections_glmm)



###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
