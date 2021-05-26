# Example structure of a statistical model
# Gergana Daskalova
# 4th April 2019
# gndaskalova@gmail.com

library(tidyr)
library(dplyr)
library(brms)

options(mc.cores = parallel::detectCores())

# Population change and cumulative driver intensity
# Note that the model structures will be similar for the remaining models as well,
# with differences in the fixed effect of interest

# See full pre-registration for the structure of each model

#	Set weakly-regularising priors
prior1 <- c(set_prior(prior = 'normal(0, 6)', class = 'b', coef = 'duration'),
                       set_prior(prior = 'normal(0, 6)', class = 'b', coef = 'cumulative'), 	# global slope
                       set_prior(prior = 'normal(0, 6)', class = 'Intercept', coef = ''), 		# global intercept
                       set_prior(prior = 'cauchy(0, 2)', class = 'sd'))								# group-level intercepts

pop_cumulative_model <- brm(bf(value ~ cumulative * latitudinal_band + 
                                 duration + (1 + taxa + study_id |taxa / study_id)), 
                   data = popbio[popbio$type == "Population",], 
                   prior = prior1, iter = 6000,
                   warmup = 2000,
                   inits = '0',
                   control = list(adapt_delta = 0.95),
                   cores = 3, chains = 3)

# Check model and save outputs
summary(pop_cumulative_model)
plot(pop_cumulative_model)
save(pop_cumulative_model, file = "data/output/pop_cumulative_model.RData")

# Note that we will also explore using more informative priors 
# so that we carry through the uncertainty around the population/biodiversity trends
# into our second stage analysis of trends ~ driver intensity

# We will also explore if it's possible to account for correlation in the random effects
# in the prior structure

# Determining whether such priors will allow model convergence with our data is only
# possible once we have started the analyses
# thus we have not done this before the pre-registration