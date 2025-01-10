rm(list = ls())
pacman::p_load(data.table, fixest)
source("Code/R/Functions.R")

# ---- define the global variables ----
var_exo <- c("hp2wt", "air", "dpm", "size") # exogenous variables
var_end <- c("p") # endogenous variables
var_x <- c(var_exo, var_end) # all the rhs variables
var_iv_instr <- c("hp2wt_instr", "air_instr", "dpm_instr", "size_instr") # instruments
var_rand <- c("const") # the variables with random coefficients


# ---- the main estimation steps ----

# load the data
dt <- readRDS("Data/Out/carpanel_q4.rds")
# draw n=500 normal random draws to approximate the integral
n <- 500
draw <- rnorm(n)

# optim
setFixest_notes(FALSE) # suppress fixest notes (otherwise, too many notes will be printed in the console)
# esimate the non-linear parameters sigma and rho
sigma_opt <- optim(c(2, 0.5), blp_moment_condition, data = dt, var_iv_new = var_iv_instr, var_rand_coef = var_rand)
sigma <- sigma_opt$par # the non-linear parameters
# get the linear parameters
blp_result <- blp_intermediate(sigma, dt, var_iv_instr, var_rand)
ivreg <- blp_result[[2]]
summary(ivreg, cluster = "modelid")
etable(ivreg,
    cluster = "modelid", fitstat = ~ n + ar2 + ivf1 + ivf1.p + sargan + sargan.p,
    file = "Results/param_linear.tex", tex = TRUE, replace = TRUE
)
beta <- ivreg$coefficients # the linear parameters
delta <- blp_result[[1]]
mu <- blp_result[[3]]
save(dt, delta, sigma, ivreg, beta, mu, var_exo, var_end, var_rand, file = "Data/Out/blp_results.rda")

# joint estimation to get the standard errors with starting values from the previous estimation
beta <- as.numeric(beta)
theta_init <- c(sigma, beta)

theta_opt <- optim(theta_init, blp_moment_condition_full, data = dt, var_iv_new = var_iv_instr, var_rand_coef = var_rand, hessian = TRUE)
theta <- theta_opt$par # all the parameters
theta

# ---- calculate the standard errors ----

# need G W Z residual
pacman::p_load(numDeriv)

f <- function(x) {
    result <- blp_intermediate_full(theta = x, data = dt, var_iv_new = var_iv_instr, var_rand_coef = var_rand)[[2]]
    return(result)
}
G <- jacobian(f, theta_init)

blp_results_full <- blp_intermediate_full(theta_init, dt, var_iv_instr, var_rand)
g <- blp_results_full$g
Z <- blp_results_full$Z
W <- blp_results_full$W
residual <- blp_results_full$residual

blp_se <- gmm_se_hetero(G, Z = Z, W = W, residual = residual)
