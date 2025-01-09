rm(list = ls())
pacman::p_load(data.table, fixest, gmm)

var_exo <- c("hp2wt", "air", "dpm", "size")
var_end <- c("p")
var_rand <- c("const") # the random coefficient can be on the const term as well.

# the main part

source("Code/R/Functions.R")
# load the data
dt <- readRDS("Data/Out/carpanel_q4.rds")
# draw 500 normal random draws to approximate the integral
n <- 500
# draw in the normal(0, 1)
draw <- rnorm(n)
# check the variables
var_x <- c(var_exo, var_end)
# var_iv_blp <- c("const_rival", "const_ownothers", "dpm_rival", "dpm_ownothers", "hp2wt_rival", "hp2wt_ownothers", "size_rival", "size_ownothers", "air_rival", "air_ownothers")
# var_iv_nl <- c("const_rival", "const_ownothers", "dpm_rival", "dpm_ownothers", "hp2wt_rival", "hp2wt_ownothers", "air_rival", "air_ownothers", "const_rival_g", "hp2wt_rival_g", "dpm_rival_g", "air_rival_g")
# var_iv_local <- c("size_local", "wb_local", "hp_local", "hp2wt_local", "disp_local")
# var_iv_quad <- c("wb_quad", "hp_quad", "size_quad", "hp2wt_local", "disp_quad")
var_iv_instr <- c("hp2wt_instr", "air_instr", "dpm_instr", "size_instr")

# optim
setFixest_notes(FALSE) # suppress fixest notes
sigma_opt <- optim(c(2, 0.5), blp_moment_condition, data = dt, var_iv_new = var_iv_instr, var_rand_coef = var_rand)
sigma <- sigma_opt$par
blp_result <- blp_intermediate(sigma, dt, var_iv_instr, var_rand)
ivreg <- blp_result[[2]]
summary(ivreg, cluster = "modelid")
etable(ivreg,
    cluster = "modelid", fitstat = ~ n + ar2 + ivf1 + ivf1.p + sargan + sargan.p,
    file = "Results/param_linear.tex", tex = TRUE, replace = TRUE
)
beta <- ivreg$coefficients
delta <- blp_result[[1]]
mu <- blp_result[[3]]
save(dt, delta, sigma, ivreg, beta, mu, var_exo, var_end, var_rand, file = "Data/Out/blp_results.rda")

# end


# joint estimation to get the standard errors with starting values from the previous estimation
beta <- as.numeric(beta)
theta_init <- c(sigma, beta)

# theta <- theta_init
# data <- dt
# var_iv_new <- var_iv_instr
# var_rand_coef <- var_rand
# theta[3:length(theta)]

theta_opt <- optim(theta_init, blp_moment_condition_full, data = dt, var_iv_new = var_iv_instr, var_rand_coef = var_rand)
theta <- theta_opt$par

# to get the standard errors
# the following doesn't work:)
# hessian <- theta_opt$hessian
# se <- sqrt(diag(solve(hessian)))
