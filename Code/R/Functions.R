# In this file, I define the functions used in the estimation

# given delta_j and simulated draws, calculate the market share
share <- function(delta, mu) {
    # J is the number of products, n is the number of draws
    # delta is of J x 1, mu is of J x n
    n <- ncol(mu) # get the number of draws
    numer <- exp(delta + mu) # for each product and each draw, calculate the numerator, dim(numer) = J x n
    denom <- 1 + colSums(numer) # for each draw, sum up the products, dim(denom) = 1 x n
    sij <- numer / denom # dim(sij) = J x n
    sj <- rowSums(sij) / n # for each product, sum up the draws and divide by n, dim(sj) = J x 1
    return(sj)
}

# given the observed market share sj, back out delta_j, and the number of iterations
blp_contraction <- function(sj, mu, delta_init) {
    # dim(sj) = J x 1, dim(mu) = J x n, dim(delta_init) = J x 1
    J <- length(delta_init) # get the number of products
    delta_before <- rep(0, J)
    delta_next <- delta_init
    iter <- 0
    while (max(abs(delta_next - delta_before)) > 1e-7 && iter < 500) {
        delta_before <- delta_next
        iter <- iter + 1
        delta_next <- delta_before + log(sj) - log(share(delta_before, mu)) # update delta
    }
    delta <- delta_next
    flag_cv <- (iter == 500 || max(is.na(delta)) == 1)
    return(list(delta, flag_cv, iter))
}

# gmm objective function (estimate the linear parameters to concentrate it out in the objective function)
blp_moment_condition <- function(theta, data, var_iv_new, var_rand_coef) {
    # theta is the parameter vector, dim(theta) = n(var_rand_coef) + rho + n(linear parameters)
    # data is the panel data, ~ name + year
    # var_iv_new is the instruments used in the 1) IV regression (for estimation of linear parameters) 2) gmm objective function (for estimation of sigma and rho)
    # var_rand_coef is the variable with random coefficient
    delta_new <- c()
    sigma <- theta[1]
    rho <- theta[2]
    # flag_cv_new <- c()
    # iter_new <- c()

    # step 1: construct the delta_j for each product
    for (i in unique(data$year)) {
        J <- data[year == i, .N]
        sj <- data[year == i, mktshr]
        mu <- sigma * outer(as.vector(as.matrix(data[year == i, ..var_rand_coef])), draw) # construct the mu matrix of J x n
        delta_init <- data[year == i, log(mktshr) - log(shr_0)]
        contraction_result <- blp_contraction(sj, mu, delta_init)
        delta_market <- contraction_result[[1]]
        delta_new <- c(delta_new, delta_market)
        # flag_cv_new <- c(flag_cv_new, contraction_result[[2]])
        # iter_new <- c(iter_new, contraction_result[[3]])
    }

    # step 2: construct the lagged difference
    # here we need to introduce dynamics
    # for each y construct y_t-\rho*y_{t-1}
    data$y <- delta_new
    data <- panel(data, ~ name + year)
    data[, y_ld := y - l(y, 1) * rho]
    for (i in c(var_exo, var_end)) {
        data[, (paste0(i, "_ld")) := get(i) - l(get(i), 1) * rho]
    }

    # step 3: esitmation of linear parameters
    iv_formula <- as.formula(paste("y_ld ~", paste(paste0(var_exo, split = "_ld"), collapse = " + "), "|", paste(paste0(var_end, split = "_ld"), collapse = " + "), "~", paste(var_iv_new, collapse = " + ")))
    iv_reg <- feols(iv_formula, data = data)
    residual <- iv_reg$residuals
    data <- unpanel(data) # unpanel the data
    data_used <- data[obs(iv_reg)]

    # step 4: construct the gmm objective function
    var_xz <- c(paste0(var_exo, split = "_ld"), var_iv_new)
    Z <- as.matrix(data_used[, ..var_xz])
    W <- solve(t(Z) %*% Z)
    g <- t(Z) %*% residual
    gmm_obj <- t(g) %*% W %*% g
    return(gmm_obj)
}

# gmm objective function (estimate all parameters simultaneously without concentrating out the linear parameters first)
blp_moment_condition_full <- function(theta, data, var_iv_new, var_rand_coef) {
    delta_new <- c()
    sigma <- theta[1]
    rho <- theta[2]
    # flag_cv_new <- c()
    # iter_new <- c()

    # step 1: construct the delta_j for each product
    for (i in unique(data$year)) {
        J <- data[year == i, .N]
        sj <- data[year == i, mktshr]
        mu <- sigma * outer(as.vector(as.matrix(data[year == i, ..var_rand_coef])), draw)
        delta_init <- data[year == i, log(mktshr) - log(shr_0)]
        contraction_result <- blp_contraction(sj, mu, delta_init)
        delta_market <- contraction_result[[1]]
        delta_new <- c(delta_new, delta_market)
        # flag_cv_new <- c(flag_cv_new, contraction_result[[2]])
        # iter_new <- c(iter_new, contraction_result[[3]])
    }

    # step 2: construct the lagged difference
    # here we need to introduce dynamics
    # for each y construct y_t-\rho*y_{t-1}

    data$y <- delta_new
    data <- panel(data, ~ name + year)
    data[, y_ld := y - l(y, 1) * rho]
    for (i in c("const", var_end, var_exo)) {
        data[, (paste0(i, "_ld")) := get(i) - l(get(i), 1) * rho]
    }
    data <- unpanel(data)

    # step 3: construct the gmm objective function
    y_lhs <- as.vector(data[, y_ld])
    obs_used <- !is.na(y_lhs)
    y_lhs <- y_lhs[obs_used]

    rhs <- paste0(c("const", var_end, var_exo), split = "_ld")
    x_rhs <- as.matrix(data[, ..rhs])
    x_rhs <- x_rhs[obs_used, ]

    residual <- y_lhs - x_rhs %*% theta[3:length(theta)] # residual
    var_xz <- c(paste0(var_exo, split = "_ld"), var_iv_new)
    Z <- as.matrix(data[obs_used, ..var_xz]) # instruments

    W <- solve(t(Z) %*% Z)
    g <- t(Z) %*% residual
    gmm_obj <- t(g) %*% W %*% g
    return(gmm_obj)
}

# for the blp_moment_condition function, the linear parameters are estimated before constructing the gmm function
# blp_intermediate is to return the linear parameters estimated within the blp_moment_condition function
blp_intermediate <- function(theta, data, var_iv_new, var_rand_coef) {
    delta_new <- c()
    sigma <- theta[1]
    rho <- theta[2]
    mu_new <- c()
    # flag_cv_new <- c()
    # iter_new <- c()
    for (i in unique(data$year)) {
        J <- data[year == i, .N]
        sj <- data[year == i, mktshr]
        mu <- sigma * outer(as.vector(as.matrix(data[year == i, ..var_rand_coef])), draw)
        delta_init <- data[year == i, log(mktshr) - log(shr_0)]
        contraction_result <- blp_contraction(sj, mu, delta_init)
        delta_market <- contraction_result[[1]]
        delta_new <- c(delta_new, delta_market)
        mu_new <- rbind(mu_new, mu)
        # flag_cv_new <- c(flag_cv_new, contraction_result[[2]])
        # iter_new <- c(iter_new, contraction_result[[3]])
    }
    data$y <- delta_new
    data <- panel(data, ~ name + year)
    data[, y_ld := y - l(y, 1) * rho]
    for (i in c(var_exo, var_end)) {
        data[, (paste0(i, "_ld")) := get(i) - l(get(i), 1) * rho]
    }
    iv_formula <- as.formula(paste("y_ld ~", paste(paste0(var_exo, split = "_ld"), collapse = " + "), "|", paste(paste0(var_end, split = "_ld"), collapse = " + "), "~", paste(var_iv_new, collapse = " + ")))
    reg_iv <- feols(iv_formula, data = data, cluster = "modelid")
    return(list(delta_new, reg_iv, mu_new))
}

# return the delta, instruments, weights, and residuals from the gmm joint estimation
blp_intermediate_full <- function(theta, data, var_iv_new, var_rand_coef) {
    delta_new <- c()
    sigma <- theta[1]
    rho <- theta[2]
    # flag_cv_new <- c()
    # iter_new <- c()

    # step 1: construct the delta_j for each product
    for (i in unique(data$year)) {
        J <- data[year == i, .N]
        sj <- data[year == i, mktshr]
        mu <- sigma * outer(as.vector(as.matrix(data[year == i, ..var_rand_coef])), draw)
        delta_init <- data[year == i, log(mktshr) - log(shr_0)]
        contraction_result <- blp_contraction(sj, mu, delta_init)
        delta_market <- contraction_result[[1]]
        delta_new <- c(delta_new, delta_market)
        # flag_cv_new <- c(flag_cv_new, contraction_result[[2]])
        # iter_new <- c(iter_new, contraction_result[[3]])
    }

    # step 2: construct the lagged difference
    # here we need to introduce dynamics
    # for each y construct y_t-\rho*y_{t-1}

    data$y <- delta_new
    data <- panel(data, ~ name + year)
    data[, y_ld := y - l(y, 1) * rho]
    for (i in c("const", var_end, var_exo)) {
        data[, (paste0(i, "_ld")) := get(i) - l(get(i), 1) * rho]
    }
    data <- unpanel(data)

    # step 3: construct the gmm objective function
    y_lhs <- as.vector(data[, y_ld])
    obs_used <- !is.na(y_lhs)
    y_lhs <- y_lhs[obs_used]

    rhs <- paste0(c("const", var_end, var_exo), split = "_ld")
    x_rhs <- as.matrix(data[, ..rhs])
    x_rhs <- x_rhs[obs_used, ]

    residual <- y_lhs - x_rhs %*% theta[3:length(theta)] # residual
    var_xz <- c(paste0(var_exo, split = "_ld"), var_iv_new)
    Z <- as.matrix(data[obs_used, ..var_xz]) # instruments

    W <- solve(t(Z) %*% Z)
    g <- t(Z) %*% residual
    gmm_obj <- t(g) %*% W %*% g
    return(list(delta = delta_new, g = g, Z = Z, W = W, residual = residual))
}

# calucate the standard errors following the heteroskedasticity robust formula (without efficient weighting matrix)
# G_tWG_inv %*% (t(G) %*% W %*% Omega %*% W %*% G) %*% G_tWG_inv
gmm_se_hetero <- function(G, Z, W, residual) {
    Omega <- t(Z) %*% diag(as.vector(residual^2)) %*% Z / nrow(Z)
    G_tWG_inv <- solve(t(G) %*% W %*% G)
    robust_var_cov_matrix <- G_tWG_inv %*% (t(G) %*% W %*% Omega %*% W %*% G) %*% G_tWG_inv
    se <- sqrt(diag(robust_var_cov_matrix))
    return(list(se, robust_var_cov_matrix))
}
