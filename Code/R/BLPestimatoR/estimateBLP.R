estimateBLP<-function (blp_data, par_theta2, solver_method = "BFGS", solver_maxit = 10000,
    solver_reltol = 1e-06, standardError = "heteroskedastic",
    extremumCheck = FALSE, printLevel = 2, ...)
{
    call_arguments <- match.call(expand.dots = TRUE)
    nobs <- blp_data$parameters$nobs
    K <- blp_data$parameters$K
    if (!is(blp_data, "blp_data"))
        stop("Input has wrong class. Call BLP_data() first.")
    Z <- blp_data$data$Z
    W <- try(solve((t(Z) %*% Z)))
    if (any(class(W) == "try-error"))
        stop("Problems with singular matrizes. This might be caused by (nearly) linear dependent regressors or weak instruments.")
    xzwz <- t(blp_data$data$X_lin) %*% Z %*% W %*% t(Z)
    xzwzx <- xzwz %*% blp_data$data$X_lin
    invxzwzx <- try(solve(xzwzx))
    if (any(class(invxzwzx) == "try-error"))
        stop("Problems with singular matrices. This might be caused by (nearly) linear dependent regressors or weak instruments.")
    blp_data$data$W <- W
    blp_data$data$xzwz <- xzwz
    blp_data$data$invxzwzx <- invxzwzx
    start_theta2 <- .prepare_theta2(par_theta2, final_col_names_par = c("unobs_sd",
        blp_data$parameters$demographic_names), final_row_names_par = colnames(blp_data$data$X_rand), 
        K = blp_data$parameters$K, M = blp_data$parameters$total_demogr)
    cat("blp_data were prepared with the following arguments:\n")
    print(blp_data$call_arguments)
    if (printLevel > 0) {
        cat("Starting a BLP demand estimation with ", blp_data$parameters$nobs, 
            " observations in ", blp_data$parameters$nmkt, " markets...\n")
        cat("[integration::method", blp_data$integration$integration_method,
            " integration::amountDraws", blp_data$integration$amountDraws, 
            "]\n")
        cat("[blp::inner_tol", blp_data$parameters$inner_tol,
            " blp::inner_maxit", blp_data$parameters$inner_maxit,
            "]\n")
        cat("[solver::method", solver_method, " solver::maxit", 
            solver_maxit, " solver::reltol", solver_reltol, "]\n")
    }
    blp_results <- new.env(parent = emptyenv())
    blp_results$deltaOld <- blp_data$data$delta
    blp_results$innerItAll <- c()
    blp_results$negShares <- FALSE
    blp_results$gradient <- rep(NA_real_, start_theta2$total_par)
    start_time <- Sys.time()
    res <- optim(par = start_theta2$par_theta2, fn = gmm_obj,
        gr = gmm_gr, method = solver_method, control = list(reltol = solver_reltol,
            maxit = solver_maxit), indices = start_theta2$indices, 
        blp_results = blp_results, blp_data = blp_data, printLevel = printLevel,
        ...)
    solverMessage <- if (res$convergence == 0)
        "Successful convergence"
    else paste("See error code (optim package)", res$convergence)
    outer_it_out <- res$counts[1]
    end_time <- Sys.time()
    time <- end_time - start_time
    cat("------------------------------------------ \n")
    cat(paste("Solver message:", solverMessage, "\n"))
    if (!(solverMessage == "Successful convergence")) 
        stop("Cannot compute post estimation results due to failed minimization routine.")
    cat("------------------------------------------ \n")
    cat("Final GMM evaluation at optimal parameters: \n")
    innerItAll_out <- blp_results$innerItAll
    blp_results$deltaOld <- rep(0, nobs)
    finalTmp <- gmm_obj(par_theta2 = res$par, indices = start_theta2$indices,
        blp_results = blp_results, blp_data = blp_data, printLevel = 3)
    delta_out <- blp_results$deltaOld
    theta_rc_out <- res$par
    theta_lin_out <- blp_results$bet
    sij_out <- blp_results$sij
    local_min_out <- finalTmp
    gradient_out <- blp_results$gradient
    jacob_out <- blp_results$jacobian
    xi_out <- blp_results$xi
    names_rc <- kronecker(start_theta2$final_col_names_par, start_theta2$final_row_names_par,
        paste, sep = "*")
    relevantRcDem_index <- start_theta2$indices[, "row"] + max(start_theta2$indices[,
        "row"]) * (start_theta2$indices[, "col"] - 1)
    names(theta_rc_out) <- names_rc[relevantRcDem_index]
    X_lin <- blp_data$data$X_lin
    Z <- blp_data$data$Z
    W <- blp_data$data$W
    a <- t(cbind(X_lin, jacob_out)) %*% Z
    tmpSE <- try(solve(a %*% W %*% t(a)))
    lin_len <- dim(X_lin)[2]
    valid_SE <- (standardError %in% c("heteroskedastic", "homoskedastic",
        "cluster")) && (length(standardError) == 1)
    if (!valid_SE) {
        message("Invalid standard error option is provided. Switching to heteroskedastic standard errors...")
        standardError <- "heteroskedastic"
    }
    if (any(class(tmpSE) == "try-error"))
        stop("Standard errors cannot be computed due to singular matrices.")
    if (standardError == "heteroskedastic") {
        cat("Using the heteroskedastic asymptotic variance-covariance matrix... \n")
        omega <- xi_out^2
        b <- as.matrix(t(Z) %*% Diagonal(length(omega), omega) %*%
            Z)
        COV <- tmpSE %*% a %*% W %*% b %*% W %*% t(a) %*% tmpSE
    }
    if (standardError == "homoskedastic") {
        cat("Using the homoskedastic asymptotic variance-covariance matrix... \n")
        COV <- c((t(xi_out) %*% xi_out)/nobs) * tmpSE
    }
    if (standardError == "cluster") {
        group_structure <- blp_data$data$group_structure
        if (any(is.na(group_structure)) || is.null(group_structure)) 
            stop("Valid group_structure is not availalbe in blp_data. Clustered standard errors require a variable that describes the cluster affiliation.")
        group_structure <- data.frame(group_structure = group_structure)
        group_matrix <- model.matrix(as.Formula("~0+group_structure"),
            group_structure)
        tmp <- c(xi_out) * group_matrix
        omega <- tmp %*% t(tmp)
        b <- t(Z) %*% omega %*% Z
        COV <- tmpSE %*% a %*% W %*% b %*% W %*% t(a) %*% tmpSE
    }
    seLinear_out <- sqrt(diag(COV))[1:lin_len]
    seRc_out <- sqrt(diag(COV))[-(1:lin_len)]
    WaldStatistic <- t(matrix(theta_rc_out)) %*% solve(COV[-(1:lin_len),
        -(1:lin_len)]) %*% matrix(theta_rc_out)
    if (extremumCheck) {
        hessian <- invisible(hessian(func = gmm_obj, x = res$par,
            indices = start_theta2$indices, blp_results = blp_results, 
            blp_data = blp_data, printLevel = 0))
        hessianEig <- eigen(hessian)$values
        isMin_out <- sum(hessianEig > 0) == start_theta2$total_par
        isMin_out <- if (isMin_out)
            "positive"
        else "negative"
        cat(paste("Extremum Check:", isMin_out))
    }
    else {
        isMin_out <- NA
    }
    output <- list(theta_rc = theta_rc_out, theta_lin = theta_lin_out, 
        se_rc = seRc_out, se_linear = seLinear_out, local_min = local_min_out,
        gradient = gradient_out, time = time, outer_it = outer_it_out,
        inner_it = innerItAll_out, delta = delta_out, xi = xi_out,
        `#demogrCoef` = blp_data$parameters$total_demogr, `#nmkt` = blp_data$parameters$nmkt,
        `#nobs` = nobs, `#exoCoef` = length(colnames(blp_data$data$X_exg)),
        indices = start_theta2$indices, rand_coef_rownames = start_theta2$final_row_names_par,
        rand_coef_colnames = start_theta2$final_col_names_par, 
        drawsRcMktShape = blp_data$integration$draws_mktShape,
        drawsDemMktShape = blp_data$integration$dD, weights = blp_data$integration$weights,
        sij = sij_out, WaldStatistic = WaldStatistic, IslocalMin = isMin_out,
        outerCrit = solver_reltol, innerCrit = blp_data$parameters$inner_tol, 
        intMethod = blp_data$integration$method, intdraws = length(blp_data$integration$weights),
        standardErrorMethod = standardError, call = call_arguments)
    class(output) <- "blp_est"
    return(output)
}