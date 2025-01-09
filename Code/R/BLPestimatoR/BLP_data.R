BLP_data<-function (model, market_identifier, product_identifier, par_delta,
    group_structure = NULL, additional_variables = NULL, productData,
    demographic_draws, integration_accuracy, integration_method,
    integration_draws, integration_weights, integration_seed,
    blp_inner_tol = 1e-09, blp_inner_maxit = 10000)
{
    call_arguments <- match.call(expand.dots = TRUE)
    formula <- Formula::as.Formula(model)
    model_rhs_length <- length(formula)[2]
    stopifnot(length(formula)[1] == 1L, model_rhs_length %in%
        1:4)
    tmp <- model.frame(formula, productData, na.action = na.fail)
    tmp <- NULL
    f1 <- formula(formula, lhs = 1, rhs = 0)
    shares <- model.response(model.frame(f1, productData), type = "numeric")
    if (any(shares < 0) || any(shares > 1))
        stop("Shares contain values out of [0,1].")
    nobs <- length(shares)
    if ((!is.character(market_identifier)) || (length(market_identifier) !=
        1))
        stop("market_identifier is not valid.")
    if (!market_identifier %in% names(productData))
        stop("market_identifier is not available in the provided data.")
    market_id_char_in <- vapply(market_identifier, function(x) as.character(get(x, 
        productData)), character(nobs))
    nmkt <- length(unique(market_id_char_in))
    market_id_numeric <- .indexing.markets(market_id_char_in)
    market_id_numeric_o <- order(market_id_numeric)
    if (!product_identifier %in% names(productData))
        stop("market_identifier is not available in the provided data.")
    if ((!is.character(product_identifier)) || (length(product_identifier) !=
        1))
        stop("product_identifier is not valid.")
    product_id_char_in <- vapply(product_identifier, function(x) as.character(get(x,
        productData)), character(nobs))
    tmp <- table(paste0(market_id_char_in, "_", product_id_char_in))
    if (any(tmp > 1))
        warning("Combination of market_identifier and product_identifier is not unique.")
    if (!missing(blp_inner_tol)) {
        if ((!is.finite(blp_inner_tol)) || (length(blp_inner_tol) !=
            1)) {
            cat("Invalid blp_inner_tol. Set to default (1e-9).\n")
            blp_inner_tol <- 1e-09
        }
    }
    if (!missing(blp_inner_maxit)) {
        if ((!is.finite(blp_inner_maxit)) || (length(blp_inner_maxit) != 
            1)) {
            cat("Invalid blp_inner_maxit. Set to default (10000).\n")
            blp_inner_maxit <- 10000
        }
    }
    missing_delta <- missing(par_delta)
    if (!missing_delta) {
        if ((!is.character(par_delta)) || (length(par_delta) !=
            1))
            stop("par_delta is not valid.")
        par_delta_var_name <- par_delta
        par_delta <- get(par_delta, productData)
        delta_error <- !all(is.finite(exp(par_delta)))
    }
    else delta_error <- FALSE
    if (missing_delta || delta_error) {
        par_delta_var_name <- "delta"
        par_delta <- rep(0, nobs)
        cat("Mean utility (variable name: `delta`) is initialized with 0 because of missing or invalid par_delta argument.\n")
    }
    f2 <- formula(formula, lhs = 0, rhs = 1)
    X_lin <- model.matrix(f2, productData)
    tmp <- apply(X_lin, 2, function(x) round(sum(abs(diff(x))), 
        3) == 0)
    if (sum(tmp) > 1)
        stop("Do not include a column of constants. Constants are used by default and can be omitted in the formula.")
    if (model_rhs_length >= 2) {
        f3 <- formula(formula, lhs = 0, rhs = 2)
        X_exg <- model.matrix(f3, productData)
        tmp <- apply(X_exg, 2, function(x) round(sum(abs(diff(x))),
            3) == 0)
        if (sum(tmp) > 1)
            stop("Do not include a column of constants. Constants are used by default and can be omitted in the formula.")
    }
    else X_exg <- NULL
    if (model_rhs_length >= 3) {
        f4 <- formula(formula, lhs = 0, rhs = 3)
        X_rand <- model.matrix(f4, productData)
        K <- dim(X_rand)[2]
        tmp <- apply(X_rand, 2, function(x) round(sum(abs(diff(x))),
            3) == 0)
        if (sum(tmp) > 1) 
            stop("Do not include a column of constants. Constants are used by default and can be omitted in the formula.")
    }
    else X_rand <- NULL
    if (model_rhs_length >= 4) {
        f5 <- formula(formula, lhs = 0, rhs = 4)
        IV <- model.matrix(f5, productData)
        tmp <- apply(IV, 2, function(x) round(sum(abs(diff(x))),
            3) == 0)
        if (sum(tmp) > 1)
            stop("Do not include a column of constants. Constants are used by default and can be omitted in the formula.")
    }
    else IV <- NULL
    has_own_int <- !missing(integration_weights) && !missing(integration_draws)
    has_int_method <- !missing(integration_accuracy) && !missing(integration_method)
    if (has_own_int == has_int_method)
        stop("Provide either the name and accuracy of a valid integration method or your own weights and draws.")
    final_order_draws <- colnames(X_rand)
    if (has_own_int) {
        integration_method <- "provided_by_user"
        weights <- na.fail(as.matrix(as.numeric(integration_weights),
            ncol = 1))
        integration_list_names <- names(integration_draws)
        if (!(length(integration_draws) == K))
            stop("Provided list of integration draws has not enough entries. Number of random coefficients determines length.")
        if (!setequal(integration_list_names, final_order_draws))
            stop("Names of list entries for draws (unobs. heterogeneity) do not match with names of random coefficients. Remember to name any constant \"(Intercept)\" .\n")
        final_order <- match(final_order_draws, integration_list_names)
        integration_draws <- integration_draws[final_order]
        draws_mktShape <- .draws_listToMatrix(drawList = integration_draws,
            amountDraws = length(weights), market_identifier_pd = get(market_identifier,
                productData), market_identifier_list_name = market_identifier,
            use = "rc")
    }
    if (has_int_method) {
        tmp <- get_integration_input(dim = K, method = integration_method,
            accuracy = integration_accuracy, nmkt = nmkt, seed = integration_seed)
        draws_mktShape <- tmp$nodesMktShape
        weights <- tmp$weights
    }
    stopifnot(all(dim(draws_mktShape) == c(nmkt, length(weights) *
        K)))
    if (!missing(demographic_draws)) {
        stopifnot(is.list(demographic_draws))
        demographic_names <- names(demographic_draws)
        M <- length(demographic_names)
        dD <- .draws_listToMatrix(drawList = demographic_draws,
            amountDraws = length(weights), market_identifier_pd = get(market_identifier,
                productData), market_identifier_list_name = market_identifier,
            use = "demographics")
        stopifnot(all(dim(dD) == c(nmkt, length(weights) * M)))
    }
    else {
        demographic_names <- NULL
        dD <- NULL
        M <- 0
    }
    shares <- shares[market_id_numeric_o]
    X_rand <- X_rand[market_id_numeric_o, , drop = FALSE]
    X_lin <- X_lin[market_id_numeric_o, , drop = FALSE]
    X_exg <- X_exg[market_id_numeric_o, , drop = FALSE]
    IV <- IV[market_id_numeric_o, , drop = FALSE]
    if (is.null(dD))
        dD <- matrix(NA)
    if (is.null(group_structure)) {
        group_structure <- NULL
    }
    else {
        if ((!is.character(group_structure)) || (length(group_structure) !=
            1))
            stop("group_structure is not valid.")
        group_structure <- as.character(get(group_structure,
            productData))
        group_structure <- group_structure[market_id_numeric_o]
    }
    market_id_numeric <- market_id_numeric[market_id_numeric_o]
    market_id_char_in <- market_id_char_in[market_id_numeric_o]
    product_id_char_in <- product_id_char_in[market_id_numeric_o]
    par_delta <- par_delta[market_id_numeric_o]
    cdindex <- as.numeric(c(0, cumsum(table(market_id_numeric))))
    Z <- cbind(X_exg, IV)
    Z <- Z[, unique(colnames(Z))]
    if (!is.null(additional_variables)) {
        additional_data <- data.frame(identifier = paste0(market_id_char_in,
            product_id_char_in))
        for (i in 1:length(additional_variables)) {
            vn_i <- additional_variables[i]
            if (!vn_i %in% names(productData))
                stop(paste0(vn_i, " is not available in the provided data."))
            additional_data[[vn_i]] <- productData[[vn_i]][market_id_numeric_o]
        }
    }
    else {
        additional_data <- NULL
    }
    integration <- list(drawsRcMktShape = draws_mktShape, drawsDemMktShape = dD,
        weights = weights, method = integration_method, amountDraws = length(weights))
    parameters <- list(inner_tol = blp_inner_tol, inner_maxit = blp_inner_maxit,
        nobs = nobs, cdindex = cdindex, market_id = market_id_numeric,
        nmkt = nmkt, K = K, total_demogr = M, market_id_numeric_o = market_id_numeric_o,
        demographic_names = demographic_names, market_id_char_in = market_id_char_in, 
        market_id_varname = market_identifier, product_id = product_id_char_in,
        product_id_varname = product_identifier, par_delta_varname = par_delta_var_name,
        share_varname = as.character(f1)[2])
    data <- list(X_lin = X_lin, X_exg = X_exg, X_rand = X_rand,
        shares = shares, Z = Z, group_structure = group_structure,
        delta = par_delta, additional_data = additional_data)
    out <- list(call_arguments = call_arguments, integration = integration,
        parameters = parameters, data = data)
    class(out) <- "blp_data"
    return(out)
}