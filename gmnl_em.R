library("gmnl")

## suml function from mlogit (Croissant, 2013) ====
suml <- function(x){
  n <- length(x)
  if (!is.null(dim(x[[1]]))) {
    d <- dim(x[[1]])
    s <- matrix(0,d[1], d[2])
    for (i in 1:n) {
      x[[i]][is.na(x[[i]])] <- 0
      s <- s + x[[i]]
    }
  }
  else{
    s <- rep(0, length(x[[n]]))
    for (i in 1:n) {
      x[[i]][is.na(x[[i]])] <- 0
      s <- s + x[[i]]
    }
  }
  s
}

## suml function for array ====
suml.array <- function(x){
  n <- length(x)
  d <- dim(x[[1]])
  s <- array(0, dim = c(d[1], d[2], d[3]))
  for (i in 1:n) {
    x[[i]][is.na(x[[i]])] <- 0
    s <- s + x[[i]]
  }
  s
}

## Multinomial Model based on mlogit (Croissant, 2013) ====
ll.mlogit <- function(theta, y, X, gradient = TRUE, 
                      hessian =  TRUE, weights = NULL, 
                      get.bi = FALSE){
  if (is.null(weights)) weights <- 1
  exb  <- lapply(X, function(x) exp(crossprod(t(x), theta)))
  sexb <- suml(exb)
  Pni  <- lapply(exb, function(x){v <- x / sexb; 
  v[is.na(v)] <- 0;
  as.vector(v)})
  Pn <- Reduce("+", mapply("*", Pni, y, SIMPLIFY = FALSE))
  ll <- sum(log(Pn) * weights) 
  
  # Gradient
  if (gradient | hessian ) {
    Px <- suml(mapply("*", X, Pni, SIMPLIFY = FALSE)) 
    yx <- suml(mapply("*", X, y, SIMPLIFY = FALSE))
    gradi <- (yx - Px) * weights ; colnames(gradi) <- names(theta)
    gradi[is.na(gradi)] <- 0
    attr(ll, "gradient") <- gradi
  }
  # Hessian
  if (hessian) {
    dxpx <- lapply(X, function(x){
      d <- x - Px ;
      d[is.na(d)] <- 0 ;
      d})
    hess <- -suml(mapply(function(x, y) crossprod(x * y, y), Pni, dxpx , SIMPLIFY = FALSE))
    attr(ll, "hessian")
  }
  if (get.bi) {
    Pni <- Reduce("cbind", Pni)
    colnames(Pni) <- names(y)
    attr(ll, "prob.alt") <- Pni
    attr(ll, "prob.ind") <- Pn
  }
  ll
}

gmnl_em <- function (formula, data, subset, weights, na.action, model = c("mnl", 
                                                                       "mixl", "smnl", "gmnl", "lc", "mm", "em"), start = NULL, ranp = NULL, 
                  R = 40, Q = 2, haltons = NA, mvar = NULL, seed = 12345, correlation = FALSE, 
                  bound.err = 2, panel = FALSE, hgamma = c("direct", "indirect"), 
                  reflevel = NULL, init.tau = 0.1, init.gamma = 0.1, notscale = NULL, 
                  print.init = FALSE, gradient = TRUE, typeR = TRUE, ...) 
{
  print("Using EM algorithm...")
  start.time <- proc.time()
  callT <- match.call(expand.dots = TRUE)
  callF <- match.call(expand.dots = FALSE)
  formula <- callF$formula <- gFormula(formula)
  nframe <- length(sys.calls())
  model <- match.arg(model)
  has.rand <- !is.null(ranp)
  if (model == "mixl" && !has.rand) 
    stop("mixl model needs ranp to be specified")
  if (model == "gmnl" && !has.rand) 
    stop("gmnl model needs ranp to be specified")
  if (model == "mm" && !has.rand) 
    stop("mn model needs ranp to be specified")
  if ((model == "lc" || model == "mn" || model == "em") && Q < 2) 
    stop("Classes cannot be lower than 2")
  if (model == "mnl") {
    if (is.null(callT$method)) 
      callT$method <- "nr"
  }
  else {
    if (is.null(callT$method)) 
      callT$method <- "bfgs"
  }
  if (!inherits(data, "mlogit.data")) 
    stop("Data must be of class mlogit.data")
  mf <- callT
  m <- match(c("formula", "data", "subset", "na.action", "weights"), 
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  if (!is.null(reflevel)) 
    attr(mf, "index")[["alt"]] <- relevel(attr(mf, "index")[["alt"]], 
                                          reflevel)
  index <- attr(mf, "index")
  alt <- index[["alt"]]
  chid <- index[["chid"]]
  alt.lev <- levels(alt)
  J <- length(alt.lev)
  if (panel) {
    if (model == "mnl") 
      stop("Panel is not relevant for mnl model")
    id <- index[["id"]]
    if (is.null(id)) 
      stop("No individual index")
    id <- split(index[["id"]], alt)[[1]]
  }
  else {
    id <- attr(mf, "index")[["id"]] <- NULL
  }
  y <- model.response(mf)
  X <- model.matrix(formula, mf)
  has.mvar <- has.othervar(formula, 4)
  if (has.mvar) {
    if (model == "mnl" || model == "smnl") 
      stop(paste("Variables for mean are not relevant for", 
                 paste(model, collapse = ": ")))
    if (is.null(mvar)) 
      stop("mvar is null")
    if (!is.list(mvar)) 
      stop("mvar is not a list")
    rvar <- names(mvar)[!(names(mvar) %in% names(ranp))]
    if (length(rvar) > 0) {
      udstr <- paste("The following variables are not specified in the argument ranp:", 
                     paste(unique(rvar), collapse = ", "))
      stop(udstr)
    }
    Z <- model.matrix(formula, mf, rhs = 4)
    for (i in 1:length(mvar)) {
      rvar <- mvar[[i]][!(mvar[[i]] %in% colnames(Z))]
      if (length(rvar) > 0) {
        udstr <- paste("The following variables are not specified in the formula: ", 
                       paste(unique(rvar), collapse = ", "))
        stop(udstr)
      }
    }
  }
  else Z <- NULL
  has.het <- has.othervar(formula, 5)
  if ((model == "lc" || model == "mm" || model == "em") && !has.het) 
    stop("lc needs variables for class probabilities")
  if (has.het) {
    if (model == "mnl" || model == "mixl") 
      stop(paste("variables for scale are not relevant for", 
                 paste(model, collapse = ": ")))
    if (model == "lc" || model == "mm" || model == "em") 
      H <- model.matrix(formula, mf, rhs = 5, Q = Q)
    else H <- model.matrix(formula, mf, rhs = 5)
  }
  else H <- NULL
  if (any(names(mf) == "(weights)")) {
    weights <- mf[["(weights)"]] <- mf[["(weights)"]]/mean(mf[["(weights)"]])
    weights <- split(weights, alt)[[1]]
  }
  else weights <- NULL
  freq <- table(alt[y])
  Xl <- vector(length = J, mode = "list")
  names(Xl) <- levels(alt)
  for (i in levels(alt)) {
    Xl[[i]] <- X[alt == i, , drop = FALSE]
  }
  yl <- split(y, alt)
  yl <- lapply(yl, function(x) {
    x[is.na(x)] <- FALSE
    x
  })
  mean.names <- colnames(X)
  names.stds <- start.stds <- c()
  if (has.rand) {
    if (!correlation) {
      ndist <- ranp[!(ranp %in% c("cn", "ln", "n", "u", 
                                  "t", "sb"))]
      if (length(ndist) > 0) {
        udstr <- paste("unknown distribution", paste(unique(ndist), 
                                                     collapse = ", "))
        stop(udstr)
      }
    }
    else {
      ndist <- ranp[!(ranp %in% c("cn", "ln", "n"))]
      if (length(ndist) > 0) {
        udstr <- paste("Correlated parameters is suitable for distribution from normal, such as cn, ln or n")
        stop(udstr)
      }
    }
    if (model == "gmnl") {
      ndist <- ranp[!(ranp %in% c("n", "u", "t", "ln"))]
      if (length(ndist) > 0) {
        udstr <- paste("Coefficients from G-MNL model can only be distributed as n, u, ln or t")
        stop(udstr)
      }
    }
    namesX <- colnames(X)
    novar <- names(ranp)[!((names(ranp) %in% namesX))]
    if (length(novar) > 0) {
      uvar <- paste("The following random variables are not in the data: ", 
                    paste(unique(novar), collapse = ", "))
      stop(uvar)
    }
    Vara <- sort(match(names(ranp), namesX))
    Varc <- (1:ncol(X))[-Vara]
    fixed <- !(length(Varc) == 0)
    Xa <- lapply(Xl, function(x) x[, Vara, drop = F])
    Xc <- lapply(Xl, function(x) x[, Varc, drop = F])
    allX <- if (fixed) 
      mapply(cbind, Xc, Xa, SIMPLIFY = FALSE)
    else Xa
    mean.names <- colnames(allX[[1]])
    nrap <- length(Vara)
    if (!correlation) {
      names.stds <- paste("sd", namesX[Vara], sep = ".")
      start.stds <- rep(0.1, nrap)
      names(start.stds) <- names.stds
    }
    else {
      names.stds <- c()
      Ka <- length(ranp)
      for (i in 1:Ka) {
        names.stds <- c(names.stds, paste("sd", names(ranp)[i], 
                                          names(ranp)[i:Ka], sep = "."))
        start.stds <- rep(0.1, 0.5 * nrap * (nrap + 1))
        names(start.stds) <- names.stds
      }
    }
  }
  else {
    mean.names <- colnames(X)
    allX <- Xl
  }
  names.het <- start.het <- c()
  if (has.het) {
    if (model == "lc" | model == "mm" || model == "em") {
      start.het <- rep(0, ncol(H))
      names.het <- colnames(H)
      classes <- attr(H, "alt")
      Hl <- vector(length = Q, mode = "list")
      names(Hl) <- levels(classes)
      for (i in levels(classes)) {
        Hl[[i]] <- H[classes == i, , drop = FALSE]
      }
    }
    else {
      start.het <- rep(0, ncol(H))
      names.het <- paste("het", colnames(H), sep = ".")
      names(start.het) <- names.het
    }
  }
  names.phi <- start.phi <- c()
  if (has.mvar) {
    names.phi <- unlist(lapply(names(mvar), function(x) outer(x, 
                                                              mvar[[x]], FUN = paste, sep = ".")))
    start.phi <- rep(0, length(names.phi))
    names(start.phi) <- names.phi
  }
  if (model == "lc" || model == "em") {
    lc.names <- c()
    for (i in 1:Q) {
      lc.names <- c(lc.names, paste("class", i, colnames(X), 
                                    sep = "."))
    }
  }
  if (model == "mm") {
    lc.names <- c()
    ls.names <- c()
    for (i in 1:Q) {
      lc.names <- c(lc.names, paste("class", i, colnames(X), 
                                    sep = "."))
      ls.names <- c(ls.names, paste("class", i, names.stds, 
                                    sep = "."))
    }
  }
  if (model == "smnl") 
    allnames <- c(mean.names, "tau", names.het)
  if (model == "mixl") 
    allnames <- c(mean.names, names.phi, names.stds)
  if (model == "gmnl") 
    allnames <- c(mean.names, names.het, names.stds, "tau", 
                  "gamma")
  if (model == "lc" || model == "em") 
    allnames <- c(lc.names, names.het)
  if (model == "mm") 
    allnames <- c(lc.names, ls.names, names.het)
  if (!is.null(start)) {
    if (length(start) != length(allnames)) 
      stop("Incorrect number of initial parameters")
    theta <- start
    names(theta) <- allnames
  }
  else {
    K <- ncol(allX[[1]])
    theta <- rep(0, K)
    names(theta) <- mean.names
    if (model != "mnl") {
      calls <- call("maxLik")
      calls$start <- theta
      calls$method <- "nr"
      calls$weights <- as.name("weights")
      calls$X <- as.name("allX")
      calls$y <- as.name("yl")
      calls$logLik <- as.name("ll.mlogit")
      mean <- coef(eval(calls, sys.frame(which = nframe)))
      if (model == "lc" || model == "mm" || model == "em") {
        lc.mean <- c()
        init.shift <- seq(-0.02, 0.02, length.out = Q)
        for (i in 1:Q) {
          lc.mean <- c(lc.mean, mean + init.shift[i])
        }
      }
      if (model == "mm") {
        ls.mean <- c()
        init.shift <- seq(-0.02, 0.02, length.out = Q)
        for (i in 1:Q) {
          ls.mean <- c(ls.mean, start.stds + init.shift[i])
        }
      }
      theta <- switch(model, smnl = c(mean, init.tau, start.het), 
                      mixl = c(mean, start.phi, start.stds), gmnl = c(mean, 
                                                                      start.phi, start.het, start.stds, init.tau, 
                                                                      init.gamma), em = , lc = c(lc.mean, start.het), mm = c(lc.mean, 
                                                                                                                      ls.mean, start.het))
      names(theta) <- allnames
    }
  }
  if (model == "smnl" || model == "gmnl") {
    if (!is.null(notscale)) {
      if (length(notscale) != ncol(X)) 
        stop("Not scaled variables vector is not the same length as initial parameters")
      names(notscale) <- mean.names
      wns <- names(notscale[notscale == 1])
      cat("\nThe following variables are not scaled:\n")
      print(wns)
    }
    else {
      notscale <- rep(0, ncol(X))
      names(notscale) <- mean.names
    }
  }
  if ((model == "mixl" || model == "gmnl" || model == "mm") & 
      is.null(start)) {
    ln <- names(ranp[ranp == "ln"])
    sb <- names(ranp[ranp == "sb"])
    if (length(ln) != 0) {
      if (model == "mixl" || model == "gmnl") {
        if (sum(theta[ln] < 0) >= 1) 
          stop("Some variables specified as ln have negative values in the clogit")
        theta[ln] <- log(theta[ln])
      }
      else {
        log.names <- c()
        for (i in 1:Q) log.names <- c(log.names, paste("class", 
                                                       i, ln, sep = "."))
        if (sum(theta[log.names] < 0) >= 1) 
          stop("Some variables specified as ln have negative values in the clogit")
        theta[log.names] <- log(theta[log.names])
      }
    }
    if (length(sb) != 0) {
      if (model == "mixl" || model == "gmnl") {
        if (sum(theta[sb] < 0) >= 1) 
          stop("Some variables specified as sb have negative values in the clogit")
        theta[sb] <- log(theta[sb])
      }
      else {
        sb.names <- c()
        for (i in 1:Q) sb.names <- c(sb.names, paste("class", 
                                                     i, sb, sep = "."))
        if (sum(theta[sb.names] < 0) >= 1) 
          stop("Some variables specified as sb have negative values in the clogit")
        theta[sb.names] <- log(theta[sb.names])
      }
    }
  }
  if (print.init) {
    cat("\nStarting Values:\n")
    print(theta)
  }
  opt <- callT
  opt$start <- theta
  m <- match(c("method", "print.level", "iterlim", "start", 
               "tol", "ftol", "steptol", "fixed", "constraints"), names(opt), 
             0L)
  opt <- opt[c(1L, m)]
  opt[[1]] <- as.name("maxLik")
  opt[c("X", "y")] <- list(as.name("Xl"), as.name("yl"))
  if (model == "mnl") 
    opt$logLik <- as.name("ll.mlogit")
  opt$gradient <- as.name("gradient")
  if (model == "smnl") {
    cat("Estimating SMNL model", "\n")
    opt$logLik <- as.name("ll.smlogit")
    opt[c("R", "seed", "bound.err", "id", "H", "notscale", 
          "typeR")] <- list(as.name("R"), as.name("seed"), 
                            as.name("bound.err"), as.name("id"), as.name("H"), 
                            as.name("notscale"), as.name("typeR"))
  }
  if (model == "mixl") {
    cat("Estimating MIXL model", "\n")
    opt$logLik <- as.name("ll.mixlog")
    opt[c("R", "seed", "ranp", "id", "Z", "correlation", 
          "haltons", "mvar")] <- list(as.name("R"), as.name("seed"), 
                                      as.name("ranp"), as.name("id"), as.name("Z"), as.name("correlation"), 
                                      as.name("haltons"), as.name("mvar"))
  }
  if (model == "gmnl") {
    cat("Estimating GMNL model", "\n")
    opt$logLik <- as.name("ll.gmlogit")
    hgamma <- match.arg(hgamma)
    opt[c("R", "seed", "ranp", "id", "H", "correlation", 
          "haltons", "bound.err", "hgamma", "notscale", "typeR")] <- list(as.name("R"), 
                                                                          as.name("seed"), as.name("ranp"), as.name("id"), 
                                                                          as.name("H"), as.name("correlation"), as.name("haltons"), 
                                                                          as.name("bound.err"), as.name("hgamma"), as.name("notscale"), 
                                                                          as.name("typeR"))
  }
  if (is.null(weights)) 
    weights <- 1
  opt$weights <- as.name("weights")
  if (model == "lc") {
    cat("Estimating LC model", "\n")
    opt$logLik <- as.name("ll.mlogitlc")
    opt[c("H", "Q", "id")] <- list(as.name("Hl"), as.name("Q"), 
                                   as.name("id"))
  }
  if (model == "em") {
    cat("Estimating EM model", "\n")
    opt[[1]] = as.name("em.lc")
    opt$logLik <- as.name("ll.em.lc")
    opt[c("H", "Q", "id")] = list(as.name("Hl"), as.name("Q"), 
                                   as.name("id"))
  }
  if (model == "mm") {
    cat("Estimating MM-MNL model", "\n")
    opt$logLik <- as.name("ll.mnlogit")
    opt[c("R", "seed", "ranp", "id", "H", "correlation", 
          "haltons", "Q")] <- list(as.name("R"), as.name("seed"), 
                                   as.name("ranp"), as.name("id"), as.name("Hl"), as.name("correlation"), 
                                   as.name("haltons"), as.name("Q"))
  }
  
  x <- eval(opt, sys.frame(which = nframe))
  actpar <- activePar(x)
  
  names(opt)[[2]] <- "theta"
  betahat <- coef(x)
  opt$gradient <- FALSE
  opt$get.bi <- TRUE
  opt$fixed <- opt$steptol <- opt$iterlim <- opt$method <- opt$print.level <- opt$tol <- opt$ftol <- opt$logLik <- opt$start <- NULL
  opt$constraints <- NULL
  if (model != "mnl") {
    opt[[1]] <- switch(model, smnl = as.name("ll.smlogit"), 
                       mixl = as.name("ll.mixlog"), gmnl = as.name("ll.gmlogit"), 
                       em = , lc = as.name("ll.mlogitlc"), mm = as.name("ll.mnlogit")) # em will use the lc log-likelihood function
    if (has.rand) {
      if (!correlation) {
        if (model != "mm") {
          diag.sd <- paste("sd", namesX[Vara], sep = ".")
          betahat <- ifelse(names(betahat) %in% diag.sd, 
                            abs(betahat), betahat)
          names(betahat) <- names(coef(x))
        }
        else {
          betahat <- ifelse(names(betahat) %in% ls.names, 
                            abs(betahat), betahat)
          names(betahat) <- names(coef(x))
        }
      }
    }
    opt[[2]] <- betahat
    again <- eval(opt, sys.frame(which = nframe))
    bi <- attr(again, "bi")
    Qir <- attr(again, "Qir")
    x$estimate <- betahat
  }
  else {
    opt$hessian <- FALSE
    opt[[1]] <- as.name("ll.mlogit")
    names(opt)[[2]] <- "theta"
    opt[[2]] <- betahat
    again <- eval(opt, sys.frame(which = nframe))
    bi <- NULL
    Qir <- NULL
  }
  prob.alt <- attr(again, "prob.alt")
  prob.ind <- attr(again, "prob.ind")
  residuals <- Reduce("cbind", yl) - prob.alt
  colnames(prob.alt) <- colnames(residuals) <- levels(alt)
  if (gradient) 
    gradientObs <- x$gradientObs[, actpar]
  else gradientObs <- NULL
  logLik <- structure(list(maximum = logLik(x), gradient = x$gradient[actpar], 
                           nobs = nrow(X)/J, gradientObs = gradientObs, hessian = hessian(x)[actpar, 
                                                                                             actpar], iterations = nIter(x), type = maximType(x), 
                           code = returnCode(x), nparam = length(x$estimate), message = returnMessage(x)), 
                      class = "logLik")
  result <- structure(list(coefficients = x$estimate[actpar], 
                           logLik = logLik, mf = mf, formula = formula, time = proc.time() - 
                             start.time, freq = freq, draws = haltons, model = model, 
                           R = R, ranp = ranp, residuals = residuals, correlation = correlation, 
                           prob.alt = prob.alt, prob.ind = prob.ind, bi = bi, Qir = Qir, 
                           notscale = notscale, Q = Q, call = callT), class = "gmnl")
  result
}


