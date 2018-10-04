## expectation maximization function for Latent class model ====
em.lc <- function(start, y, X, H, Q, method, logLik, id = NULL, weights = NULL,
                      gradient = TRUE, get.bi =  FALSE){
  
  # Parameters and globals
  K     <- ncol(X[[1]])
  N     <- nrow(X[[1]])
  J     <- length(X)
  panel <- !is.null(id)
  if (panel) if (length(weights) == 1) weights <- rep(weights, N)
  
  maxIter = 1000
  ftol_abs = 1.0e-3
  
  theta = start
  
  for (ii in 1:maxIter) {
    
    print(paste0("ii = ", ii))
    print(theta)
    
    beta  <- matrix(theta[1L:(K * Q)], nrow = K, ncol = Q)
    rownames(beta) <- colnames(X[[1]]); colnames(beta) <- paste("class", 1:Q, sep = ':')
    gamma <- theta[-c(1L:(K * Q))]
    
    # Make weigths
    ew <- lapply(H, function(x) exp(crossprod(t(x), gamma)))
    sew <- suml(ew)
    Wnq <- lapply(ew, function(x){v <- x / sew; 
                  v[is.na(v)] <- 0;
                  as.vector(v)})
    Wnq <- Reduce(cbind, Wnq) # N*Q matrix: Probability for individual n in segement q
    
    # Make Multinomial Probability
    ep <- vector(mode = "list", length = J)
    for (j in 1:J) ep[[j]] <- matrix(NA, N, Q)
    for (j in 1:J) ep[[j]] <- exp(tcrossprod(X[[j]], t(beta))) 
    
    sep  <- suml(ep)
    Pnjq <- lapply(ep, function(x) x / sep) # list of N * Q
    Pnq  <- suml(mapply("*", Pnjq, y, SIMPLIFY = FALSE)) # Selected Probability
    if (panel) Pnq <- apply(Pnq, 2, tapply, id, prod)
    WPnq <- Wnq * Pnq # N * Q
    Ln   <- apply(WPnq, 1, sum)
    
    # E-step
    Wnq_post = WPnq/Ln
    
    # M-step
    out = maxLik(method = "bfgs", start =theta, logLik = logLik, X = X, y = y, Wnq_post = Wnq_post, 
                 H = H, Q = Q, id = id, gradient = gradient, weights = weights)
    theta_hat = coef(out)
    
    if (all(abs(theta_hat-theta) <= ftol_abs)) {
      break
    } else {
      theta = theta_hat
    }
  }
  
  out
  
  # ## Gradient
  # if (gradient) {
  #   lambda <- mapply(function(y, p) y - p, y, Pnjq, SIMPLIFY =  FALSE) # J list of N * Q
  #   Qnr  <-  WPnq / Ln # N * Q
  #   if (panel) Qnr <- Qnr[id, ] 
  #   eta  <- lapply(lambda, function(x) x * Qnr) # J list of N * Q
  #   etar <- lapply(eta,  function(x) x[, rep(1:Q, each = K)])
  #   Xg   <- lapply(X,  function(x) x[, rep(1:K, Q)]) # J list of N * (K*Q)
  #   grad.beta <- suml(mapply("*", Xg, etar, SIMPLIFY = FALSE))
  #   
  #   if (panel) {
  #     Wnq <- Wnq[id, ]
  #     H   <- lapply(H, function(x) x[id, ])
  #   }
  #   Wg <- vector(mode = "list", length = Q)
  #   IQ <- diag(Q)
  #   #for(q in 1:Q) Wg[[q]] <- matrix(NA, N, 1)
  #   for (q in 1:Q) Wg[[q]] <- rowSums(Qnr * (repRows(IQ[q, ], N) - repCols(Wnq[, q], Q)))
  #   grad.gamma <- suml(mapply("*", H, Wg, SIMPLIFY = FALSE)) 
  #   gari <- cbind(grad.beta, grad.gamma)
  #   colnames(gari) <- names(theta)
  #   attr(lnL, "gradient") <- gari * weights
  # }
  # if (get.bi) {
  #   attr(lnL, 'Wnq') <- Wnq
  #   if (panel) Wnq <- Wnq[id, ]
  #   Pw <- lapply(Pnjq, function(x) x * Wnq)
  #   attr(lnL, "prob.alt") <- sapply(Pw, function(x) apply(x, 1, sum))
  #   attr(lnL, "prob.ind") <- Ln
  #   attr(lnL, "bi") <- bi
  #   attr(lnL, 'Qir') <- Qnr # WPnq / Ln
  # }
  # lnL 
}

ll.em.lc = function(theta, y, X, H, Q, Wnq_post, id = NULL, weights = NULL,
                                   gradient = FALSE, get.bi =  FALSE){
  # Parameters and globals
  K     <- ncol(X[[1]])
  N     <- nrow(X[[1]])
  J     <- length(X)
  panel <- !is.null(id)
  if (panel) if (length(weights) == 1) weights <- rep(weights, N)
  
  beta  <- matrix(theta[1L:(K * Q)], nrow = K, ncol = Q)
  rownames(beta) <- colnames(X[[1]]); colnames(beta) <- paste("class", 1:Q, sep = ':')
  gamma <- theta[-c(1L:(K * Q))]
  
  # if (get.bi) bi <- t(beta)
  
  # Make weigths
  ew <- lapply(H, function(x) exp(crossprod(t(x), gamma)))
  sew <- suml(ew)
  Wnq <- lapply(ew, function(x){v <- x / sew; 
                v[is.na(v)] <- 0;
                as.vector(v)})
  Wnq <- Reduce(cbind, Wnq) # N*Q matrix: Probability for individual n in segement q
  
  # Make Multinomial Probability
  ep <- vector(mode = "list", length = J)
  for (j in 1:J) ep[[j]] <- matrix(NA, N, Q)
  for (j in 1:J) ep[[j]] <- exp(tcrossprod(X[[j]], t(beta))) 
  
  sep  <- suml(ep)
  Pnjq <- lapply(ep, function(x) x / sep) # list of N * Q
  Pnq  <- suml(mapply("*", Pnjq, y, SIMPLIFY = FALSE)) # Selected Probability
  if (panel) Pnq <- apply(Pnq, 2, tapply, id, prod)
  
  WPnq <- Wnq * Pnq # N * Q
  Ln   <- apply(Wnq_post * log(WPnq), 1, sum)
  # Ln = apply((Wnq_post * (log(Pnq) + log(Wnq))), 1, sum)
  
  lnL <- if (panel) sum(Ln * weights[!duplicated(id)]) else sum(Ln * weights)
  

  
#  if (get.bi)  Qnr  <-  WPnq / Ln
#  lnL <- if (panel) sum(log(Ln) * weights[!duplicated(id)]) else sum(log(Ln) * weights)
  
  
  # ## Gradient
  # if (gradient) {
  #   lambda <- mapply(function(y, p) y - p, y, Pnjq, SIMPLIFY =  FALSE) # J list of N * Q
  #   Qnr  <-  WPnq / Ln # N * Q
  #   if (panel) Qnr <- Qnr[id, ] 
  #   eta  <- lapply(lambda, function(x) x * Qnr) # J list of N * Q
  #   etar <- lapply(eta,  function(x) x[, rep(1:Q, each = K)])
  #   Xg   <- lapply(X,  function(x) x[, rep(1:K, Q)]) # J list of N * (K*Q)
  #   grad.beta <- suml(mapply("*", Xg, etar, SIMPLIFY = FALSE))
  #   
  #   if (panel) {
  #     Wnq <- Wnq[id, ]
  #     H   <- lapply(H, function(x) x[id, ])
  #   }
  #   Wg <- vector(mode = "list", length = Q)
  #   IQ <- diag(Q)
  #   #for(q in 1:Q) Wg[[q]] <- matrix(NA, N, 1)
  #   for (q in 1:Q) Wg[[q]] <- rowSums(Qnr * (repRows(IQ[q, ], N) - repCols(Wnq[, q], Q)))
  #   grad.gamma <- suml(mapply("*", H, Wg, SIMPLIFY = FALSE)) 
  #   gari <- cbind(grad.beta, grad.gamma)
  #   colnames(gari) <- names(theta)
  #   attr(lnL, "gradient") <- gari * weights
  # }
  # if (get.bi) {
  #   attr(lnL, 'Wnq') <- Wnq
  #   if (panel) Wnq <- Wnq[id, ]
  #   Pw <- lapply(Pnjq, function(x) x * Wnq)
  #   attr(lnL, "prob.alt") <- sapply(Pw, function(x) apply(x, 1, sum))
  #   attr(lnL, "prob.ind") <- Ln
  #   attr(lnL, "bi") <- bi
  #   attr(lnL, 'Qir') <- Qnr # WPnq / Ln
  # }
  
  lnL 
}
