
glmqm <- function(formula,
                  data,
                  subset,
                  weights = NULL,
                  weights.x = NULL,
                  weights.var = NULL,
                  na.action,
                  family = gaussian,
                  contrasts = NULL,
                  tau = 0.5,
                  start.fit,
                  offset,
                  control = glmqmControl(),                        
                  rlm.control = rlmControl(),
                  glmrob.control = glmrobControl(),
                  ...) 
{
  
  #######################################################
  ## quantiles verification
  if (any(tau <= 0) | any(tau >= 1)) 
    stop("Quantile index out of range")
  
  
  #######################################################
  
  #######################################################
  ## matchig data and other elements
  call <- match.call()
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  fami <- family$family
  if (is.null(fami)) 
    stop(gettextf("'%s' is not a valid family (see ?family)", 
                  as.character(call[["family"]])), domain = NA)
  if (!(fami %in% c("binomial", "poisson", "gaussian"))) {
    stop(gettextf("GLM M-Quantile fitting not yet implemented for family %s", 
                  fami), domain = NA)
  }
  if (missing(data)) 
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  
  if (length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm)) 
      names(Y) <- nm
  }
  
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0)
  #######################################################
  
  #######################################################
  ## verification of weights
  weights <- model.weights(mf)
  
  if (!is.null(weights) && any(weights < 0)) 
    stop("'weights' must be non-negative")
  
  ########################################################
  ### glmrob specific
  if (fami %in% c('poisson','binomial')) {
    
    ## settings
    maxit <- control$maxit 
    acc <- control$acc
    k <- control$k
    method <- glmrob.control$method
    method.control <- glmrob.control$method.control
    weights.on.x <- glmrob.control$weights.on.x
    trace.lev <-  glmrob.control$trace.lev
    fami <- family$family
    start.glmrob <- glmrob.control$start
    method.control$tcc <- k
    ## end settings
    
    ## checking offset
    offset <- model.offset(mf)  
    
    if (identical(method, "model.frame")) 
      meth. <- if (method == "WBY") "BY"
    else method
    if (is.null(control)) 
      #control <- get(paste0("glmrob", meth., ".control"))(...)
      control <- get(paste0("glmrob", meth., ".control"))
    
    #if (missing(weights.on.x) || is.character(weights.on.x)) 
    #  weights.on.x <- match.arg(weights.on.x)
    
    # if (!is.null(offset) && length(offset) != NROW(Y)) 
    #   stop(gettextf("Number of offsets is %d, should rather equal %d (number of observations)", 
    #                 length(offset), NROW(Y)), domain = NA)
    # 
    # else if (!(is.function(weights.on.x) || is.list(weights.on.x) || 
    #            (is.numeric(weights.on.x) && length(weights.on.x) == 
    #             NROW(Y)))) 
    #   stop("'weights.on.x' must be a string, function, list or numeric n-vector")
    # 
    # if (!is.null(start) && !is.numeric(start)) {
    #   if (!is.character(start)) 
    #     stop("'start' must be a numeric vector, NULL, or a character string")
    # }
  }
  
  #######################################################
  
  result <- switch(fami,
                   gaussian = lmqm.fit(y = Y, X = X, 
                                       weights = weights, 
                                       weights.x = weights,
                                       weights.var = weights, 
                                       q = tau, 
                                       control = control,
                                       rlm.control = rlm.control),
                   poisson = glmqm.fit(y = Y, X = X,
                                       family = family, 
                                       weights = weights,
                                       offset = offset,
                                       q = tau,
                                       control = control,
                                       glmrob.contr = glmrob.control),
                   binomial = glmqm.fit(y = Y, X = X,
                                        family = family, 
                                        weights = weights,
                                        offset = NULL,
                                        q = tau,
                                        control = control,
                                        glmrob.contr = glmrob.control))
  
  ########################################################
  fit <- list()
  fit$call <- call
  fit$family <- family
  fit$na.action <- attr(mf, "na.action")
  fit$contrasts <- attr(X, "contrasts")
  fit$nobs <- length(Y)
  fit$tau <- tau
  fit$fitted.values <- result$fitted.values
  fit$residuals <- result$residuals
  fit$weights <- result$q.weights
  fit$coefficients <- result$coefficients
  fit$control <- control
  fit$rlm.control <- rlm.control
  #fit$glmrob.control <- glmrob.control
  fit$iter <- result$iter
  
  class(fit) <- c('glmqm', 'glm', 'lm')
  fit
}



###

lmqm.fit <- function(y,
                     X, 
                     weights, ## case.weights
                     weights.x, ## var.weights (weight on x)
                     weights.var,
                     q,
                     control,
                     rlm.control,
                     ...) {
  
  init <- rlm.control$init
  scale.est <- rlm.control$scale.est
  method <- rlm.control$method
  
  k <- control$k
  psi <- control$psi
  
  maxit <- control$maxit
  acc <- control$acc
  test.vec <- rlm.control$test.vec
  
  temp.rlm <- MASS:::rlm.default(y = y,
                                 x = X,
                                 w = weights.var,
                                 weights = weights,
                                 init = init,
                                 psi = psi,
                                 scale.est = scale.est,
                                 k2 = k,
                                 method = method,
                                 maxit = maxit,
                                 acc = acc,
                                 test.vec = test.vec)
  
  ## mq-estimates
  done <- FALSE
  conv <- NULL
  resid <- temp.rlm$residuals
  n1 <- nrow(X) - ncol(X)
  if (scale.est != "MM")
    scale <- mad(resid / sqrt(weights.var), 0)
  theta <- 2 * pnorm(k) - 1
  gamma <- theta + k ^ 2 * (1 - theta) - 2 * k * dnorm(k)
  qest <- matrix(0, nrow = ncol(X), ncol = length(q))
  qwt <- matrix(0, nrow = nrow(X), ncol = length(q))
  qfit <- matrix(0, nrow = nrow(X), ncol = length(q))
  qres <- matrix(0, nrow = nrow(X), ncol = length(q))
  qiter <- numeric(length = length(q))
  
  irls.delta <- function(old, new) {
    sqrt(sum((old - new) ^ 2) / max(1e-20, sum(old ^ 2)))}
  
  irls.rrxwr <- function(x, w, r) {
    w <- sqrt(w)
    max(abs((matrix( r * w, 1, length(r) ) %*% x) / 
              sqrt(matrix(w, 1, length(r)) %*% (x ^ 2)))) / sqrt(sum(w * r ^ 2))
  }
  
  for (i in 1:length(q)) {
    for (iiter in 1:maxit) {
      if (!is.null(test.vec))
        testpv <- get(test.vec)
      if (method != "MM") {
        if (scale.est == "MAD")
          scale <- median(abs(resid / sqrt(weights.var))) / 0.6745
        else
          scale <-
            sqrt(sum(pmin(resid^2 / weights.var, (k * scale) ^ 2)) / (n1 * gamma))
        if (scale == 0) {
          done <- TRUE
          break
        }
      }
      w <- psi(resid / (scale * sqrt(weights.var))) * weights
      ww <- 2 * (1 - q[i]) * w
      ww[resid > 0] <- 2 * q[i] * w[resid > 0]
      w <- ww
      temp <- lm.wfit(X, y, w, method = "qr")
      coef <- temp$coef
      resid <- temp$residuals
      if (!is.null(test.vec))
        convi <- irls.delta(testpv, get(test.vec))
      else
        convi <- irls.rrxwr(X, wmod, resid)
      conv <- c(conv, convi)
      done <- (convi <= acc)
      if (done)
        break
    }
    if (!done)
      warning(paste("rlm failed to converge in", maxit, "steps at q = ", q[i]))
    qest[, i] <- coef
    qwt[, i] <- w
    qfit[, i] <- temp$fitted.values
    qres[, i] <- resid
    qiter[i] <- iiter
  }
  
  
  return(list(
    fitted.values = qfit,
    residuals = qres,
    q.values = q,
    q.weights = qwt,
    coefficients = qest,
    iter = qiter
  ))
}


### 
glmqm.fit <- function(y,
                      X,
                      family, 
                      weights, 
                      offset,
                      q,
                      control,
                      glmrob.contr,
                      ...) {
  
  maxit <- control$maxit 
  acc <- control$acc
  k <- control$k
  method <- glmrob.contr$method
  method.control <- glmrob.contr$method.control
  wox <- glmrob.contr$weights.on.x
  trace.lev <-  glmrob.contr$trace.lev
  fami <- family$family
  start.glmrob <- glmrob.contr$start
  method.control$tcc <- k
  
  temp.glmrob <- switch(
    method,
    cubif = stop("For method 'cubif', use glmRob() from package 'robust'"),
    Mqle = robustbase:::glmrobMqle(X = X,
                                   y = y,
                                   weights = weights,
                                   start = start.glmrob,
                                   offset = offset,
                                   family = family,
                                   weights.on.x = wox,
                                   control = method.control,
                                   intercept = attr(X, "intercept") > 0,
                                   trace = trace.lev),
    BY = ,
    WBY = {
      if (fami != "binomial")
        stop(gettextf("method='%s' is only applicable for binomial family, but family=\"\"", 
                      method, family$family), domain = NA)
      robustbase:::glmrobBY(X = X,
                            y = y,
                            weights = weights,
                            start = start.glmrob,
                            method = method,
                            weights.on.x = wox,
                            control = method.control,
                            intercept = attr(X, "intercept") > 0,
                            trace.lev = trace.lev)
    },
    MT = {
      robustbase:::glmrobMT(x = X,
                            y = y,
                            weights = weights,
                            start = start.glmrob,
                            offset = offset,
                            family = family,
                            weights.on.x = wox,
                            control = method.control,
                            intercept = attr(X, "intercept") > 0,
                            trace.lev = trace.lev)
    },
    stop("invalid 'method': ", method)
  )
  
  fit.init <- temp.glmrob$fitted.values
  resid.init <- y - temp.glmrob$fitted.values
  coef <- temp.glmrob$coef
  w.x <- temp.glmrob$w.x 
  phi.init <- 1
  
  
  ########################################################
  ## this will be rewritten to c++
  
  irls.delta <- function(old, new) {
    abs(max(old - new)) / abs(max(old))
  }
  
  psi_q <- function(resid, scale, k, weights, q) {
    tmp <- psi.huber(resid / scale, k = k) * weights * (resid / scale)
    tmp1 <- 2 * (1 - q) * tmp
    tmp1[resid > 0] <- 2 * q * tmp[resid > 0]
    tmp1
  }
  
  done <- FALSE
  conv <- NULL
  qest <- matrix(0, nrow = ncol(X), ncol = length(q))
  qfit <- matrix(0, nrow = nrow(X), ncol = length(q))
  qres <- matrix(0, nrow = nrow(X), ncol = length(q))
  qwt  <- matrix(0, nrow = nrow(X), ncol = length(q))
  qiter <- numeric(length = length(q))
  
  ## switch function 
  
  if (fami == 'binomial')
  {
    for (i in 1:length(q)) {
      
      #We define the starting values
      resid <- resid.init
      fit <- fit.init
      a.j <- weights
      w <- weights
      coef <- coef
      n <- length(resid.init)
      
      for (iiter in 1:maxit) {
        
        resid.old <- resid
        coef.old <- coef
        
        
        # We define the probability mu=exp(xb)/(1+exp(xb))
        probab <- fit
        mu <- probab
        
        #We define the variance
        V <- probab * (1 - probab)
        
        #We define the scale
        scale <- drop(sqrt(V))
        
        #We standardize the residuals
        r.stand <- (y - mu) / scale
        
        #we compute i1 and i2
        jinf <- floor(mu - k * scale)
        jsup <- floor(mu + k * scale)
        ni <-  rep(1, n)
        
        #We compute the values of a_j(b)
        if (k == Inf)
        {
          a.j <- numeric(ni)
        }
        if (k != Inf)
        {
          indic <- ifelse(jinf + 1 <= 1 & jsup >= 1, 1, 0)
          
          a.j <- -k * pbinom(jinf, 1, probab) +
            k * (1 - pbinom(pmin(jsup, ni), 1, probab)) +
            1 / sqrt(V) * ifelse(ni == 1, probab * indic,
                                 mu * (pbinom(pmin(jsup - 1, ni - 1), 
                                              pmax(0, 1), probab) - 
                                         pbinom(jinf - 1, pmax(0, 1), probab))) - 
            mu / sqrt(V) * (pbinom(pmin(jsup, 1), 1, probab) - pbinom(jinf, 1, probab))
        }
        
        a.j <- 2 * a.j * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))
        
        #we define a part of w_j
        w <- drop((mu * (1 - mu) / scale) * c(w.x))
        
        tmp <- psi_q(resid, scale, k, weights, q = q[i])
        
        #we compute psi_q(r )-E(psi_q(r ))
        A <- (tmp - a.j)
        
        #We compute the -E(psi_prime)
        res1 <- drop(1 - mu)
        res0 <- drop(mu)
        
        B.tmp <- psi_q(res1, scale, k, weights, q = q[i])
        
        tmp <- psi_q(res0, scale, k, weights, q = q[i])
        
        B <- drop(drop(V) * (B.tmp + tmp))
        
        #We estimate betas
        temp <- coef + solve(t(X * w * B) %*% X) %*% t(X * w) %*% A
        
        coef <- temp
        eta <- X %*% coef
        fit <- exp(eta) / (1 + exp(eta))
        resid <- y - fit
        convi <- irls.delta(coef.old, coef)
        done <- (convi <= acc)
        
        if (done)
          break
      }
      
      if (!done)
        warning(paste("MQlogit failed to converge in", maxit, "steps at q = ", q[i]))
      
      qest[, i] <- coef
      qwt[, i]  <- w
      qfit[, i] <- fit
      qres[, i] <- resid
      qiter[i] <- iiter
    }
  }  
  else  ## poisson
  {
    for (i in 1:length(q)) {
      #We define the starting values
      resid <- resid.init
      fit <- fit.init
      phi <- phi.init
      a.j <- weights
      w <- weights
      coef <- temp.glmrob$coef
      
      for (iiter in 1:maxit) {
        resid.old <- resid 
        coef.old <- coef
        
        # We define the probability mu=t*exp(xb)
        probab <- fit
        mu <- probab
        deriv.mu <- mu
        #We define the variance
        V <- phi * probab
        #We define the scale
        scale <- c(sqrt(V))
        #We standardize the residuals
        r.stand <- (y - mu) / sqrt(V)
        #we compute i1 and i2
        jinf <- floor(mu - k * sqrt(V))
        jsup <- floor(mu + k * sqrt(V))
        
        #We compute the values of a_j(b)
        if (k == Inf)
        {
          a.j <- rep(1, n)
        }
        if (k != Inf)
        {
          a.j <- (-k) * ppois(jinf, mu) + k * (1 - ppois(jsup, mu)) +
            mu / sqrt(V) * (ppois(jinf, mu) - ppois(jinf - 1, mu) - 
                              (ppois(jsup, mu) - ppois(jsup - 1, mu)))
        }
        a.j <- 2 * a.j * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))
        
        #we define a part of w_j
        w <- diag(c(mu) / scale) * diag(c(w.x))
        
        #we compute psi_q(res)
        
        tmp <- psi_q(resid, scale, k, weights, q = q[i])
        
        #we compute psi_q(r )-E(psi_q(r ))
        A <- (tmp - a.j) 
        if (k == Inf)
        {
          esp.carre.cond <- rep(1, n)
        }
        if (k != Inf)
        {
          esp.carre.cond <- k * (ppois(jinf, mu) - 
                                   ppois(jinf - 1, mu) + 
                                   (ppois(jsup, mu) - ppois(jsup - 1, mu))) + 
            (mu^2 /  V^(3 / 2)) * 
            (ppois(jinf - 1, mu) - ppois(jinf - 2, mu) - 
               (ppois(jinf, mu) - ppois(jinf - 1, mu)) - 
               (ppois(jsup - 1, mu) - ppois(jsup - 2, mu)) + 
               (ppois(jsup, mu) -  ppois(jsup -  1, mu))) + 
            (mu / V^(3 / 2)) * (ppois(jsup - 1, mu) - ppois(jinf, mu))
        }
        b.j <- 2 * esp.carre.cond * (q[i] * (r.stand > 0) + (1 - q[i]) * (r.stand <= 0))
        B <- diag(c(V * b.j))
        
        #We estimate betas
        temp <- coef + solve(t(X) %*% w %*% B %*% X) %*% t(X) %*% w %*% A
        coef <- temp
        eta <- X %*% coef
        
        if (is.null(offset)) {
          fit <- exp(eta)
        } else {
          fit <- offset * exp(eta)
        }
        
        
        resid <- y - fit
        convi <- irls.delta(coef.old, coef)
        conv <- c(conv, convi)
        done <- (convi <= acc)
        
        if (done)
          break
      }
      
      if (!done)
        warning(paste("MQPoisson failed to converge in", maxit,
                      "steps at q = ", q[i]))
      
      qest[, i] <- coef
      qwt[, i]  <- diag(w)
      qfit[, i] <- fit
      qres[, i] <- resid
      qiter[i] <- iiter
    }
    
  }
  
  ########################################################
  return(list(
    fitted.values = qfit,
    residuals = qres,
    q.values = q,
    q.weights = qwt,
    coefficients = qest,
    qiter = qiter
  ))
  
}

