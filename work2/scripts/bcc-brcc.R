FF <- function(x,Digits=4,Width=4){(formatC(x,digits=Digits,width=Width,format="f"))}

# Two-Parameter Beta Distribution -----------------------------------------
beta.bc <- function(n, mle)
{
  p <- mle[1]; q <- mle[2]
  bias.p <- -(0.5e0 * trigamma(p + q) ^ 2 * psigamma(p, 2) - trigamma(q) * trigamma(p + q) * psigamma(p, 2) - 0.5e0 * trigamma(q) ^ 2 * psigamma(p + q, 2) + 0.5e0 * trigamma(q) ^ 2 * psigamma(p, 2) - 0.5e0 * psigamma(p + q, 2) * trigamma(p) * trigamma(q) - 0.5e0 * trigamma(p + q) ^ 2 * psigamma(q, 2) + 0.5e0 * trigamma(p + q) * trigamma(p) * psigamma(q, 2)) * gamma(q) * gamma(p) / n / gamma(p + q) / ((trigamma(p) + trigamma(q)) * trigamma(p + q) - trigamma(p) * trigamma(q)) ^ 2 / beta(p, q)
  bias.q <- (0.5e0 * trigamma(p + q) ^ 2 * psigamma(p, 2) - 0.5e0 * trigamma(q) * trigamma(p + q) * psigamma(p, 2) + 0.5e0 * psigamma(p + q, 2) * trigamma(p) * trigamma(q) - 0.5e0 * trigamma(p + q) ^ 2 * psigamma(q, 2) + trigamma(p + q) * trigamma(p) * psigamma(q, 2) + 0.5e0 * trigamma(p) ^ 2 * psigamma(p + q, 2) - 0.5e0 * trigamma(p) ^ 2 * psigamma(q, 2)) * gamma(q) * gamma(p) / n / gamma(p + q) / ((trigamma(p) + trigamma(q)) * trigamma(p + q) - trigamma(p) * trigamma(q)) ^ 2 / beta(p, q)
  c(alpha = bias.p, beta = bias.q)
}


# Beta Control Chart with Bias Correction ---------------------------------

bcc <- function(y, alpha = 0.005, bias = FALSE)
{
  fit     <- fitdist(data = y, distr = "beta")
  thetas  <- coef(fit)
  if(bias == TRUE)
  {
    pdf    <- quote(gamma(theta1 + theta2) / (gamma(theta1) * gamma(theta2)) * x ^ (theta1 - 1) * (1 - x) ^ (theta2 - 1))
    lpdf   <- quote(lgamma(theta1 + theta2) - lgamma(theta1) - lgamma(theta2) + theta1 * log(x) + theta2 * log(1 - x))
    cs     <- coxsnell.bc(density = pdf, logdensity = lpdf, n = length(y), parms = c("theta1", "theta2"), mle = thetas, lower = 0, upper = 1)
    # thetas <- thetas - beta.bc(n = length(y), mle = as.vector(thetas))
  }
  LCL     <- qbeta(p = alpha / 2, shape1 = thetas[1], shape2 = thetas[2])
  UCL     <- qbeta(p = 1 - (alpha / 2), shape1 = thetas[1], shape2 = thetas[2])
  CL      <- as.vector(thetas[1] / (thetas[1] + thetas[2]))
  limites <- c(LCL, CL, UCL)
  return(limites)
}


# Beta Regression Control Chart -------------------------------------------

brcc <- function(y, X, Z, alpha = 0.005, diagnostic = FALSE, npdf)
{
  y <- as.vector(y)
  X <- as.matrix(cbind(1, X))
  Z <- as.matrix(cbind(1, Z))
  
  n <- length(y)
  r <- ncol(X)
  s <- ncol(Z)
  
  loglike <- function(parms, y)
  {
    beta   <- parms[1:r]
    gamma  <- parms[(r + 1):(r + s)]
    
    eta1   <- X %*% beta
    eta2   <- Z %*% gamma
    
    mu     <- exp(eta1) / (1 + exp(eta1))
    sigma  <- exp(eta2) / (1 + exp(eta2))
    aux1   <- (1 - sigma ^ 2) / sigma ^ 2
    
    theta1 <- mu * aux1
    theta2 <- (1 - mu) * aux1
    
    ll    <- sum(dbeta(x = y, shape1 = theta1, shape2 = theta2, log = TRUE))
    
    return(ll)
  }
  
  
  ystar    <- log(y / (1 - y))
  betaols  <- lm.fit(X, ystar)$coefficients
  etaols   <- X %*% betaols
  muols    <- exp(etaols) / (1 + exp(etaols))
  varols   <- as.numeric(t(ystar - etaols) %*% (ystar-etaols)) / ( (n - r) * ((1 / (muols * (1 - muols))) ^ 2))
  phiols   <- (muols * (1 - muols) / varols ) - 1
  sigmaini <- sqrt(1 / (1 + phiols))
  gamaols  <- lm.fit(Z, log(sigmaini / (1 - sigmaini)))$coefficients
  ini      <- c(betaols, gamaols)
  
  opt <- optim(ini, loglike, method = "BFGS", control = list(fnscale = -1, maxit = 500, reltol = 1e-9), hessian = TRUE, y = y)
  
  if(opt$conv != 0) return(warning("FUNCTION DID NOT CONVERGE!"))
  
  beta  <- opt$par[1:r]
  gamma <- opt$par[(r + 1):(r + s)]
  eta1  <- X %*% beta
  eta2  <- Z %*% gamma
  mu    <- as.vector(exp(eta1) / (1 + exp(eta1)))
  sigma <- as.vector(exp(eta2) / (1 + exp(eta2)))
  aux1  <- (1 - sigma ^ 2) / sigma ^ 2
  
  theta1 <- mu * aux1
  theta2 <- (1 - mu) * aux1
  
  LCL     <- qbeta(p = alpha / 2, shape1 = theta1, shape2 = theta2)
  UCL     <- qbeta(p = 1 - (alpha / 2), shape1 = theta1, shape2 = theta2)
  limites <- data.frame(LCL = LCL, UCL = UCL)
  
  final <- list()
  
  ## Resumo dos parâmetros estimados
  coef     <- opt$par
  vcov     <- solve(-opt$hessian)
  stderror <- sqrt(diag(vcov))
  lowerci  <- coef - stderror * qnorm(0.975)
  upperci  <- coef + stderror * qnorm(0.975)
  zstat    <- abs(coef / stderror)
  pvalues  <- 2 * (1 - pnorm(zstat))
  model_summary <- cbind(coef, stderror, lowerci, upperci, pvalues)
  colnames(model_summary) <- c("Estimate","Std. Error","Lower CI", "Upper CI", "Pr(>|z|)")
  
  ## Critérios de discriminação
  loglik <- opt$value
  p      <- length(opt$par)
  aic    <- -2 * loglik + 2 * p
  aicc   <- -2 * loglik + (2 * n * p) / (n - p - 1)
  bic    <- -2 * loglik + p * log(n)
  model_measures <- c(loglik, aicc, bic)
  
  ## Exibição dos resultados
  print(model_summary)
  print(" ",quote=F)
  print(c("Log-likelihood:", round(loglik, 4)), quote = F)
  print(c("AIC:", round(aic, 4), " AICc:", round(aicc, 4), " BIC:", round(bic, 4)), quote = F)
  
  ## Resíduos do modelo
  my.residuals <- function(y, mu, sigma, type = "standard")
  {
    if(type == "standard")
    {
      vary <- mu * (1 - mu) * sigma
      res  <- (y - mu) / vary
    }
    if(type == "rpp2")
    {
      phi    <- (1 - sigma ^ 2)/ sigma ^ 2
      PHI    <- diag(phi)
      ystar  <- log(y / (1 - y))
      mustar <- digamma(mu * phi) - digamma((1 - mu) * phi)
      glmu   <- (1 / (1 - mu) + mu / (1 - mu) ^ 2) / (mu /(1 - mu))
      nu     <- trigamma(mu * phi) + trigamma((1 - mu) * phi)
      W      <- diag(phi * nu * (1 / glmu ^ 2))
      H      <- sqrt(W %*% PHI) %*% X %*% solve(t(X) %*% PHI %*% W %*% X) %*% t(X) %*% sqrt(PHI %*% W)
      h      <- diag(H)
      res   <- (ystar - mustar) / (nu * (1 - h))
    }
    return(res)
  }
  
  
  if(diagnostic == TRUE)
  {
    M   <- 99
    i   <- 1
    mat <- matrix(nrow = n, ncol = M)
    while (i <= M) 
    {
      ysim <- rbeta(n = n, shape1 = theta1, shape2 = theta2)
      out  <- optim(par = as.vector(coef), method = "BFGS", fn = loglike, control = list(fnscale = -1, maxit = 500, reltol = 1e-9), y = ysim)
      if(out$conv == 0) 
      {
        beta     <- out$par[1:r]
        gamma    <- out$par[(r + 1):(r + s)]
        eta1     <- X %*% beta
        eta2     <- Z %*% gamma
        mut      <- as.vector(exp(eta1) / (1 + exp(eta1)))
        sigmat   <- as.vector(exp(eta2) / (1 + exp(eta2)))
        mat[, i] <- sort(abs(my.residuals(y = ysim, mu = mut, sigma = sigmat, type = "rpp2"))) 
        i <- i + 1 
      }
      # cat(i, "\n")
    }
    l           <- 1:n
    rpp2.obs    <- sort(abs(my.residuals(y = y, mu = mu, sigma = sigma, type = "rpp2")))
    rpp2.teo    <- qnorm((l + n - 1/8) / (2 * n + 0.5))
    rpp2.min    <- apply(mat, 1, quantile, probs = 0.025)
    rpp2.max    <- apply(mat, 1, quantile, probs = 0.975)
    rpp2.median <- apply(mat, 1, median)  
    
    Ry <- c(min(rpp2.min), max(rpp2.max))
    Rx <- range(rpp2.teo)
    
    pdf(paste0("res-", npdf, ".pdf"))
    par(mar = c(3.2, 3.2, 1.0, 1.0), cex = 1.8)
    plot(x = rpp2.teo, y = rpp2.obs, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', xlim = Rx, ylim = Ry, cex = 0.4, 
         bty = 'o', pch = 19)
    lines(x = rpp2.teo, y = rpp2.min)
    lines(x = rpp2.teo, y = rpp2.max)
    lines(x = rpp2.teo, y = rpp2.median, lty = 2)
    axis(1, seq(Rx[1], Rx[2], l = 5), FF(seq(Rx[1], Rx[2], l = 5), 2))
    axis(2, seq(Ry[1], Ry[2], l = 5), FF(seq(Ry[1], Ry[2], l = 5), 2))
    mtext("Percentil da Normal(0, 1)", side = 1.1, line = 2.2, cex = 1.8)
    mtext("Resíduos padronizados ponderados 2",   side = 2.1, line = 2.2, cex = 1.8)
    graphics.off()
  }
  
  final$crtlim   <- limites
  final$model    <- model_summary
  # final$measures <- model_measures
  
  return(final)
}

# Regression control chart ------------------------------------------------

rcc <- function(model, data, alpha = 0.005)
{
  n       <- nrow(data)
  fit     <- lm(formula = model, data = data)
  mut     <- fit$fitted.values
  k       <- length(fit$coefficients)
  sigma   <- sqrt(1 / (n - k - 1) * sum(fit$residuals^2))
  LCL     <- qnorm(p = alpha / 2, mean = mut, sd = sigma)
  UCL     <- qnorm(p = 1 - alpha / 2, mean = mut, sd = sigma)
  limites <- data.frame(LCL = LCL, UCL = UCL)
  return(limites)
}


# p-Control Chart with corrections ----------------------------------------

pcontrol <- function(prop, n, correction = c("none", "seno", "ryan", "chen", "joekes"))
{
  correction  <- match.arg(correction)
  p           <- mean(prop)
  CL.shewhart <- p + c(-1, 1) * 3 * sqrt(p * (1 - p) / n)
  CL.seno     <- asin(sqrt(p)) + c(-1, 1) * 3 * sqrt(1 / (4 * n))
  CL.ryan     <- CL.shewhart + 1.25 / n
  CL.chen     <- CL.shewhart + 4 / (3 * n) * (1 - 2 * p)
  CL.joekes   <- CL.chen - (p * (1 - p) + 2) / (6 * n ^ 2 * sqrt((p * (1 - p) / n)))
  
  CL <- switch(correction,
               none   = CL.shewhart,
               seno   = CL.seno,
               ryan   = CL.ryan,
               chen   = CL.chen,
               joekes = CL.joekes)
  CL <- c(CL[1], p, CL[2])
  
  return(CL)
}

