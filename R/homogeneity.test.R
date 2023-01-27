## ID: homogeneity.test.R, last updated 2021-04-09, F.Osorio

homogeneity.test <-
function(object, test = "LRT")
{
  ## local functions
  restrictedScore <-
  function(z, weights, center, Scatter) {
    n <- nrow(z)
    p <- ncol(z)
    y <- sweep(z, 2, center)
    y <- sqrt(weights) * y
    S <- crossprod(y)
    s <- solve(Scatter, S) - n * diag(p)
    s <- s %*% solve(Scatter)
    s <- .5 * dupl.prod(n = p, x = as.vector(s), transposed = TRUE, side = "left")
    s <- as.vector(s)
    s
  }

  ## extract object info
  fit <- object
  if (!inherits(fit, "studentFit"))
    stop("Use only with 'studentFit' objects")
  covType <- fit$covariance$type
  if (covType != "UN")
    stop("Only implemented for objects with 'unstructured' Scatter.")
  n <- fit$dims[1]
  p <- fit$dims[2]
  eta <- fit$eta

  ## restricted parameter estimation
  f0 <- studentFit(fit$x, family = Student(eta = eta), covStruct = "HOMO", control = fit$control)
  null.fit <- list(call = fit$call, center = f0$center, Scatter = f0$Scatter, sigma2 = f0$sigma2,
                   Phi = f0$Phi, eta = f0$eta, distances = f0$distances, weights = f0$weights,
                   logLik = f0$logLik, numIter = f0$numIter)

  switch(test,
          LRT = {
            stat <- 2. * (fit$logLik - f0$logLik)
            names(stat) <- "LRT"
            method <- "Likelihood ratio test"
          },
          Wald = {
            phi  <- fit$Scatter[lower.tri(fit$Scatter, diag = TRUE)]
            phi0 <- f0$Scatter[lower.tri(f0$Scatter, diag = TRUE)]
            dif  <- phi - phi0
            aCov <- fisher.matrix(fit)[-(1:p),-(1:p)]
            rows <- cols <- seq.int(from = 1, length.out = p * (p + 1) / 2)
            acol <- ncol(aCov)
            if (fit$eta != 0.0)
              aCov <- aCov[rows,cols] - outer(aCov[rows,acol], aCov[rows,acol]) / aCov[acol,acol]
            else
              aCov <- aCov[rows,cols]
            h <- aCov %*% dif
            stat <- n * sum(dif * h)
            names(stat) <- "Wald"
            method <- "Wald test"
          },
          score = {
            s <- restrictedScore(fit$x, f0$weights, f0$center, f0$Scatter)
            aCov <- fisher.matrix(f0)[-(1:p),-(1:p)]
            rows <- cols <- seq.int(from = 1, length.out = p * (p + 1) / 2)
            acol <- ncol(aCov)
            if (f0$eta != 0.0)
              aCov <- aCov[rows,cols] - outer(aCov[rows,acol], aCov[rows,acol]) / aCov[acol,acol]
            else
              aCov <- aCov[rows,cols]
            h <- solve(aCov, s)
            stat <- sum(s * h) / n
            names(stat) <- "Score"
            method <- "Score test"
          },
          gradient = {
            phi  <- fit$Scatter[lower.tri(fit$Scatter, diag = TRUE)]
            phi0 <- f0$Scatter[lower.tri(f0$Scatter, diag = TRUE)]
            dif <- phi - phi0
            s <- restrictedScore(fit$x, f0$weights, f0$center, f0$Scatter)
            stat <- sum(s * dif)
            names(stat) <- "Gradient"
            method <- "Gradient test"
          },
          stop(paste("unimplemented test:", test))
          )
  df <- p - 1
  pval <- 1 - pchisq(stat, df = df)

  ## output object
  dimnames(f0$Scatter) <- dimnames(fit$Scatter)
  z <- list(statistic = stat, parameter = df, p.value = pval, estimate = fit$Scatter,
            null.value = f0$Scatter, method = method, null.fit = null.fit)
  if (is.null(fit$call$data))
    z$data <- fit$call$x
  else
    z$data <- fit$call$data
  z$family <- fit$family
  class(z) <- "homogeneity.test"
  z
}

print.homogeneity.test <- function(x, digits = 4, ...)
{
  ## local functions
  print.symmetric <-
  function(z, digits = digits, ...)
  {
    ll <- lower.tri(z, diag = TRUE)
    z[ll] <- format(z[ll], ...)
    z[!ll] <- ""
    print(z, ..., quote = F)
  }

  cat("\n")
  cat(paste(x$method, "for equality of variances", sep = " "), "\n")
  cat("\n")
  cat("data:", x$data, "\n")
  out <- character()
  out <- c(out, paste(names(x$statistic), "statistic =",
                      format(round(x$statistic, digits = digits))))
  out <- c(out, paste("df =", x$parameter))
  out <- c(out, paste("p-value =", format(round(x$p.value, digits = digits))))
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("alternative hypothesis: true variances are not equal.\n")
  #print(x$null.value, ...)
  cat("\n")
  cat("sample estimate:\n")
  print.symmetric(x$estimate, digits = digits)
  invisible(x)
}
