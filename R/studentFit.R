## ID: studentFit.R, last updated 2022-08-21, F.Osorio

studentFit <-
function(x, data, family = Student(eta = .25), covStruct = "UN", subset, na.action, control = MVT.control())
{ ## multivariate-t fitter
  Call <- match.call()
  if (missing(x))
    stop("'x' is not supplied")
  if (inherits(x, "formula")) {
    mt <- terms(x, data = data)
    if (attr(mt, "response") > 0)
      stop("response not allowed in formula")
    attr(mt, "intercept") <- 0
    mf <- match.call(expand.dots = FALSE)
    names(mf)[names(mf) == "x"] <- "formula"
    mf$family <- mf$covStruct <- mf$control <- NULL
    mf[[1L]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    na.act <- attr(mf, "na.action")
    z <- model.matrix(mt, mf)
  } else {
    z <- as.matrix(x)
    if (!missing(subset))
      z <- z[subset, , drop = FALSE]
    if (!missing(na.action))
      z <- na.omit(z)
    else
      z <- na.fail(z)
  }
  if (!is.numeric(z))
    stop("studentFit applies only to numerical variables")
  znames <- dimnames(z)[[2]]
  dz <- dim(z)
  n <- dz[1]
  p <- dz[2]
  storage.mode(z) <- "double"
  switch(covStruct,
         "UN" = {
           covType <- "UN"
           covName <- "Unstructured"
           covStruct <- 0
         },
         "CS" = {
           covType <- "CS"
           covName <- "Equicorrelation"
           covStruct <- 1
         },
         "DIAG" = {
           covType <- "DIAG"
           covName <- "Diagonal"
           covStruct <- 2
         },
         "HOMO" = {
           covType <- "HOMO"
           covName <- "Variance homogeneous"
           covStruct <- 3
         },
         stop("structure not implemented.")
  )

  ## initial estimates
  o <- cov.weighted(z)
  center <- o$mean
  Scatter <- o$cov
  distances <- Mahalanobis(z, center, Scatter)

  ## extract family info
  if (!inherits(family, "Student.family"))
    stop("Use only with 'Student.family' objects")
  if (is.null(family$family))
    stop("'family' not recognized")
  kind <- family$kind
  if ((kind < 0) || (kind > 1))
    stop("not valid 'family' object")
  settings <- c(kind, unlist(family$pars))

  ## set control values
  if (missing(control))
    control <- MVT.control()
  if (kind == 0)
    control$maxiter <- 0
  ctrl <- unlist(control)
  ctrl <- c(ctrl, 0)

  ## initialization of CS and HOMO cases
  sigma2 <- 0
  rho <- 0
  Phi <- double(p * p)

  if (covStruct == 1) {
    ## initial estimates: equicorrelation case
    sigma2 <- sum(diag(Scatter)) / p
    rho <-  2 * sum(Scatter[lower.tri(Scatter)]) / (sigma2 * p * (p - 1))
    Scatter <- corCS(rho, p)
    Scatter <- sigma2 * Scatter
  }

  if (covStruct == 3) {
    ## initial estimates: variance homogeneity case
    Phi <- cov2cor(Scatter)
    sigma2 <- sum(diag(solve(Phi, Scatter))) / p
    Scatter <- sigma2 * Phi
    ctrl <- c(ctrl, 0)
  }

  ## Call fitter
  now <- proc.time()
  fit <- .C("fitter_dispatcher",
            z = z,
            dims = as.integer(dz),
            settings = as.double(settings),
            covStruct = as.integer(covStruct),
            center = as.double(center),
            Scatter = as.double(Scatter),
            sigma2 = as.double(sigma2),
            rho = as.double(rho),
            Phi = as.double(Phi),
            distances = as.double(distances),
            weights = as.double(rep(1, n)),
            logLik = double(1),
            control = as.double(ctrl))
  speed <- proc.time() - now

  ## creating the output object
  out <- list(call = Call,
              x = z,
              dims = dz,
              family = family,
              settings = fit$settings,
              start = list(center = center, Scatter = Scatter),
              center = fit$center,
              Scatter = matrix(fit$Scatter, ncol = p),
              sigma2 = NULL,
              rho = NULL,
              Phi = NULL,
              logLik = fit$logLik,
              numIter = fit$control[4],
              control = control,
              weights = fit$weights,
              distances = fit$distances,
              covariance = list(type = covType, name = covName),
              speed = speed,
              converged = FALSE)
  names(out$center) <- znames
  dimnames(out$Scatter) <- list(znames, znames)
  eta <- out$settings[2]
  if (!control$fix.shape && kind != 0) {
    eta <- signif(eta, 4)
    out$family$call <- call(out$family$family, eta = eta)
  }
  out$eta <- out$settings[2]
  if (kind == 0)
    out$converged <- TRUE
  else if (out$numIter < control$maxiter)
    out$converged <- TRUE
  if (covStruct == 1) {
    out$sigma2 <- fit$sigma2
    out$rho <- fit$rho
  }
  if (covStruct == 3) {
    out$sigma2 <- fit$sigma2
    out$Phi <- matrix(fit$Phi, ncol = p)
  }

  class(out) <- "studentFit"
  out
}

print.studentFit <-
function(x, digits = 4, ...)
{
  ## local functions
  print.symmetric <- function(z, digits = digits, ...) {
    # print upper triangular part of covariance matrix
    ll <- lower.tri(z, diag = TRUE)
    z[ll] <- format(z[ll], ...)
    z[!ll] <- ""
    print(z, ..., quote = F)
  }
  cat("Call:\n")
  x$call$family <- x$family$call
  dput(x$call, control = NULL)
  if (x$converged)
    cat("Converged in", x$numIter, "iterations\n")
  else
    cat("Maximum number of iterations exceeded")
  cat("\nCenter:\n ")
  print(format(round(x$center, digits = digits)), quote = F, ...)
  cat("\nScatter matrix estimate:\n")
  if (x$dims[2] <= 5)
    print.symmetric(x$Scatter, digits = digits)
  else {
    print.symmetric(x$Scatter[1:5,1:5], digits = digits)
    cat("...")
  }
  nobs <- x$dims[1]
  cat("\nNumber of Observations:", nobs, "\n")
  invisible(x)
}
