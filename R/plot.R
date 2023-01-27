## ID: plot.R, last updated 2021-03-05, F.Osorio

envelope.student <- function(object, reps = 50, conf = 0.95, plot.it = TRUE)
{ ## simulated envelope
  envel <- function(n, mean, Sigma, covType, eta, reps, conf) {
    conf <- 1 - conf
    # initialize progress bar
    cat(" Progress:\n")
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
    elims <- matrix(0, nrow = n, ncol = reps)
    for (i in 1:reps) {
      x <- rmt(n, mean = mean, Sigma = Sigma, eta = eta)
      fit <- studentFit(x, family = Student(eta = eta), covStruct = covType)
      z <- WH.student(x, center = fit$center, cov = fit$Scatter, eta = fit$eta)
      elims[,i] <- sort(z)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
    band <- matrix(0, nrow = n, ncol = 2)
    for (i in 1:n)
      band[i,] <- quantile(elims[i,], probs = c(conf / 2, 1 - conf / 2))
    band
  }

  n <- object$dims[1]
  covType <- object$covariance$type
  z <- WH.student(object$x, center = object$center, cov = object$Scatter, eta = object$eta)

  if (plot.it) {
    band  <- envel(n, object$center, object$Scatter, covType, object$eta, reps, conf)
    ylim <- range(z, band)
    qqnorm(z, ylim = ylim, main = "Transformed distances Q-Q plot")
    par(new = TRUE)
    qqnorm(band[,1], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
    par(new = TRUE)
    qqnorm(band[,2], axes = F, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
  }
  invisible(list(transformed = z, envelope = band))
}
