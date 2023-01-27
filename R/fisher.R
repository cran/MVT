## ID: fisher.R, last updated 2021-04-09, F.Osorio

## this is an internal function that performs the computation of the
## Fisher information matrix for the multivariate-t model 
fisher.matrix <-
function(object) {
  ## local functions
  c.eta <- function(eta) eta / (1 - 2 * eta)
  c.phi <- function(eta, p) (1 + p * eta) / (1 + (p + 2) * eta)
  c.mu <- function(eta, p) c.phi(eta, p) / (1 - 2 * eta)
  beta.dot <- function(eta, p) {
    dif <- trigamma(.5 * (1 + p * eta) / eta) - trigamma(.5 / eta)
    -.5 * dif / eta^2
  }

  ## extract elements from 'object'
  Scatter <- object$Scatter
  eta <- object$eta

  ## checking elements
  p <- ncol(Scatter)
  if (nrow(Scatter) != p)
    stop("'Scatter' must be a square dispersion matrix")
  if (!isSymmetric(Scatter))
    Scatter <- asSymmetric(Scatter)

  ## invert and vectorizing Scatter matrix
  inv <- solve(Scatter)
  if (!isSymmetric(inv))
    inv <- asSymmetric(inv)
  vec.Scatter <- as.vector(inv)

  ## Fisher information matrix about 'center'
  fisher.center <- c.mu(eta, p) * inv
  if (!isSymmetric(fisher.center))
    fisher.center <- asSymmetric(fisher.center)

  ## Fisher information matrix about 'Scatter'
  fisher.Scatter <- 2 * c.phi(eta, p) * symm.prod(n = p, kronecker.prod(inv, inv), side = "right")
  fisher.Scatter <- fisher.Scatter + (c.phi(eta, p) - 1) * outer(vec.Scatter, vec.Scatter)
  fisher.Scatter <- .25 * dupl.cross(n = p, k = p, fisher.Scatter)
  if (!isSymmetric(fisher.Scatter))
    fisher.Scatter <- asSymmetric(fisher.Scatter)

  ## crossed Fisher information about 'Scatter' and 'eta'
  fisher.cross <- -(c.eta(eta) * (p + 2)) / ((1 + p * eta) * (1 + (p + 2) * eta))
  fisher.cross <- fisher.cross * dupl.prod(n = p, vec.Scatter, transposed = TRUE, side = "left")

  ## Fisher information about 'eta'
  fisher.eta <- 1 + p * eta * (1 - 4 * eta) - 8 * eta^2
  fisher.eta <- fisher.eta / ((1 + p * eta) * (1 + (p + 2) * eta))
  fisher.eta <- fisher.eta * p / (1 - 2 * eta)^2
  fisher.eta <- -.5 * (fisher.eta - beta.dot(eta, p)) / eta^2

  ## forming the Fisher information matrix
  Dp.cols <- p * (p + 1) / 2
  fisher <- matrix(0, nrow = p + Dp.cols + 1, ncol = p + Dp.cols + 1)
  fisher[1:p, 1:p] <- fisher.center
  rows <- cols <- seq.int(from = p + 1, length.out = Dp.cols)
  fisher[rows, cols] <- fisher.Scatter
  cols <- ncol(fisher)
  fisher[rows, cols] <- fisher[cols, rows] <- fisher.cross
  fisher[cols, cols] <- fisher.eta

  fisher
}
