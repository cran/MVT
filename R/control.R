## ID: control.R, last updated 2021-03-05, F.Osorio

MVT.control <-
function(maxiter = 2000, tolerance = 1e-6, fix.shape = FALSE)
{ ## control parameters for EM algorithm
  list(maxiter = maxiter, tolerance = tolerance, fix.shape = fix.shape)
}
