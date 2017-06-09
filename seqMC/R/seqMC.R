#' Sequential Monte Carlo
#'
#' @param f function, when called with parameter k (time point) and 
#' 		x_k (state vector at time k), it would return x_{k+1}
#' @param prob_y_given_x function, when called with parameter k (time point), 
#' 		y_k (observation vector at time k) and x_k (state vector at time k), it would return
#' 		the conditional probability/density: Prob(y_k | x_k )
#' @param x0 matrix, sample of state vector at time 0, one sample vector per col.
#' @param y matrix, observations, col 1 is observation at time 1, col 2 is
#'		observation at time 2, ... etc.
#' @return sample from posterior distribution of state vectors, a 3D array, with
#' 		dimension of d x N x K, where d is the length of a state vector, N is the number of samples,
#'		K is the number of time steps.
#'
#' @examples
#'
#' f <- function(k, x) {
#'	 0.5 * x + 25 * x / (1 + x * x) + 8.0 * cos(1.2 * (k-1)) + rnorm(length(x), sd=sqrt(10.0))
#' }
#'
#' prob_y_given_x <- function(k, y, x) {
#' 	 as.numeric(dnorm(y - x * x / 20.0))
#' }
#'
#' ### simulate true path ###
#' K = 50
#' x = rep(0.0, K+1)
#' for (k in 1:K) {
#' 	 x[k+1] = f(k, x[k])
#' }
#' x = x[-1]
#' y = x * x / 20 + rnorm(length(x))
#'
#' ### estimate the posterior of state vector ####
#' N = 4000
#' x0 = matrix(rnorm(N, sd=sqrt(2)), nrow=1, ncol=N)
#' xhat = seqMC(f, prob_y_given_x, x0, matrix(y, nrow=1))
#' xhat_mean = apply(xhat, 3, mean)
#'
#' alpha = 0.05
#' xhat_ci_lower = apply(xhat, 3, quantile, probs=alpha/2)
#' xhat_ci_upper = apply(xhat, 3, quantile, probs=1-alpha/2)
#' 
#' plot(x[-1], ylim=c(-40, 40), pch='*')
#' lines(xhat_ci_lower, lty='dotted')
#' lines(xhat_mean)
#' lines(xhat_ci_upper, lty='dotted')
#'
#' @export
seqMC <- function(f, prob_y_given_x, x0, y) {
	stopifnot(is.matrix(x0), is.matrix(y))
	N = ncol(x0)
	K = ncol(y)
	
	normalize <- function(x) {
		x / sum(x)
	}
	
	bootstrap <- function(x, prob) {
		stopifnot(is.matrix(x), is.vector(prob), ncol(x) == length(prob))
		N = ncol(x)
		x[, sample.int(N, size=N, replace=TRUE, prob=prob), drop=FALSE]
	}
	
	systematic_resample <- function(x, prob) {
		stopifnot(is.matrix(x), is.vector(prob), ncol(x) == length(prob))
		N = ncol(x)
		t = floor(N * cumsum(c(0.0, prob)) + runif(1))
		m = t[-1] - t[-(N+1)]
		x[, rep.int(1:N, m), drop=FALSE]
	}
	
	x = list(K, mode="list")
	for (k in 1:K) {
		xk = f(k, x0)
		q =  normalize(prob_y_given_x(k, y[,k], xk))
		# xk = bootstrap(xk, q)
		xk = systematic_resample(xk, q)
		x[[k]] = xk 
		x0 = xk
	}
	x = unlist(x)
	dim(x) = c(dim(x0), K)
	invisible(x)
}





