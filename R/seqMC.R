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
#' @param sample_method character, specify sample method in the resample stage.
#'      Default systematic, means "systematic resampling".
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
seqMC <- function(f, prob_y_given_x, x0, y,
	              sample_method=c("systematic", "residual", "bootstrap")) {
	stopifnot(is.matrix(x0), is.matrix(y))
	sample_method = match.arg(sample_method)
	N = ncol(x0)
	K = ncol(y)
	
	normalize <- function(x) {
		x / sum(x)
	}
	
	bootstrap_sample <- function(prob, n=length(prob)) {
		sample.int(length(prob), size=n, replace=TRUE, prob=prob)
	}
	
	residual_sample <- function(prob, n=length(prob)) {
		size = length(prob)
		count = n * prob
		count_int = floor(count)
		n_int = sum(count_int)
		ans = rep(1:size, count_int)
		if (n_int < n) {
			n_left = n - n_int
			prob_left = (count - count_int) / n_left
			ans = c(ans, bootstrap_sample(prob_left, n_left))
		}
		ans
	}
	
	systematic_sample <- function(prob, n=length(prob)) {
		zeta = floor(n * cumsum(c(0.0, prob)) + runif(1))
		m = zeta[-1] - zeta[-length(zeta)]
		rep.int(1:length(prob), m)
	}
	
	resampler = get(sprintf("%s_sample", sample_method))
	
	x = list(K, mode="list")
	for (k in 1:K) {
		xk = f(k, x0)
		q = normalize(prob_y_given_x(k, y[,k], xk))
		stopifnot(is.matrix(xk), is.vector(q), ncol(xk) == length(q))
		# resample cols of xk with probability q #
		xk = xk[, resampler(q), drop=FALSE]
		x[[k]] = xk 
		x0 = xk
	}
	x = unlist(x)
	dim(x) = c(dim(x0), K)
	invisible(x)
}
