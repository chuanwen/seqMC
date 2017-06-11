#' Sequential Monte Carlo
#'
#' @param f function, when called with parameter t (time point) and 
#' 		x_t (state vector at time t), it would return x_(t+1)
#' @param prob_y_given_x function, when called with parameter t (time point), 
#' 		y_t (observation vector at time t) and x_t (state vector at time t), it would return
#' 		the conditional probability/density: Prob(y_t | x_t )
#' @param x0 matrix of N cols, sample of state vector at time 0, one sample vector per col,
#'      N is the number of samples.
#' @param y matrix of T cols, observations, col 1 is observation at time 1, col 2 is
#'		observation at time 2, ... etc. T is the number of time points.
#' @param sample_method character, specify sample method in the resample stage.
#'      Default systematic, means "systematic resampling".
#' @return sample from posterior distribution of state vectors, a 3D array, with
#' 		dimension of d x N x T, where d is the length of a state vector, N is the number of samples,
#'		T is the number of time steps.
#'
#' @examples
#'
#' f <- function(t, x) {
#'	 0.5 * x + 25 * x / (1 + x * x) + 8.0 * cos(1.2 * (t-1)) + rnorm(length(x), sd=sqrt(10.0))
#' }
#'
#' prob_y_given_x <- function(t, y, x) {
#' 	 as.numeric(dnorm(y - x * x / 20.0))
#' }
#'
#' ### simulate true path ###
#' T = 50
#' x = rep(0.0, T+1)
#' for (t in 1:T) {
#' 	 x[t+1] = f(t, x[t])
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
	N = ncol(x0) # x0 is matrix of N cols, with N = number of samples (particles)
	T = ncol(y)  # y is matrix of T cols, with T = number of time points.
	
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
	
	x = list(T, mode="list")
	for (t in 1:T) {
		xt = f(t, x0)
		w = normalize(prob_y_given_x(t, y[,t], xt))
		stopifnot(is.matrix(xt), is.vector(w), ncol(xt) == length(w))
		# resample cols of xt with probability w #
		xt = xt[, resampler(w), drop=FALSE]
		x[[t]] = xt # save posterior distribution/sample at time t.
		x0 = xt
	}
	x = unlist(x)
	dim(x) = c(dim(x0), T)
	invisible(x)
}
