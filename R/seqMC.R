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

#' Sequential Monte Carlo
#'
#' @inheritParams seqMC
#' @param y matrix of T cols, observations, col 1 is observation at time 1, col 2 is
#'		observation at time 2, ... etc. T is the number of time points.
#' @return sample from posterior distribution of state vectors, a 3D array, with
#' 		dimension of d x N x T, where d is the length of a state vector, N is the number of samples,
#'		T is the number of time steps.
#'
#' @export
batchSeqMC <- function(f, prob_y_given_x, x0, y,
	              sample_method=c("systematic", "residual", "bootstrap")) {
	stopifnot(is.matrix(x0), is.matrix(y))
	sample_method = match.arg(sample_method)
	mod = seqMC(f, prob_y_given_x, x0, sample_method)
	
	T = ncol(y)  # y is matrix of T cols, with T = number of time points.
	x = list(T, mode="list")
	for (t in 1:T) {
		mod = update(mod, y[, t])
		x[[t]] = mod$x
	}
	x = unlist(x)
	dim(x) = c(dim(x0), T)
	invisible(x)
}

#' seqMC creates a seqMC object
#'
#' @param f function, when called with parameter t (time point) and 
#' 		x_t (state vector at time t), it would return x_(t+1)
#' @param prob_y_given_x function, when called with parameter t (time point), 
#' 		y_t (observation vector at time t) and x_t (state vector at time t), it would return
#' 		the conditional probability: Prob(y_t | x_t )
#' @param x0 matrix, sample of state vector at time 0, each col is a sample of state at time 0.
#' @param sample_method character, specify sample method in the resample stage.
#'      Default systematic, means "systematic resampling".
#' @return a seqMC object, which can be updated at each time point. e.g. \code{obj = update(ojb, y)}
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
#' x0 = matrix(rnorm(4000, sd=2), nrow=1, ncol=4000)
#'
#' mod = seqMC(f, prob_y_given_x, x0)
#'
#' ### simulate true path ###
#' x0 = 0.1
#' T = 50
#' x = rep(0.0, T)
#' x[1] = f(0, x0)
#' for (t in 1:(T-1)) {
#' 	 x[t+1] = f(t, x[t])
#' }
#' y = x * x / 20 + rnorm(length(x))
#'
#' ### estimate the posterior of state vector given y[t] ####
#' xhat = sapply(1:T, function(t) {
#'    mod <<- update(mod, y[t])
#'    c(mean(mod$x), quantile(mod$x, c(0.025, 0.975)))
#' })
#' 
#' plot(x, ylim=c(-40, 40), pch='*')
#' lines(xhat[1,])
#' lines(xhat[2,], lty='dotted')
#' lines(xhat[3,], lty='dotted')
#'
#' @export
seqMC <- function(f, prob_y_given_x, x0, sample_method=c("systematic", "residual", "bootstrap")) {
	stopifnot(is.matrix(x0)) # x0 is matrix of N cols, with N = number of samples (particles)
	sample_method = match.arg(sample_method)
	resampler = get(sprintf("%s_sample", sample_method))
	structure(list(resampler=resampler, f=f, prob_y_given_x=prob_y_given_x, N=ncol(x0), t=0, x=x0), class="seqMC")
}


#' @export
print.seqMC <- function(obj, ...) {
	obj$x = NULL
	print.default(obj, ...)
}

#' Update seqMC object after an observation (if y is not missing) or simply update the object to next time point
#' (if y is missing).
#'
#' @param obj seqMC object.
#' @param y observation at time obj$t + 1.
#' @return seqMC object updated to next time point.
#' @export
update.seqMC <- function(obj, y) {
	obj$x = obj$f(obj$t, obj$x)
	obj$t = obj$t + 1
	if(!missing(y)) {
		w = normalize(obj$prob_y_given_x(obj$t, y, obj$x))
		obj$x = obj$x[, obj$resampler(w), drop=FALSE]
	}
	obj
}


