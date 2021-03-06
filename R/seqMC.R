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
batchSeqMC <- function(f, logprob_y_given_x, x0, y,
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
#' @param logprob_y_given_x function, when called with parameter t (time point), 
#' 		y_t (observation vector at time t) and x_t (state vector at time t), it would return
#' 		the conditional log_probability: log(Prob(y_t | x_t ))
#' @param x0 matrix, sample of state vector at time 0, each col is a sample of state at time 0.
#' @param y0 observation at time 0 (can be missing).
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
#' logprob_y_given_x <- function(t, y, x) {
#' 	 as.numeric(-(y - x * x / 20.0)**2/2.0)
#' }
#'
#' x0 = matrix(rnorm(4000, sd=2), nrow=1, ncol=4000)
#'
#' mod = seqMC(f, logprob_y_given_x, x0)
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
#'    estimate.seqMC(mod)
#' })
#' 
#' plot(x, ylim=c(-40, 40), pch='*')
#' lines(xhat[1,])
#' lines(xhat[2,], lty='dotted')
#' lines(xhat[3,], lty='dotted')
#'
#' @export
seqMC <- function(f, logprob_y_given_x, x0, y0, sample_method=c("systematic", "residual", "bootstrap")) {
	stopifnot(is.matrix(x0)) # x0 is matrix of N cols, with N = number of samples (particles)
	sample_method = match.arg(sample_method)
	resampler = get(sprintf("%s_sample", sample_method))
	ans = structure(list(resampler=resampler,
		                 f=f,
						 logprob_y_given_x=logprob_y_given_x,
						 N=ncol(x0),
						 t=0,
						 x=x0,
						 w=rep(1, ncol(x0))), class="seqMC")
	if (!missing(y0)) {
		update_weights.seqMC(ans, y0)
	} else {
		ans
	}
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
#' @param y vector or matrix, if y is a vector, it represents an observation at a single time step,
#' if y is a matrix, then each col is an observation at a time point, number of cols in y equal to number of 
#' time steps.
#' @return seqMC object updated after given observation(s).
#' @export
update.seqMC <- function(obj, y) {
	if (missing(y)) {
		return(next_time.seqMC(obj))
	}
	y = if (!is.matrix(y)) matrix(y, ncol=1) else y
	for (col in 1:ncol(y)) {
		obj = resample.seqMC(obj)
		obj = next_time.seqMC(obj)
		obj = update_weights.seqMC(obj, y[, col])
	}
	obj
}

effective_sample_size <- function(w) {
	sum(w)**2 / sum(w**2)
}

resample.seqMC <- function(obj) {
	obj$x = obj$x[, obj$resampler(obj$w/obj$N), drop=FALSE]
	obj$w = rep(1, ncol(obj$x))		
	obj
}

next_time.seqMC <- function(obj) {
	obj$x = obj$f(obj$t, obj$x)
	obj$t = obj$t + 1
	obj
}

update_weights.seqMC <- function(obj, y) {
	logprob = obj$logprob_y_given_x(obj$t, y, obj$x)
	obj$w = obj$w * exp(logprob - max(logprob))
	# sum of obj$w should be equal to it's length.
	obj$w = obj$w * (length(obj$w) / sum(obj$w))
	obj
}

weightedMean <- function(x, weights=NULL) {
	if (!is.null(weights)) {
		sum(x*weights)/sum(weights)
	} else {
		mean(x)	
	}
}

#' @export
estimate.seqMC <- function(obj) {
	mean = apply(obj$x, 1, weightedMean, weight=obj$w)
	ci = apply(obj$x, 1, Hmisc::wtd.quantile, weights=obj$w, probs=c(0.025, 0.975))
	cbind(mean, t(ci))
}


