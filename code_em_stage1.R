# Packages
library(ggplot2)
library(MASS)
library(dplyr)
library(knitr)
library(mclust)
library(INLA)
library(inlabru)
library(parallel)
library(spatstat)

#######################################################################
# Parameter transformations
#######################################################################

# gamma copula transformation
gamma.t <- function(x, a, b) {
  bru_forward_transformation(qgamma, x, a, b)
}
# uniform copula transformation
unif.t <- function(x, a, b) {
  bru_forward_transformation(qunif, x, min = a, max = b)
}
# log-gaussian copula transformation
loggaus.t <- function(x, m, s) {
  bru_forward_transformation(qlnorm, x, meanlog = m, sdlog = s)
}

# exponential copula transformation
exp.t <- function(x, e) {
  bru_forward_transformation(qexp, x, rate = e)
}

# inversions

# gamma copula transformation
gamma.t.inv <- function(x, a, b) {
  bru_inverse_transformation(pgamma, x, a, b)
}
# uniform copula transformation
unif.t.inv <- function(x, a, b) {
  bru_inverse_transformation(punif, x, min = a, max = b)
}
# log-gaussian copula transformation
loggaus.t.inv <- function(x, m, s) {
  bru_inverse_transformation(plnorm, x, meanlog = m, sdlog = s)
}

# exponential copula transformation
exp.t.inv <- function(x, e) {
  bru_inverse_transformation(pexp, x, rate = e)
}

#######################################################################
# Scale change
#######################################################################

internal_to_natural <- function(param, link.functions, inverse = FALSE) {
  if (inverse) {
    values <- c(
      link.functions$mu.inv(param[1]),
      link.functions$K0.inv(param[2]),
      link.functions$w.inv(param[3]),
      link.functions$sig.inv(param[4])
    )
  } else {
    values <- c(
      link.functions$mu(param[1]),
      link.functions$K0(param[2]),
      link.functions$w(param[3]),
      link.functions$sig(param[4])
    )
  }
  if (inverse) {
    names(values) <- c("theta_mu", "theta_k0", "theta_w", "theta_sigma")
  } else {
    names(values) <- c("mu", "k0", "w", "sigma")
  }
  values
}

#######################################################################
# Functions to create ciruclar bands and time intervals
#######################################################################

# Check the points that are within the study region.
within_region <- function(row, poly) {
  xx <- row["xx"]
  yy <- row["yy"]
  radio <- row["xy_end"]
  
  circle_center <- st_point(c(xx, yy)) 
  circle <- st_buffer(st_sfc(circle_center), dist = radio)
  st_crs(circle) <- st_crs(bdy)
  is_inside <- st_within(circle, poly, sparse = FALSE)
  
  if (!is_inside) {
    return("No")
  } else {
    return("Yes")
  }
}
# Count the points that meet both conditions (inside band and polygon)
points_inside_polygon_and_band <- function(x_coords, y_coords,x_centre, y_centre, r , polygon_sf=bdy ) {
  polygon_owin <- as.owin(polygon_sf)
  # Calculate the size of the square window that includes the entire band.
  extra_space <- r[2]  # Outer radius defines the additional size.
  xrange <- c(x_centre - extra_space, x_centre + extra_space)
  yrange <- c(y_centre - extra_space, y_centre + extra_space)
  # Create the extended window.
  extended_window <- owin(xrange, yrange)
  points_ppp <- ppp(x_coords, y_coords, window = extended_window)
  inside_polygon <- inside.owin(points_ppp$x, points_ppp$y, polygon_owin)
  # Calculate the distances to the center and verify if they are within the band.
  sq_norm <- (x_coords - x_centre)^2 + (y_coords - y_centre)^2
  inside_band <- sq_norm >= r[1]^2 & sq_norm <= r[2]^2
  # Count the points that meet both conditions
  count <- sum(inside_polygon & inside_band)
  return(count)
}

# Weighted the integral
weight_integral2 <- function(row) {
  # Number of points on the grid (higher number for greater accuracy)
  n_points <- 10000
  r <- c(as.numeric(row["xy_start"]), as.numeric(row["xy_end"]))
  x_centre <- as.numeric(row["xx"])
  y_centre <- as.numeric(row["yy"])
  angulos = runif(n_points, 0, 2*pi)
  radiis = sqrt(runif(n_points, r[1]^2, r[2]^2))
  # Generate random points within the circular band.
  x_coords <- x_centre + radiis * cos(angulos)
  y_coords <- y_centre + radiis * sin(angulos)
  # Count how many points are inside the polygon boundary.
  points_inside_band <-
    sum(points_inside_polygon_and_band(x_coords, y_coords, x_centre, y_centre, r , polygon_sf=bdy))
  # Calculate the proportion of points that fall inside the square.
  area_band_approx <- points_inside_band / n_points
  return(area_band_approx)
}

# Find points defining the bins for an observed point (temporal)
breaks_exp <- function(tt_, T2_, coef_, delta_, N_exp_ = 10) {
  tt_breaks <- tt_ + delta_ * ((1 + coef_)^(0:N_exp_))
  tt_breaks <- tt_breaks[tt_breaks < T2]
  if (T2_ - tt_ < delta_) {
    return(c(tt_, T2_))
  }
  if (T2 - tt_breaks[length(tt_breaks)] < delta_) {
    tt_breaks[length(tt_breaks)] <- T2_
  }
  if (tt_breaks[length(tt_breaks)] < T2_) {
    tt_breaks <- c(tt_breaks, T2_)
  }
  return(c(tt_, tt_breaks))
}

# Create the time grid
time.grid <- function(data.point, coef.t, delta.t,
                      T2., displaygrid = FALSE, N.exp.) {
  tt. <- data.point$ts
  idx.p <- data.point$idx.p
  # time bins
  # find bins break points
  t_b <- breaks_exp(tt., T2., coef_ = coef.t, delta_ = delta.t, N_exp_ = N.exp.)
  
  time.bins <- data.frame(
    t.start = t_b[-length(t_b)],
    t.end = t_b[-1]
  ) %>%
    mutate(t.bin.name = paste0(round(t.start, 3), "-", round(t.end, 3)))
  
  if (nrow(time.bins) - 1 == 0) {
    time.bins$t.ref_layer <- paste0("last-", idx.p)
  } else {
    time.bins$t.ref_layer <- c(1:(nrow(time.bins) - 1), paste0("last-", idx.p))
  }
  time.bins <- cbind(time.bins, data.point, row.names = NULL)
  time.bins
}

#######################################################################
# Functions to compute the integral
#######################################################################

# Temporal integral
It_df <- function(param_, time.df) {
  tth <- as.numeric(time.df$ts)
  T1b <- as.numeric(time.df$t.start)
  T2b <- as.numeric(time.df$t.end)
  param_w <- param_[3]
  T.l <- pmax(tth, T1b)
  1 / param_w * (exp(-param_w * (T.l - tth)) - exp(-param_w * (T2b - tth)))
}

# Spatial integral
Is_df <- function(param_, space.df) {
  r1 <- as.numeric(space.df$xy_start)
  r2 <- as.numeric(space.df$xy_end)
  param_sig <- param_[4]
  int <- NULL
  for (i in seq_along(r1)) {
    int[i] <- 2 * pi * param_sig^2 * (exp(-r1[i]^2 / (2 * param_sig^2)) - exp(-r2[i]^2 / (2 * param_sig^2)))
  }
  int
}

# Compute the temporal integral
compute.grid <- function(param., list.input_) {
  It.vec <- It_df(param_ = param., time.df = list.input_$time.sel)
  It.vec[list.input_$Imapping]
}

# Compute the weighted spatial integral
compute.grid_s <- function(param., list.input_) {
  Is.vec <- Is_df(param_ = param., space.df = list.input_$space.sel)
  Is.vec_weights <- Is.vec[list.input_$Imapping_s] * list.input_$weights
  Is.vec_weights
}

# Function to compute the whole integral
logLambda.i.inla <- function(theta_, list.input_, link.functions) {
  comp. <- compute.grid(param. = theta_, list.input_ = list.input_)
  comp.s <- compute.grid_s(param. = theta_, list.input_ = list.input_)
  
  out <- log(theta_[2]) + log(theta_[3]) - log(2 * pi * theta_[4]^2) + log(comp. + 1e-10) + log(comp.s + 1e-10)
  out
}

#######################################################################
# EXPECTATION-MAXIMIZATION
#######################################################################

# Functions for the Maxim likelihood

#######################################################################

# E-step
updatep <- function(t, x, y, p, k0, w, mu, sigma, cutoff) {
  lam <- c()
  N <- length(t)
  p[1, 1] <- mu
  lam[1] <- sum(p[1, 1])
  p[1, 1] <- p[1, 1] / sum(p[1, 1])
  
  for (i in 2:N) {
    j <- seq_len(i-1)
    j <- j[i - j < cutoff] # j subset
    # probability i triggered by j is proportional to triggering
    # kernel evaluated at inter-point times and distances
    p[i, j] <- k0 * w * exp(-w * (t[i] - t[j])) * exp((-(x[i] - x[j])^2 - (y[i] - y[j])^2) / (2 * sigma^2)) / (2 * pi * sigma^2)
    # probablity i is background event proportional to mu background rate
    p[i, i] <- mu
    # save intensity at each event for analysis
    lam[i] <- sum(p[i, 1:i])
    # normalize probabilities to sum to 1
    p[i, 1:i] <- p[i, 1:i] / sum(p[i, 1:i])
  }
  return(list(p = p, lam = lam))
}

# M-Step (update parameters)
updatepar_posterior <- function(param,t_diff,x_diff,y_diff,  p, T, domain_area, link.functions, list.input, N) {
  
  p_lower <- p - diag(diag(p), nrow = N, ncol = N)
  p_diag <- diag(p)
  sum_p_diag <- sum(p_diag)
  sum_p_lower <- sum(p_lower)
  sum_p_lower_t_diff <- sum(p_lower * t_diff)
  sum_p_lower_x_norm_2 <- sum(p_lower * (x_diff^2 + y_diff^2))
  # Objective function to be maximized
  target <- function(param, T, domain_area, link.functions, list.input) {
    mu <- link.functions$mu(param[1])
    k0 <- link.functions$K0(param[2])
    w <- link.functions$w(param[3])
    sigma <- link.functions$sig(param[4])
    theta_ <- c(mu, k0, w, sigma)
    
    sum_triggering_integral <- sum(
      exp(
        logLambda.i.inla(theta_,
                         list.input_ = list.input,
                         link.functions = link.functions
        )
      )
    ) 
    log(mu) * sum_p_diag - mu * T * domain_area +
      log(k0 * w / (2 * pi * sigma^2)) * sum_p_lower -
      w * sum_p_lower_t_diff -
      sum_p_lower_x_norm_2 / (2 * sigma^2) -
      sum_triggering_integral+
      sum(dnorm(param, 0, 1, log = TRUE)) #funtion to be maximized
  }
  
  opt <-
    optim(
      par = param,
      target,
      T = T,
      domain_area = domain_area,
      link.functions = link.functions,
      list.input=list.input,
      method = "BFGS",
      control = list(fnscale = -1),
    )
  message(paste0("Convergence: ", opt$convergence))

  return(list(
    mu_total = link.functions$mu(opt$par[1]),
    k0 = link.functions$K0(opt$par[2]),
    w = link.functions$w(opt$par[3]),
    sigma = link.functions$sig(opt$par[4]),
    param = opt$par
  ))
}

# Function to check if matrix p has stopped changing significantly.
# We will use this function for our EM process stopping criterion.

matrix_stopped_changing <- function(p_new, p_old, tol = 1e-6) {
  max_norm_diff <- norm(p_new - p_old, type = "M") # Calculate the maximum norm of the difference.
  has_stopped <- max_norm_diff < tol
  return(list(has_stopped = has_stopped, max_norm_diff = max_norm_diff))
}

#######################################################################
# EM
#######################################################################

#' Title
#'
#' @param data.bru data.frame of observations - columns must have names x, y, ts 
#' @param coef.t.  scalar - time-bin param
#' @param delta.t. scalar - time-bin param
#' @param N.max. scalar - time-bin param
#' @param coef_ scalar - space-bin param
#' @param delta_ scalar - space-bin param
#' @param N_exp_ scalar - space-bin param
#' @param initial list of initial values for the parameters - must have names k0, w, mu, sigma
#' @param Tend scalar - end of the time interval (0, Tend)
#' @param link.functions list of link functions representing the priors
#' @param cutoff scalar - parameter for EM
#' @param max_iter scalar - parameter for EM
#'
#' @return list with four elements: 
#' result -  data.frame with posterior mode of the parameters
#' show_convergence - data.frame of parameters value per iteration
#' list,input - list of inputs containing the space-time grid and the mapping
#' p - matrix number of obs x number of obs (not clear what it represents)

run_EM <- function(data.bru,
                   bdy, 
                   T2,
                   initial = list(k0 = .5, w = 1, mu = 10, sigma = .01),
                   link.functions,
                   coef.t., 
                   delta.t., 
                   N.max.,
                   coef_,
                   delta_,
                   N_exp_,
                   cutoff = 1000,
                   max_iter = 100 # number of maximum iterations for EM
                   
){
  
  t <- data.bru$ts
  x <- data.bru$x 
  y <- data.bru$y
  domain_area = as.numeric(st_area(bdy))
  
  df.time <- NULL
  for (idx in 1:nrow(data.bru)) {
    # print(idx)
    result <- time.grid(
      data.point = data.bru[idx, ],
      coef.t = coef.t.,
      delta.t = delta.t.,
      T2. = T2,
      N.exp. = N.max.,
    )
    df.time <- rbind(df.time, result)
  }
  
  radio <- delta_ * ((1 + coef_)^(0:N_exp_))
  radio
  xy_start <- c(0, radio[-(N_exp_ + 1)])
  xy_end <- radio
  print("Creating circular bins... This may take several minutes...")
  df.space <- list()
  for (i in 1:dim(data.bru)[1]) {
    #print(paste0(c("id:",i)))
    xx <- data.bru$x[i]
    yy <- data.bru$y[i]
    df.space[[i]] <- as.data.frame(cbind(
      idx.p = i, xx, yy, xy_start, xy_end,
      s.ref_layer = seq_along(xy_start)
    ))
    df.space[[i]]$within <- apply(df.space[[i]], 1, function(row) within_region(row, bdy))
    weight_values <- rep(1, nrow(df.space[[i]])) #Create a vector of NA with the length of df
    if (any(df.space[[i]]$within == "No")) {
      # If there are such rows, assign the values calculated by apply() only to those rows
      weight_values[df.space[[i]]$within == "No"] <- apply(df.space[[i]][df.space[[i]]$within == "No", ], 1, weight_integral2)
    }
    df.space[[i]]$weight <- weight_values
  }
  df.space <- do.call(rbind, df.space)
  
  df.j <- merge(df.time, df.space, by = "idx.p")
  df.j <- df.j %>% dplyr::select(t.start, t.end, t.bin.name, t.ref_layer, idx.p, x, y, ts, xy_start, xy_end, weight, s.ref_layer)

  # input that will be stored in list.input and needed for the efficient calculation of the expected number of triggered events
  t.names <- unique(df.j$t.ref_layer)
  time.sel <-
    df.j[vapply(t.names, \(bname) match(TRUE, df.j$t.ref_layer == bname), 0L), , drop = FALSE]
  Imapping <- match(df.j$t.ref_layer, t.names)
  
  s.names <- unique(df.j$s.ref_layer)
  space.sel <-
    df.j[vapply(s.names, \(bname) match(TRUE, df.j$s.ref_layer == bname), 0L), , drop = FALSE]
  Imapping_s <- match(df.j$s.ref_layer, s.names)

  list.input <- list(
    df_grid = df.j,
    Imapping = Imapping,
    Imapping_s = Imapping_s,
    time.sel = time.sel,
    space.sel = space.sel,
    weights = df.j$weight
  )
  
  N <- length(t)
  t_diff <- matrix(t, N, N) - matrix(t, N, N, byrow = TRUE)
  x_diff <- matrix(x, N, N) - matrix(x, N, N, byrow = TRUE)
  y_diff <- matrix(y, N, N) - matrix(y, N, N, byrow = TRUE)
  T <- T2
  # p is a matrix storing branching probabilities
  p <- matrix(0, N, N)
  
  k0 <- initial$k0 
  w <- initial$w 
  mu <- initial$mu 
  sigma <- initial$sigma
  
  param <-
    internal_to_natural(c(mu, k0, w, sigma),
                        link.functions = link.functions,
                        inverse = TRUE)
  
  print('Starting EM')
  show_convergence <- c()
  medir_t <- proc.time()
  # number of EM iterations
  for (k in 1:max_iter) {
    # E-step
    t_iter <-proc.time()
    p_old <- p
    res_updatep <- updatep(t = t, x = x, y = y, p = p, k0 = k0, w = w, mu = mu, sigma = sigma, cutoff = cutoff)
    p <- res_updatep$p
    lam <- res_updatep$lam
    # M-step
    res_updatepar <- updatepar_posterior(
      param,
      t_diff=t_diff, x_diff= x_diff, y_diff=y_diff, p = p, T = T, domain_area = domain_area,
      link.functions = link.functions,
      list.input=list.input, N
    )
    k0 <- res_updatepar$k0
    w <- res_updatepar$w
    sigma <- res_updatepar$sigma
    mu <- res_updatepar$mu_total
    param <- res_updatepar$param
    # result is the outcome, that is, the output, the estimated parameters in each iteration
    resultado <- data.frame(k0 = k0, w = w, sigma = sigma, mu = mu)
    stop_criterion <- matrix_stopped_changing(p, p_old, tol = 0.001)
    show_convergence <- rbind(show_convergence, cbind(n_iter = k, resultado, diff_p_norm=round(stop_criterion$max_norm_diff, 4), time_iteration=round((proc.time()-t_iter)[3], 3)))
    print(unlist(show_convergence[k, ]))
    
    if (stop_criterion$has_stopped) { #stopping criterion
      print(paste("Stop criterion reached at iteration", k))
      break
    }
  }
  total_t <- proc.time() - medir_t 
  result_em <- resultado #outcome from EM algorithm 
  return(list(result = result_em, show_convergence = show_convergence,
              list.input = list.input, p = p, param=param, total_t=total_t) )
  
}


