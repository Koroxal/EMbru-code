# Triggering function
gt <- function(th, t, ti, x, xi, y, yi) {
  output <- rep(0, length(ti))
  t.diff <- t - ti
  x.diff <- x - xi
  y.diff <- y - yi
  neg <- t.diff <= 0
  if (sum(!neg) > 0) {
    log.out <- log(th[2]) + log(th[3]) - th[3] * t.diff[!neg] - log(2 * pi * th[4]^2) - (x.diff^2 + y.diff^2) / (2 * th[4]^2)
    output[!neg] <- exp(log.out)
  } else {
    output
  }
  output
}

# Hawkes process conditional intensity
lambda_ <- function(th, t, ti.v, x, xi.v, y, yi.v, k) {
  if (is.null(ti.v) | all(ti.v > t)) {
    a <- th[1]
  } else {
    a <- th[1] + sum(gt(th, t, ti.v, x, xi.v, y, yi.v))
  }
  a
}
# Compute the integral of the triggering functions
logLambda.i.inla <- function(th.K0, th.w, th.sig, list.input_, link.functions) {
  theta_ <- c(
    0,
    link.functions$K0(th.K0[1]),
    link.functions$w(th.w[1]),
    link.functions$sig(th.sig[1])
  )
  # compute the integral efficiently for each bin of each observation
  comp. <- compute.grid(param. = theta_, list.input_ = list.input_) # temporal integral
  comp.s <- compute.grid_s(param. = theta_, list.input_ = list.input_) #spatial integral
  
  out <- log(theta_[2]) + log(theta_[3]) - log(2 * pi * theta_[4]^2) + log(comp. + 1e-10) + log(comp.s + 1e-10)
  out
}

# Function to calculate the Hawkes process conditional log-intensity for a set of observations (tt, xx, yy) given the history of the process. It requires also the parameters in the internal scale.
loglambda.inla <- function(th.mu, th.K0, th.w, th.sig, tt, xx, yy, th, link.functions) {
  # if no link.functions are provided
  if (is.null(link.functions)) {
    th.p <- c(th.mu[1], th.K0[1], th.w[1], th.sig[1])
  } else {
    th.p <- c(
      link.functions$mu(th.mu[1]),
      link.functions$K0(th.K0[1]),
      link.functions$w(th.w[1]),
      link.functions$sig(th.sig[1])
    )
  }
  out <- mean(unlist(mclapply(tt, function(x) {
    th_x <- th < x
    ids <- which(x == tt)  
    
    sapply(ids, function(id) {
      log(lambda_(th = th.p, t = x, ti.v = th[th_x], x = xx[id], xi.v = xx[th_x], y = yy[id], yi.v = yy[th_x]))
    })
  }, mc.cores = 5)))
  
  out
}

# Predictor function for the surrogate Poisson model
predictor.fun <- function(th.mu, th.K0, th.w, th.sig,
                          list.input, T1, T2, domain_area,
                          link.functions = NULL) {
  out <- rep(0, list.input$n)
  out[list.input$idx.bkg] <-
    log(link.functions$mu(th.mu[1])) + log(T2 - T1) + log(domain_area)
  out[list.input$idx.trig] <- logLambda.i.inla(
    th.K0 = th.K0, th.w = th.w,
    th.sig = th.sig,
    list.input_ = list.input,
    link.functions = link.functions
  )
  
  out[list.input$idx.sl] <- loglambda.inla(
    th.mu = th.mu, th.K0 = th.K0,
    th.w = th.w,
    th.sig = th.sig,
    tt = list.input$sample.s$ts,
    th=list.input$sample.s$ts,
    xx = list.input$sample.s$x,
    yy = list.input$sample.s$y,
    link.functions = link.functions
  )
  out
}

#' Title
#'
#' @param init_param vector of initial values of parameters in the internal scale (theta_mu, theta_k0, theta_w, theta_sigma)
#' @param sample.s data.frame of observations with columns idx.p, x, y, ts 
#' @param link.functions list of link functions representing the priors and their inverse.
#' @param df.j.grid data.frame representing the spatio-temporal grid used to approximate the integral of the triggering function
#' @param T_start numeric representing the start of the time interval
#' @param T_end numeric representing the end of the time interval
#' @param X_start numeric representing the start of the space interval (x-axis), 
#' @param X_end numeric representing the end of the space interval (x-axis), 
#' @param Y_start numeric representing the start of the space interval (y-axis), 
#' @param Y_end numeric representing the end of the space interval (y-axis)
#' @param bru_verb numeric setting the visual output of inlabru (default 4)
#' @param bru_iter numeric setting the maximum number of inlabru iterations (default 100)
#' @param bru_max_step numeric setting the maximum step for the internal routine (default 1.5)
#' @param bru_rel_tol numeric setting the relative tolerance for convergence (default 0.1)
#' @return
#' @export
#'
#' @examples
running_inla <- function(init_param, sample.s, link.functions, df.j.grid, 
                         T1, T2, bdy,
                         bru_verb = 4, bru_iter = 100, bru_max_step = 1.5, bru_rel_tol = 0.1){
  domain_area = as.numeric(st_area(bdy))
  # Initial values
  th.init <- list(
    th.mu = init_param[1],
    th.K0 = init_param[2],
    th.w = init_param[3],
    th.sig = init_param[4]
  )
  # inlabru options list
  bru.opt.list <- list(
    bru_verbose = bru_verb, # type of visual output
    bru_max_iter = bru_iter, # maximum number of iterations
    bru_method = list(max_step = bru_max_step, 
                      rel_tol = bru_rel_tol), # options
    bru_initial = th.init, # parameters initial values
    control.mode = list(x = init_param) # not sure why this is needed.
  ) 
  
  # create data.frame representing background part of the integrated intensity
  df.0 <- data.frame(counts = 0, exposures = 1, part = "background")
  
  # modify data.frame representing the triggered part of the integrated intensity
  df.j.grid$counts <- 0
  df.j.grid$exposures <- 1
  df.j.grid$part <- "triggered"
  # input that will be stored in list.input and needed for the efficient calculation of the expected number of triggered events
  
  t.names <- unique(df.j.grid$t.ref_layer)
  time.sel <-
    df.j.grid[vapply(t.names, \(bname) match(TRUE, df.j.grid$t.ref_layer == bname), 0L), , drop = FALSE]
  Imapping <- match(df.j.grid$t.ref_layer, t.names)
  
  s.names <- unique(df.j.grid$s.ref_layer)
  space.sel <-
    df.j.grid[vapply(s.names, \(bname) match(TRUE, df.j.grid$s.ref_layer == bname), 0L), , drop = FALSE]
  Imapping_s <- match(df.j.grid$s.ref_layer, s.names)
  
  
  # create data.frame representing the sum of the log intensities
  df.s <- data.frame(counts = nrow(sample.s), exposures = 0, part = "SL")
  
  # bind the data.frames together
  data.input <- bind_rows(df.0, df.s, df.j.grid)
  
  # create list.input
  list.input <- list(
    n = nrow(data.input),
    df_grid = df.j.grid,
    Imapping = Imapping,
    Imapping_s = Imapping_s,
    time.sel = time.sel,
    space.sel = space.sel,
    sample.s = sample.s,
    weights = df.j.grid$weight,
    idx.bkg = data.input$part == "background",
    idx.trig = data.input$part == "triggered",
    idx.sl = data.input$part == "SL"
  )
  
  # create formula representing the logintensity of the surrogate Poisson counts model
  merged.form <- counts ~ predictor.fun(
    th.mu = th.mu, th.K0 = th.K0,
    th.w = th.w,
    th.sig = th.sig,
    list.input = list.input,
    T1 = T1, T2 = T2,domain_area=domain_area,
    link.functions = link.functions
  )
  # create components representing the parameters in the internal scale
  cmp.part <- counts ~ -1 +
    th.mu(1, model = "linear", mean.linear = 0, prec.linear = 1) +
    th.K0(1, model = "linear", mean.linear = 0, prec.linear = 1) +
    th.w(1, model = "linear", mean.linear = 0, prec.linear = 1) +
    th.sig(1, model = "linear", mean.linear = 0, prec.linear = 1)
  
  fit_bru <- bru(
    components = cmp.part,
    like(
      formula = merged.form,
      data = data.input,
      family = "poisson",
      E = exposures
    ),
    options = bru.opt.list
  )
  return(fit_bru)
}


