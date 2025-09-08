################################################################################

# Compute the Hessian from Bayesian EM (stage 1 of EMbru)

################################################################################

# Target function
target_hessian <- function(param, p, T, domain_area, link.functions, list.input, N, data.bru) {
  
  t_diff <- matrix(data.bru$t, N, N) - matrix(data.bru$t, N, N, byrow = TRUE)
  x_diff <- matrix(data.bru$x, N, N) - matrix(data.bru$x, N, N, byrow = TRUE)
  y_diff <- matrix(data.bru$y, N, N) - matrix(data.bru$y, N, N, byrow = TRUE)
  
  p_lower <- p - diag(diag(p), nrow = N, ncol = N)
  p_diag <- diag(p)
  sum_p_diag <- sum(p_diag)
  sum_p_lower <- sum(p_lower)
  sum_p_lower_t_diff <- sum(p_lower * t_diff)
  sum_p_lower_x_norm_2 <- sum(p_lower * (x_diff^2 + y_diff^2))
  
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
    sum(dnorm(param, 0, 1, log = TRUE))
}

# Function to compute the second partial derivative with respect to specific parameters
segunda_derivada_parametro <- function(f, params, p, T, domain_area, link.f.be, list.input, N, data.bru, epsilon = 1e-5) {
  n <- length(params)
  segunda_derivada <- matrix(0, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      # Function to calculate the fist derivate respect i
      primera_derivada_i <- function(params, i, epsilon) {
        params_plus <- params
        params_plus[i] <- params_plus[i] + epsilon
        params_minus <- params
        params_minus[i] <- params_minus[i] - epsilon
        return((f(params_plus, p, T, domain_area, link.f.be, list.input, N, data.bru) - f(params_minus, p, T, domain_area, link.f.be, list.input, N, data.bru)) / (2 * epsilon))
      }
      
      # Function to calculate the second derivate respect j
      segunda_derivada_j <- function(params, j, epsilon) {
        params_plus <- params
        params_plus[j] <- params_plus[j] + epsilon
        params_minus <- params
        params_minus[j] <- params_minus[j] - epsilon
        return((primera_derivada_i(params_plus, i, epsilon) - primera_derivada_i(params_minus, i, epsilon)) / (2 * epsilon))
      }
      
      # Calculate second partial derivate
      segunda_derivada[i, j] <- segunda_derivada_j(params, j, epsilon)
    }
  }
  return(segunda_derivada)
}  

# Load data
load("synthetic_data.RData")

# Result from Bayes EM (stage 1 of EMbru)
load("result_em_stage1EMbru.RData")
result_em <- resultado_em$result
param <- resultado_em$param

# Define necessary objects
source("code_em_stage1.R")
bdy <- st_sfc(st_polygon(list(matrix(c(0, 0, 1, 0, 1, 1, 0, 1, 0, 0), ncol = 2, byrow = TRUE))))
T2 <- ceiling(max(data.bru$ts))

# Priors
link.functions <- list(
  mu = \(x) gamma.t(x, 900, 30),
  K0 = \(x) loggaus.t(x, -0.3, 0.3),
  w = \(x)  gamma.t(x, 9, 3),
  sig = \(x) loggaus.t(x, -5, 2),
  mu.inv = \(x) gamma.t.inv(x, 900, 30),
  K0.inv = \(x) loggaus.t.inv(x, -0.3, 0.3),
  w.inv = \(x)  gamma.t.inv(x, 9, 3),
  sig.inv = \(x) loggaus.t.inv(x, -5, 2)
)
link.f.be <- link.functions
# Compute Hessian

#Compute the second derivative (Hessian)
segunda_derivada <- segunda_derivada_parametro(target_hessian, param, p=resultado_em$p, T=T2, domain_area = as.numeric(st_area(bdy)) , link.f.be=link.f.be, list.input = resultado_em$list.input, N=nrow(data.bru), data.bru=data.bru)
print(-segunda_derivada)

#Take the negative of the Hessian (Fisher information matrix)
Hessian_EM <- -segunda_derivada
Hessian_EM[lower.tri(Hessian_EM)] <- 0

# Invert the Hessian
inv_Hessian_EM <- solve(Hessian_EM)

# standard errors:
sqrt(diag(inv_Hessian_EM))


################################################################################

# Compute the Hessian from inlabru (stage 2 of EMbru)

################################################################################

load("fit_inlabru_stage2EMbru.RData")
Hessian_EMbru <- fit_inlabru$misc$configs$config[[1]]$Q
inv_Hessian_EMbru <- fit_inlabru$misc$configs$config[[1]]$Qinv
sqrt(diag(inv_Hessian_EMbru))



