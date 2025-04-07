library(sf)
library(dplyr)
library(ggplot2)
library(inlabru)
library(INLA)
library(tidyr)

#######################################################################
# Load data
#######################################################################

load("synthetic_data.RData")
summary(data.bru)

#######################################################################
# Create polygon: Unit square
#######################################################################

library(sf)
square <- st_sfc(st_polygon(list(matrix(c(0, 0,  
                                          1, 0,  
                                          1, 1,  
                                          0, 1,  
                                          0, 0),  
                                        ncol = 2, byrow = TRUE))))

bdy <- square

#######################################################################
# Plot data
#######################################################################

ggplot() +
  # Dibujar la ventana de observación con un borde negro más grueso
  geom_sf(data = square, fill = "lightgray", color = "black", linewidth = 1) +
  geom_point(data = data.bru[data.bru$ts>90,], aes(x = x, y = y), color = "black", size = 1) +
  # Expandir los límites para ver más espacio alrededor
  scale_x_continuous(limits = c(-0.2, 1.2)) +
  scale_y_continuous(limits = c(-0.2, 1.2)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),  # Tamaño de los números en el eje X
    axis.text.y = element_text(size = 14),  # Tamaño de los números en el eje Y
    axis.title.x = element_text(size = 16),  # Tamaño de la etiqueta del eje X
    axis.title.y = element_text(size = 16)   # Tamaño de la etiqueta del eje Y
  )

#######################################################################
# Settings
#######################################################################

T1 = 0
T2 = ceiling(max(data.bru$ts))
data.bru <- data.bru %>% mutate(idx.p = 1:nrow(data.bru)) %>% dplyr::select(idx.p, x, y, ts)

#######################################################################
# Prior distributions
#######################################################################

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

#######################################################################
# Temporal bins (parameters)
#######################################################################

coef.t. <- 1 # binning parameter (delta)
delta.t. <- .5 # binning parameter (delta)
N.max. <- 7 # binning parameter (n.max)
delta.t. * ((1 + coef.t.)^(0:N.max.))

#######################################################################
# Circular spatial bins (parameters)
#######################################################################

coef_ <- 0.31
delta_ <- 0.02
N_exp_ <- 16
radio <- delta_ * ((1 + coef_)^(0:N_exp_))
radio
xy_start <- c(0, radio[-(N_exp_ + 1)])
xy_end <- radio

#######################################################################
# Initial values
#######################################################################

k0 <- .5
w <- 1
mu <- 10
sigma <- .01

initial = list(k0=k0, w=w, mu=mu, sigma=sigma)

#######################################################################
# EXPECTATION-MAXIMIZATION (Stage 1 of EMbru)
#######################################################################

source("code_em_stage1.R")
set.seed(1234)
resultado_em <- run_EM(data.bru,
                       bdy, 
                       T2,
                       initial = initial,
                       link.functions,
                       coef.t., 
                       delta.t., 
                       N.max.,
                       coef_,
                       delta_,
                       N_exp_,
                       cutoff = 1000,
                       max_iter = 100) #outcome from EM algorithm 

show_convergence <- resultado_em$show_convergence
result_em <- resultado_em$result
param <- resultado_em$param

# Parameter estimates:

print(result_em)

# Plot convergence:

ggplot(show_convergence %>%
         pivot_longer(
           cols = c(mu, k0, w, sigma),
           names_to = "param",
           values_to = "value"
         )) +
  geom_line(aes(x = n_iter, y = value)) +
  facet_wrap(~param, scales = "free_y")

#######################################################################
# INLABRU (Stage 2 of EMbru)
#######################################################################

param_from_em <- param # starting values from EM
sample.s <- data.bru
source("code_inla_stage2.R")
set.seed(1234)
medir_t_inla <- proc.time()
fit_inlabru <- running_inla(init_param = param_from_em, 
                            sample.s = data.bru,
                            link.functions = link.functions, 
                            df.j.grid = resultado_em$list.input$df_grid,
                            T1 = T1, T2 = T2, 
                            bdy=bdy,
                            bru_verb = 4, bru_iter = 100, bru_max_step = 1.5, bru_rel_tol = 0.1)
total_t_inla <- proc.time()-medir_t_inla
inlabru_param = internal_to_natural(param = fit_inlabru$summary.fixed$mean, link.functions = link.functions)

# Posterior means of the parameters:

print(inlabru_param)

# Plot convergence:

bru_convergence_plot(fit_inlabru)
inlabru:::make_track_plots(fit_inlabru)$default

#######################################################################
# Results EMbru
#######################################################################

rbind(
  EM = internal_to_natural(param = param, link.functions = link.functions),
  inlabru_0.025 = internal_to_natural(param = fit_inlabru$summary.fixed$"0.025quant", link.functions = link.functions),
  inlabru = internal_to_natural(param = fit_inlabru$summary.fixed$mean, link.functions = link.functions),
  inlabru_0.975 = internal_to_natural(param = fit_inlabru$summary.fixed$"0.975quant", link.functions = link.functions)
) 

#######################################################################
# Plot posteriors
#######################################################################

# posterior mean
post.mu <- data.frame(inla.tmarginal(link.functions$mu, fit_inlabru$marginals.fixed$th.mu),param = 'mu', model = 'EMbru')
post.K0 <- data.frame(inla.tmarginal(link.functions$K0, fit_inlabru$marginals.fixed$th.K0), param = 'k0', model = 'EMbru')
post.w <- data.frame(inla.tmarginal(link.functions$w, fit_inlabru$marginals.fixed$th.w), param = 'w', model = 'EMbru')
post.sig <- data.frame(inla.tmarginal(link.functions$sig, fit_inlabru$marginals.fixed$th.sig), param = 'sigma', model = 'EMbru')

data_gg <- rbind(post.mu, post.K0, post.w, post.sig)

# EM point estimates
line_data <- data.frame(
  param = c("mu", "k0", "w", "sigma"),
  xintercept = unlist(c(resultado_em$result[[4]],resultado_em$result[[1]],resultado_em$result[[2]],resultado_em$result[3])),
  model="EM"
)
line_data$param <- factor(line_data$param, levels = unique(data_gg$param))
line_data$model <- factor(line_data$model)

ggplot(data_gg, 
       aes(x = x, y = y, color = model, linetype = model)) + 
  geom_line(size = 1.2) + 
  geom_vline(data = line_data, aes(xintercept = xintercept, linetype = model, color = model), 
             size = 1.2, show.legend = FALSE) +  # Añadido show.legend = FALSE
  facet_wrap(~ param, scales = 'free', labeller = label_parsed, ncol = 4) +  
  xlab("value") + 
  ylab("pdf") + 
  theme_bw() +
  scale_linetype_manual(
    values = c("EMbru" = "solid", "EM" = "dotted")
  ) +
  scale_color_manual(
    values = c("EMbru" = "#F8766D", "EM" = "#619CFF")
  ) +
  guides(color = guide_legend(title = ""), 
         linetype = guide_legend(title = ""))

#######################################################################
## We approximate the number of points by calculating the approximate value of the integral to assess the goodness of our estimates.
#######################################################################

inlabru_param = internal_to_natural(param = fit_inlabru$summary.fixed$mean, link.functions = link.functions)
domain_area = as.numeric(st_area(bdy))

th.mu <- internal_to_natural(param = fit_inlabru$summary.fixed$mean, link.functions = link.functions)[1]
th.K0 <- internal_to_natural(param = fit_inlabru$summary.fixed$mean, link.functions = link.functions)[2]
th.w <- internal_to_natural(param = fit_inlabru$summary.fixed$mean, link.functions = link.functions)[3]
th.sig <- internal_to_natural(param = fit_inlabru$summary.fixed$mean, link.functions = link.functions)[4]

lambda.N <- function(th.mu, th.K0, th.w, th.sig, T1, T2, domain_area, list.input, link.functions){
  theta_etas <- c(link.functions$mu(th.mu[1]),
                  link.functions$K0(th.K0[1]),
                  link.functions$w(th.w[1]),
                  link.functions$sig(th.sig[1]))
  theta_etas[1]*(T2 - T1)*domain_area + sum(exp(logLambda.i.inla(th = theta_etas,
                                                                 list.input_ = list.input, link.functions = link.functions)))
}

source("code_em_stage1.R")
lambda.N.post <- predict(fit_inlabru,
                         data.frame(), 
                         ~ lambda.N(th.mu, th.K0, th.w, th.sig,
                                    T1, T2, domain_area,
                                    resultado_em$list.input,
                                    link.functions))

c(lambda.N.post[1:5], true = nrow(data.bru))

