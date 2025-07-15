#library(simDAG)
#library(ggforce)
#set.seed(1234)

# Documentation : https://cran.r-project.org/web/packages/simDAG/simDAG.pdf
# Documentation : https://cran.r-project.org/web/packages/simDAG/vignettes/v_sim_from_dag.html

#dag <- empty_dag()
#dag <- dag + 
#  node("S", type = "rbernoulli", p = 0.5, output="numeric") +
#  node("X1", type="gaussian", parents=c("S"), 
#       betas=c(1.1), intercept=1, error=.5) +
#  node("X2", type="gaussian", parents=c("S", "X1"), betas=c(0.1, 0.3), intercept=-2, error=.5) #+
#node("X3", type="multinomial", parents=c("S", "X1", "X2"),
#     betas=matrix(c(0.2, 0.4, 0.1), ncol=3), intercepts=1,
#     output = "factor") # ncol for categ nb

#plot(dag) # to get the causal graph
#summary(dag) # to get the equations

#sim_dat <- sim_from_dag(dag=dag, n_sim=1000, keep_simguide = TRUE)
#sim_dat$

#dag2 <- do(dag, names="S", values=1)
#plot(dag2)
#sim_dat <- sim_from_dag(dag=dag, n_sim=5)
#sim_dat2 <- sim_from_dag(dag=dag2, n_sim=5)

########################################################################################################
library(R6causal)

# A multinomial sampling function returning category 1, 2, or 3
# probs depend on S, X1, and X2 through a simple softmax linear model
multinomial_sample <- function(S, X1, X2, U) {
  # Linear predictors for each category (logits)
  logits <- matrix(
    data = c(
      0.1 + 0.3 * S + 0.2 * X1 + 0.1 * X2 + U,
      -0.2 + 0.1 * S - 0.1 * X1 + 0.4 * X2 + U,
      rep(0, length(S))
    ),  # reference category baseline,
    nrow = length(S),
    ncol = 3
  )
  exp_logits <- exp(logits)
  probs <- exp_logits / apply(exp_logits, 1, sum)
  
  # Sample category from multinomial with these probs
  sample_row <- function(p) sample(1:3, size = 1, prob = p)
  apply(probs, 1, sample_row)
}

# Background variables have gaussian dist. --> exogeneous variables
model <- SCM$new(
  name = "simple SCM",
  uflist = list(
    us  = function(n) {return(rnorm(n))},
    ux1 = function(n) {return(rnorm(n))},
    ux2 = function(n) {return(rnorm(n))},
    ux3 = function(n) {return(rnorm(n))}
  ),
  vflist = list(
    s  = function(us) {return(as.numeric(us < 0.4))},
    x1 = function(ux1, s) {return(as.numeric(-2 + 0.5*s + ux1))},
    x2 = function(ux2, s, x1) {return(as.numeric(1 + 0.4*s + 0.4*x1 + ux2))},
    x3 = function(ux3, s, x1, x2) {return(multinomial_sample(s, x1, x2, ux3))}
  )
)

model$plot()
model$simulate(n = 100, seed = 1234)
model$simdata

# A simple intervention
model_cf <- model$clone()  # making a copy
model_cf$intervene("s", 1) # applying the intervention
model_cf$simulate(n = 100, seed = 1234)
model_cf$simdata
model_cf$plot()
