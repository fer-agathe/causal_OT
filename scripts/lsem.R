# Load libraries
library(mediation)

set.seed(2024)
n <- 10000

# Binary treatment
A <- rbinom(n, 1, 0.5)

# X1 = a1 * A + noise
a1 <- 0.6
X1 <- a1 * A + rnorm(n, mean = 0, sd = .25)

# X2 = a2 * A + b2 * X1 + noise
a2 <- 0.3
b2 <- 0.4
X2 <- a2 * A + b2 * X1 + rnorm(n, mean = 0, sd = .25)

# X3 = a3 * A + b3 * X1 + noise
a3 <- 0.7
b3 <- 0.2
X3 <- a3 * A + b3 * X1 + rnorm(n, mean = 0, sd = .25)

# Y = a4 * A + b4 * X1 + c4 * X2 + d4 * X3 + noise
a4 <- 0.2
b4 <- 0.5
c4 <- 0.4
d4 <- 0.1
Y <- a4 * A + b4 * X1 + c4 * X2 + d4 * X3 + rnorm(n, mean = 0, sd = .25)

data <- data.frame(A, X1, X2, X3, Y)

# Regressions
fit_X1 <- lm(X1 ~ A, data = data)
a1_hat <- fit_X1$coefficients[2] # True value: 0.6
fit_X2 <- lm(X2 ~ A + X1, data = data)
a2_hat <- fit_X2$coefficients[2] # True value: 0.3
b2_hat <- fit_X2$coefficients[3] # True value: 0.4
fit_X3 <- lm(X3 ~ A + X1, data = data)
a3_hat <- fit_X3$coefficients[2] 
b3_hat <- fit_X3$coefficients[3]
fit_Y  <- lm(Y ~ A + X1 + X2 + X3, data = data)
a4_hat <- fit_Y$coefficients[2] # True value: 0.2
b4_hat <- fit_Y$coefficients[3] # True value: 0.5
c4_hat <- fit_Y$coefficients[4] # True value: 0.4
d4_hat <- fit_Y$coefficients[5]

# Direct effect
direct <- a4
direct_hat <- a4_hat
abs(direct - direct_hat)

# Indirect effect A -> X1 -> Y
indirect_X1 <- b4 * a1 
indirect_X1_hat <- b4_hat * a1_hat
abs(indirect_X1 - indirect_X1_hat)

# Indirect effect A -> X2 -> Y
indirect_X2 <- c4 * a2 
indirect_X2_hat <- c4_hat * a2_hat
abs(indirect_X2 - indirect_X2_hat)

# Indirect effect A -> X3 -> Y
indirect_X3 <- d4 * a3
indirect_X3_hat <- d4_hat * a3_hat
abs(indirect_X3 - indirect_X3_hat)

# Indirect effect A -> X1 -> X2 -> Y
indirect_X12 <- a1 * b2 * c4
indirect_X12_hat <- a1_hat * b2_hat * c4_hat
abs(indirect_X12 - indirect_X12_hat)

# Indirect effect A -> X1 -> X3 -> Y
indirect_X13 <- a1 * b3 * d4
indirect_X13_hat <- a1_hat * b3_hat * d4_hat
abs(indirect_X13 - indirect_X13_hat)

# Total indirect effect
indirect <- indirect_X1 + indirect_X2 + indirect_X3 + indirect_X12 + indirect_X13
indirect_hat <- indirect_X1_hat + indirect_X2_hat + indirect_X3_hat + indirect_X12_hat + indirect_X13_hat
abs(indirect - indirect_hat)

# Total effect 
fit_Y_total <- lm(Y ~ A, data = data)
total_hat_fitY <- fit_Y_total$coefficients[2]
total_hat <- direct_hat + indirect_hat
total <- direct + indirect
abs(total - total_hat)

# With mediation package
med_mod_X1 <- multimed(
  outcome = "Y", 
  med.main = "X1", 
  med.alt = c("X2", "X3"), 
  treat = "A", 
  data = data
)
# Indirect effect for X1: A -> X1 -> Y
delta_0_medX1 <- mean((med_mod_X1$d0.lb + med_mod_X1$d0.ub) / 2)
# Direct + Other indirect effects: A -> Y, A -> X2 -> Y, A -> X3 -> Y, A -> X1 -> X2 -> Y, A -> X1 -> X3 -> Y
zeta_1_medX1 <- mean((med_mod_X1$z1.lb + med_mod_X1$z1.ub) / 2)
# Total effect for X1
tot_effect_medX1 <- delta_0_medX1 + zeta_1_medX1

med_mod_X2 <- multimed(
  outcome = "Y", 
  med.main = "X2", 
  med.alt = c("X1", "X3"), 
  treat = "A", 
  data = data
)
# Indirect effect for X2: A -> X2 -> Y, A -> X1 -> X2 -> Y
delta_0_medX2 <- mean((med_mod_X2$d0.lb + med_mod_X2$d0.ub) / 2)
# Direct + Other indirect effects: A -> Y, A -> X1 -> Y, A -> X3 -> Y, A -> X1 -> X3 -> Y
zeta_1_medX2 <- mean((med_mod_X2$z1.lb + med_mod_X2$z1.ub) / 2)
# Total effect for X2
tot_effect_medX2 <- delta_0_medX2 + zeta_1_medX2

med_mod_X3 <- multimed(
  outcome = "Y", 
  med.main = "X3", 
  med.alt = c("X1", "X2"), 
  treat = "A", 
  data = data
)
# Indirect effect for X3: A -> X3 -> Y, A -> X1 -> X3 -> Y
delta_0_medX3 <- mean((med_mod_X3$d0.lb + med_mod_X3$d0.ub) / 2)
# Direct + Other indirect effects: A -> Y, A -> X1 -> Y, A -> X2 -> Y, A -> X1 -> X2 -> Y
zeta_1_medX3 <- mean((med_mod_X3$z1.lb + med_mod_X3$z1.ub) / 2)
# Total effect for X3
tot_effect_medX3 <- delta_0_medX3 + zeta_1_medX3


