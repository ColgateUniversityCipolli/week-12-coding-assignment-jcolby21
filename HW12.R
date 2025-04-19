#1 a
# (a) Critical value for T20 at alpha = 0.05
alpha <- 0.05
df_20 <- 20 - 1
t_crit_20 <- qt(1 - alpha, df = df_20)

#1 b 
# (b) Critical value for T30 at alpha = 0.05
alpha <- 0.05
df_30 <- 30 - 1
t_crit_30 <- qt(1 - alpha, df = df_30)


#1 c
library(VGAM)

n_sim <- 10000
n1 <- 20
n2 <- 30
b <- 4

# Initialize counter
type1_error <- 0

# Start simulations
for (i in 1:n_sim) {
  # Generate Laplace data under H0: mu = 0
  x <- rlaplace(n2, location = 0, scale = b)
  
  # Month 20 test
  x1 <- x[1:n1]
  t1 <- (mean(x1)) / (sd(x1) / sqrt(n1))
  crit1 <- qt(0.95, df = n1 - 1)
  
  if (t1 > crit1) {
    type1_error <- type1_error + 1  # Stopped early: false positive
  } else {
    # Month 30 test
    x2 <- x[1:n2]
    t2 <- (mean(x2)) / (sd(x2) / sqrt(n2))
    crit2 <- qt(0.95, df = n2 - 1)
    
    if (t2 > crit2) {
      type1_error <- type1_error + 1  # Final test: false positive
    }
  }
}

# Empirical Type I error rate
type1_error_rate <- type1_error / n_sim

#1 d
library(VGAM)

#Given Info
n_sim <- 10000
n1 <- 20
n2 <- 30
b <- 4

# Search bounds for alpha
lower_alpha <- 0.01
upper_alpha <- 0.1
target_error <- 0.05
tolerance <- 0.001

# Binary search function
find_alpha <- function() {
  while ((upper_alpha - lower_alpha) > tolerance) {
    alpha_mid <- (lower_alpha + upper_alpha) / 2
    
    # Get corresponding critical t-values
    t_crit_20 <- qt(1 - alpha_mid, df = n1 - 1)
    t_crit_30 <- qt(1 - alpha_mid, df = n2 - 1)
    
    # Simulate data under H0
    type1_error <- 0
    for (i in 1:n_sim) {
      x <- rlaplace(n2, location = 0, scale = b)
      
      # Test at n = 20
      x1 <- x[1:n1]
      t1 <- mean(x1) / (sd(x1) / sqrt(n1))
      
      if (t1 > t_crit_20) {
        type1_error <- type1_error + 1
      } else {
        # Test at n = 30
        x2 <- x[1:n2]
        t2 <- mean(x2) / (sd(x2) / sqrt(n2))
        
        if (t2 > t_crit_30) {
          type1_error <- type1_error + 1
        }
      }
    }
    
    type1_rate <- type1_error / n_sim
    
    # Narrow bounds
    if (type1_rate > target_error) {
      upper_alpha <- alpha_mid
    } else {
      lower_alpha <- alpha_mid
    }
  }
  return((lower_alpha + upper_alpha) / 2)
}

# Run the search
adjusted_alpha <- find_alpha()

#2 a
simulate_left_tail <- function(beta_params, n = 15, n_sim = 10000, alpha = 0.05) {
  a <- beta_params[1]
  b <- beta_params[2]
  mu <- a / (a + b)
  errors <- 0
  
  for (i in 1:n_sim) {
    x <- rbeta(n, a, b)
    x <- x - mu
    t_stat <- mean(x) / (sd(x) / sqrt(n))
    crit_val <- qt(alpha, df = n - 1)
    if (t_stat < crit_val) errors <- errors + 1
  }
  
  return(errors / n_sim)
}

left_tail_errors <- c(
  Beta_10_2 = simulate_left_tail(c(10, 2)),
  Beta_2_10 = simulate_left_tail(c(2, 10)),
  Beta_10_10 = simulate_left_tail(c(10, 10))
)

print(left_tail_errors)

#2 b 
simulate_right_tail <- function(beta_params, n = 15, n_sim = 10000, alpha = 0.05) {
  a <- beta_params[1]
  b <- beta_params[2]
  mu <- a / (a + b)
  errors <- 0
  
  for (i in 1:n_sim) {
    x <- rbeta(n, a, b)
    x <- x - mu
    t_stat <- mean(x) / (sd(x) / sqrt(n))
    crit_val <- qt(1 - alpha, df = n - 1)
    if (t_stat > crit_val) errors <- errors + 1
  }
  
  return(errors / n_sim)
}

right_tail_errors <- c(
  Beta_10_2 = simulate_right_tail(c(10, 2)),
  Beta_2_10 = simulate_right_tail(c(2, 10)),
  Beta_10_10 = simulate_right_tail(c(10, 10))
)

print(right_tail_errors)

#2 c
simulate_two_tail <- function(beta_params, n = 15, n_sim = 10000, alpha = 0.05) {
  a <- beta_params[1]
  b <- beta_params[2]
  mu <- a / (a + b)
  errors <- 0
  
  for (i in 1:n_sim) {
    x <- rbeta(n, a, b)
    x <- x - mu
    t_stat <- mean(x) / (sd(x) / sqrt(n))
    crit_val <- qt(1 - alpha/2, df = n - 1)
    if (abs(t_stat) > crit_val) errors <- errors + 1
  }
  
  return(errors / n_sim)
}

two_tail_errors <- c(
  Beta_10_2 = simulate_two_tail(c(10, 2)),
  Beta_2_10 = simulate_two_tail(c(2, 10)),
  Beta_10_10 = simulate_two_tail(c(10, 10))
)

print(two_tail_errors)

#2 d

# Create a data frame taking in data from parts a-c
type1_errors <- data.frame(
  Distribution = rep(c("Beta(10,2)", "Beta(2,10)", "Beta(10,10)"), each = 3),
  Test_Type = rep(c("Left-tailed", "Right-tailed", "Two-tailed"), times = 3),
  Error_Rate = c(0.0315, 0.0779, 0.0643,
                 0.0843, 0.0272, 0.0593,
                 0.0525, 0.0469, 0.0519)
)

library(ggplot2)

# Plot
ggplot(type1_errors, aes(x = Distribution, y = Error_Rate, fill = Test_Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  labs(title = "Effect of Skewness on Type I Error Rate",
       y = "Type I Error Rate", x = "Distribution",
       fill = "Test Type") +
  theme_minimal(base_size = 14) +
  ylim(0, 0.12) +
  geom_text(aes(label = sprintf("%.3f", Error_Rate)),
            position = position_dodge(0.9),
            vjust = -0.5, size = 4)
