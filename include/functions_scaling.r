
# Scaling notes from documentation ----
# https://projects.coin-or.org/Ipopt/wiki/HintsAndTricks

# When you formulate your optimizaton problem, you should try to make it "well-scaled" to make it easier for
# Ipopt (or any nonlinear optimization package) to solve it. For this, you should try to make the "typical" 
# values of the non-zero first partial derivatives of the objective and constraint functions to be on the order
# of, say, 0.01 to 100. For example, if you multiply a problem function by a number K, then the first partial
# derivatives for this function are also multiplied by K. On the other hand, if you replace a variable xi by a 
# number K, then the partial derivatives with respect to this variable are divided by K.

# By default, Ipopt performs some very simple scaling of the problem functions, by looking at the gradients
# of each function evaluated at the user-provided starting point. If any entry for a given function is larger
# than 100, then this function is scaled down to make the largest entry in the gradient 100 (see the Ipopt
# options nlp_scaling_method and nlp_scaling_max_gradient). Of course, if some of the gradient elements are
# huge and some are very small, the variables corresponding to the small entries are almost ignored. If you set 
# the print_level to at least 5, then you can see by how much the functions are scaled. For sufficiently 
# large print_level you can also see the individual scaling factors for each constraint, and also the values 
# of the first derivatives, if you want to find out which derivatives are very large or very small.

# If you have trouble finding a solution with Ipopt, it might sometimes help to increase to accuracy of the 
# computation of the search directions, which is done by solving a linear system. There are a number of ways
# to do this: First, you can tell Ipopt to scale the linear system before it is sent to the linear solver; 
# this assumes that you compiled your code with the Harwell routine MC19. By default this option is only 
# activated if it seems necessary (see linear_system_scaling and linear_scaling_on_demand). Also, you can 
# increase the pivot tolerance of the linear solver (if that is supported by the linear solver), for example, 
# by increasing the value of ma27_pivtol if you are using the linear solver MA27.

# Other options you might want to play with if you have trouble solving a problem are: mu_strategy, 
# obj_scaling_factor, mu_init, bound_relax_factor...

# My scaling comments ----
# ...try to make the "typical" values of the non-zero first partial derivatives of the objective and constraint
# functions to be on the order of, say, 0.01 to 100.

scale_inputs <- function(inputs, scale_goal=1){
  
  ccsum <- inputs$cc_sparse %>%
    group_by(i, cname) %>%
    summarise(nzmin=min(nzcc), nzmdn=median(nzcc), nzmax=max(nzcc),
              .groups="drop")
  
  scale_factors <- tibble(cname=inputs$constraint_names, cvalue=inputs$constraints) %>%
    left_join(ccsum, by = "cname") %>%
    mutate(scale=abs(nzmdn / scale_goal)) %>%
    mutate_at(vars(cvalue, starts_with("nz")), list(~ . / scale))
  
  inputs_scaled <- inputs
  
  inputs_scaled$cc_sparse <- inputs$cc_sparse %>%
    left_join(scale_factors %>% select(cname, scale), by = "cname") %>%
    mutate(nzcc = nzcc / scale)
  
  inputs_scaled$constraints <- scale_factors$cvalue
  
  inputs_scaled$clb <- inputs$clb / scale_factors$scale
  inputs_scaled$cub <- inputs$cub / scale_factors$scale
  
  return(inputs_scaled)
}

