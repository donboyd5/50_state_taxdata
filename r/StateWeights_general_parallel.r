
# GOAL: create a more-flexible program that relies on functions as much as possible 
# possible precursor to writing a package

# code folding ----
# alt-o, shift-alt-o
# alt-l, shift-alt-l
# alt-r

# notes ----


# libraries ----
source(here::here("include", "libraries.r"))
library(nloptr)
# remotes::install_github("tidyverse/dplyr") if needed
devtools::session_info()
(.packages()) %>% sort

# globals ----
dbox <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/"
(fns <- paste0(c("acs_10krecs_5states", "acs_100krecs_20states", "acs_200krecs_50states", "acs_400krecs_50states"), ".rds"))

# functions ----
source(here::here("include", "functions_scaling.r"))
source(here::here("include", "functions_prep_dev.r")) # soon we will replace functions_prep.r with the dev version
source(here::here("include", "functions_check_inputs.r"))
source(here::here("include", "functions_ipopt.r"))

# GET AND PREPARE SAVED DATA ----

# to create data used in this program, run the program below:
#
#    create_ACS_subset_from_extract.r
#
# to create a subset suitable for targeting

# choose which file to use 
samp1 <- readRDS(here::here("data", fns[4])) %>% 
  select(-nrecs, -pop) # note that we no longer need nrecs; pop ordinarily would not be in the data so drop here and create later
glimpse(samp1)
summary(samp1)
count(samp1, mar)
# djb note that we have nrecs and pop variables -- I want to make sure they are not needed for anything ----
# if you want to target the total number of weighted records we need a variable that is 1 for all records ----

# PREPARE DATA  ----
#.. modify the sample (don't think we need a function for this) ----
# - define income groups
# - create an indicator for each income variable as to whether it is nonzero
# _ expand categoricals into dummies as needed
# if we don't have a variable such as pop where all values are 1, we should create it as it makes it easy to get weighted record counts
samp2 <- samp1 %>%
  mutate(pid=row_number(), # pid -- an id variable for each person in the file
         incgroup=ntile(pincp, 10), # divide the data into 10 income ranges
         pop=1, # it's useful to have a variable that is 1 on every record
         # convert categoricals to dummies if we will base targets upon them
         mar1=ifelse(mar==1, 1, 0), # married
         mar5=ifelse(mar==5, 1, 0), # single
         marx15=ifelse(mar %nin% c(1, 5), 1, 0)
  )
summary(samp2)
ht(samp2)

#.. define the kinds of (weighted) targets we want and prepare the file accordingly ----
# sum:    sum of values
# nnz:    number of nonzero values
# sumneg: sum of negative values
# nneg:   number of zero values
# sumpos: sum of positive value
# npos:   number of positive values

# For the PUF the SOI data provide only the first two kinds of targets, but for the ACS we could have any of them.

# TRY TO AVOID DEPENDENT CONSTRAINTS - redundancy - as they can make the problem very hard to solve.
# For example, suppose there are 3 kinds of income (wages, interest, retirement) plus a total (sum of the 3)
#   -- don't create targets for each of the 3 kinds plus a target for the total -- leave one item out
# Another, less obvious example: don't target the total number of returns plus the number for each marital status - leave one out.

nnz_vars <- c("pop", "mar1", "mar5", "pincp", "wagp") # note that I leave the 3rd marital status out -- marx15
sum_vars <- c("pincp", "wagp", "intp", "pap", "retp", "ssip", "ssp") # DO NOT include total plus all sums - leave one out (otherincp)
sumneg_vars <- "otherincp"

# define a vector of variable names for "possible" targets (a superset) -- we may not target all
possible_target_vars <- make_target_names(
  list(nnz_vars, sum_vars, sumneg_vars),
  c("nnz", "sum", "sumneg"))
possible_target_vars

# prepare data by creating variables with those names:
#   nnz, nneg, and npos will have 1 for rows where the respective variable is nz, neg, or pos, respectively, and 0 otherwise
#   sum will have its respective variable's value
#   sumneg and sumpos will have the variable's value if negative or positive, respectively, and 0 otherwise
samp <- prep_data(samp2, possible_target_vars)
glimpse(samp)

ht(samp)
summary(samp) # make sure there are no NAs -- there shouldn't be
count(samp, stabbr)
count(samp, incgroup)

# SUMMARIZE the data ----
#.. count records ----
# the following requires the dev version of dplyr
#   remotes::install_github("tidyverse/dplyr")
# may not be worth installing just for this
count_recs(samp) %>% kable(format="rst", format.args = list(big.mark=",")) # this is just for information, not needed to move forward

#.. get weighted summary values needed for developing constraints ----
summary_vals <- get_summary_vals(samp, .weight=pwgtp, .sum_vars=possible_target_vars, stabbr, incgroup)
summary_vals

# look at the summary values
tlist <- get_tabs(summary_vals, .popvar=pop_nnz)
names(tlist)
tlist$wsum
tlist$wcount
tlist$wmean
tlist$pctnz

# Create a data frame with all targets for all states and income groups ----
# for the PUF, we will create this using information from Historical Table 2
# for the ACS, we construct the targets from the ACS data itself
all_target_values <- summary_vals

# now we have:
#  1) properly-prepared data that includes all income groups, 
#  2) a set of targets for all income groups and all states

# the next step is to prepare data and targets for a single income group
# so that we can solve for optimal values
# we will wrap that in functions so that we can loop through income groups


# wrap everything we need for a single income group into a function that returns a list ----

# SINGLE INCOME GROUP ----
#.. define target incgroup, target variable names, and target values for each state in the income group ----
target_incgroup <- 2 # define target income group

possible_target_vars
(target_vars <- possible_target_vars[c(1, 3, 6, 7)])
target_vars <- possible_target_vars
target_vars <- setdiff(possible_target_vars, c("pap_sum", "ssip_sum", "intp_sum", "otherincp_sumneg"))

# define target values and states, for this income group
targets_wide <- all_target_values %>%
  filter(incgroup==target_incgroup) %>%
  select(stabbr, incgroup, nrecs, all_of(target_vars)) # a small list of variables to target; we have nrecs because we created it in summary_vals
targets_wide # these are the targets we want to hit
target_states <- targets_wide$stabbr

#.. get data for this income group: will include pid, incgroup, weight_total (national weight), and the target variables ----
incgroup_data <- prep_incgroup(.data=samp, .incgroup=target_incgroup, .weight=pwgtp, .target_vars=target_vars)
# djb now go to wfs ----






targets_df <- get_targets(targets_wide, target_vars, .pid=incgroup_data$pid, .weight_total=incgroup_data$weight_total)
tail(targets_df)


# CONSTRUCT initial weights ----
#.. prepare the data we will reweight to hit (or come close to) the targets ----
# let's stack the data so that each person appears 1 time for each state that is targeted
# create a stub with a record for every person and every targeted state

# define initial weight for each person-state combination. This is important as the
# optimization will try to keep the final weights near these initial weights. The
# more we can improve these initial weights, the better the results are likely to be
iweights <- get_initial_weights(targets_wide, incgroup_data, .popvar=pop_nnz)
ht(iweights)

# DEFINE constraint coefficients (cc) ----
cc_sparse <- get_constraint_coefficients(incgroup_data, target_vars, iweights, targets_df) # can take a while on big problems

# DEFINE inputs for ipopt ----
inputs <- get_inputs(.targets_df=targets_df, .iweights=iweights, .cc_sparse=cc_sparse, .targtol=.02, .xub=50, .conscaling=TRUE, scale_goal=1)


#.. examine constraints, revise if needed ----
names(inputs)
inputs$n_constraints
inputs$n_targets
ht(inputs$cc_sparse)
check <- check_constraints(inputs)
summary(check)
probs <- c(0, .5, .95, .96, .97, .98, .99, .995, .999, .9995, 1)
quantile(abs(check$pdiff), na.rm=TRUE, probs)
check %>% arrange(desc(abs(pdiff)))

#.. OPTIONAL revise constraints ----

#.. VERIFY (if desired) that problem is set up properly ----
# eval_f_xm1sq(inputs$x0, inputs)
# eval_grad_f_xm1sq(inputs$x0, inputs)
# eval_g(inputs$x0, inputs)
# eval_jac_g(inputs$x0, inputs)


# OPTIMIZE ----
# CAUTION: If you use ma77, be sure to delete temporary files it created in the project
# folder before rerunning or RStudio may bomb.
# For LARGE problems use linear_solver=ma77, obj_scaling=1, and mehrotra_algorithm=yes
opts <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "max_iter"= 50,
             "linear_solver" = "ma57", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             # "ma57_automatic_scaling" = "yes", # if using ma57
             # "ma57_pre_alloc" = 3, # 1.05 is default; even changed, cannot allocate enough memory, however
             # "ma77_order" = "amd",  # metis; amd -- not clear which is faster
             "mehrotra_algorithm" = "yes",
             "obj_scaling_factor" = 1, # 1e-3, # default 1; 1e-1 pretty fast to feasible but not to optimal
             # "nlp_scaling_method" = "equilibration-based", # NO - use default gradient_based
             "nlp_scaling_max_gradient" = 100, # default is 100 - seems good 
             # "jac_c_constant" = "yes", # does not improve on moderate problems
             # "jac_d_constant" = "yes", # does not improve on  moderate problems
             # "hessian_constant" = "yes", # KEEP default NO - if yes Ipopt asks for Hessian of Lagrangian function only once and reuses; default "no"
             # "hessian_approximation" = "limited-memory", # KEEP default of exact
             # "derivative_test" = "first-order",
             # "derivative_test_print_all" = "yes",
             # "derivative_test_perturbation" = 1e-6,
             "output_file" = here::here("out", "test_os100cs1.out"))

#setwd(here::here("temp3"))
#getwd()

inputs$objscale <- 100 # this is the right direction -- eliminates grad_f deriv errors

# inputs$objscale <- 1
# inputs$objscale <- .01

resos100cs1 <- ipoptr(x0 = inputs$x0,
                 lb = inputs$xlb,
                 ub = inputs$xub,
                 eval_f = eval_f_xm1sq, # arguments: x, inputs; eval_f_xtop eval_f_xm1sq
                 eval_grad_f = eval_grad_f_xm1sq, # eval_grad_f_xm1sq, # eval_grad_f_xtop eval_grad_f_xm1sq
                 eval_g = eval_g, # constraints LHS - a vector of values
                 eval_jac_g = eval_jac_g,
                 eval_jac_g_structure = inputs$eval_jac_g_structure,
                 eval_h = eval_h_xm1sq, # the hessian is essential for this problem eval_h_xtop eval_h_xm1sq
                 eval_h_structure = inputs$eval_h_structure,
                 constraint_lb = inputs$clb,
                 constraint_ub = inputs$cub,
                 opts = opts,
                 inputs = inputs)

res1$solution
res100$solution
resp01$solution
resos100cs1$solution

a <- res1$solution
b <- resos100cs1$solution

sum((a - b)^2)
cor(a, b)



d2 <- result$solution
d2
d2g <- eval_grad_f_xm1sq(d2, inputs)
d2[1]
d2g[1]

eval_f_xm1sq(d2[1], inputs)
sum(inputs$iweight * (d2-1)^2)


library(numDeriv)
grad2 <- function(x, inputs){
  grad(func=eval_f_xm1sq, x, inputs=inputs)
}

x <- runif(length(inputs$x0))
a <- eval_grad_f_xm1sq(x, inputs)
b <- grad2(x, inputs)
a; b
cbind(a[1:10], b[1:10])


eval_f_xm1sq2(d2, inputs)
diff <- 1e-8
d2a <- d2
d2a[1] <- d2[1] + diff
f1 <- eval_f_xm1sq2(d2, inputs)
f2 <- eval_f_xm1sq2(d2a, inputs)
(f2 - f1) / diff

eval_grad_f_xm1sq2(d2, inputs)

a <- sum(inputs$iweight * (d2-1)^2)
b <- sum(inputs$iweight * (d2-1 + .000001)^2)
b - a

eval_f_xm1sq2 <- function(x, inputs) {
  # objective function - evaluates to a single number
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  # here are the objective function, the 1st deriv, and the 2nd deriv
  # http://www.derivative-calculator.net/
  # w*(x-1)^2                  objective function
  # 2*w*(x-1)                  first deriv
  # 2*w                        second deriv
  
  # make it easier to read:
  w <- inputs$iweight / inputs$objscale
  print(w[1])
  print(x[1])
  
  obj <- sum(w * (x-1)^2)
  
  return(obj)
}

eval_grad_f_xm1sq2 <- function(x, inputs){
  # gradient of objective function - a vector length x 
  # giving the partial derivatives of obj wrt each x[i]
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  
  # http://www.derivative-calculator.net/
  # w*(x-1)^2                  objective function
  # 2*w*(x-1)                  first deriv
  # 2*w                        second deriv
  
  # make it easier to read:
  w <- inputs$iweight / inputs$objscale
  print(w[1])
  print(x[1])
  
  gradf <- 2 * w * (x-1)
  
  return(gradf)
}

# EVALUATE RESULTS ----
# saveRDS(result, here::here("results", "result.rds"))
# result <- readRDS(here::here("results", "result.rds"))

str(result)
quantile(result$solution, probs=c(0, .001, .01, .05, .1, .25, .5, .75, .9, .95, .99, .999, 1))
result$solution %>% sort() %>% ht(25)
result$solution %>% sort() %>% head(10000)
start <- 26e3
result$solution %>% sort() %>% .[start:(start + 1e3)]
result$constraints

# now show the actual constraints and results
# slicen <- 15
# constraints_start %>%
#   mutate(confinal=eval_g(result$solution, inputs)) %>%
#   slice(1:slicen, (n() - slicen + 1):n()) %>%
#   kable(digits=0, format="rst", format.args=list(big.mark=","))
str(inputs)
stub <- iweights %>%
  mutate(x=result$solution,
         weight=iweight_state * x)
stub

probs <- c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99)
quantile(stub$weight, probs) %>% round(3)
sum(stub$weight)

stub %>%
  group_by(pid) %>%
  summarise(weight_total=first(weight_total),
            iweight=sum(iweight_state),
            weight=sum(weight)) %>%
  mutate(diff=weight - weight_total) %>%
  arrange(-abs(diff))

comp <- base_data %>%
  pivot_longer()
  mutate(x=result$solution,
         weight=iweight * x) %>%
  group_by(stabbr) %>%
  summarise_at(vars(all_of(target_vars)), list(initial = ~sum(. * iweight), solution = ~sum(. * weight))) %>%
  pivot_longer(-stabbr) %>%
  separate(name, into=c("name1", "name2", "type"), fill="left") %>%
  unite(name, name1, name2, na.rm=TRUE) %>%
  bind_rows(targets %>% select(stabbr, all_of(target_vars)) %>% pivot_longer(-stabbr) %>% mutate(type="target")) 

slicen <- 15
comp %>%
  pivot_wider(names_from = type) %>%
  select(stabbr, name, target, initial, solution) %>%
  mutate(idiff=initial - target,
         diff=solution - target,
         pdiff=diff / target * 100) %>%
  arrange(-abs(pdiff)) %>%
  slice(1:slicen, (n() - slicen + 1):n()) %>%
  kable(digits=c(rep(0, 7), 2), format="rst", format.args=list(big.mark=","))

# what states did the record weights come from?
stub %>%
  mutate(x=result$solution,
         weight=iweight * x) %>%
  group_by(stabbr_true, stabbr) %>%
  summarise_at(vars(pwgtp, weight), sum) %>%
  pivot_wider(names_from = stabbr, values_from = weight) %>%
  select(1:15) %>%
  kable(digits=1, format="rst", format.args=list(big.mark=","))

stub %>%
  mutate(x=result$solution,
         weight=iweight * x) %>%
  group_by(stabbr_true, stabbr) %>%
  summarise_at(vars(pwgtp, weight), sum) %>%
  group_by(stabbr) %>%
  mutate(pct_true=pwgtp / sum(pwgtp) * 100) %>%
  group_by(stabbr_true) %>%
  mutate(pct_est=weight / sum(weight) * 100) %>%
  select(-weight) %>%
  pivot_wider(names_from = stabbr, values_from = pct_est) %>%
  select(c(1:12, MA, MI, NY, OH, PA, TX)) %>%
  kable(digits=c(0, 0, rep(2, 50)), format="rst", format.args=list(big.mark=","))


# how did we do against the adding-up goal?
ratios <- stub %>%
  mutate(x=result$solution,
         weight=iweight * x) %>%
  group_by(p) %>%
  summarise(pwgtp=first(pwgtp), weight=sum(weight)) %>%
  mutate(ratio=weight / pwgtp)
quantile(ratios$ratio, probs = c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1))


# ma77 notes
#  !!! DELETE ANY ma77* FILES IN WORKING DIRECTORY BEFORE RUNNING ma77 !!!!----
# when we get to here we should have puf4 and targets

# maybe adjust threads for faster ma86?? I haven't found anything that works
# print(Sys.unsetenv("OMP_NUM_THREADS")) # unset the threads environment variable

# Sys.getenv(c("R_HOME", "OMP_NUM_THREADS"))
# print(Sys.setenv(OMP_NUM_THREADS = 6))
# Sys.getenv(c("R_HOME", "OMP_NUM_THREADS"))


