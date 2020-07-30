
# code folding ----
# alt-o, shift-alt-o
# alt-l, shift-alt-l
# alt-r

# notes ----

# Use subsets of the ACS file to test different approaches to weighting state income tax files.

# important variables ----
# serialno -- unique identifier
# sporder -- person number within household (in case we ever want to link to ACS household data)
# st -- state FIPS code, we can drop it
# stabbr -- state abbreviation
# pwgtp -- person weight
# adjinc -- multiplicative factor to convert all 5 ACS years to a common year (already applied)
# mar -- marital status (1 married, 2 widowed, 3 divorced, 4 separated, 5 never married or < 15 years old)
# sex -- 1 male, 2 female
# agep -- 0-99 where 0 is < 1, and 99 includes older ages as well
# pincp -- personal income
# wagp -- wage income
# intp -- interest income 
# pap -- public assistance income
# retp -- retirement income (does not include Social Security)
# ssip -- Supplemental Security Income
# ssp Social Security income


# libraries ----
source(here::here("include", "libraries.r"))

# globals ----
dbox <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/"
(fns <- paste0(c("acs_10krecs_5states", "acs_100krecs_20states", "acs_200krecs_50states", "acs_400krecs_50states"), ".rds"))

# functions ----
source(here::here("include", "functions_ipopt.r"))
source(here::here("include", "functions_scaling.r"))

# GET SAVED DATA ----

# to create data used in this program, as needed, run:
#
#    create_ACS_subset_from_extract.r
#
# to create a subset suitable for targeting

samp <- readRDS(here::here("data", fns[1]))  # choose which file to use 

# modify the sample:
# - define income groups
# - create an indicator for each income variable as to whether it is nonzero
samp <- samp %>%
  mutate(pid=row_number(), # pid is person id
         incgroup=ntile(pincp, 10))  %>% # divide the data into 10 income ranges
  mutate_at(vars(pincp, wagp, intp, pap, retp, ssip, ssp, otherincp), list(n = ~ as.numeric(. !=0)))
glimpse(samp)
ht(samp)
summary(samp) # make sure there are no NAs -- there shouldn't be
count(samp, stabbr)
count(samp, incgroup)
count(samp, incgroup, stabbr)
count(samp, sex)
quantile(samp$agep)

# CALCULATE SUMMARIES for each income group and state ----
# we will use these summaries for developing constraints
summary_vals <- samp %>%
  group_by(stabbr, incgroup) %>%
  summarise_at(vars(nrecs, pop, pincp, intp, pap, retp, ssip, ssp, wagp, otherincp, ends_with("_n")),
               list(~ sum(. * pwgtp))) %>%
  ungroup

summary_vals

# sums of values:
summary_vals %>%
  select(-ends_with(("_n"))) %>%
  kable(digits=0, format.args = list(big.mark=","), format="rst")

# sums of weighted counts:
summary_vals %>%
  select(stabbr, incgroup, nrecs, pop, ends_with(("_n"))) %>%
  kable(digits=0, format.args = list(big.mark=","), format="rst")

# weighted means
summary_vals %>%
  select(-ends_with(("_n"))) %>%
  mutate_at(vars(pincp:otherincp), list(~ . / pop)) %>%
  kable(digits=0, format.args = list(big.mark=","), format="rst")

# percentages of population that have nonzero values
summary_vals %>%
  select(stabbr, incgroup, nrecs, pop, ends_with(("_n"))) %>%
  mutate_at(vars(pincp_n:otherincp_n), list(~ . / pop * 100)) %>%
  kable(digits=1, format.args = list(big.mark=","), format="rst")

# MULTIPLE-STATES with ADDING-UP constraints ----

# define a state and income group that we will target; also, pick only a few variables to target
target_incgroup <- 2
target_vars <- c("pop", "pincp", "wagp", "intp", "ssip_n") # keep it small for our test
# CAUTION: in code below we'll keep variables in this order so that we can index into matrices by column
# we can put the variables above in any order we want, but once we choose that order we want to be sure to
# respect it below; this will be clearer later

data_true <- samp %>%
  filter(incgroup==target_incgroup) %>%
  select(pid, stabbr, incgroup, pwgtp, nrecs, all_of(target_vars))
data_true  # all states, 1 income range

# define the targets we want to hit for a given state and income, using data for all states in that income group
targets <- summary_vals %>%
  filter(incgroup==target_incgroup) %>%
  select(stabbr, incgroup, nrecs, all_of(target_vars)) # a small list of variables to target
targets # these are the targets we want to hit
targets_long <- targets %>%
  pivot_longer(cols=all_of(target_vars)) %>%
  mutate(cname=paste0(name, "_", stabbr)) %>%
  arrange(cname)
targets_long

# now create an add-on data frame that has the adding-up constraint for each person
targets_adding_up <- data_true %>%
  mutate(cname=paste0("p", str_pad(pid, width=8, pad="0")), value=pwgtp) %>%
  select(cname, pid, cname, value) %>%
  arrange(cname)
targets_adding_up %>% ht
summary(targets_adding_up)

# we also need these as a vector -- make sure they are in the same order
(constraint_names <- c(targets_long$cname, targets_adding_up$cname))
(constraints <- c(targets_long$value, targets_adding_up$value))
constraints_df <- tibble(cname=constraint_names, cvalue=constraints)


#.. prepare the data we will reweight to hit (or come close to) the targets ----
# let's stack the data so that each person appears 1 time for each state that is targeted
# create a stub with a record for every person and every targeted state
base_data <- samp %>%
  filter(incgroup==target_incgroup) %>%
  select(pid, stabbr_true=stabbr, incgroup, pwgtp, nrecs, all_of(target_vars))
base_data %>% ht

stub <- expand_grid(pid=unique(base_data$pid), stabbr=unique(base_data$stabbr_true)) %>%
  arrange(pid, stabbr) %>%
  mutate(j=row_number()) %>%
  left_join(base_data) %>%
  left_join(targets %>% select(stabbr, state_pop=pop)) %>%
  group_by(stabbr) %>%
  mutate(iweight=pwgtp / sum(pwgtp) * state_pop) %>%
  ungroup

stub %>%
  group_by(stabbr) %>%
  summarise(n=n(), pwgtp_sum=sum(pwgtp), iweight_sum=sum(iweight), pincp=sum(pincp * iweight))
targets

#.. create constraint coefficients ----
# before we create a dense matrix that had 1 column for every constraint and then we made it long to get
# our sparse matrix
# Now that we have 1,000 additional constraints let's find another way so that we don't have to have 1,000 columns
# In a real-world problem we might have 20k or 30k adding-up constraints, or more

# let's do it in 2 steps: a dense matrix for the regular constraints, and go directly to sparse for the adding up
# constraints
cc_dense1 <- stub %>%
  mutate_at(vars(all_of(target_vars)), list(~ . * iweight))

cc_sparse1 <- cc_dense1 %>%
  pivot_longer(cols = all_of(target_vars), values_to = "nzcc") %>% # nzcc for nonzero constraint coefficient
  filter(nzcc!=0) %>%
  mutate(cname=paste0(name, "_", stabbr),
         i=match(cname, constraints_df$cname)) %>%
  select(pid, j, i, cname, nzcc, iweight, name, stabbr, stabbr_true, pwgtp, incgroup) %>%
  arrange(i, j) # this ordering is crucial for the Jacobian
ht(cc_sparse1)
count(cc_sparse1, i)
# We now have the first part of the sparse matrix, for the "true" constraints (i.e., other than adding-up constraints)

# Let's build the second part, for the adding up constraints. It will have 1 row per person per state,
#   the same as the number of rows in cc_dense1 above
# Each row will have the following variables:
#   pid  the identifier of the person in the original data
#   i  the row number for the constraint in the imaginary dense matrix, which is its position in the constraints vector
#   j  the column number for the x variable in the imaginary dense matrix, which is its row number in the dense1 matrix
#   cname the constraint name
#   nzcc  the nonzero constraint coefficient, which is the amount by which the constraint value will change for a
#      unit change in the x value, which will simply be the initial weight value
# we can keep some other variables if we want
cc_sparse2 <- cc_dense1 %>%
  select(pid, j, nzcc=iweight, stabbr_true, stabbr) %>%
  # let's get cname from the targets_adding_up data frame, although we could simply create it
  left_join(targets_adding_up %>% select(pid, cname), by = "pid") %>%
  # and get i, the index for the constraint column
  mutate(i=match(cname, constraints_df$cname))
glimpse(cc_sparse2)

# make sure that the constraint numbering of cc_sparse1 and cc_sparse2 form a continuous series of integers
# as they will be used as row indexes into the imaginary dense matrix
icheck <- c(cc_sparse1$i, cc_sparse2$i) %>% unique %>% sort
itrue <- 1:nrow(constraints_df)
all.equal(icheck, itrue)

cc_sparse <- bind_rows(cc_sparse1, cc_sparse2) %>%
  arrange(i, j) # this ordering is crucial for the Jacobian
ht(cc_sparse)

#.. set up the optimization ----
inputs_unscaled <- list()
inputs_unscaled$p <- 2
inputs_unscaled$iweight <- cc_dense1$iweight # the initial weight
inputs_unscaled$cc_sparse <- cc_sparse
inputs_unscaled$constraints_df <- constraints_df
inputs_unscaled$n_variables <- length(inputs_unscaled$iweight)
inputs_unscaled$n_constraints <- nrow(inputs_unscaled$constraints_df)
inputs_unscaled$objscale <- 1
# define constraint lower and upper bounds -- for now they will be the same
tol <- 0 # tolerance
tol <- rep(0, inputs_unscaled$n_constraints)
tol[1:nrow(targets_long)] <- 0.01
tol

bounds <- tibble(cvalue=inputs_unscaled$constraints_df$cvalue) %>%
  mutate(clb=ifelse(cvalue==0, -Inf, cvalue - abs(cvalue) * tol),
         cub=ifelse(cvalue==0, +Inf, cvalue + abs(cvalue) * tol))

inputs_unscaled$clb <- bounds$clb
inputs_unscaled$cub <- bounds$cub


# do one of the following
inputs <- inputs_unscaled

str(inputs)

xlb <- rep(0, inputs$n_variables)
xub <- rep(20, inputs$n_variables)
x0 <- rep(1, inputs$n_variables)
eval_jac_g_structure <- define_jac_g_structure_sparse(inputs$cc_sparse, ivar="i", jvar="j")
eval_jac_g_structure[5]
length(eval_jac_g_structure)

eval_h_structure <- lapply(1:inputs$n_variables, function(x) x) # diagonal elements of our Hessian

eval_f_xtop(x0, inputs)# ; eval_f_absapprox(xlb, inputs); eval_f_absapprox(xub, inputs)
eval_grad_f_xtop(x0, inputs)

constraints_start <- inputs$constraints_df %>%
  mutate(clb=inputs$clb, cub=inputs$cub, conx0=eval_g(x0, inputs),
         pdiff=conx0 / cub * 100 - 100)

constraints_start %>%
  arrange(-abs(pdiff))

end1 <- min(10, trunc((inputs$n_constraints)/2))
start2 <- max(inputs$n_constraints - 10, 1)
rownums <- c(1:end1, start2:inputs$n_constraints)
constraints_start %>%
  filter(row_number() %in% rownums) %>%
  kable(digits=0, format="rst", format.args=list(big.mark=","))

eval_jac_g(x0, inputs)[1:10]
eval_h_xtop(x0, obj_factor=1, hessian_lambda=rep(1, inputs$n_constraints), inputs)[1:10] # length is n_variables
# length(unlist(eval_h_structure)) # length should be n_variables


# CAUTION: If you use ma77, be sure to delete temporary files it created in the project
# folder before rerunning or RStudio may bomb.
opts <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "linear_solver" = "ma57", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             "max_iter"= 100,
             "obj_scaling_factor" = 1, # 1e7, # default 1; 1e-1 pretty fast to feasible but not to optimal
             # "mehrotra_algorithm" = "yes",
             #"nlp_scaling_max_gradient" = 1e-3, # 1e-3, # default 100
             # "derivative_test" = "first-order",
             #"derivative_test_print_all" = "yes",
             "output_file" = here::here("out", "multiple_states_add.out"))

# Approx timings for a problem close in size and structure to what we'll be doing. This problem has:
#   1 income range (#2)
#   50 states
#   20k people
#   5 targets per state (pop, pincp, wagp, intp, ssip_n) each with a 1% tolerance around the target,
#   plus and adding-up requirement for each person
# Thus, it has:
#   1 million variables to solve for (20k people x 50 states),
#   250 inequality constraints (5 targets x 50 states), as the 1% tolerance creates lower and upper bounds
#   20k equality constraints, one per person, to ensure that the 
#   2.6 million nonzero constraint coefficients out of of 20.25 billion possible (1m vars x 20.25k constraints)
#   Hence the problem is very sparse (0.01% of imaginary matrix has nonzero values)
# Timings (AMD FX-8350 w/8 cores and 32gb RAM) and notes:
#   default ma27 1,880 secs (31.3 mins) appears to use about 8 gb of memory on top of other system needs
#   ma57  200 secs (3.3 mins) appears to use about 12gb of memory on top of other system needs, maybe less
#   ma77  460 secs (7.7 mins) out-of-memory solver, lowest memory usage; usage will not scale up greatly a problem size does; requires care
#   ma86  1,338 secs (22.3 mins), appears to use about 8gb of memory on top of other system needs
#   ma97 517 secs (8.6 mins)
# conclusion:  use ma57 if problem fits in memory, and ma77 if it doesn't
# even if our real-world problems take 2x as long, 10 income ranges in parallel might not be possible in 10-20 mins

# note that:
# - feasible solution found at iteration 6 out of 38
# - obj function improved 9.6% from there to optimal solution at iter 38
# - 70% of that improvement came in next 10 iterations -- that is, 6.8% improvement by iteration 16
# in other words, if we need more speed during test runs, we might try to stop it shortly after feasibility
# (we'd have to figure out an approximate way to do this)


result <- ipoptr(x0 = x0,
                 lb = xlb,
                 ub = xub,
                 eval_f = eval_f_xtop, # arguments: x, inputs
                 eval_grad_f = eval_grad_f_xtop,
                 eval_g = eval_g, # constraints LHS - a vector of values
                 eval_jac_g = eval_jac_g,
                 eval_jac_g_structure = eval_jac_g_structure,
                 eval_h = eval_h_xtop, # the hessian is essential for this problem
                 eval_h_structure = eval_h_structure,
                 constraint_lb = inputs$clb,
                 constraint_ub = inputs$cub,
                 opts = opts,
                 inputs = inputs)

#.. examine results ----
# saveRDS(result, here::here("results", "bigresult.rds"))
# result <- readRDS(here::here("results", "bigresult.rds"))

str(result)
quantile(result$solution, probs=c(0, .001, .01, .05, .1, .25, .5, .75, .9, .95, .99, .999, 1))
result$solution %>% sort() %>% tail(50)
result$constraints

# show the scaled constraints and results
slicen <- 15
constraints_start %>%
  mutate(confinal=eval_g(result$solution, inputs)) %>%
  slice(1:slicen, (n() - slicen + 1):n()) %>%
  kable(digits=0, format="rst", format.args=list(big.mark=","))

# now show the actual constraints and results
comp <- stub %>%
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
  select(1:min(ncol(.), 15)) %>%
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
  # select(c(1:12, MA, MI, NY, OH, PA, TX)) %>%
  kable(digits=c(0, 0, rep(2, 50)), format="rst", format.args=list(big.mark=","))


# how did we do against the adding-up goal?
ratios <- stub %>%
  mutate(x=result$solution,
         weight=iweight * x) %>%
  group_by(pid) %>%
  summarise(pwgtp=first(pwgtp), weight=sum(weight)) %>%
  mutate(ratio=weight / pwgtp)
quantile(ratios$ratio, probs = c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1))


# ma77 notes
#  !!! DELETE ANY ma77* FILES IN WORKING DIRECTORY BEFORE RUNNING ma77 !!!!----
# when we get to here we should have puf4 and targets

# print(Sys.unsetenv("OMP_NUM_THREADS")) # unset the threads environment variable

# Sys.getenv(c("R_HOME", "OMP_NUM_THREADS"))
# print(Sys.setenv(OMP_NUM_THREADS = 6))
# Sys.getenv(c("R_HOME", "OMP_NUM_THREADS"))