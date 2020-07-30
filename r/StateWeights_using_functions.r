
# GOAL: create a more-flexible program that relies on functions as much as possible 
# possible precursor to writing a package

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
# remotes::install_github("tidyverse/dplyr") if needed
devtools::session_info()
(.packages()) %>% sort

# globals ----
dbox <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/"
(fns <- paste0(c("acs_10krecs_5states", "acs_100krecs_20states", "acs_200krecs_50states", "acs_400krecs_50states"), ".rds"))

# functions ----
source(here::here("include", "functions_scaling.r"))
source(here::here("include", "functions_prep.r"))
source(here::here("include", "functions_check_inputs.r"))
source(here::here("include", "functions_ipopt.r"))

# GET AND PREPARE SAVED DATA ----

# to create data used in this program, run the program below:
#
#    create_ACS_subset_from_extract.r
#
# to create a subset suitable for targeting

samp1 <- readRDS(here::here("data", fns[4]))  # choose which file to use 
glimpse(samp1)
summary(samp1)
count(samp1, mar)

# PREPARE DATA  ----
#.. modify the sample (don't think we need a function for this) ----
# - define income groups
# - create an indicator for each income variable as to whether it is nonzero
# _ expand categoricals into dummies as needed
samp2 <- samp1 %>%
  mutate(pid=row_number(), # pid -- an id variable for each person in the file
         incgroup=ntile(pincp, 10), # divide the data into 10 income ranges
         # convert categoricals to dummies
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
# for the PUF we only have the first two kinds of targets, but for the ACS we could have any of them
# TRY TO AVOID DEPENDENT CONSTRAINTS - redundancy
nnz_vars <- c("pop", "mar1", "mar5", "pincp", "wagp")
sum_vars <- c("pincp", "wagp", "intp", "pap", "retp", "ssip", "ssp") # DO NOT include total plus all sums - leave one out (otherincp)
sumneg_vars <- "otherincp"

# define variable names for possible targets -- we may not target all
possible_targets <- make_target_names(
  list(nnz_vars, sum_vars, sumneg_vars),
  c("nnz", "sum", "sumneg"))

possible_targets

# prepare data by creating variables with those names:
#   nnz, nneg, and npos will have 1 for rows where the respective variable is nz, neg, or pos, respectively, and 0 otherwise
#   sum will have its respective variable's value
#   sumneg and sumpos will have the variable's value if negative or positive, respectively, and 0 otherwise
samp <- prep_data(samp2, possible_targets)
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
count_recs(samp) # this is just for information, not needed to move forward

#.. get weighted summary values needed for developing constraints ----
summary_vals <- get_summary_vals(samp, .weight=pwgtp, .sum_vars=possible_targets, stabbr, incgroup)
summary_vals

# look at the summary values
tlist <- get_tabs(summary_vals)
names(tlist)
tlist$wsum
tlist$wcount
tlist$wmean
tlist$pctnz


# DEFINE targets and constraints ----
# define a state and income group that we will target; also, pick only a few variables to target
target_incgroup <- 2

# create a set of target variables drawn from all_targets
possible_targets
# target_vars <- c("pop_nnz", "mar1_nnz", "wagp_nnz", 
#                  "pincp_sum", "wagp_sum", "retp_sum", "ssip_sum",
#                  "otherincp_sumneg") # keep it small for our test
target_vars <- possible_targets
target_vars
setdiff(target_vars, possible_targets) # make sure all desired targets are in the possible set

# CAUTION: in code below we'll keep variables in target_vars order so that we can index into matrices by column
# we can put the variables above in any order we want, but once we choose that order we want to be sure to
# respect it below; this will be clearer later


# make a dataset that has only what we need
base_data <- samp %>%
  filter(incgroup==target_incgroup) %>%
  select(pid, stabbr_true=stabbr, incgroup, weight_total=pwgtp, all_of(target_vars))
base_data %>% ht

# define the targets we want to hit for a given state and income, using data for all states in that income group
targets_wide <- summary_vals %>%
  filter(incgroup==target_incgroup) %>%
  select(stabbr, incgroup, nrecs, all_of(target_vars)) # a small list of variables to target
targets_wide # these are the targets we want to hit
target_states <- targets_wide$stabbr
anyDuplicated(target_states) # better not be any duplicates


targlist <- get_targets(targets_wide, target_vars, .pid=base_data$pid, .weight_total=base_data$weight_total)
names(targlist)

# targlist$targets_long %>% ht
# targlist$targets_adding_up %>% ht
# summary(targlist$targets_adding_up)
# summary(targets_adding_up)
targlist$constraint_names %>% ht
targlist$constraints %>% ht

# CONSTRUCT initial weights ----
#.. prepare the data we will reweight to hit (or come close to) the targets ----
# let's stack the data so that each person appears 1 time for each state that is targeted
# create a stub with a record for every person and every targeted state

# define initial weight for each person-state combination. This is important as the
# optimization will try to keep the final weights near these initial weights. The
# more we can improve these initial weights, the better the results are likely to be
iweights <- get_initial_weights(targets_wide, base_data)
ht(iweights)

# DEFINE constraint coefficients (cc) ----
cc_sparse <- get_constraint_coefficients(base_data, target_vars, iweights, targlist) # can take a while on big problems

# BUILD inputs list ----
inputs <- get_inputs(.targlist=targlist, .iweights=iweights, .cc_sparse=cc_sparse, .targtol=.02, .xub=50, .conscaling=FALSE, scale_goal=100)

inputs_save <- inputs
# inputs <- inputs_save

#.. examine inputs ----
names(inputs)
inputs$n_constraints
inputs$n_targets
ht(inputs$cc_sparse)
check <- check_constraints(inputs)
check %>% arrange(desc(abs(pdiff)))

#.. OPTIONAL revise inputs ----
check <- check %>%
  mutate(clb2=ifelse(i==7, clb, clb),
         cub2=ifelse(i==7, 10e3, cub))

# inputs$clb <- check$clb2
# inputs$cub <- check$cub2
#.. END OPTIONAL revise inputs ----


# CAUTION: If you use ma77, be sure to delete temporary files it created in the project
# folder before rerunning or RStudio may bomb.
# For LARGE problems use linear_solver=ma77, obj_scaling=1, and mehrotra_algorithm=yes
opts <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "max_iter"= 100,
             "linear_solver" = "ma77", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             # "ma57_automatic_scaling" = "yes", # if using ma57
             # "ma57_pre_alloc" = 3, # 1.05 is default; even changed, cannot allocate enough memory, however
             # "ma77_order" = "amd",  # metis; amd -- not clear which is faster
             "mehrotra_algorithm" = "yes",
             "obj_scaling_factor" = 1, # 1e-3, # default 1; 1e-1 pretty fast to feasible but not to optimal
             # "nlp_scaling_method" = "equilibration-based", # NO - use default gradient_based
             # "nlp_scaling_max_gradient" = 100, # default is 100 - seems good
             # "jac_c_constant" = "yes", # does not improve on moderate problems
             # "jac_d_constant" = "yes", # does not improve on  moderate problems
             # "hessian_constant" = "yes", # KEEP default NO - if yes Ipopt asks for Hessian of Lagrangian function only once and reuses; default "no"
             # "hessian_approximation" = "limited-memory", # KEEP default of exact
             # "derivative_test" = "first-order",
             # "derivative_test_print_all" = "yes",
             "output_file" = here::here("out", "test.out"))

result <- ipoptr(x0 = inputs$x0,
                 lb = inputs$xlb,
                 ub = inputs$xub,
                 eval_f = eval_f_xtop, # arguments: x, inputs
                 eval_grad_f = eval_grad_f_xtop,
                 eval_g = eval_g, # constraints LHS - a vector of values
                 eval_jac_g = eval_jac_g,
                 eval_jac_g_structure = inputs$eval_jac_g_structure,
                 eval_h = eval_h_xtop, # the hessian is essential for this problem
                 eval_h_structure = inputs$eval_h_structure,
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

# print(Sys.unsetenv("OMP_NUM_THREADS")) # unset the threads environment variable

# Sys.getenv(c("R_HOME", "OMP_NUM_THREADS"))
# print(Sys.setenv(OMP_NUM_THREADS = 6))
# Sys.getenv(c("R_HOME", "OMP_NUM_THREADS"))