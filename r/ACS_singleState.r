
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
(fns <- paste0(c("acs_10krecs_5states", "acs_100krecs_20states", "acs_200krecs_50states"), ".rds"))

# functions ----
source(here::here("include", "functions_ipopt.r"))

# GET SAVED DATA ----

# to create data used in this program, as needed, run:
#
#    create_ACS_subset_from_extract.r
#
# to create a subset suitable for targeting

samp <- readRDS(here::here("data", fns[2])) # choose which file to use 

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


# SINGLE-STATE optimization - try to hit the targets for a single & income range ----
#.. define (1) a stabbr-incgroup to target, (2) the variables to target, and (3) the target values ----
# for now just use the actual values (we could introduce noise if we want)

# define a state and income group that we will target; also, pick only a few variables to target
target_stabbr <- "CA"
target_incgroup <- 2
target_vars <- c("pop", "pincp", "wagp", "intp", "ssip_n") # keep it small for our test
# CAUTION: in code below we'll keep variables in this order so that we can index into matrices by column
# we can put the variables above in any order we want, but once we choose that order we want to be sure to
# respect it below; this will be clearer later

# since we have it, let's keep the "true" data for comparison to the results of targeting
# this, of course, is unknown when we are working with the PUF - we only know the summary values
data_true <- samp %>%
  filter(stabbr==target_stabbr, incgroup==target_incgroup) %>%
  select(pid, stabbr, incgroup, pwgtp, nrecs, all_of(target_vars))
data_true 

# define the targets we want to hit for a given state and income, using data for all states in that income group
targets <- summary_vals %>%
  filter(stabbr==target_stabbr, incgroup==target_incgroup) %>%
  select(stabbr, incgroup, nrecs, all_of(target_vars)) # a small list of variables to target
targets # these are the targets we want to hit
# we also need these as a vector -- make sure they are in the same order
(constraints <- targets %>% select(all_of(target_vars)) %>% unlist(., use.names=FALSE))
# we want to pull all data in this income group and try to make it look like the targeted state in this group

#.. prepare the data we will reweight to hit (or come close to) the targets ----
# determine which variables we need counts for
data_targ <- samp %>%
  filter(incgroup==targets$incgroup) %>% # we want records from the targeted income group, but for all states
  mutate(j=row_number(), # j because each row here will be a column in a constraint-coefficients
         weight=pwgtp / sum(pwgtp) * targets$pop) %>% # create initial weight that will sum to targeted weight sum
  select(pid, j, stabbr, incgroup, pwgtp, weight, all_of(target_vars)) # note that this keeps the variables in desired order
data_targ
sum(data_targ$weight); targets$pop; sum(data_true$pwgtp) # verify that initial weight sums to targeted population
# the idea is that we will adjust these initial weights so that not only do they hit the targeted
# population, but they also hit the desired values
# We will do this by multiplying each weight by a ratio
# we'll penalize distances in this ratio from 1

# do a quick check to see how close we come to the targets with this initial weight
check <- data_targ %>%
  summarise_at(vars(all_of(target_vars)),
               list(~ sum(. * weight)))
bind_rows(targets, check) # later we may write a function for comparisons but this is good enough for now

#.. define constraint coefficients: change in constraint for 1-unit change in ratio of new weight to initial weight ----
# imagine a matrix where:
#   we have a row for each constraint (target), and
#   a row for each x variable (each weight to be found), and
#   each cell is the amount by which the constraint changes if the weight for that record changes by 1
# for example suppose 
#   row 3 is the target for total wages, and
#   row 7 is the constraint for the # of people who have SSI income (i.e., where ssip > 0)
#   
# construct the matrix in dense form; cc stands for constraint coefficient
cc_dense <- data_targ %>%
  mutate_at(vars(all_of(target_vars)), list(~ . * weight))

cc_sparse <- cc_dense %>%
  pivot_longer(cols = all_of(target_vars), values_to = "nzcc") %>% # nzcc for nonzero constraint coefficient
  filter(nzcc!=0) %>%
  mutate(i=match(name, target_vars)) %>%
  select(j, i, name, nzcc, weight, stabbr, pwgtp, incgroup) %>%
  arrange(i, j) # this ordering is crucial for the Jacobian

#.. set up the optimization problem ----
inputs <- list()
inputs$p <- 2
inputs$iweight <- data_targ$weight # the initial weight
inputs$cc_sparse <- cc_sparse
inputs$n_variables <- length(inputs$iweight)
inputs$n_constraints <- length(target_vars)
inputs$objscale <- 1
# inputs$constraint_scales <- constraint_scales # later we may need to scale constraints
str(inputs)

xlb <- rep(0, inputs$n_variables)
xub <- rep(10, inputs$n_variables)
x0 <- rep(1, inputs$n_variables)

# define constraint lower and upper bounds -- for now they will be the same
tol <- 0 # tolerance
clb <- constraints - abs(constraints) * tol
cub <- constraints + abs(constraints) * tol
cbind(clb, constraints, cub) %>% round(0)

eval_jac_g_structure <- define_jac_g_structure_sparse(inputs$cc_sparse, ivar="i", jvar="j")
eval_jac_g_structure[5]
length(eval_jac_g_structure)

eval_h_structure <- lapply(1:inputs$n_variables, function(x) x) # diagonal elements of our Hessian

eval_f_xtop(x0, inputs)# ; eval_f_absapprox(xlb, inputs); eval_f_absapprox(xub, inputs)
eval_grad_f_xtop(x0, inputs)
eval_g(x0, inputs)
check
eval_jac_g(x0, inputs)
eval_h_xtop(x0, obj_factor=1, hessian_lambda=rep(1, inputs$n_constraints), inputs) # length is n_variables
# length(unlist(eval_h_structure)) # length should be n_variables


# ma86 was faster in one test I did
opts <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "linear_solver" = "ma86", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             "max_iter"= 100,
             #"obj_scaling_factor" = 10, # 1e7, # default 1
             #"nlp_scaling_max_gradient" = 1e-3, # 1e-3, # default 100
             # "derivative_test" = "first-order",
             #"derivative_test_print_all" = "yes",
             "output_file" = here::here("out", "single_state.out"))

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
                 constraint_lb = clb,
                 constraint_ub = cub,
                 opts = opts,
                 inputs = inputs)

str(result)
quantile(result$solution)
result$constraints
eval_g(x0, inputs)
eval_g(result$solution, inputs)

# calculate the new weights using the solution
# stack the 3 files

data_new <- data_targ %>%
  mutate(x=result$solution,
         weight_new=weight * x)

check2 <- data_new %>%
  summarise_at(vars(all_of(target_vars)),
               list(~ sum(. * weight_new)))

bind_rows(targets, check, check2)

# how much of CA incgroup2 came from each state?
data_new %>%
  group_by(stabbr) %>%
  summarise(n=n(), pwgtp=sum(pwgtp), weight_new=sum(weight_new)) %>%
  mutate_at(vars(-stabbr), list(pct= ~ . / sum(.) * 100))

# at this point, given that all is working, we can add constraints if we want
