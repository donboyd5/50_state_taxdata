
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

# GET SAVED DATA ----

# to create data used in this program, as needed, run:
#
#    create_ACS_subset_from_extract.r
#
# to create a subset suitable for targeting

samp <- readRDS(here::here("data", fns[1])) # choose which file to use 

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


# MULTIPLE-STATE OPTIMIZATION, 1 income range, no adding-up constraint ----

# define a state and income group that we will target; also, pick only a few variables to target
target_incgroup <- 2
target_vars <- c("pop", "pincp", "wagp", "intp", "ssip_n") # keep it small for our test
# CAUTION: in code below we'll keep variables in this order so that we can index into matrices by column
# we can put the variables above in any order we want, but once we choose that order we want to be sure to
# respect it below; this will be clearer later

# since we have it, let's keep the "true" data for comparison to the results of targeting
# this, of course, is unknown when we are working with the PUF - we only know the summary values
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

# we also need these as a vector -- make sure they are in the same order
(constraint_names <- targets_long$cname)
(constraints <- targets_long$value)
# constraints vector order is: 1 variable all states, then next variable all states, then ...
# we want to pull all data in this income group and try to make it look like the targeted state in this group


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

cc_dense <- stub %>%
  mutate_at(vars(all_of(target_vars)), list(~ . * iweight))

cc_sparse <- cc_dense %>%
  pivot_longer(cols = all_of(target_vars), values_to = "nzcc") %>% # nzcc for nonzero constraint coefficient
  filter(nzcc!=0) %>%
  mutate(cname=paste0(name, "_", stabbr),
         i=match(cname, constraint_names)) %>%
  select(pid, j, i, cname, nzcc, iweight, name, stabbr, stabbr_true, pwgtp, incgroup) %>%
  arrange(i, j) # this ordering is crucial for the Jacobian
ht(cc_sparse)
count(cc_sparse, i)

#.. set up the optimization ----
inputs <- list()
inputs$p <- 2
inputs$iweight <- cc_dense$iweight # the initial weight
inputs$cc_sparse <- cc_sparse
inputs$n_variables <- length(inputs$iweight)
inputs$n_constraints <- length(constraints)
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
tibble(constraint_names, clb, constraints, cub) %>%
  kable(digits=0, format="rst", format.args=list(big.mark=","))

eval_jac_g_structure <- define_jac_g_structure_sparse(inputs$cc_sparse, ivar="i", jvar="j")
eval_jac_g_structure[5]
length(eval_jac_g_structure)

eval_h_structure <- lapply(1:inputs$n_variables, function(x) x) # diagonal elements of our Hessian

eval_f_xtop(x0, inputs)# ; eval_f_absapprox(xlb, inputs); eval_f_absapprox(xub, inputs)
eval_grad_f_xtop(x0, inputs)

tibble(constraint_names, clb, constraints, conx0=eval_g(x0, inputs), cub) %>%
  kable(digits=0, format="rst", format.args=list(big.mark=","))

eval_jac_g(x0, inputs)
eval_h_xtop(x0, obj_factor=1, hessian_lambda=rep(1, inputs$n_constraints), inputs) # length is n_variables
# length(unlist(eval_h_structure)) # length should be n_variables


opts <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "linear_solver" = "ma57", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             "max_iter"= 100,
             #"obj_scaling_factor" = 10, # 1e7, # default 1
             #"nlp_scaling_max_gradient" = 1e-3, # 1e-3, # default 100
             # "derivative_test" = "first-order",
             #"derivative_test_print_all" = "yes",
             "output_file" = here::here("out", "multiple_states.out"))

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

#.. examine results ----
str(result)
quantile(result$solution)
result$constraints
tibble(constraint_names, clb, constraints, cub, 
       conx0=eval_g(x0, inputs), confinal=eval_g(result$solution, inputs)) %>%
  kable(digits=0, format="rst", format.args=list(big.mark=","))


stub %>%
  mutate(x=result$solution,
         weight=iweight * x) %>%
  group_by(stabbr) %>%
  summarise_at(vars(all_of(target_vars)), list(~ sum(. * weight)))
targets

# what states did the records come from?
stub %>%
  mutate(x=result$solution,
         weight=iweight * x) %>%
  group_by(stabbr_true, stabbr) %>%
  summarise_at(vars(pwgtp, weight), sum) %>%
  pivot_wider(names_from = stabbr, values_from = weight)

# how did we do against the adding-up goal?
ratios <- stub %>%
  mutate(x=result$solution,
         weight=iweight * x) %>%
  group_by(p) %>%
  summarise(pwgtp=first(pwgtp), weight=sum(weight)) %>%
  mutate(ratio=weight / pwgtp)
quantile(ratios$ratio, probs = c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1))

ratios %>% filter(ratio > 1.05)
stub %>%
  mutate(x=result$solution,
         weight=iweight * x) %>%
  filter(p==22)
summary(samp %>% filter(incgroup==2))

