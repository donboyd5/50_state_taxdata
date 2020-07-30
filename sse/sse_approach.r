

# libraries ----
library(magrittr)
library(plyr) # needed for ldply; must be loaded BEFORE dplyr
library(tidyverse)
options(tibble.print_max = 65, tibble.print_min = 65) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats

library(btools) # has a few useful functions -- devtools::install_github("donboyd5/btools")

library(optimx)
library(ipoptr)

devtools::session_info()


# data location ----
dbox <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/"

# functions ----
source(here::here("sse", "sse_functions_prep.r"))
source(here::here("sse", "sse_functions_optim.r"))
source(here::here("sse", "sse_functions_eval.r"))

# Prepare data ----
# choose file
(fnames <- paste0(c("acs_10krecs_5states", "acs_100krecs_20states", "acs_200krecs_50states", "acs_400krecs_50states"), ".rds"))
samp1 <- readRDS(here::here("data", fnames[1]))
glimpse(samp1)

# prepare file
samp2 <- samp1 %>% 
  select(-nrecs, -pop) %>%
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

# define the kinds of (weighted) targets we want and prepare the file accordingly
# sum:    sum of values
# nnz:    number of nonzero values
# sumneg: sum of negative values
# nneg:   number of zero values
# sumpos: sum of positive value
# npos:   number of positive values

nnz_vars <- c("pop", "mar1", "mar5", "pincp", "wagp") # note that I leave the 3rd marital status out -- marx15
sum_vars <- c("pincp", "wagp", "intp", "pap", "retp", "ssip", "ssp") # DO NOT include total plus all sums - leave one out (otherincp)
sumneg_vars <- "otherincp"

# define a vector of variable names for "possible" targets (a superset) -- we may not target all
possible_target_vars <- make_target_names(
  list(nnz_vars, sum_vars, sumneg_vars),
  c("nnz", "sum", "sumneg"))
possible_target_vars

# prepare data by creating variables with those names:
samp <- prep_data(samp2, possible_target_vars)
glimpse(samp)


# Prepare targets for the full file ----
summary_vals <- samp %>%
  group_by(stabbr, incgroup) %>%
  summarise(nrecs=n(),
            across(.cols=(all_of(possible_target_vars)), ~ sum(.x * pwgtp)),
            .groups="drop")

# Prepare a single income group ----
target_incgroup <- 2 # define target income group

possible_target_vars

# use a subset of possible targets, or all (run only 1 of the next 2 lines)
(target_vars <- possible_target_vars[c(1, 3, 6, 7)])
target_vars <- possible_target_vars

incgroup_data <- samp %>%
  filter(incgroup==target_incgroup) %>%
  arrange(pid) %>%
  mutate(i=row_number()) %>%
  select(i, pid, incgroup, weight_total=pwgtp, all_of(target_vars))

targets_wide <- summary_vals %>%
  filter(incgroup==target_incgroup) %>%
  select(stabbr, incgroup, nrecs, all_of(target_vars))
targets_wide

targets_df <- targets_wide %>%
  pivot_longer(cols=all_of(target_vars), names_to = "vname_calctype", values_to="target") %>%
  separate(vname_calctype, c("vname", "calctype"), sep="_", remove=FALSE) %>%
  mutate(targname=get_targname(vname_calctype, stabbr)) %>%
  arrange(targname) %>%
  mutate(targnum=row_number()) %>%
  select(targnum, targname, stabbr, vname, calctype, target)
ht(targets_df) # 1 row per target per state

# Construct initial weights using national ratios -- 1 row per person per state ----
iweights <- expand_grid(pid=unique(incgroup_data$pid), stabbr=targets_wide$stabbr) %>%
  left_join(incgroup_data %>% select(i, pid, weight_total), by = "pid") %>%
  left_join(targets_wide %>% select(stabbr, pop_target=pop_nnz), by="stabbr") %>%
  group_by(stabbr) %>%
  mutate(iweight_state=weight_total / sum(weight_total) * pop_target) %>%
  ungroup %>%
  arrange(i, stabbr) %>%
  mutate(j=row_number(),
         xname=paste0("p", pid, "_", stabbr)) %>% # index for variable we are solving for
  select(j, xname, i, pid, stabbr, pop_target, weight_total, iweight_state)
iweights
ht(iweights)

# Define objective function elements (ofe) ----
# dense 1 row per variable (person-state)
ofe_dense <- iweights %>%
  left_join(incgroup_data %>% select(pid, all_of(target_vars)), by = "pid") %>%
  mutate_at(vars(all_of(target_vars)), list(~ . * iweight_state)) %>% # coefficients
  arrange(pid, stabbr)

# objective function elements
# sparse format, 1 row per variable (person-state) per target, if nonzero coefficient
# rows are individuals pid, columns are j xname (pid-stabbr) variables 
# cells are iweight_state
# get totals for each person, and ipid
ofe_sparse <- ofe_dense %>% # has j and xname already
  pivot_longer(cols = all_of(target_vars),
               names_to="vname_calctype",
               values_to = "nz_coef") %>%
  filter(nz_coef!=0) %>%
  mutate(targname=get_targname(vname_calctype, stabbr)) %>%
  left_join(targets_df %>% select(targnum, targname), by = "targname") %>% # get targnum
  select(i, pid, stabbr, targnum, j, targname, xname, weight_total, iweight_state, nz_coef)
summary(ofe_sparse) # make sure there are no NA values 
ht(ofe_sparse)


# Define inputs ----

inputs <- get_inputs_sse(targets_df, iweights, ofe_sparse, target_scaling=TRUE, scale_goal=100)

# step 1 optim without constraints ----
eval_f_sse(inputs$x0, inputs) # obj function at starting point

a <- proc.time()
sse_opt <- optimx::optimr(par=inputs$x0,
                          fn = eval_f_sse,
                          gr = eval_grad_f_sse,
                          lower = inputs$xlb,
                          upper = inputs$xub,
                          method = "L-BFGS-B", # L-BFGS-B
                          control = list(trace=1, 
                                          maxit=1000, 
                                          kkt=FALSE,
                                          starttests=FALSE),
                          inputs=inputs)
b <- proc.time()
b - a
str(sse_opt)

# how did we do against targets?
names(inputs)
inputs$targets_df
check_targ <- check_targets(sse_opt$par, inputs)

check_targ %>%
  select(targnum, targname, target_unscaled, calctarg, diff, pdiff) %>%
  arrange(-abs(pdiff))

# how bad were they at the starting point?
check_targets(inputs$x0, inputs) %>%
  select(targnum, targname, target_unscaled, calctarg, diff, pdiff) %>%
  arrange(-abs(pdiff))


# how close are the sums of state weights to total weights?
cw <- check_weights(sse_opt$par, inputs)
cw %>% ht
sse_weights(cw)
cw %>%
  ggplot(aes(weight_total, weight)) +
  geom_point() +
  geom_abline(slope= 1 , intercept = 0) +
  labs(x="National weight in the data file", 
       y="Sum of state weights, calculated")


# step 2 prepare adjusted starting point for ipopt ----
# calculate starting point that satisfies adding-up requirement
wadjust <- inputs$iweights_df %>%
  mutate(xdf=sse_opt$par,
         wstate=iweight_state * xdf) %>%
  group_by(i, pid) %>%
  mutate(wstate_adj=wstate / sum(wstate) * weight_total,
         xadj=wstate_adj / iweight_state) %>%
  ungroup

# how much worse is the SSE of targets at this starting point
eval_f_sse(sse_opt$par, inputs)
eval_f_sse(wadjust$xadj, inputs)
eval_f_sse(inputs$x0, inputs)

# verify that weights are good at this starting point
cw <- check_weights(wadjust$xadj, inputs)
cw %>% ht
sse_weights(cw)
cw %>%
  ggplot(aes(weight_total, weight)) +
  geom_point() +
  geom_abline(slope= 1 , intercept = 0) +
  labs(x="National weight in the data file", 
       y="Sum of state weights, calculated")



# step 3 ipopt optimization with adding up constraints ----

opts_ipopt <- list("print_level" = 0,
                   "file_print_level" = 5, # integer
                   "max_iter"= 200,
                   "linear_solver" = "ma86", # mumps pardiso ma27 ma57 ma77 ma86 ma97
                   #"mehrotra_algorithm" = "yes", # default no
                   # "jac_c_constant" = "yes", # default no
                   #"hessian_approximation" = "limited-memory", # default exact; limited-memory
                   #"limited_memory_update_type" = "sr1", # default bfgs; sr1 docs say "not working well" but this really helps
                   #"hessian_constant" = "yes", # default no; yes
                   "obj_scaling_factor" = 1e1, # 1e-3, # default 1; 1e-1 pretty fast to feasible but not to optimal
                   # "nlp_scaling_max_gradient" = 100, # default is 100 - seems good
                   "output_file" = here::here("out", "test.out"))

# setwd(here::here("temp1"))
#getwd()
ip1 <- ipoptr(x0 = wadjust$xadj, # inputs$x0 wadjust$xadj sse_opt$par
              lb = inputs$xlb,
              ub = inputs$xub,
              eval_f = eval_f_sse, # arguments: x, inputs; eval_f_xtop eval_f_xm1sq eval_f_wfs
              eval_grad_f = eval_grad_f_sse, # eval_grad_f_xm1sq, # eval_grad_f_xtop eval_grad_f_xm1sq
              eval_g = eval_g_addup, # constraints LHS - a vector of values
              eval_jac_g = eval_jac_g_addup,
              eval_jac_g_structure = inputs$eval_jac_g_structure,
              # eval_h = eval_h_xm1sq, # the hessian is essential for this problem eval_h_xtop eval_h_xm1sq
              # eval_h_structure = inputs$eval_h_structure,
              constraint_lb = inputs$constraints,
              constraint_ub = inputs$constraints,
              opts = opts_ipopt,
              inputs = inputs)


# how did we do against targets?
names(inputs)
inputs$targets_df
check_targ <- check_targets(ip1$solution, inputs)

check_targ %>%
  select(targnum, targname, target_unscaled, calctarg, diff, pdiff) %>%
  arrange(-abs(pdiff))

# verify that sums of state weights equal total weights
cw <- check_weights(ip1$solution, inputs)
cw %>% ht
sse_weights(cw)
cw %>%
  ggplot(aes(weight_total, weight)) +
  geom_point() +
  geom_abline(slope= 1 , intercept = 0) +
  labs(x="National weight in the data file", 
       y="Sum of state weights, calculated")
