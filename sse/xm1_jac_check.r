# Create small problem for testing hessian

# TODO: construct without any dense info
#       create structure
#       test in ipopt


# libraries ----
library(magrittr)
library(plyr) # needed for ldply; must be loaded BEFORE dplyr
library(tidyverse)
options(tibble.print_max = 65, tibble.print_min = 65) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats

library(btools) # has a few useful functions -- devtools::install_github("donboyd5/btools")

library(optimx)
library(ipoptr)

library(numDeriv)

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
glimpse(samp)

# prepare file
samp2_base <- samp1 %>% 
  select(-nrecs, -pop) %>%
  mutate(pid=row_number(), # pid -- an id variable for each person in the file
         incgroup=ntile(pincp, 10), # divide the data into 10 income ranges
         pop=1, # it's useful to have a variable that is 1 on every record
         # convert categoricals to dummies if we will base targets upon them
         mar1=ifelse(mar==1, 1, 0), # married
         mar5=ifelse(mar==5, 1, 0), # single
         marx15=ifelse(mar %nin% c(1, 5), 1, 0)
  )
count(samp2_base, stabbr, incgroup)

# create subset
set.seed(123)
samp2 <- samp2_base %>%
  group_by(stabbr, incgroup) %>%
  sample_frac(size=.02) %>%
  ungroup %>%
  filter(stabbr %in% c("CA", "FL", "NY"))
count(samp2, stabbr, incgroup)  
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
ns(samp)


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
(target_vars <- possible_target_vars[c(1, 6)])
# target_vars <- possible_target_vars

# djb adjust here ----

samp_puf <- prep_data(samp, possible_target_vars)
ns(samp_puf)
incgroup_data <- prep_incgroup(.data=samp_puf, .incgroup=target_incgroup, .weight=pwgtp, .target_vars=target_vars)
head(incgroup_data)

#incgroup_data %>% filter(pid==197)

t1 <- incgroup_data %>%
  summarise(across(all_of(target_vars), ~ sum(.x * weight_total)))

t2 <- targets_wide %>%
  summarise(across(all_of(target_vars), sum))

bind_rows(t1 %>% mutate(type="data"), t2 %>% mutate(type="target")) %>%
  pivot_longer(cols=-type) %>%
  pivot_wider(names_from = type) %>%
  mutate(pdiff=target / data * 100 - 100)

# djb revise the targets by the ratio of data to target for N1_nnz ?? ----

targets_df <- get_targets(targets_wide, target_vars, .pid=incgroup_data$pid, .weight_total=incgroup_data$weight_total)
head(targets_df)
tail(targets_df)
targets_df$i

targets_df %>%
  filter(targtype != "addup")
targets_df %>%
  filter(vname_calctype=="A00200_sum") %>%
  summarise(value=sum(value))

targets_df %>%
  filter(targtype!="aggregate")

# CONSTRUCT initial weights ----
#.. prepare the data we will reweight to hit (or come close to) the targets ----
# let's stack the data so that each person appears 1 time for each state that is targeted
# create a stub with a record for every person and every targeted state

# define initial weight for each person-state combination. This is important as the
# optimization will try to keep the final weights near these initial weights. The
# more we can improve these initial weights, the better the results are likely to be
iweights <- get_initial_weights(targets_wide, incgroup_data, .popvar=pop_nnz)
head(iweights)

# djb -- check whether state weights add to total weights --
# iweights %>%
#   filter(pid <= 2000) %>%
#   group_by(pid) %>%
#   summarise(wtotal=first(weight_total), wstates=sum(iweight_state))
# end ---

# DEFINE constraint coefficients (cc) ----
# cc_sparse <- get_constraint_coefficients(incgroup_data, target_vars, iweights, targets_df) # can take a while on big problems

cc_sparse <- get_constraint_coefficients(.incgroup_data=incgroup_data,
                                         .target_vars=target_vars,
                                         .iweights=iweights,
                                         .targets_df=targets_df) # can take a while on big problems
cc_sparse %>%
  filter(targtype=="aggregate")

cc_sparse %>%
  filter(pid==741, str_detect(cname, "pincp_sum")) %>%
  summarise(nzcc=sum(nzcc))

# samp_puf %>% filter(pid==197)
incgroup_data %>% filter(pid==741) %>% as.data.frame()

# count(cc_sparse, targtype)
# # djb ?? ----
# # cc_sparse <- cc_sparse %>%
# #   filter(targtype=="aggregate")
# # djb end ----
# quantile(cc_sparse$i)
# quantile(cc_sparse$j)


# djb look at constraint coefficients ----
# summary(cc_sparse)
# cc_sparse %>%
#   filter(j <= 10) %>%
#   group_by(j, cname, targtype, STATE) %>%
#   summarise(nzcc=sum(nzcc))

# end ----

# DEFINE inputs for ipopt ----
# inputs_us <- get_inputs(.targets_df=targets_df, .iweights=iweights, .cc_sparse=cc_sparse,
#                      .targtol=.02, .adduptol=0.2, .xub=50, .conscaling=FALSE, scale_goal=100)

inputs <- get_inputs(.targets_df=targets_df, .iweights=iweights, .cc_sparse=cc_sparse,
                     .targtol=.05, .adduptol=0, .xub=50, .conscaling=FALSE, scale_goal=100)
names(inputs)
# djb end adjust ----


# investigate jacobian ----
library(numDeriv)
jacfn <- eval_jac_g(inputs$x0, inputs)
system.time(jacfd <- jacobian(eval_g, inputs$x0, inputs=inputs))
dim(jacfd)
jacfd[1, 1:10]
jacfn[1:10]
sum(jacfd)
sum(jacfn)
as.vector(jacfd[jacfd!=0]) 
jacfn

jacfn_str <- define_jac_g_structure_sparse(inputs$cc_sparse)
jacfd_str <- make.sparse(jacfd)

identical(jacfn_str, jacfd_str)
# end jac check -----

opts <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "max_iter"= 20,
             "linear_solver" = "ma86", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             #"linear_system_scaling" = "mc19",
             #"linear_scaling_on_demand" = "no", # default is yes -- no means ALWAYS scale
             # "ma57_automatic_scaling" = "yes", # if using ma57
             # "ma57_pre_alloc" = 3, # 1.05 is default; even changed, cannot allocate enough memory, however
             # "ma77_order" = "amd",  # metis; amd -- not clear which is faster
             "mehrotra_algorithm" = "yes",
             #"obj_scaling_factor" = 1, # 1e-3, # default 1; 1e-1 pretty fast to feasible but not to optimal
             #"nlp_scaling_method" = "none", # NO - use default gradient_based, none, equilibration-based
             # "nlp_scaling_max_gradient" = 100, # default is 100 - seems good
             #"jac_c_constant" = "yes", # does not improve on moderate problems equality constraints
             # "jac_d_constant" = "yes", # does not improve on  moderate problems inequality constraints
             #"hessian_constant" = "yes", # KEEP default NO - if yes Ipopt asks for Hessian of Lagrangian function only once and reuses; default "no"
             # "hessian_approximation" = "limited-memory", # KEEP default of exact
             # "limited_memory_update_type" = "bfgs", # default bfgs; sr1 docs say "not working well" but this really helps
             # "derivative_test" = "first-order",
             #"derivative_test_print_all" = "yes",
             "output_file" = "C:/RPrograms PC/OSPC/50_state_taxdata/out/check2.out")

# inputs$objscale <- 1

# setwd(here::here("temp1"))
# getwd()
result <- ipoptr(x0 = inputs$x0,
                 lb = inputs$xlb,
                 ub = inputs$xub,
                 eval_f = eval_f_xm1sq, # arguments: x, inputs; eval_f_xtop eval_f_xm1sq
                 eval_grad_f = eval_grad_f_xm1sq, # eval_grad_f_xtop eval_grad_f_xm1sq
                 eval_g = eval_g, # constraints LHS - a vector of values
                 eval_jac_g = eval_jac_g,
                 eval_jac_g_structure = inputs$eval_jac_g_structure,
                 eval_h = eval_h_xm1sq, # the hessian is essential for this problem eval_h_xtop eval_h_xm1sq
                 eval_h_structure = inputs$eval_h_structure,
                 constraint_lb = inputs$clb,
                 constraint_ub = inputs$cub,
                 opts = opts,
                 inputs = inputs)


