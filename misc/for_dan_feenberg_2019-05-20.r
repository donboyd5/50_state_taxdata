
# For Dan Feenberg
# 5/20/2019

# Based on prior versions:
    # Don Boyd
    # 7/13/2017

# Possible extended Ipopt options??
# http://www.pserc.cornell.edu/Matpower/docs/ref/matpower5.1/ipopt_options.html#_top


#****************************************************************************************************
#                Load libraries ####
#****************************************************************************************************
# url <- "http://users.nber.org/~taxsim/R/taxsim9_0.1.tar.gz"
# install.packages(url, repos = NULL, type = "source")

library("magrittr")
library("plyr") # needed for ldply; must be loaded BEFORE dplyr
library("tidyverse")
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr

library("scales")
library("hms") # hms, for times.
library("stringr") # stringr, for strings.
library("lubridate") # lubridate, for date/times.
library("forcats") # forcats, for factors.
library("readxl") # readxl, for .xls and .xlsx files.
library("haven") # haven, for SPSS, SAS and Stata files.
library("vctrs")

# library("precis")

# library("grDevices")
# library("knitr")
# 
# library("zoo") # for rollapply

# devtools::install_github("donboyd5/btools")
# library("btools")
# devtools::install_github("donboyd5/bdata")
# library("bdata")

library("ipoptr")


#****************************************************************************************************
#                Load functions ####
#****************************************************************************************************
# source("./r/functions_general.r")
# source("./r/functions_problem_setup(6).r")
# source("./r/functions_ipopt.r")
# source("./r/functions_get_test_problems.r")


#****************************************************************************************************
#                Get SOI targets data ####
#****************************************************************************************************
#.. Define a few minimal targets for a few states, based upon the 2006 Historical Table 2 AGI >= $200k ----
# get soi state summaries
soi2006 <- readRDS("./data/soi_states_2006.rds") # amounts are in $k
write_csv(soi2006, "./data/soi_states_2006.csv")
head(soi2006)
count(soi2006, incrange)
count(soi2006, stabbr)
count(soi2006, line, vname) %>% as.data.frame()

# create some variable names to map to puf variable names
soimap <- read_csv(
"line, soivar, note
1, wt, weight
2, MARS2, n joint returns
4, XTOT, total exemptions
5, E00100, agi
6, E00200_n, n wages
7, E00200, $ wages
15, netcg_n, n net capital gain
16, netcg, $ net cap gain
29, E18425_n, n state and local income tax
30, E18425, $ state and local income tax
33, E18500_n, n state and local real estate taxes paid
34, E18500, $ state and local real estate taxes paid
42, E04800, $ taxable income
58, E06500, $ total income tax
55, E09600_n, # of AMT returns
56, E09600, $ AMT")
soimap


# make a small set of targets
soilines <- c(1, 2, 4, 5, 30)
targets1 <- left_join(soimap %>% filter(line %in% soilines),
                      soi2006 %>% filter(incrange=="ge200k", stabbr %in% c("AL", "CA", "NY", "US"))) %>%
  mutate(value=ifelse(line %in% c(5, 7), value * 1000, value))
targets1

# we need to create an "other" area for each target
targets_other <- targets1 %>%
  filter(stabbr!="US") %>%
  group_by(line, soivar, note, vname, incrange) %>%
  summarise(sumval=sum(value)) %>%
  left_join(targets1 %>% filter(stabbr=="US") %>% rename(usval=value) %>% select(-stabbr)) %>%
  mutate(value=usval - sumval, stabbr="other") %>%
  select(-sumval, -usval)

targets <- bind_rows(targets1 %>% filter(stabbr != "US"),
                     targets_other) %>%
  arrange(line, stabbr)
targets


#****************************************************************************************************
#                Get puf high income returns ####
#****************************************************************************************************
# .. first, get the puf ----
# ..ONETIME save puf as rds file
# soid <- "d:/Data/SOI PUF/" # directory in which the SAS puf is located
# fn.sas <- "puf_2006.sas7bdat" # This is the sas file originally provided by SOI -- has 145,898 records
# puf <- read_sas(paste0(soid, fn.sas))
# glimpse(puf)
# saveRDS(puf, "./data/puf2006.rds")
# ..END ONETIME save

puf <- readRDS("./data/puf2006.rds")
glimpse(puf) # 145,898 recs
count(puf, STATE)
# quantile(puf$RECID)
# quantile(puf$S006)


# preliminary processing - get JUST the uncoded returns with agi >= 200k
puf2 <- puf %>% 
  mutate(wt=S006 / 100,
         RECID=as.integer(RECID)) %>%
  filter(STATE==0, E00100 >= 200e3)


#****************************************************************************************************
#                Prepare PUF for constraint evaluation ####
#****************************************************************************************************
puf3 <- puf2 %>%
  # define any variables or one-zero indicators needed for constraints, to be multiplied by wt
  mutate(netcg=E01000 + E01100,
         wt_one=1, # when multiplied by wt will give the weighted value
         MARS2=ifelse(MARS==2, 1, 0)) %>%
  select(RECID, wt_one, wt, MARS2, XTOT, E00100, E18425) # same variables as targets
glimpse(puf3)


#****************************************************************************************************
#    Preprocessing: calibrate puf weights so that they sum to Historical Table 2 ####
#****************************************************************************************************
(soi_wtsum <- targets %>% filter(soivar=="wt") %>% summarise(value=sum(value)) %>% .[["value"]])
(puf_wtsum <- puf3 %>% summarise(wt=sum(wt)) %>% .[["wt"]])

puf4 <- puf3 %>% 
  rename(wtold=wt) %>%
  mutate(wt=wtold * soi_wtsum / puf_wtsum)
glimpse(puf4)
puf4 %>% summarise(wt=sum(wt))

pufdf <- puf4

#****************************************************************************************************
#                functions for putting constraints into desired form ####
#****************************************************************************************************

pufsum_f <- function(pufdf){
  # summarise puf in steps
  pufsum.a <- pufdf %>% 
    summarise_at(vars(wt), funs(sum(.) / 1e3)) %>%
    gather(soivar, pufval)
  
  pufsum.x <- pufdf %>%
    summarise_at(vars(XTOT, MARS2), funs(sum(. * wt) / 1e3)) %>%
    gather(soivar, pufval)
  
  pufsum.b <- pufdf %>% 
    summarise_at(vars(E00100, E18425),
                 funs(val=sum(. * wt) / 1e6)) %>%
    gather(soivar, pufval) %>%
    mutate(soivar=ifelse(str_detect(soivar, "_val"),
                         str_replace(soivar, "_val", ""), 
                         soivar)) 
  
  pufsum <- bind_rows(pufsum.a, pufsum.x, pufsum.b)
  
  return(pufsum)
}


get_ccmatdf <- function(pufdf, cnames, sts){
  # Create sparse constraint coefficients matrix, in data frame form
  # create needed vectors
  ids <- pufdf$RECID %>% sort
  
  # adding-up constraints
  print("....getting the adding-up constraints for each return (sum of shares allocated to states must equal 1)")
  # xij sum for each person i across states j must be 1 -- this is dense
  ccmat.ids.df <- expand.grid(RECID=ids, 
                              stabbr=sts) %>%
    mutate(stabbr=as.character(stabbr)) %>%
    arrange(RECID, stabbr) %>%
    mutate(cname=paste0("sum.", RECID),
           j=row_number(),
           value=1)
  
  state_constraints <- get_all_state_constraints(cnames, pufdf, sts, ids)
  
  # put the constraints together
  ccmat.df <- bind_rows(ccmat.ids.df, 
                        state_constraints) %>%
    mutate(i=match(cname, unique(cname))) %>%
    select(RECID, cname, i, j, value)
  
  return(ccmat.df)
}

nrecs <- problem$nrecs
pufsumdf <- problem$pufsum

get_rhs <- function(nrecs, pufsumdf, cnames, rhs_shares=NULL){
  #  Create RHS of constraints
  if(is.null(rhs_shares)) rhs.sharetotals <- rep(1, nrecs) else
    rhs.sharetotals <- rhs_shares
  
  # rhsvars <- ifelse(cnames=="wt_one", "wt", cnames)
  rhsvars <- cnames
  
  targets <- c(rhs.sharetotals,
               stack(pufsumdf[, rhsvars])$values)
  return(targets)
}

scale_median <- function(ccmat.df, scale, mdnval=50){
  # scale so that the median of the nonzero constraint coefficients for each constraint is 50
  # get the median value of each raw constraint coefficient
  mdn <- ccmat.df %>% 
    group_by(i) %>% # group by each constraint
    filter(value!=0) %>% # value is the constraint coefficient for each j (variable), for this constraint (i)
    summarise(value=median(value)) %>%
    .[["value"]] # result is a vector
  
  scale <- ifelse(mdn!=0, (1 / mdn) * mdnval, scale)
  return(scale)
}

scale_2norm <- function(ccmat.df, scale){
  # Modeling and Optimization in Space Engineering, edited by Giorgio Fasano, János D. Pintér, pp.51-52
  # scale so that the 2-norm  of each constraint is 1, where 2-norm is sqrt of sum of squared coefficients
  norm_vec <- function(vec) sqrt(sum(vec^2))
  
  norm2 <- ccmat.df %>% 
    group_by(i) %>% # group by each constraint
    filter(value!=0) %>% # value is the constraint coefficient for each j (variable), for this constraint (i)
    summarise(value=norm_vec(value)) %>%
    .[["value"]] # result is a vector
  
  scale <- ifelse(norm2!=0, (1 / norm2), scale)
  return(scale)
}


scale_constraints <- function(ccmat.df, targets, scale_type="none"){
  # define a scale value for each constraint and then apply
  
  scale <- rep(1, length(targets)) # default scaling is 1
  
  if(scale_type=="none"){
    # do nothing
  } else if(scale_type=="median50") {
    scale <- scale_median(ccmat.df, scale)
    
  } else if(scale_type=="2norm") {
    scale <- scale_2norm(ccmat.df, scale)
  }
  
  # print(ht(scale))
  
  targets.scaled <- targets * scale
  
  ccmat.scaled <- ccmat.df %>%
    mutate(value = value * scale[i]) # use i in ccmat.df as index into scale
  
  constraints.scaled <- list()
  constraints.scaled$targets.scaled <- targets.scaled
  constraints.scaled$ccmat.scaled <- ccmat.scaled
  
  return(constraints.scaled)
}


# get_rhs <- function(pufdf, pufsumdf, rhs_shares=NULL){
#   #  Create RHS of constraints
#   if(is.null(rhs_shares)) rhs.sharetotals <- rep(1, nrow(pufdf)) else
#     rhs.sharetotals <- rhs_shares
#   
#   targets <- c(rhs.sharetotals,
#                pufsumdf$wt.sum,
#                pufsumdf$MARS1,
#                pufsumdf$E00100,
#                pufsumdf$E18425)
#   return(targets)
# }

get_problem <- function(pufdf, cnames, rhs_targets, rhs_shares=NULL){
  print("Defining problem list")
  problem <- list()
  problem$sts <- rhs_targets$stabbr
  problem$cnames <- cnames
  problem$nrecs <- nrow(pufdf)
  problem$pufsum <- rhs_targets
  
  print("..getting unscaled RHS targets")
  problem$targets.unscaled <- get_rhs(problem$nrecs, problem$pufsum, problem$cnames, rhs_shares)
  
  print("..getting unscaled constraint coefficients, may take some time")
  problem$ccmat.df.unscaled <- get_ccmatdf(pufdf, problem$cnames, problem$sts)
  
  problem$num.vars <- max(problem$ccmat.df.unscaled$j)
  problem$num.constraints <- max(problem$ccmat.df.unscaled$i)
  problem$nnz <- nrow(problem$ccmat.df.unscaled)
  
  print("..getting indexes of nonzero constraint coefficients, to define sparsity structure")
  problem$ccmat.indexes <- dlply(problem$ccmat.df.unscaled, "i", function(x) return(x$j))
  
  print("Done defining problem")
  cat("\n")
  
  return(problem)
}


get_inputs_list <- function(problem){
  inputs <- list()
  inputs$eps <- .0001
  inputs$xij_zeros <- function(nrows, sts) matrix(data=0, nrow=nrows, ncol=length(sts), byrow=TRUE, dimnames=list(NULL, sts))
  inputs$x_to_xij <- function(x, sts) matrix(data=x, ncol=length(sts), byrow=TRUE, dimnames=list(NULL, sts))
  inputs$xij_to_x <- function(xij) c(t(xij)) # unflattens in the right order
  inputs$sts <- problem$sts
  inputs$nrecs <- problem$nrecs
  inputs$ccmat.df <- problem$constraints.scaled
  inputs$num.vars <- problem$num.vars
  inputs$num.constraints <- problem$num.constraints
  inputs$nnz <- problem$nnz
  inputs$targets <- problem$targets.scaled
  inputs$eval_jac_g_structure <- problem$ccmat.indexes
  inputs$eval_h_structure <- lapply(1:problem$num.vars, function(x) x) # diagonal elements of our Hessian
  return(inputs)
}


#****************************************************************************************************
#                functions for ipoptr ####
#****************************************************************************************************
eval_f_entmax <- function(x, inputs) {
  # objective function - evaluates to a single number
  #   x * log(x) {entmax} entropy maximization functions for ipopt  
  
  # ipoptr requires that ALL functions receive the same arguments, so I pass the inputs list to ALL functions
  
  obj <- sum((x + inputs$eps) * log(x + inputs$eps ))
  return(obj)
}


eval_grad_f_entmax <- function(x, inputs){
  # gradient of objective function - a vector length x 
  # giving the partial derivatives of obj wrt each x[i]
  
  # ipoptr requires that ALL functions receive the same arguments, so I pass the inputs list to ALL functions
  
  # http://www.derivative-calculator.net/
  gradf <- 1 + log(x + inputs$eps)
  return(gradf)
}


eval_h_entmax <- function(x, obj_factor, hessian_lambda, inputs){
  # The Hessian matrix has many zero elements and so we set it up as a sparse matrix
  # We only keep the (potentially) non-zero values that run along the diagonal.
  
  # http://www.derivative-calculator.net/
  hess <- obj_factor * 1 / (x + inputs$eps)
  
  return(hess)
}


#****************************************************************************************************
#                Set up and solve the total problem ####
#****************************************************************************************************
#  !!! DELETE ANY ma77* FILES IN WORKING DIRECTORY BEFORE RUNNING ma77 !!!!----
# when we get to here we should have puf4 and targets

# print(Sys.unsetenv("OMP_NUM_THREADS")) # unset the threads environment variable

Sys.getenv(c("R_HOME", "OMP_NUM_THREADS"))
print(Sys.setenv(OMP_NUM_THREADS = 6))
Sys.getenv(c("R_HOME", "OMP_NUM_THREADS"))

# constraint names
cnames.all <- c("wt_one", "XTOT", "MARS2", "E00100", "E18425")

# c7 <- c(1:5, 7:8)
# c.test <- c(1, 2, 4, 7)

glimpse(puf4)
problem0 <- get_problem(puf4, cnames, rhs_targets=targets)
names(problem0)

inputs0 <- get_inputs_list(problem0, scale_type="median50") # scale_type: none, median50, 2norm

# Ipopt options: https://www.coin-or.org/Ipopt/documentation/node40.html
# ..ma77: https://www.coin-or.org/Ipopt/documentation/node57.html

# review problem size:
inputs0$nrecs
length(inputs0$sts)
inputs0$num.vars
inputs0$num.constraints
inputs0$nnz

# names(inputs0)

# look at values at x0
# x0 <- x0_nstates(inputs0)
# x0 <- x0_wtshare(inputs0)
# x0 <- rep(0, inputs0$num.vars)

# eval_f(x0, inputs0)
# eval_f_lown(x0, inputs0)

# dftarg <- get_constraints_df(inputs0, x0)
# sum(dftarg$diff)
# sum(abs(dftarg$diff))
# glimpse(dftarg)

# obj_scaling_factor: USED to be:  100 is VERY bad, .01 is worse than 1 but better than 100; 1000 better?

opts <- list("print_level" = 5,
             "file_print_level" = 5,
             "jac_c_constant" = "yes",
             "jac_d_constant" = "yes",
             "hessian_constant" = "yes", # if yes Ipopt asks for Hessian of Lagrangian function only once and reuses; default "no"
             # "hessian_approximation" = "limited-memory", # exact is default
             # "tol" = 1e-4, # 1e-8
             "linear_solver" = "ma77", # mumps; pardiso ma27 ma57 ma77 ma86 ma97
             # "nlp_scaling_method" = "gradient-based", # gradient-based; equilibration-based, none
             # "linear_system_scaling" = "mc19", # mc19, none, slack-based
             "obj_scaling_factor" = 1, # 1; -inf, inf; 100-1000 best w/[ma77, amd], so far; 500; oddly, 0.1 not awful
             "mehrotra_algorithm" = "yes", # mehrotra causes substantial reduction in # iterations (for our max entropy problem)
             # "mu_strategy" = "adaptive", # monotone, adaptive; adaptive does NOT improve it on its own
             # "mu_oracle" = "quality-function", # quality-function, probing, loqo (loqo CANNOT USE WITH MEHROTRA)
             # "ma77_maxstore" = 0, # 0; 400e6 and 300e6 make it crash R (I think); 200e6 is ok
             "ma77_order" = "amd",  # metis; amd -- amd seems about 10-15% faster
             # "pardiso_iterative" = "yes",
             # "ma77_print_level" = 2,
             # "ma77_file_size" = 2^22, # 2^21=2,097,152; limits unclear to me; measured in the number of entries (int or double), NOT bytes
             # "ma77_nemin" = 100,  # 8 [1, inf); 4 no impact, 100 no impact
             # "accept_every_trial_step" = "yes", # no impact
             #"derivative_test" = "only-second-order", # none first-order second-order only-second-order
             "max_iter" = 200,
             "output_file" = "full14_low_n77.out")


