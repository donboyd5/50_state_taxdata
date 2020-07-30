
# TODO ----
# step 7 create new state weights
# step 8 test results


# Problem targets ----
# N09600_nnz and especially A09600_sum bad for most states in AGI_STUB 1
# also bad in other stubs, especially small states


# code folding ----
# alt-o, shift-alt-o
# alt-l, shift-alt-l
# alt-r


# libraries ----
library(magrittr)
library(plyr) # needed for ldply; must be loaded BEFORE dplyr
library(tidyverse)
options(tibble.print_max = 65, tibble.print_min = 65) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats
library(scales)
library(hms) # hms, for times
library(lubridate) # lubridate, for date/times
library(vctrs)

library(purrr)
# library(repurrrsive) ??
library(multidplyr)

library(knitr)
library(btools) # has a few useful functions -- devtools::install_github("donboyd5/btools")
library(ipoptr)

library(microweight) # devtools::install_github("donboyd5/microweight")

devtools::session_info()

# functions ----
source(here::here("puf_tpc_approach", "functions_puf_tpc.r"))
source(here::here("puf_tpc_approach", "functions_puf_ipoptr.r"))
source(here::here("puf_tpc_approach", "functions_ipopt.r"))
source(here::here("puf_tpc_approach", "functions_targets.r"))


# globals ----
dbdir <- "C:/Users/donbo/Dropbox (Personal)/"

# useful cut points for quantiles
probs1_99 <- c(0, .01, .05, .10, .25, .5, .75, .9, .95, .99, 1)

# here is Peter's list of desired targets
desired_target_vars <- c("N1_nnz", "MARS1_nnz", "MARS2_nnz",
                         "A00100_sum", "posAGI_nnz",
                         "A00200_sum", "N00200_nnz",
                         "A01000_sum", "N01000_nnz",
                         "A04470_sum", "N04470_nnz",
                         "A17000_sum", "N17000_nnz",
                         "A04800_sum", "N04800_nnz",
                         "A05800_sum", "N05800_nnz",
                         "A09600_sum", "N09600_nnz",
                         "A00700_sum", "N00700_nnz")


# 1. get targets, examine ----
puf_targ <- read_csv(file = paste0(dbdir, 'puf_2017_targets.csv')) %>%
  select(-X1)
ht(puf_targ)
(available_target_vars <- setdiff(names(puf_targ), c("STATE", "AGI_STUB")))

setdiff(desired_target_vars, available_target_vars)
setdiff(available_target_vars, desired_target_vars)

#.. check the targets ----
# get US totals for each income range and in total, for each target
count(puf_targ, STATE)
ustargtot <- puf_targ %>%
  select(AGI_STUB, all_of(available_target_vars)) %>%
  group_by(AGI_STUB) %>%
  summarise(across(all_of(available_target_vars), function(x) sum(x) / 1e3), .groups = "drop") %>%
  janitor::adorn_totals()

# useful for comparison vs puf values
ustargtot[, 1:7] %>%
  kable(digits=0, format.args = list(big.mark = ','))

# get each state as % of all states
shares <- puf_targ %>%
  select(AGI_STUB, STATE, all_of(available_target_vars)) %>%
  group_by(AGI_STUB) %>%
  mutate(across(all_of(available_target_vars), ~ . / sum(.) * 100))
shares %>% filter(AGI_STUB==2)

# which have a state share of total that is far from their pop share?
shareslong <- puf_targ %>%
  select(AGI_STUB, STATE, all_of(available_target_vars)) %>%
  pivot_longer(-c(AGI_STUB, STATE), names_to="targname", values_to="target") %>%
  group_by(AGI_STUB, targname) %>%
  mutate(target_pct=target / sum(target, na.rm=TRUE) * 100) %>%
  group_by(AGI_STUB, STATE) %>%
  mutate(pop_pct=target_pct[targname=="N1_nnz"],
         diff=target_pct - pop_pct,
         pdiff=diff / pop_pct * 100,
         apdiff=abs(pdiff)) %>%
  ungroup

shareslong %>%
  arrange(AGI_STUB, desc(apdiff)) %>%
  group_by(AGI_STUB) %>%
  filter(row_number() < 5)

# bad stub 10: NH, MS, RI, AR, IA WY
shareslong %>%
  filter(AGI_STUB==10, STATE=="WY")

# SD 2 (9600) WV 3 (9600)
shareslong %>%
  filter(AGI_STUB==2, STATE=="MI")


# worst 20 -- all are in incgroup 1
shareslong %>%
  filter(AGI_STUB == 1) %>%
  filter(STATE == "NY") %>%
  arrange(AGI_STUB, desc(apdiff)) %>%
  filter(row_number() < 20)

# check n returns
mars <- puf_targ %>%
  select(AGI_STUB, STATE, N1_nnz, MARS1_nnz, MARS2_nnz) %>%
  mutate(resid=N1_nnz - MARS1_nnz - MARS2_nnz,
         pct=resid / N1_nnz * 100)

mars %>% arrange(pct) # DC group 1 is negative

quantile(mars$pct)

mars %>% arrange(-pct)


# 2. get puf ----
puf_2017 <- read_csv(file = paste0(dbdir, 'puf_2017_filter.csv'))
glimpse(puf_2017)
ht(puf_2017)

# put an id on the file that allows us to identify original records, then sort
# and put a sortid on the file
puf_2017b <- puf_2017 %>%
  mutate(recid=row_number()) %>%
  arrange(AGI_STUB, recid) %>%
  mutate(sortid=row_number())
glimpse(puf_2017b)


# 3. prepare puf ----
available_target_vars

# get vector of variables for each kind of target, and remove suffix from name
nnz_vars <-remove_suffix("_nnz", available_target_vars) # variables for which we want the weighted number of nonzero values
sum_vars <- remove_suffix("_sum", available_target_vars) # variables for which we want the weighted sum
sumneg_vars <- remove_suffix("_sumneg", available_target_vars) # variables for which we want the weighted sum

pufprep <- puf_2017b %>%
  mutate_at(nnz_vars,
            list(nnz = ~ 1 * (. != 0))) %>%
  mutate_at(sum_vars,
            list(sum = ~ . * (. != 0))) %>%
  mutate_at(sumneg_vars,
            list(sumneg = ~ . * (. < 0)))
pufprep

# prepare a stripped-down version - an important file
pufstrip <- pufprep %>%
  select(sortid, AGI_STUB, s006, contains("_nnz"), contains("_sum"), contains("_sumneg"))
glimpse(pufstrip)

(puf_target_vars <- setdiff(names(pufstrip), c("sortid", "AGI_STUB", "s006"))) # 21 target vars

# get puf totals and compare them to selected target totals, created earlier
puftots <- pufstrip %>%
  group_by(AGI_STUB) %>%
  summarise(across(all_of(puf_target_vars), function(x) sum(x * s006) / 1e3), .groups="drop") %>%
  janitor::adorn_totals()

names(ustargtot)
cols <- 1:8
cols <- c(1, 9:14)
cols <- c(1, 15:21)
(vars <- names(ustargtot)[cols])
puftots[, vars] %>%
  kable(digits=0, format.args = list(big.mark = ','))

ustargtot[, vars] %>%
  kable(digits=0, format.args = list(big.mark = ','))
# good, they are the close or identical
# note that there are some NA values in the targets and the corresponding
# puf values are not always zero -- esp. A04470_sum and A17000_sum


# 4. define targets ----
# pick a subset of possible_target_vars
(possible_target_vars <- intersect(puf_target_vars, available_target_vars)) # 21 targets

# These are the possible targets until I go in and get more
# "N1_nnz"
# "MARS1_nnz", "MARS2_nnz"
# "posAGI_nnz",
# "A00100_sum"
# "N00200_nnz", "A00200_sum" # Wages
# "N00700_nnz", "A00700_sum" # State income tax refunds
# "N01000_nnz", "A01000_sum" # Net capital gain or loss
# "N04470_nnz", "A04470_sum" # Total deductions (standard or itemized)
# "N04800_nnz", "A04800_sum" # Taxable income
# "N05800_nnz", "A05800_sum" # Income tax before credits
# "N09600_nnz", "A09600_sum" # Alternative minimum tax
# "N17000_nnz", "A17000_sum" # Total medical and dental expense deduction

# A09600_sum N09600_nnz have zeroes -- must adjust

# define target_vars_df which will have different target vars for different
# STATE-AGI_STUB combinations
# we need to define fewer targets for income group 1, I think
possible_target_vars

# drop the 5800 and 9600 targets as they give a LOT of problems
default_targets <- c("N1_nnz", "MARS2_nnz",
                     "posAGI_nnz", "A00100_sum",
                     "N00200_nnz", "A00200_sum",
                     "N00700_nnz", "A00700_sum",
                     "N01000_nnz", "A01000_sum",
                     "N04470_nnz", "A04470_sum",
                     "N04800_nnz", "A04800_sum",
                     "N17000_nnz", "A17000_sum")



#.. adjust targets to get rid of zero and NA values ----
# calculate an adjusted target for any zero or NA target so that its per-return
# value is equal to 50% of the lowest nonzero per-return value for other states for that target
# N1_nnz is the return variable
target_vars <- default_targets

good_values <- puf_targ %>%
  select(AGI_STUB, STATE, all_of(target_vars)) %>%
  pivot_longer(-c(AGI_STUB, STATE, N1_nnz)) %>%
  group_by(AGI_STUB, name) %>%
  mutate(allbad=sum(is.na(value))==length(name)) %>%
  ungroup %>%
  filter(!allbad)
# %>% filter(AGI_STUB==1, name=="A17000_sum")filter(AGI_STUB==1, name=="A17000_sum")
# now do per-return adjustments
fraction_of_lowest <- .5
good_values_adj = good_values %>%
  group_by(AGI_STUB, name) %>%
  mutate(per_return=value / N1_nnz,
         isgood=!is.na(per_return) & per_return!=0,
         abs_per_return=ifelse(isgood, abs(per_return), NA),
         min_isgood=min(abs_per_return, na.rm=TRUE),
         whichmin=which.min(abs_per_return),
         min_nzval_per_ret=min_isgood * sign(value[whichmin]),
         value_adj=ifelse(is.na(value) | value==0,
                          N1_nnz * min_nzval_per_ret * fraction_of_lowest,
                          value)) %>%
  ungroup

good_values_adj %>%
  filter(is.na(value) | value==0) %>%
  arrange(AGI_STUB, name, STATE)


good_values_adj %>%
  filter(value != value_adj)

targets_adj <- good_values_adj %>%
  select(AGI_STUB, STATE, N1_nnz, name, value_adj) %>%
  pivot_wider(values_from = value_adj) %>%
  select(AGI_STUB, STATE, all_of(target_vars)) # ensure that columns are in the original order

targets_adj %>% filter(AGI_STUB==10)


# create lists of each kind of target
# targs_default <- list(default_targets)
# targs_ig1 <- list(incgroup1_targets)

# data frame, each AGI_STUB has a vector of target vars
# get bad targets

target_vars_df <- expand_grid(STATE=unique(puf_targ$STATE),
                              AGI_STUB=unique(puf_targ$AGI_STUB)) %>%
  left_join(targets_adj, by = c("STATE", "AGI_STUB")) %>%
  pivot_longer(-c(STATE, AGI_STUB)) %>%
  filter(!is.na(value)) %>%
  group_by(STATE, AGI_STUB) %>%
  summarise(target_vec=list(name), .groups="drop")

df <- target_vars_df %>% filter(STATE=="NY", AGI_STUB==2)
df$target_vec[[1]]


# 5. prepare initial state weights one state at a time ----
# I found in prior investigation that LM is the best method
# notes on parallel:
# https://jennybc.github.io/purrr-tutorial/ls03_map-function-syntax.html
# https://cfss.uchicago.edu/notes/split-apply-combine/

#.. test results for a single state and income group ----
# OA 1 RI 1 DC 1 OA 5 AL 2 NY 1 WY 2
st <- "WY"; ig <- 2 # NH MS RI
st <- "SD"; ig <- 2 #
st <- "WV"; ig <- 3 #
st <- "OA"; ig <- 5 #
st <- "NY"; ig <- 1 #
st <- "KY"; ig <- 10 #

df <- target_vars_df %>% filter(STATE==st, AGI_STUB==ig)
df$target_vec[[1]]

# # SD 2 (9600) WV 3 (9600)
# tvdf2 <- target_vars_df %>%
#   mutate(target_vec=ifelse(STATE==st & AGI_STUB==ig,
#                            list(setdiff(target_vec, "posAGI_nnz")),
#                            target_vec))
# tvdf2 %>% filter(STATE==st & AGI_STUB==ig)
# tvdf2 <- target_vars_df

good_values_adj %>%
  filter(AGI_STUB==ig, STATE==st)

tvdf2 <- df %>%
  mutate(target_vec=list(setdiff(target_vec[[1]], c("posAGI_nnz", "N17000_nnz", "A17000_sum"))))
tvdf2 <- df %>%
  mutate(target_vec=list(setdiff(target_vec[[1]], c("posAGI_nnz"))))
tvdf2$target_vec[[1]]

res <- one_state(st=st, incgroup=ig, target_vars_df=tvdf2, quiet=FALSE, maxiter=50)

# res <- one_state(st=st, incgroup=ig, target_vars_df=tvdf2, method='Newton', quiet=FALSE)
# res <- one_state(st=st, incgroup=ig, target_vars_df=tvdf2, method='Broyden', quiet=FALSE)

names(res)
res$h; res$s; res$k
res$etime
res$solver_message
res$sse_unweighted
res$sse_weighted
# res$targets %>% kable(digits=0, format.args=list(big.mark = ','))
# res$targets_calc %>% kable(digits=0, format.args=list(big.mark = ','))
# res$targets_diff %>% kable(digits=0, format.args=list(big.mark = ','))
# res$targets_pctdiff %>% round(2)
# cbind(res$sortid, res$wh)

# calc vs ALL targets
res$whs
comp <- pufstrip %>%
  filter(AGI_STUB==ig) %>%
  left_join(tibble(STATE=st, sortid=res$sortid, whs=res$whs[, 1]),
            by=c("sortid")) %>%
  pivot_longer(-c(sortid, AGI_STUB, STATE, s006, whs)) %>%
  mutate(calc=whs * value) %>%
  group_by(AGI_STUB, STATE, name) %>%
  summarise(calc=sum(calc, na.rm=TRUE), .groups="drop") %>%
  left_join(puf_targ %>%
              filter(AGI_STUB==ig, STATE==st) %>%
              pivot_longer(-c(AGI_STUB, STATE),
                           values_to="target"),
            by = c("AGI_STUB", "STATE", "name")) %>%
  left_join(good_values_adj %>%
              select(AGI_STUB, STATE, name, target_adj=value_adj),
            by = c("AGI_STUB", "STATE", "name")) %>%
  select(AGI_STUB, STATE, name, target, target_adj, calc) %>%
  mutate(target_adj=ifelse(is.na(target_adj), target, target_adj),
         diff=calc - target,
         pdiff=diff / target * 100,
         diff_adj = calc - target_adj,
         pdiff_adj = diff_adj / target_adj * 100,
         targeted=name %in% tvdf2$target_vec[[1]]) %>%
  arrange(-targeted, -abs(pdiff_adj))
comp



#.. loop through many states and income ranges ----
cluster <- new_cluster(6)

# CAUTION: unfortunately no progress reporting when run in parallel
parallel <- TRUE
# parallel <- FALSE

if(parallel){
  # set the latest versions of functions, etc. up for the run
  cluster_copy(cluster, c('one_state', 'check_xmat', 'puf_targ',
                          'pufstrip', 'target_vars_df')) # functions and data not in a library
  cluster_library(cluster, c("dplyr", "microweight", "tidyr", "purrr"))
}

a <- proc.time()
res_df <- puf_targ %>%
  select(AGI_STUB, STATE) %>%
  # filter(AGI_STUB %in% 2:3, STATE %in% c("AL", "CA", "IL")) %>%
  # filter(AGI_STUB == 1) %>%
  # filter(STATE == "OA") %>%
  group_by(AGI_STUB, STATE) %>%
  nest()  %>%
  {if (parallel) partition(., cluster) else .} %>%
  mutate(res=map(STATE, one_state, incgroup=AGI_STUB, target_vars_df=target_vars_df, quiet=TRUE)) %>%
  {if (parallel) collect(.) else .} %>%
  ungroup %>%
  arrange(AGI_STUB, STATE) # must sort if parallel
b <- proc.time()
b - a # seconds
(b - a) / 60 # minutes 35 minutes

# save results as the above can take a long time to run
# system.time(saveRDS(res_df, here::here("ignore", "single_states.rds"))) # 33 secs
# res_df <- readRDS(here::here("ignore", "single_states.rds"))

names(res_df$res[[1]])

res_df %>%
  unnest_wider(res)

# examine solver messages and sse_weighted
tmp <- res_df %>%
  ungroup %>%
  hoist(res, "h", "solver_message", "sse_weighted")
tmp
count(tmp, h)
count(tmp, solver_message)
count(tmp, AGI_STUB, solver_message)
tmp %>%
  group_by(AGI_STUB) %>%
  summarise(sum(is.na(sse_weighted)))
quantile(tmp$sse_weighted, na.rm=TRUE)

tmp %>% filter(AGI_STUB==10, STATE=="RI")

tmp %>%
  arrange(-sse_weighted)

rec <- res_df %>%
  filter(AGI_STUB==10, STATE=="KY") %>%
  .$res
rec[[1]]$targets_pctdiff %>% round(2)


#.. retrieve original national weights, one set per AGI_STUB ----
wts_wh <- res_df %>%
  select(AGI_STUB, res) %>%
  group_by(AGI_STUB) %>%
  slice(n=1) %>%
  ungroup %>% # ungroup is needed but I don't know why
  hoist(res, "sortid", "wh") %>%
  select(AGI_STUB, sortid, wh) %>%
  unnest(c(sortid, wh))

wts_wh %>%
  group_by(AGI_STUB) %>%
  slice_head(n=4)

#.. retrieve initial state weights and calculate national sum for each record ----
# res_df$res[[1]]$whs[1:10, ] # examine the whs matrix; AGI_STUB 1 is bad
wts_whs <- res_df %>%
  unnest_wider(res) %>%
  select(AGI_STUB, STATE, sortid, whs) %>%
  mutate(whs=map(whs, function(mat) mat[, 1])) %>% # get whs vector for each state
  unnest(c(sortid, whs)) %>%
  pivot_wider(names_from = STATE, values_from=whs)

wts_all <- wts_wh %>%
  left_join(wts_whs, by = c("AGI_STUB", "sortid")) %>%
  mutate(wh_sum=rowSums(across(-c(AGI_STUB, sortid, wh)))) %>%
  select(AGI_STUB, sortid, wh, wh_sum, everything())
glimpse(wts_all)


#.. calculate & compare NATIONAL targets with new initial state weights ----
wts_all %>%
  select(AGI_STUB, sortid, wh, wh_sum) %>%
  mutate(ratio = wh_sum / wh) %>%
  group_by(AGI_STUB) %>%
  slice_head(n=4)

wts_all %>%
  mutate(ratio=wh_sum / wh) %>%
  group_by(AGI_STUB) %>%
  do(qtiledf(.$ratio, probs=probs1_99))

calctargs <- pufstrip %>%
  right_join(wts_all %>% select(AGI_STUB, sortid, wh_sum),
             by = c("sortid", "AGI_STUB")) %>%
  pivot_longer(-c(sortid, AGI_STUB, s006, wh_sum)) %>%
  mutate(target=s006 * value, calc=wh_sum * value) %>%
  group_by(AGI_STUB, name) %>%
  summarise(target=sum(target, na.rm=TRUE),
            calc=sum(calc, na.rm=TRUE),
            .groups="drop") %>%
  mutate(diff=calc - target,
         pdiff=diff / target * 100)
calctargs

calctargs %>%
  arrange(-abs(pdiff))

#.. calculate & compare STATE targets with new initial state weights ----
calc_targ <- calc_targets(wtsdf=wts_all, pufdf=pufstrip)

target_comp <- bind_rows(puf_targ %>% mutate(type="target"),
                         calc_targ %>% mutate(type="calc")) %>%
  pivot_longer(cols=-c(AGI_STUB, STATE, type),
               names_to = "targname") %>%
  pivot_wider(names_from=type) %>%
  mutate(diff=calc - target,
         pdiff=diff / target * 100) %>%
  left_join(target_vars_df %>%
              select(AGI_STUB, STATE, targname=target_vec) %>%
              unnest(targname) %>%
              mutate(targeted=TRUE),
            by = c("AGI_STUB", "STATE", "targname")) %>%
  mutate(targeted=ifelse(is.na(targeted), FALSE, TRUE))

target_comp %>%
  filter(STATE=="CT") %>%
  arrange(-abs(pdiff))

target_comp %>%
  filter(STATE=="NY", AGI_STUB==10) %>%
  arrange(-abs(pdiff))


# 6. reweight the national file to hit national targets ----
# get the national targets
ustargs <- puf_targ %>%
  group_by(AGI_STUB) %>%
  select(-STATE) %>%
  summarise(across(.fns=sum), .groups="drop")
ustargs

# test my reweight function ----
# source(here::here("R", "reweight.r"))

#.. test the reweighting function on a single AGI_STUB ----
d <- one_stub(1)
names(d)
d$solver_message
quantile(d$result$solution)
d$targets_df

d$targets_df %>%
  mutate(across(c(target, targinit, targcalc, targtol, targinit_diff, targcalc_diff),
                function(x) x = ifelse(str_detect(targname, "_sum"), x / 1000, x))) %>%
  kable(digits=c(rep(0, 8), 2, 2, 2), format.args = list(big.mark = ","))


#.. loop through the stubs in serial ----
# CAUTION: DO NOT DO THIS IN PARALLEL with multidplyr
# apparently multidplyr cannot handle all needed functions
# but it is fast enough in serial

a <- proc.time()
opt_df <- tibble(AGI_STUB=1:10) %>%
  group_by(AGI_STUB) %>%
  nest()  %>%
  mutate(output=map(AGI_STUB, one_stub)) %>%
  ungroup
b <- proc.time()
b - a # seconds
(b - a) / 60 # minutes -- < 1 minute in serial

# saveRDS(opt_df, here::here("ignore", "opt_df.rds"))
# opt_df <- readRDS(here::here("ignore", "opt_df.rds"))


names(opt_df$output[[1]])

tmp <- opt_df %>%
  hoist(output, "solver_message")
tmp


# 7. retrieve adjusted national weights and scale initial state weights ----
names(opt_df$output[[1]])
# opt_df$output[[1]]$weights_df # AGI_STUB sortid iweight     x whs_adj
# "solver_message" "etime"          "objective"      "weights"        "targets_df"     "result"

opt_wts <- opt_df %>%
  hoist(output, "weights") %>%
  select(AGI_STUB, weights) %>%
  unnest(col=weights)

wts_adj <- wts_all %>%
  mutate(wh_adj=opt_wts$weights,
         x=wh_adj / wh_sum) %>%
  mutate(across(-c(AGI_STUB, sortid, x, wh, wh_sum, wh_adj),
                function(wt) wt * x)) %>%
  select(-x) %>%
  mutate(wh_adjsum=rowSums(across(all_of(c(state.abb, "DC", "OA"))))) %>% # this is just a check to verify sums
  select(AGI_STUB, sortid, wh, wh_sum, wh_adj, wh_adjsum, everything())

# make sure all looks good
wts_adj %>%
  select(AGI_STUB, sortid, wh, wh_sum, wh_adj, wh_adjsum) %>%
  ht

wts_adj %>%
  mutate(diff=wh_adjsum - wh_adj) %>%
  do(qtiledf(.$diff))

# looks good so drop the check variable
wts_adj <- wts_adj %>%
  select(-wh_adjsum)
wts_adj

#.. calculate state-AGI_STUB targets using scaled initial state weights ----
calc_targ_adj <- calc_targets(wtsdf=wts_adj, pufdf=pufstrip)

target_comp_adj <- bind_rows(puf_targ %>% mutate(type="target"),
                             calc_targ_adj %>% mutate(type="calc")) %>%
  pivot_longer(cols=-c(AGI_STUB, STATE, type),
               names_to = "targname") %>%
  pivot_wider(names_from=type) %>%
  mutate(diff=calc - target,
         pdiff=diff / target * 100) %>%
  left_join(target_vars_df %>%
              select(AGI_STUB, STATE, targname=target_vec) %>%
              unnest(targname) %>%
              mutate(targeted=TRUE),
            by = c("AGI_STUB", "STATE", "targname")) %>%
  mutate(targeted=ifelse(is.na(targeted), FALSE, TRUE))

# pick a state and income group, compare target differences with
# initial state weights and with adjusted initial state weights
st <- "NY"; ig <- 10

# initial first, then adjusted
target_comp %>%
  filter(STATE==st, AGI_STUB==ig) %>%
  arrange(-abs(pdiff))

target_comp_adj %>%
  filter(STATE==st, AGI_STUB==ig) %>%
  arrange(-abs(pdiff))

#.. question: are the adjusted initial state weights good enough or do we need final step? ----
# use adjusted targets from below
check <- good_values_adj %>%
  select(AGI_STUB, STATE, name, target=value_adj) %>%
  left_join(target_comp_adj %>%
              select(AGI_STUB, STATE, name=targname, calc, targeted),
            by = c("AGI_STUB", "STATE", "name")) %>%
  mutate(diff=calc - target,
         pdiff=diff / target * 100) %>%
  ungroup
summary(check)

check %>%
  filter(AGI_STUB==2, name %in% target_vars) %>% # target vars from below
  group_by(name) %>%
  summarise(sse=sum(pdiff^2), .groups="drop") %>%
  arrange(sse) %>%
  mutate(sse_sum=cumsum(sse),
         pct=sse / max(sse_sum) * 100,
         pct_sum=cumsum(pct))


# 8. construct final state weights using poisson method ----
#.. test a single income group ----

count(pufstrip, AGI_STUB)

ig <- 9

# CAUTION: choice for wh below is important: use one of the following as national weights
# create wh:
# - wh for original file weights or
# - whsum for national weights constructed as sum of state weights from the initial pass in step 5
# - wh_adj for the national weights from step 5, adjusted to hit/approximate national targets
wh <- wts_adj %>%
  filter(AGI_STUB==ig) %>%
  .$wh # the new national weight -- use either wh_sum or wh_adj -- or wh for original weight

# target_vars <- default_targets
target_vars <- target_vars_df %>%
  filter(STATE=="AL", AGI_STUB==ig) %>%
  .$target_vec %>%
  unlist
# target_vars <- setdiff(target_vars, "posAGI_nnz") # seems to be a problem

xmat <- pufstrip %>%
  filter(AGI_STUB==ig) %>%
  select(all_of(target_vars)) %>%
  as.matrix

check <- check_xmat(xmat)
check
# print(check)
if(check$remove > 0) {
  print(paste0("WARNING: Removing linearly dependent column: ", target_vars[check$remove]))
  print("CAUTION: Not checking for further linear dependence...")
  xmat <- xmat[, -check$remove]
  target_vars <- target_vars[-check$remove]
}

targets <- puf_targ %>%
  filter(AGI_STUB==ig) %>%
  select(AGI_STUB, STATE, all_of(target_vars))

targets_adj <- good_values_adj %>%
  filter(AGI_STUB==ig) %>%
  select(AGI_STUB, STATE, N1_nnz, name, value_adj) %>%
  pivot_wider(values_from = value_adj) %>%
  select(AGI_STUB, STATE, all_of(target_vars)) # ensure that columns are in the original order

# after reviewing, use targets_adj instead of targets
targmat <- targets %>%
  select(all_of(target_vars)) %>%
  as.matrix
rownames(targmat) <- targets$STATE

targmat_adj <- targets_adj %>%
  select(all_of(target_vars)) %>%
  as.matrix
rownames(targmat_adj) <- targets$STATE

# compare the two
# targmat
# targmat_adj
# targmat_adj - targmat



# # note targmat_adj, also note that matrix was made invertible first
res.lm <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'LM',
                 maxiter=50, opts=list(factor=50)) # ssew 140644.5 38 mins for 10 iterations; ssew  28113, ~120 mins for 30 iterations

res.lm2 <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'LM',
                 betavec=as.vector(res.lm$beta_opt_mat), maxiter=10) # ssew  101986.4 15.5 mins, 4 iterations [on top of prior 10]

res.n <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'Newton',
                   maxiter=75)  # ssew 92065.48, 38 mins with 10 iterations, ssew 84.1 mins with 30 iterations

res.n2 <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'Newton',
                    betavec=as.vector(res.n$beta_opt_mat), maxiter=10) # ssew 108.7 39 mins

res.n3 <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'Newton',
                    betavec=as.vector(res.n2$beta_opt_mat), maxiter=10) # ssew 84.1 40 mins

res.n4 <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'Newton',
                    betavec=as.vector(res.n3$beta_opt_mat), maxiter=10) # 75439 51 mins
res.n5 <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'Newton',
                    betavec=as.vector(res.n4$beta_opt_mat), maxiter=10) # 73567 51 mins
res.n6 <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'Newton',
                    betavec=as.vector(res.n5$beta_opt_mat), maxiter=10) # 72156 51 mins

res.b <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'Broyden'
                   , maxiter=500) # 640176.7,  44 mins

res.n6lm <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'LM',
                     betavec=as.vector(res.n6$beta_opt_mat), maxiter=10) # 69556 51 mins

res.n6lm2 <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'LM',
                      betavec=as.vector(res.n6lm$beta_opt_mat), maxiter=10) # 68932 35.6 mins

res.n6lm2n <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'Newton',
                    betavec=as.vector(res.n6lm2$beta_opt_mat), maxiter=10)
res.n6lm2n <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'Newton',
                        betavec=as.vector(res.n6lm2n$beta_opt_mat), maxiter=10)

res.n6lm2nlm <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'LM',
                        betavec=as.vector(res.n6lm2n$beta_opt_mat), maxiter=10) # 68274.76 51 mins

# saveRDS(res.n6lm2nlm, here::here("ignore", paste0("final_ig", ig, "res.n6lm2nlm.rds")))

# res <- readRDS(here::here("ignore", paste0("final_ig", ig, "res.n6lm2nlm.rds")))

memory()
res <- res.lm
res <- res.n
res$sse_weighted
res$etime / 60
res$targets_pctdiff %>% round(2)

check2 <- as_tibble(res$targets_pctdiff, rownames = "STATE") %>%
  pivot_longer(-STATE) %>%
  mutate(sse=value^2)

# get top differences
check2 %>%
  slice_max(sse, n=30)

# get sse by variable
check2 %>%
  group_by(name) %>%
  summarise(maxapdiff=max(abs(value)),
            sse=sum(sse), .groups="drop") %>%
  arrange(sse) %>%
  mutate(sse_sum=cumsum(sse),
         pct=sse / max(sse_sum) * 100,
         pct_sum=cumsum(pct))

# get sse by state
check2 %>%
  group_by(STATE) %>%
  summarise(maxapdiff=max(abs(value)),
            sse=sum(sse), .groups="drop") %>%
  arrange(sse) %>%
  mutate(sse_sum=cumsum(sse),
         pct=sse / max(sse_sum) * 100,
         pct_sum=cumsum(pct))


check %>%
  filter(AGI_STUB==2, name %in% target_vars) %>% # target vars from below
  group_by(name) %>%
  summarise(sse=sum(pdiff^2), .groups="drop") %>%
  arrange(sse) %>%
  mutate(sse_sum=cumsum(sse),
         pct=sse / max(sse_sum) * 100,
         pct_sum=cumsum(pct))


res.n6$sse_weighted
res.n6$etime / 60
res.n6$targets_pctdiff %>% round(2)



# saveRDS(res.n6lm2, here::here("ignore", paste0("final_ig", ig, ".rds")))
res_old <- readRDS(here::here("ignore", paste0("final_ig", ig, ".rds")))
names(res_old)
names(res_old$output)
res_old$output$rsstrace
res_old$sse_weighted
res_old$etime / 60
# 2.779721e+09 2.583346e+09 2.247818e+09 1.684735e+09 9.052827e+08 2.149220e+08 1.550671e+07 3.344061e+06 6.783272e+05 3.742754e+05 3.288578e+05
# at 97,  9.892193e+04

res$sse_weighted
res$etime
res$targets_pctdiff %>% round(2)
res$output$rsstrace

res.n$sse_weighted
res.n$etime / 60
res.n$targets_pctdiff %>% round(2)

res.b$sse_weighted
res.b$etime
res.b$targets_pctdiff %>% round(2)


whs_adj_pctdiff <- target_comp_adj %>%
  filter(AGI_STUB==ig) %>%
  select(STATE, targname, pdiff) %>%
  pivot_wider(names_from = targname, values_from = pdiff) %>%
  select(STATE, all_of(target_vars))

target_comp_adj %>%
  filter(AGI_STUB==ig) %>%
  filter(abs(pdiff) > 100)

whs_adj_mat <- whs_adj_pctdiff %>%
  select(all_of(target_vars)) %>%
  as.matrix
rownames(whs_adj_mat) <- whs_adj_pctdiff$STATE

res$targets_pctdiff[1:26, 1:10] %>% round(2)
whs_adj_mat[1:26, 1:10] %>% round(2)

res$targets_pctdiff[1:26, 11:20] %>% round(2)
whs_adj_mat[1:26, 11:20] %>% round(2)

whs_adj_mat <- whs_adj_pctdiff %>%
  select(all_of(target_vars)) %>%
  as.matrix

sum(res$targets_pctdiff ^ 2) #  98921
sum(whs_adj_mat ^ 2) # Inf

target_comp_adj %>%
  filter(AGI_STUB==ig) %>%
  filter(!is.infinite(abs(pdiff))) %>%
  summarise(sse=sum(pdiff^2)) # 514866

target_vars

res$h; res$s; res$k
res$opts_used
res$solver_message
res$etime
res$sse_unweighted
res$sse_weighted

#.. is it any easier if we try to establish a better starting point? ----
p <- list()
p$wh <- wh
p$xmat <- xmat
p$targets <- targmat_adj
p$h <- nrow(xmat)
p$s <- nrow(targmat_adj)
p$k <- ncol(xmat)

spoint <- get_starting_point(p) # newton
targs <- targets_mat(spoint$spoint, wh, xmat, s)
(targs - targmat_adj) %>% round(2)

spoint.lm <- get_starting_point(p)
targs <- targets_mat(spoint.lm$spoint, wh, xmat, s=nrow(targmat_adj))
(targs - targmat_adj) %>% round(2)
pdiff <- ((targs - targmat_adj) / targmat_adj * 100)
pdiff %>% round(2)
sum(pdiff^2)
spoint.lm$result$etime

res2 <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'LM', betavec=spoint.lm$spoint, maxiter = 20) # note targmat_adj
res2$solver_message
res2$sse_weighted # did not work well

res3 <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'Broyden', maxiter = 2000) # note targmat_adj
res3$etime
res3$targets_pctdiff %>% round(2)
res3$sse_weighted
as.vector(res3$beta_opt_mat)

res4 <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'Broyden', betavec= as.vector(res3$beta_opt_mat), maxiter = 1000) # note targmat_adj
res4$etime
# saveRDS(res3, here::here("ignore", paste0("final_ig", ig, "_Broyden.rds")))

#.. loop over all income groups -- maybe run overnight ----

final_stub <- function(ig, target_vars, maxiter=50, opts=NULL, weight_name="wh", fnbase="final"){
  # weight_name:
  #   wh original weight from the PUF
  #   wh_sum sum of weights after estimating model for each state-stub combination
  #   wh_adj wh_sum after adjustment to hit national targets
  wh <- wts_adj %>%
    filter(AGI_STUB==ig) %>%
    .[[weight_name]] # the new national weight -- use either wh_sum or wh_adj -- or wh for original weight

  # target_vars <- default_targets
  target_vars <- target_vars_df %>%
    filter(STATE=="AL", AGI_STUB==ig) %>%
    .$target_vec %>%
    unlist

  sortid <- pufstrip %>%
    filter(AGI_STUB==ig) %>%
    .$sortid

  xmat <- pufstrip %>%
    filter(AGI_STUB==ig) %>%
    select(all_of(target_vars)) %>%
    as.matrix

  check <- check_xmat(xmat)

  if(check$remove > 0) {
    print(paste0("WARNING: Removing linearly dependent column: ", target_vars[check$remove]))
    print("CAUTION: Not checking for further linear dependence...")
    xmat <- xmat[, -check$remove]
    target_vars <- target_vars[-check$remove]
  }

  targets <- puf_targ %>%
    filter(AGI_STUB==ig) %>%
    select(AGI_STUB, STATE, all_of(target_vars))

  targets_adj <- good_values_adj %>%
    filter(AGI_STUB==ig) %>%
    select(AGI_STUB, STATE, N1_nnz, name, value_adj) %>%
    pivot_wider(values_from = value_adj) %>%
    select(AGI_STUB, STATE, all_of(target_vars)) # ensure that columns are in the original order

  targmat_adj <- targets_adj %>%
    select(all_of(target_vars)) %>%
    as.matrix
  rownames(targmat_adj) <- targets$STATE

  output <- geoweight(wh = wh, xmat = xmat, targets = targmat_adj, method = 'LM',
                      maxiter=maxiter, opts=list(factor=50))
  output$sortid <- sortid

  fname <- paste0(fnbase, "_", "ig", ig, ".rds")
  saveRDS(output, here::here("ignore", fname))
  return(output)
}




res <- final_stub(1, maxiter=5)
res$result$etime


library(furrr) # parallel purrr
future::plan(multiprocess)
a <- proc.time()
messages <- map(1:2, final_stub, maxiter=5)
b <- proc.time()
(b - a) / 60


a2 <- proc.time()
pmessages <- future_map(1:2, final_stub, maxiter=5, fnbase="pfinal", .progress=TRUE)
b2 <- proc.time()
(b2 - a2) / 60






cluster <- new_cluster(6)

# CAUTION: unfortunately no progress reporting when run in parallel
parallel <- TRUE
# parallel <- FALSE

if(parallel){
  # set the latest versions of functions, etc. up for the run
  # functions and data not in a library
  cluster_copy(cluster, c('final_stub', 'check_xmat', 'good_values_adj',
                          'puf_targ', 'pufstrip', 'target_vars_df', 'wts_adj'))
  cluster_library(cluster, c("dplyr", "microweight", "tidyr", "purrr"))
}

a <- proc.time()
final_df <- puf_targ %>%
  select(AGI_STUB, STATE) %>%
  # filter(AGI_STUB %in% 1:2) %>%
  # filter(AGI_STUB == 1) %>%
  group_by(AGI_STUB) %>%
  nest()  %>%
  {if (parallel) partition(., cluster) else .} %>%
  mutate(res=map(AGI_STUB, final_stub, maxiter=75, fnbase="pfinal")) %>%
  {if (parallel) collect(.) else .} %>%
  ungroup %>%
  arrange(AGI_STUB) # must sort if parallel
b <- proc.time()
b - a # seconds
(b - a) / 60 # 9.3 hours

# save results as the above can take a long time to run
# system.time(saveRDS(final_df, here::here("ignore", "final_all.rds"))) # 33 secs
# final_df <- readRDS(here::here("ignore", "final_all.rds"))
memory()


# 9. check results ----
final_df
names(final_df$res[[1]])

tmp <- final_df %>%
  hoist(res, "solver_message", "sse_weighted")
tmp

tmp <- final_df %>%
  hoist(res, "targets_pctdiff")
tmp$targets_pctdiff[[10]] %>% round(2)

check <- final_df %>%
  hoist(res, "whs") %>%
  select(AGI_STUB, whs)
d <- check[1, ]
d %>% unnest %>% as_tibble
d %>%
  unnest_auto(col=whs)

wts_wh <- final_df %>%
  select(AGI_STUB, res) %>%
  group_by(AGI_STUB) %>%
  slice(n=1) %>%
  ungroup %>% # ungroup is needed but I don't know why
  hoist(res,"wh") %>%
  select(AGI_STUB, wh) %>%
  unnest(c(wh))

wts_wh <- wts_wh %>%
  mutate(sortid=pufstrip$sortid,
         s006=pufstrip$s006) %>%
  select(sortid, AGI_STUB, s006, wh)

wts_wh %>%
  group_by(AGI_STUB) %>%
  slice_head(n=4)

#.. retrieve initial state weights and calculate national sum for each record ----
# final_df$res[[1]]$whs[1:10, ] # examine the whs matrix; AGI_STUB 1 is bad
f <- function(mat){
  df <- as_tibble(mat) %>%
    mutate(rn=row_number()) %>%
    pivot_longer(-rn)
}

wts_whsfinal <- final_df %>%
  hoist(res, "whs") %>%
  select(AGI_STUB, whs) %>%
  mutate(whs=map(whs, as_tibble)) %>%
  unnest(cols=whs) %>%
  mutate(sortid=pufstrip$sortid) %>%
  select(sortid, AGI_STUB, everything())
wts_whsfinal


wts_final <- wts_wh %>%
  left_join(wts_whsfinal, by = c("sortid", "AGI_STUB"))
wts_final
saveRDS(wts_final, here::here("ignore", "wts_final.rds"))
write_csv(wts_final, here::here("ignore", "wts_final.csv"))

# check the sums
check <- wts_final %>%
  mutate(wh_sum=rowSums(across(-c(sortid, AGI_STUB, sortid, s006, wh)))) %>%
  select(sortid, AGI_STUB, s006, wh, wh_sum)
check %>%
  mutate(diff=wh_sum - wh) %>%
  group_by(AGI_STUB) %>%
  do(qtiledf(.$diff))

# check targets
calc_targ <- calc_targets(wtsdf=wts_final, pufdf=pufstrip)

target_comp <- bind_rows(puf_targ %>% mutate(type="target"),
                         calc_targ %>% mutate(type="calc")) %>%
  pivot_longer(cols=-c(AGI_STUB, STATE, type),
               names_to = "targname") %>%
  pivot_wider(names_from=type) %>%
  mutate(diff=calc - target,
         pdiff=diff / target * 100) %>%
  left_join(target_vars_df %>%
              select(AGI_STUB, STATE, targname=target_vec) %>%
              unnest(targname) %>%
              mutate(targeted=TRUE),
            by = c("AGI_STUB", "STATE", "targname")) %>%
  mutate(targeted=ifelse(is.na(targeted), FALSE, TRUE))

target_comp %>%
  filter(STATE=="NY") %>%
  arrange(-abs(pdiff))

target_comp %>%
  filter(STATE=="NY", targeted) %>%
  arrange(-abs(pdiff))

target_comp %>%
  filter(STATE=="NY", AGI_STUB==10) %>%
  arrange(-abs(pdiff))


# BELOW HERE IS APPENDIX ----
# APPENDIX: step 8 final weights using reweighting method ----
# get the initial state weights from step x
# rescale so that they hit the new national weights from step y
# use ipopt with state constraints and adding-up constraints
# let's keep states in alpha order throughout, even OA

# prepare initial weights
glimpse(wts_all2)
iweights3 <- wts_all_natladj %>%
  filter(AGI_STUB==2) %>%
  pivot_longer(-c(AGI_STUB, sortid, wh, wh_sum, wh_sumadj),
               names_to = "STATE",
               values_to = "whs_raw") %>%
  group_by(AGI_STUB, sortid) %>%
  mutate(iweight=whs_raw / sum(whs_raw) * wh_sumadj) %>%
  ungroup %>%
  arrange(AGI_STUB, sortid, STATE) %>%
  mutate(j = row_number())
ht(iweights3)


# prepare targets
target_vars <- default_targets

targagg <- puf_targ %>%
  filter(AGI_STUB==2) %>%
  select(AGI_STUB, STATE, all_of(target_vars)) %>%
  arrange(AGI_STUB, STATE) %>%
  pivot_longer(cols=-c(AGI_STUB, STATE), names_to = "var_calc", values_to = "value") %>%
  mutate(cname=paste0(STATE, "_", var_calc),
         targtype="aggregate") %>%
  separate(var_calc, c('var', 'calc')) %>%
  select(cname, value, AGI_STUB, STATE, var, calc, targtype)
ht(targagg)

targadd <- wts_all_natladj %>%
  filter(AGI_STUB == 2) %>%
  arrange(sortid) %>%
  mutate(cname=paste0("p", str_pad(sortid, width=8, side="left", pad="0")),
         targtype = "addup",
         var=as.character(sortid),
         STATE="ZZ") %>%
  select(cname, value=wh_sumadj, AGI_STUB, STATE, var, targtype)
ht(targadd)

targdf <- bind_rows(targagg, targadd) %>%
  arrange(desc(targtype), cname) %>%
  mutate(i = row_number())
ht(targdf)


# prepare income-group data
glimpse(pufstrip)
igdata <- pufstrip %>%
  filter(AGI_STUB==2)


cc_sparse <- get_cc_states(.incgroup_data = igdata,
                           .target_vars = target_vars,
                           .iweights = iweights3,
                           .targets_df = targdf)

inputs <- get_inputs(.targets_df = targdf,
                     .iweights = iweights3,
                     .cc_sparse = cc_sparse,
                     .targtol=.15, .xub=20, .conscaling=FALSE, scale_goal=1)
names(inputs)
inputs$n_variables; inputs$n_constraints; inputs$n_targets

condf <- targdf %>%
  mutate(tol=ifelse(targtype == "addup", 0, .005),
         tol=ifelse(STATE == "OA", .06, tol),
         tol=ifelse(cname == "OA_A00200_sum", .25, tol)) %>%
  mutate(clb = value - abs(tol * value),
         cub = value + abs(tol * value))
ht(condf)
count(condf, tol)
condf %>% filter(tol > .01)
# A00200_sum

inputs$clb <- condf$clb
inputs$cub <- condf$cub

# For LARGE problems use linear_solver=ma77, obj_scaling=1, and mehrotra_algorithm=yes
opts <- list("print_level" = 0,
             "file_print_level" = 5, # integer
             "max_iter"= 50,
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
             "jac_c_constant" = "yes", # does not improve on moderate problems equality constraints
             "jac_d_constant" = "yes", # does not improve on  moderate problems inequality constraints
             "hessian_constant" = "yes", # KEEP default NO - if yes Ipopt asks for Hessian of Lagrangian function only once and reuses; default "no"
             # "hessian_approximation" = "limited-memory", # KEEP default of exact
             # "limited_memory_update_type" = "bfgs", # default bfgs; sr1 docs say "not working well" but this really helps
             # "derivative_test" = "first-order",
             #"derivative_test_print_all" = "yes",
             "output_file" = here::here("test7.out"))


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
names(result)
quantile(result$solution, probs=c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1))

ht(inputs$cc_sparse)

# old bad targets ----
# get bad targets
targets_bad <- shareslong %>%
  filter(is.na(target) | target==0) %>%
  select(AGI_STUB, STATE, targname) %>%
  group_by(AGI_STUB, STATE) %>%
  summarise(badtargets=list(as.vector(targname)), .groups="drop")

remove_bad <- function(target_vec, badtargets) {
  diff <- setdiff(unlist(target_vec), unlist(badtargets))
  list(diff)
}
target_vars_df <- expand_grid(STATE=unique(puf_targ$STATE),
                              AGI_STUB=unique(puf_targ$AGI_STUB)) %>%
  left_join(targets_bad, by = c("STATE", "AGI_STUB")) %>%
  rowwise() %>%
  mutate(target_vec=case_when(AGI_STUB == 1 ~ list(basic_targets),
                              STATE == "OA" ~ list(basic_targets),
                              TRUE ~ list(default_targets)),
         target_vec=remove_bad(target_vec, badtargets))
target_vars_df %>% filter(STATE=="OA")

