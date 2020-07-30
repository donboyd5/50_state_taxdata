

library(magrittr)
library(plyr) # needed for ldply; must be loaded BEFORE dplyr
library(tidyverse)
options(tibble.print_max = 65, tibble.print_min = 65) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats
library(scales)
library(hms) # hms, for times
library(lubridate) # lubridate, for date/times
library(vctrs)

library(btools) # has a few useful functions -- devtools::install_github("donboyd5/btools")

library(microweight)

devtools::session_info()


dbdir <- "C:/Users/donbo/Dropbox (Personal)/"
puf_2017 <- read_csv(file = paste0(dbdir, 'puf_2017_filter.csv'))
ht(puf_2017)

puf_targ <- read_csv(file = paste0(dbdir, 'puf_2017_targets.csv')) %>%
  select(-X1)
ht(puf_targ)
names(puf_targ)

puf_targ %>%
  select(AGI_STUB, STATE, N1_nnz, MARS1_nnz, MARS2_nnz) %>%
  mutate(mars12=MARS1_nnz + MARS2_nnz, mdiff=N1_nnz - mars12,
         mdpct=mdiff / N1_nnz * 100,
         amdpct=abs(mdpct)) %>%
  arrange(mdpct)

possible_target_vars <- c("N1_nnz", "MARS1_nnz", "MARS2_nnz", "A00100_sum", "pos_AGI_nnz", "A00200_sum", "N00200_nnz", "A01000_sum",
                          "N01000_nnz", "A04470_sum", "N04470_nnz", "A17000_sum", "N17000_nnz", "A04800_sum", "N04800_nnz",
                          "A05800_sum", "N05800_nnz", "A09600_sum", "N09600_nnz",
                          "A00700_sum", "N00700_nnz")
feasible_target_vars <- intersect(possible_target_vars, names(puf_targ))
setdiff(possible_target_vars, names(puf_targ))


# check the targets
checkwide <- puf_targ %>%
  select(AGI_STUB, STATE, all_of(feasible_target_vars)) %>%
  group_by(AGI_STUB) %>%
  mutate(across(all_of(feasible_target_vars), ~ . / sum(.) * 100))
checkwide %>% filter(AGI_STUB==2)

check <- puf_targ %>%
  select(AGI_STUB, STATE, all_of(feasible_target_vars)) %>%
  group_by(AGI_STUB) %>%
  mutate(across(all_of(feasible_target_vars), ~ . / sum(.) * 100)) %>%
  pivot_longer(-c(AGI_STUB, STATE), names_to="target", values_to="state_pct")
ht(check)
summary(check)

check %>%
  group_by(AGI_STUB, STATE) %>%
  mutate(pop_pct=state_pct[target=="N1_nnz"],
         pdiff=state_pct / pop_pct * 100 - 100,
         apdiff=abs(pdiff)) %>%
  ungroup %>%
  arrange(desc(apdiff))



# make mypuf ----
puf_targ_djb <- puf_targ %>%
  select(AGI_STUB, STATE, all_of(feasible_target_vars)) %>%
  select(-c(A04470_sum, N04470_nnz,
            A17000_sum, N17000_nnz,
            A04800_sum, N04800_nnz))
summary(puf_targ_djb) # no NA's
count(puf_targ_djb, STATE)

# compute new state targets, leaving out the current OA
# we need national sums to redistribute across states
natsums <- puf_targ_djb %>%
  select(-STATE) %>%
  group_by(AGI_STUB) %>%
  summarise(across(everything(), sum), .groups="drop")

puf_targ_revised <- puf_targ_djb %>%
  filter(STATE != "OA") %>%
  pivot_longer(-c(AGI_STUB, STATE)) %>%
  left_join(natsums %>% pivot_longer(-AGI_STUB, values_to="value_US"), by = c("AGI_STUB", "name")) %>%
  group_by(AGI_STUB, name) %>%
  mutate(target_revised=value / sum(value) * value_US) %>%
  select(AGI_STUB, STATE, name, value=target_revised) %>%
  pivot_wider() %>%
  ungroup
summary(puf_targ_revised)
count(puf_targ_revised, STATE)

puf_targ_revised %>%
  group_by(AGI_STUB) %>%
  mutate(across(-c(AGI_STUB, STATE), function(x) x / sum(x) * 100)) %>%
  filter(AGI_STUB==2)

# now revise the puf, and then check wtd sums against target values
v_target_vars <- setdiff(names(puf_targ_revised), c("AGI_STUB", "STATE"))

ivars <- c(1, 4, 5, 6)
(target_vars <- v_target_vars[ivars])
target_vars <- setdiff(v_target_vars, c("MARS1_nnz", "MARS2_nnz", "A09600_sum", "N09600_nnz")) # A09600_sum N09600_nnz
target_vars <- setdiff(v_target_vars, c("MARS1_nnz", "MARS2_nnz")) # A09600_sum N09600_nnz
# A09600_sum N09600_nnz have zeroes -- must adjust


getbase <- function(suffix) {
  var_backend <- str_extract(target_vars, "_.*$") # gets the ending part of each variable name
  base_vars <- target_vars[which(var_backend==suffix)] %>% str_remove(suffix)
  names(base_vars) <- base_vars
  base_vars
}

nnz_vars <-getbase("_nnz") # variables for which we want the weighted number of nonzero values
sum_vars <- getbase("_sum") # variables for which we want the weighted sum
sumneg_vars <- getbase("_sumneg") # variables for which we want the weighted sum

pufdjb <- puf_2017 %>%
  mutate_at(nnz_vars,
            list(nnz = ~ 1 * (. != 0))) %>%
  mutate_at(sum_vars,
            list(sum = ~ . * (. != 0))) %>%
  mutate_at(sumneg_vars,
            list(sumneg = ~ . * (. < 0)))

pufdjb %>%
  select(AGI_STUB, s006, all_of(target_vars)) %>%
  group_by(AGI_STUB) %>%
  summarise(across(all_of(target_vars), ~ sum(.x * s006)))

puf_targ_revised %>%
  select(AGI_STUB, all_of(target_vars)) %>%
  group_by(AGI_STUB) %>%
  summarise(across(all_of(target_vars), ~ sum(.x)))

incstub <- 2
targs <- puf_targ_revised %>%
  filter(AGI_STUB==incstub) %>%
  select(all_of(target_vars)) %>%
  as.matrix

# djb alternative ----
# mutate(across(all_of(target_vars), ~ ifelse(.==0, quantile(., .25), .)))
targs1 <- puf_targ_revised %>%
  filter(AGI_STUB==incstub) %>%
  select(all_of(target_vars))
popshare <- targs1$N1_nnz / sum(targs1$N1_nnz)
targs <- targs1 %>%
  mutate(across(all_of(target_vars), ~ ifelse(.==0, sum(.) * popshare * .1, .))) %>%
  as.matrix

xmat <- pufdjb %>%
  filter(AGI_STUB==incstub) %>%
  select(all_of(target_vars)) %>%
  as.matrix

wh <- pufdjb %>%
  filter(AGI_STUB==incstub) %>%
  .$s006

p <- list()
p$s <- nrow(targs)
p$k <- ncol(targs)
p$h <- length(wh)
p$wh <- wh
p$targets <- targs
p$xmat <- xmat

puf_broy <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Broyden')
puf_broy$etime
puf_broy$targets_pctdiff %>% round(3)

puf_newt <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Newton')
puf_newt$etime
puf_newt$solver_message
puf_newt$targets
puf_newt$targets_diff %>% round(3)
puf_newt$targets_pctdiff %>% round(3)
puf_newt$beta_opt_mat %>% round(3)
quantile(puf_newt$whs, 0:10/10)

puf_lm <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'LM')
puf_lm$etime
puf_lm$solver_message
puf_lm$sse_weighted
puf_lm$sse_unweighted
puf_lm$targets_diff %>% round(3)
puf_lm$targets_pctdiff %>% round(3)


spoint <- get_starting_point(p)
spoint$result$solver_message
spoint$result$etime
spoint$result$targets_pctdiff %>% round(3)


puf_broy_2step <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Broyden', betavec=spoint$spoint)
puf_broy_2step$etime
puf_broy_2step$solver_message
puf_broy_2step$targets_pctdiff %>% round(3)


puf_newt_2step <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Newton', betavec=spoint$spoint)
puf_newt_2step$etime
puf_newt_2step$solver_message
puf_newt_2step$targets
puf_newt_2step$targets_diff %>% round(3)
puf_newt_2step$targets_pctdiff %>% round(3)

puf_lm$etime
puf_broy$etime
puf_newt$etime

puf_newt_2step$etime
puf_broy_2step$etime

puf_lm$sse_unweighted
puf_newt_2step$sse_unweighted
puf_newt$sse_unweighted
puf_broy_2step$sse_unweighted

puf_lm$sse_weighted
puf_newt_2step$sse_weighted
puf_newt$sse_weighted
puf_broy_2step$sse_weighted


# do the individual-state weights following TPC pp.6-7 method ----
getbase <- function(suffix) {
  var_backend <- str_extract(feasible_target_vars, "_.*$") # gets the ending part of each variable name
  base_vars <- feasible_target_vars[which(var_backend==suffix)] %>% str_remove(suffix)
  names(base_vars) <- base_vars
  base_vars
}

nnz_vars <-getbase("_nnz") # variables for which we want the weighted number of nonzero values
sum_vars <- getbase("_sum") # variables for which we want the weighted sum
sumneg_vars <- getbase("_sumneg") # variables for which we want the weighted sum

pufdjb <- puf_2017 %>%
  mutate_at(nnz_vars,
            list(nnz = ~ 1 * (. != 0))) %>%
  mutate_at(sum_vars,
            list(sum = ~ . * (. != 0))) %>%
  mutate_at(sumneg_vars,
            list(sumneg = ~ . * (. < 0)))


puf_targ
tvars <- setdiff(names(puf_targ), c("AGI_STUB", "STATE"))

possible_target_vars <- c("N1_nnz", "MARS1_nnz", "MARS2_nnz", "A00100_sum", "pos_AGI_nnz", "A00200_sum", "N00200_nnz", "A01000_sum",
                          "N01000_nnz", "A04470_sum", "N04470_nnz", "A17000_sum", "N17000_nnz", "A04800_sum", "N04800_nnz",
                          "A05800_sum", "N05800_nnz", "A09600_sum", "N09600_nnz",
                          "A00700_sum", "N00700_nnz")
feasible_target_vars <- intersect(possible_target_vars, names(puf_targ))
setdiff(possible_target_vars, names(puf_targ))

v_target_vars <- setdiff(names(puf_targ), c("AGI_STUB", "STATE"))

incstub <- 2

st <- "NY"

puf_targst <- puf_targ %>%
  mutate(STATE=ifelse(STATE==st, STATE, "OX")) %>%
  group_by(AGI_STUB, STATE) %>%
  summarise(across(all_of(tvars), sum), .groups="drop")

v_target_vars

ivars <- c(1, 4, 6)
ivars <- c(1:4, 6:21) # 21 is top
(target_vars <- v_target_vars[ivars])
# target_vars <- setdiff(v_target_vars, c("MARS1_nnz", "MARS2_nnz", "A09600_sum", "N09600_nnz")) # A09600_sum N09600_nnz
# target_vars <- setdiff(v_target_vars, c("MARS1_nnz", "MARS2_nnz")) # A09600_sum N09600_nnz
# A09600_sum N09600_nnz have zeroes -- must adjust

targs <- puf_targst %>%
  filter(AGI_STUB==incstub) %>%
  select(all_of(target_vars)) %>%
  as.matrix
targs

xmat <- pufdjb %>%
  filter(AGI_STUB==incstub) %>%
  select(all_of(target_vars)) %>%
  as.matrix

wh <- pufdjb %>%
  filter(AGI_STUB==incstub) %>%
  .$s006

p <- list()
p$s <- nrow(targs)
p$k <- ncol(targs)
p$h <- length(wh)
p$wh <- wh
p$targets <- targs
p$xmat <- xmat

puf_lm <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'LM')
puf_lm$etime
puf_lm$solver_message
puf_lm$targets_pctdiff %>% round(3)
puf_lm$sse_weighted
puf_lm$sse_unweighted

puf_newt <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Newton')
puf_newt$etime
puf_newt$solver_message
puf_newt$targets_pctdiff %>% round(3)
puf_newt$sse_weighted

puf_broy <- geoweight(wh=p$wh, xmat=p$xmat, targets=p$targets, method = 'Broyden')
puf_broy$etime
puf_broy$solver_message
puf_broy$targets_pctdiff %>% round(3)
puf_broy$sse_weighted



# looks like LM might be the best for this; now loop through states and get weights ----

