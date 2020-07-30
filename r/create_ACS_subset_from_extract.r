
# CAUTION: carefully set number of records to get, and set file name accordingly
# create files with:
#    10k records and 5 states
#   200k records and 50 states



# libraries ----
source(here::here("include", "libraries.r"))

# globals ----
dbox <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/"

# ONETIME: get ACS data and create subset ----
system.time(acs <- readRDS(paste0(dbox, "data/", "acs_persons.rds")))
glimpse(acs)
ns(acs)
summary(acs)

#.. define desired states and number of records, then create extract
nrecs <- c(10e3, 100e3, 200e3, 400e3) # total number of records to get (not # per state)
sts5 <- c("CA", "NY", "TX", "IL", "FL") # 5 large very different states
sts <- list(sts5, state.abb[1:20], state.abb, state.abb)
(fns <- paste0(c("acs_10krecs_5states", "acs_100krecs_20states", "acs_200krecs_50states", "acs_400krecs_50states"), ".rds"))


group <- 4
set.seed(1234)
samp1 <- acs %>%
  filter(stabbr %in% sts[[group]], !is.na(pincp), agep >= 18) %>%
  sample_n(nrecs[group]) %>%
  select(-st)
summary(samp1)

#.. clean the sample data slightly and save 
samp2 <- samp1 %>%
  mutate(otherincp=pincp - (wagp + intp + pap + retp + ssip + ssp),
         nrecs=1 / pwgtp, # the sampling ratio -- when multiplied by weight, will yield 1
         pop=1) # useful to have a variable that has a value of one for all records

glimpse(samp2)
summary(samp2)
count(samp2, stabbr)
# save this here
# saveRDS(samp2, here::here("data", ))
saveRDS(samp2, here::here("data", fns[group]))

# samp2 %>%
#   group_by(incgroup) %>%
#   summarise(imin=min(pincp), imax=max(pincp)) %>%
#   mutate(ilow=ifelse(row_number()==1, imin, (imin + lag(imax)) / 2),
#          ihigh=ifelse(row_number()==max(row_number()), imax, (imax + lead(imin)) / 2),
#          group=comma(ihigh))
