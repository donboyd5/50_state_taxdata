
# code folding ----
# alt-o, shift-alt-o
# alt-l, shift-alt-l
# alt-r

# notes ----

# Create extract from the 2017 5-year ACS that has selected income and demographic information for all
# person records (~15 million).

# Subsets of this file can be used to test different approaches to weighting state income tax files.


# libraries ----
library(magrittr)
library(plyr) # needed for ldply; must be loaded BEFORE dplyr
library(tidyverse)
options(tibble.print_max = 65, tibble.print_min = 65) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats

library(readxl)
# library(writexl) # NO
# library(openxlsx)

library(scales)
library(hms) # hms, for times
library(lubridate) # lubridate, for date/times
library(vctrs)

library(grDevices)
library(knitr)
library(kableExtra)

library(janitor)
library(btools)

# globals ----
dbox <- "C:/Users/donbo/Dropbox (Personal)/50state_taxdata/"



#****************************************************************************************************
#               ONETIME: get ACS person records with selected info and save ####
#****************************************************************************************************
# library("DBI")

# define locations
dbdir <- "D:/Data/CensusACS/20175year/RSQLITE/"
hcsvdir <- "D:/Data/CensusACS/20175year/csv_hus/"
hfnbase <- "psam_hus"

dbf <- paste0(dbdir, "acs.sqlite")

#.. get ENTIRE population so that we can match them with housholds!
acsdb <- DBI::dbConnect(RSQLite::SQLite(), dbf)
tname <- "acs2017_5year"
DBI::dbListTables(acsdb)
DBI::dbListFields(acsdb, tname)
getall <- tbl(acsdb, tname) # dplyr does lazy loading, so this does not really grab full table
str(getall)
glimpse(getall)

# intp interest income pap public assistance income retp retirement income
persons <- getall %>%
  select(serialno, sporder, st, pwgtp, adjinc, mar, sex, agep, pincp, wagp, intp, pap, retp, ssip, ssp)
# ht(persons) # tail not supported
# DON'T USE glimpse in this situation - it takes a long time
system.time(persons2 <- collect(persons, n=Inf)) # ~ 1 min
glimpse(persons2)
count(persons2, sporder)

# data(package="datasets")
codes <- unique(tigris::fips_codes %>% select(state, state_code) %>% mutate(state_code=as.numeric(state_code)))
persons3 <- persons2 %>%
  mutate(stabbr=codes$state[match(st, codes$state_code)],
         adjinc=adjinc / 1e6) %>%
  select(serialno, sporder, st, stabbr, pwgtp, adjinc, everything()) %>%
  mutate_at(vars(intp, pap, retp, ssip, ssp, wagp), ~ . * adjinc)
count(persons3, st, stabbr) # 1 min, 56 max
glimpse(persons3)
summary(persons3)

saveRDS(persons3, paste0(dbox, "data/", "acs_persons.rds"))

# rm(persons, persons2, persons3, codes)
DBI::dbDisconnect(acsdb)
