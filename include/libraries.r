
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
library(btools) # has a few useful functions -- devtools::install_github("donboyd5/btools")

library(ipoptr)
