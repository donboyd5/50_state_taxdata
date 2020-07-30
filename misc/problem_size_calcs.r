
xvars <- tibble(xvar=c(letters, LETTERS), xvals=c(26:1, 1:26))
system.time(df <- expand_grid(p=1:40e3, s=state.abb, xvars)) # 9 secs
nrow(df) / 1e6 # 104 million records
# calculate constraints 
system.time(tmp <- df %>%
              group_by(s, xvar) %>%
              summarise(x=sum(xvals))) # 16 seconds for 2600
# secs x iter x calcs-per-iter / 60 = # minutes --> 40 minutes if it gets this big
+16 * 30 * 5 / 60