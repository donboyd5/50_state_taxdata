
library(multidplyr)
cluster <- new_cluster(6)
# parallel::stopCluster(cluster)

# conditional execution in a pipeline - for parallel or not
# https://community.rstudio.com/t/conditional-pipelines/6076/2
# mtcars %>% 
#   {if (FALSE) filter(., hp == 245) else .} %>% 
#   {if (FALSE) select(., cyl) else .} %>%
#   head(10)

parallel <- TRUE
# parallel <- FALSE

# djb ----
a <- proc.time()
if(parallel){
  # set the latest versions of functions, etc. up for the run
  cluster_copy(cluster, c("getwts_ipopt", "ipopt_reweight",
                          "define_jac_g_structure_dense", "eval_f_xtop", "eval_grad_f_xtop", "eval_g_dense",
                          "eval_jac_g_dense", "eval_h_xtop",
                          "calc_constraints",
                          "bounds"))
  cluster_library(cluster, c("dplyr", "tidyr", "ipoptr", "readr"))
}
opt_pre <- ccoef %>%
  filter(ftype==ftype_opt) %>%
  # filter(ugroup %in% 6) %>%
  group_by(ugroup) %>%
  {if (parallel) partition(., cluster) else .}

opt <- opt_pre %>%
  do(getwts_ipopt(., all_constraint_vals, bounds, scale=TRUE)) %>% # try map and purrr
  {if (parallel) collect(.) else .} %>%
  ungroup
b <- proc.time()
b - a # seconds
(b - a) / 60 # minutes


https://r.789695.n4.nabble.com/parallel-processing-with-multiple-directories-td4565798.html
library(parallel)
> dirs<- list("out1","out2","out3")   # these directories are located within the current working directory
> temp<- 1:3
> testF<- function(x) {
  >    setwd(x)
  >    saveRDS(temp,"temp.drs")
  >    }
> mclapply(dirs, testF)
>
  > Any help would be appreciated!
  Suppose that jobs 1 and 3 are executed on processor 2. Then after the
first iteration the directory is ./out1 and on the second iteration
setwd() tries to change to ./out1/out3. Full path names might help.

  

