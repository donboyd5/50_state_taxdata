

# library(parDist)
library(parallelDist)

mydata<-as.data.frame(matrix(c(1,1,1,1,0,1,1,1,1,0),nrow=2))
mydata
ref<-rep(1,5)
apply(mydata,1,function(x)sqrt(sum((x-ref)^2)))

d1v <- bds[2, target_vars] %>% as.numeric
d2m <- tws[, target_vars] %>% as.matrix

# p1 az 242.209465
# p1 ak 7.737976
# p2 az

bds2 <- bds[, target_vars]

dfun <- function(vec, mat){
  # vec <- as
  # sqrt(sum(vec - mat)^2)
  rowSums((vec - mat)^2, na.rm=TRUE)
}
dfun(d1v, d2m)
# person 1 and some states
# dfun(d1v, d2m)
# [1] 61.387372  9.336565 20.967724  9.628972 10.065193  7.528750 16.681597 26.565868 13.421510  5.142451
# person 2
# [1] 123.722728  75.881051  76.328266  65.798952  67.801974  47.214718  96.879741 130.435355

tmp <- apply(bds2, 1, dfun, d2m)
str(tmp)
tmp[1:10, 1:2]
# rows are states, cols are people

probs1 <- tmp[ , 1] / sum(tmp[, 1])
sum(probs1)




target_vars

summary(tws)


glimpse(bds)
summary(bds)


get_dist <- function(targets_wide, base_data, target_vars, .popvar){
  # targets: get per-return averages for each state and scale to 0, 1
  targets_scaled <- targets_wide %>%
    mutate(across(.cols=(all_of(target_vars)), ~ .x / {{.popvar}})) %>%
    mutate(across(.cols=(all_of(target_vars)), ~ scale(.x)[, 1]))
  
  # base data: scale to 0, 1
  base_scaled <- base_data %>%
    mutate(across(.cols=(all_of(target_vars)), ~ scale(.x)[, 1]))
  
  # function to compute euclidean distance between each row of scaled base data and every row of the state targets
  edist <- function(vec, mat){
    rowSums((vec - mat)^2, na.rm=TRUE)
  }
  
  # get a matrix of distances, where every row is a state and every column is a person
  mdists <- apply(base_scaled[, target_vars] %>% as.matrix,
                  1,
                  edist,
                  targets_scaled[, target_vars] %>% as.matrix) %>%
    t() # transpose so that states are columns
  
  colnames(mdists) <- targets_scaled$stabbr
  distdf <- cbind(pid=base_data$pid, mdists) %>% as_tibble()
  distdf
}


d <- get_dist(targets_wide, base_data, target_vars, .popvar=pop_nnz)
d[1:10, 1:15]
glimpse(d)

# now we can compute an initial weight for each person-state combo
# pid, stabbr, incgroup, weight_total, iweight_state
summary_vals
stshares <- targets_wide %>%
  select(stabbr, pop_nnz) %>%
  mutate(stshare=pop_nnz / sum(pop_nnz)) %>%
  select(stabbr, stshare)


iweights2 <- base_data %>%
  select(pid, incgroup, weight_total) %>%
  left_join(d, by = "pid") %>%
  pivot_longer(cols=-c(pid, incgroup, weight_total), names_to = "stabbr", values_to = "dist") %>%
  left_join(stshares, by = "stabbr") %>%
  group_by(pid) %>%
  mutate(wdist=dist * stshare,
         prob=wdist / sum(wdist),
         iweight_state=prob * weight_total) %>%
  ungroup %>%
  arrange(pid, stabbr)

iweights3 <- base_data %>%
  select(pid, incgroup, weight_total) %>%
  left_join(d, by = "pid") %>%
  pivot_longer(cols=-c(pid, incgroup, weight_total), names_to = "stabbr", values_to = "dist") %>%
  left_join(stshares, by = "stabbr") %>%
  group_by(pid) %>%
  mutate(wdist=dist * stshare,
         prob=dist / sum(dist),
         iweight_state=prob * weight_total) %>%
  ungroup %>%
  arrange(pid, stabbr)

