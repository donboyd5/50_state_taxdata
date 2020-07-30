

nl.jacobian(inputs2$x0, heq, inputs=inputs2)

eval_grad_f_wfs(x=inputs2$x0, inputs=inputs2)
# x is 5000, evgradf is 5000 so jac of this wb 5k x 5k


hess <- nl.jacobian(inputs2$x0, eval_grad_f_wfs, inputs=inputs2)

set.seed(1234)
x2 <- runif(length(inputs2$x0))
system.time(hess2 <- nl.jacobian(x2, eval_grad_f_wfs, inputs=inputs2))

set.seed(5678)
x3 <- runif(length(inputs2$x0))
            system.time(hess3 <- nl.jacobian(x3, eval_grad_f_wfs, inputs=inputs2))


hess[1:10, 1:12]
hess2[1:10, 1:12]
hess3[1:10, 1:12]
sum(hess==0) / length(hess) * 100



library(Deriv)
obj <- expression((x1*w1*a1 + x2*w2*a2 + x3*w3*a3 - t1)^2 +
                    (x1*w1*b1 + x2*w2*b2 + x3*w3*b3 - t2)^2)
Deriv(obj, "x1", nderiv=2)

gx1 <- Deriv(obj, "x1")
Deriv(gx1, "x1")
Deriv(gx1, "x2")
Deriv(gx1, "x3")

names(inputs2)
inputs2$cc_sparse

incgroup_data
target_vars
iweights
targets_df2

hessdf <- iweights %>%
  select(pid, stabbr, iweight_state) %>%
  left_join(incgroup_data %>% select(pid, all_of(target_vars)), by = "pid") %>%
  mutate_at(vars(all_of(target_vars)), list(~ . * iweight_state)) %>% # constraint coefficients
  arrange(pid, stabbr) %>%
  mutate(j=row_number()) %>% # j is an index for x, the variables we will solve for
  select(j, pid, stabbr, iweight_state, all_of(.target_vars))


cc_dense1 <- .iweights %>%
  select(pid, stabbr, iweight_state) %>%
  left_join(.incgroup_data %>% select(pid, all_of(.target_vars)), by = "pid") %>%
  mutate_at(vars(all_of(.target_vars)), list(~ . * iweight_state)) %>% # constraint coefficients
  arrange(pid, stabbr) %>%
  mutate(j=row_number()) %>% # j is an index for x, the variables we will solve for
  select(j, pid, stabbr, iweight_state, all_of(.target_vars))

#.. sparse matrix of constraint coefficients, not including adding-up constraints
cc_sparse1 <- cc_dense1 %>%
  select(j, pid, stabbr, all_of(.target_vars)) %>%
  pivot_longer(cols = all_of(.target_vars),
               names_to="vname_calctype",
               values_to = "nzcc") %>% # nzcc for nonzero constraint coefficient
  filter(nzcc!=0) %>%
  mutate(cname=paste0(vname_calctype, "_", stabbr)) %>%
  select(j, nzcc, pid, stabbr, cname)


