
check_constraints <- function(.inputs){
  df <- .inputs$cc_sparse %>%
    group_by(i, targtype, cname) %>%
    summarise(conx0=sum(nzcc * .inputs$x0[j])) %>%
    ungroup %>%
    mutate(target=.inputs$constraints, clb=.inputs$clb, cub=.inputs$cub) %>%
    select(i, targtype, cname, clb, target, cub, conx0) %>%
    separate(cname, c("vname", "ftype", "stabbr"), remove=FALSE, fill="right") %>%
    mutate(diff=conx0 - target, pdiff=diff / target * 100)
  df
}
