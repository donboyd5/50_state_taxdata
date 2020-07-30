
check_targets <- function(x, inputs){
  inputs$targets_df %>%
    mutate(calctarg_scaled = eval_f_sse_elements(x, inputs),
           calctarg=calctarg_scaled * scale,
           target_unscaled = target * scale,
           diff=calctarg - target_unscaled,
           pdiff=diff / target_unscaled * 100)
}


check_weights <- function(x, inputs){
  wcheck <- tibble(i=1:length(inputs$constraints), weight_total=inputs$constraints) %>%
    mutate(weight=eval_g_addup(x, inputs))
  wcheck
}


sse_weights <- function(wcheck){
  wcheck %>%
    summarise(sse=sum((weight - weight_total)^2))
}
