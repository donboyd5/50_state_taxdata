
#*****************************************************************************************
#       objective function, helper function, and gradient for sse approach            ####
#*****************************************************************************************

eval_f_sse <- function(x, inputs){
  target_vec <- eval_f_sse_elements(x, inputs)
  target_diff <- target_vec - inputs$targets
  obj <- sum(target_diff^2) / inputs$objscale # * inputs$priority_weights)
  return(obj)
}


eval_grad_f_sse <- function(x, inputs){
  # http://www.derivative-calculator.net/
  # https://www.symbolab.com/solver/partial-derivative-calculator
  # Example where we only have one target element in the objective function, and only have
  #   3 records and therefore 3 weights, and need to return a vector of 3 partial derivatives
  # Notation: w1, w2, and w3 are the first 3 weights in the vector of weights
  #           a1, a2, and a3 are the 3 constants they are multiplied by (e.g., wages1, wages2, wages3)
  #           t is the sum or target we are comparing to
  #           p is the priority weight of this element of the objective function
  #           calc is the calculated value of the target with the new weights
  # The objective function to minimize is p*(w1*a1 + w2*a2 + w3*a3 - s)^2
  # The first partial derivatives of the objective function are:
  #   wrt w1:          2*a1*p*(w1*a1 + w2*a2 + w3*a3 - t) = 2*a1*p*(calc - t)
  #   wrt w2:          2*a2*p*(w1*a1 + w2*a2 + w3*a3 - t) = 2*a2*p*(calc - t)
  #   wrt w3:          2*a3*p*(w1*a1 + w2*a2 + w3*a3 - t) = 2*a3*p*(calc - t)
  # If we have multiple target elements, each vector element will have a more complex sum
  
  target_vec <- eval_f_sse_elements(x, inputs)
  target_diff <- target_vec - inputs$targets
  
  grad <- inputs$ofe_sparse %>%
    mutate(iweight=inputs$iweight[j],
           xdf=x[j],
           diff=target_diff[targnum],
           grad= 2 * nz_coef * xdf * diff) %>%
    group_by(j) %>%
    summarise(grad=sum(grad), .groups="drop")
  
  return(grad$grad / (inputs$objscale))
}


eval_f_sse_elements <- function(x, inputs) {
  # targets we would like to have hold in the solution
  # return a vector where each element evaluates a target
  target_tbl <- inputs$ofe_sparse %>%
    group_by(targnum) %>%
    summarise(target_value=sum(nz_coef * x[j]),
              .groups="keep")
  # the column constraint_value is a vector, in the order we want
  
  return(target_tbl$target_value)
}


#*****************************************************************************************
#   constraint evaluation and coefficient functions for ipoptr SPARSE -- sse approach ####
#*****************************************************************************************

eval_g_addup <- function(x, inputs){
  # return the vector of constraints (the sum of state weights) at point x
  tot_weights <- inputs$cc_sparse %>%
    mutate(weight_state=iweight_state * x[j]) %>%
    group_by(i) %>%
    summarise(weight_state=sum(weight_state), .groups = "drop")
  tot_weights$weight_state
}


eval_jac_g_addup <- function(x, inputs){
  # the Jacobian is the matrix of first partial derivatives of constraints (these derivatives may be constants)
  # this function evaluates the Jacobian at point x
  
  # return: a vector where each element gives a NONZERO partial derivative of constraints wrt change in x
  # so that the first m items are the derivs with respect to each element of first column of ccmat
  # and next m items are derivs with respect to 2nd column of ccmat, and so on
  # so that it returns a vector with length=nrows x ncolumns in ccmat
  
  # because constraints in this problem are linear, the derivatives are all constants
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  
  return(inputs$cc_sparse$iweight_state)
}


define_jac_g_structure_sparse <- function(cc_sparse, ivar="i", jvar="j"){
  # the jacobian 
  # return a list that defines the non-zero structure of the "virtual" constraints coefficient matrix
  # the list has 1 element per constraint
  #   each element of the list has a vector of indexes indicating which x variables have nonzero constraint coefficents
  
  # cc_sparse is a nonzero constraints coefficients data frame
  # ivar gives the variable name for the integer index indicating each CONSTRAINT
  # jvar gives the variable name (character) for the integer index indicating the nonzero x variables for that constraint
  
  jac_sparse <- dlply(cc_sparse, ivar, function(x) return(x[[jvar]]))
  attributes(jac_sparse) <- NULL # to be safe
  
  return(jac_sparse)
}


