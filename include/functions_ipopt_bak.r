

#****************************************************************************************************
#                constraint evaluation and coefficient functions for ipoptr SPARSE -- same for all ####
#****************************************************************************************************
eval_g <- function(x, inputs) {
  # constraints that must hold in the solution - just give the LHS of the expression
  # return a vector where each element evaluates a constraint (i.e., sum of (x * a ccmat column), for each column)
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  
  # inputs$cc_sparse has the fixed constraint coefficients in sparse form in a dataframe that has:
  #   i -- the constraint number (based on those constraints we send to ipoptr, which could be a subset of all)
  #   j -- index into x (i.e., the variable number)
  #   nzcc
  
  constraint_tbl <- inputs$cc_sparse %>%
    group_by(i) %>%
    summarise(constraint_value=sum(nzcc * x[j]))
  # the column constraint_value is a vector, in the order we want
  
  return(constraint_tbl$constraint_value)
}




eval_jac_g <- function(x, inputs){
  # the Jacobian is the matrix of first partial derivatives of constraints (these derivatives may be constants)
  # this function evaluates the Jacobian at point x
  
  # return: a vector where each element gives a NONZERO partial derivative of constraints wrt change in x
  # so that the first m items are the derivs with respect to each element of first column of ccmat
  # and next m items are derivs with respect to 2nd column of ccmat, and so on
  # so that it returns a vector with length=nrows x ncolumns in ccmat
  
  # because constraints in this problem are linear, the derivatives are all constants
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  
  return(inputs$cc_sparse$nzcc)
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
  
  return(jac_sparse)
}




#****************************************************************************************************
#                x^p + x^-p {xtop} -- functions for ipoptr ####
#****************************************************************************************************
eval_f_xtop <- function(x, inputs) {
  # objective function - evaluates to a single number
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  # here are the objective function, the 1st deriv, and the 2nd deriv
  # http://www.derivative-calculator.net/
  # w{x^p + x^(-p) - 2}                                 objective function
  # w{px^(p-1) - px^(-p-1)}                             first deriv
  # p*w*x^(-p-2)*((p-1)*x^(2*p)+p+1)                    second deriv
  
  # make it easier to read:
  p <- inputs$p
  w <- inputs$wt / inputs$objscale
  
  obj <- sum(w * {x^p + x^(-p) -2})
  
  return(obj)
}


eval_grad_f_xtop <- function(x, inputs){
  # gradient of objective function - a vector length x 
  # giving the partial derivatives of obj wrt each x[i]
  
  # ipoptr requires that ALL functions receive the same arguments, so I pass the inputs list to ALL functions
  
  # http://www.derivative-calculator.net/
  # w{x^p + x^(-p) - 2}                                 objective function
  # w{px^(p-1) - px^(-p-1)}                             first deriv
  # p*w*x^(-p-2)*((p-1)*x^(2*p)+p+1)                    second deriv
  
  # make it easier to read:
  p <- inputs$p
  w <- inputs$wt / inputs$objscale
  
  gradf <- w * (p * x^(p-1) - p * x^(-p-1))
  
  return(gradf)
}


eval_h_xtop <- function(x, obj_factor, hessian_lambda, inputs){
  # The Hessian matrix has many zero elements and so we set it up as a sparse matrix
  # We only keep the (potentially) non-zero values that run along the diagonal.
  
  # http://www.derivative-calculator.net/
  # w{x^p + x^(-p) - 2}                                 objective function
  # w{px^(p-1) - px^(-p-1)}                             first deriv
  # p*w*x^(-p-2)*((p-1)*x^(2*p)+p+1)                    second deriv
  
  # make it easier to read:
  p <- inputs$p
  w <- inputs$wt / inputs$objscale
  
  hess <- obj_factor * 
    { p*w*x^(-p-2) * ((p-1)*x^(2*p)+p+1) }
  
  return(hess)
}



