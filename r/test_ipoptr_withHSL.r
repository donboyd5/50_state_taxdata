
library(ipoptr)

# ?ipoptr
# ?'ipoptr-package'


#****************************************************************************************************
#                Do a simple check on the banana problem, with alternative linear solvers ####
#****************************************************************************************************

## Rosenbrock Banana function
eval_f <- function(x) {   
  return( 100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2 )
}

## Gradient of Rosenbrock Banana function
eval_grad_f <- function(x) { 
  return( c( -400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
             200 * (x[2] - x[1] * x[1])) )
}

# The Hessian for this problem is actually dense, 
# This is a symmetric matrix, fill the lower left triangle only.
eval_h_structure <- list( c(1), c(1,2) )

eval_h <- function( x, obj_factor, hessian_lambda ) {
  return( obj_factor*c( 2 - 400*(x[2] - x[1]^2) + 800*x[1]^2,      # 1,1
                        -400*x[1],                                 # 2,1
                        200 ) )                                    # 2,2
}

# initial values
x0 <- c( -1.2, 1 )

opts.default <- list("print_level"=5,
             "file_print_level"=12,
             "output_file"="banana_defaultLS.out",
             "tol"=1.0e-8)

# solve Rosenbrock Banana function with analytic hessian and default solver
print( ipoptr( x0=x0, 
               eval_f=eval_f, 
               eval_grad_f=eval_grad_f, 
               eval_h=eval_h,
               eval_h_structure=eval_h_structure,
               opts=opts.default))
# do you see:
#    EXIT: Optimal Solution Found.
# if so, all good


# Repeat with an alternative linear solver
# define one of the HSL solvers as the linear solver: ma57 ma77 ma86 ma97
opts.alt <- list("print_level"=5,
                 "file_print_level"=12,
                 "output_file"="banana_altLS.out",
                 "linear_solver" = "ma57", # mumps pardiso ma27 ma57 ma77 ma86 ma97
                 "tol"=1.0e-8)

print( ipoptr( x0=x0, 
               eval_f=eval_f, 
               eval_grad_f=eval_grad_f, 
               eval_h=eval_h,
               eval_h_structure=eval_h_structure,
               opts=opts.alt))

# do you see:
#    EXIT: Optimal Solution Found.
# if so, all good

