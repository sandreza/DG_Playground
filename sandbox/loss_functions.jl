"""
closure_loss1(∇², G, x⁰, b; maximum_iterations = 5)

# Description
- create a function closure for convenient loss functions

# Arguments
- ∇²:(array) The matrix
- G: (array) The inverse of the matrix
- x⁰: (vector) starting value
- b: (vector) right hand side

# Keyword Arguments
- maximum_iterations: the maximum number of iterations

# Return
- loss: (function). Input arguments y, a tridiagonal matrix, output the residuation after maximum_iterations (default 5)
"""
function closure_loss1(∇², G, x⁰, b; maximum_iterations = 5)
    ∇²_function(x) = ∇² * x
    function loss1(y)
        preconditioner(x) = y * x
        r = conjugate_gradient!(∇²_function, x⁰, b, track_residual = true, P = preconditioner, maximum_iterations = maximum_iterations)
        return r[end]
    end
    return loss1
end
