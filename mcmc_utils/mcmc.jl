using JLD2, Distributions
"""
optimize(initial_𝑪, nll; nt = 1000, restart = 0, proposal = [], scale = 1, filename = [], rescale = true, freq = 1001)
# Description
- A generic optimizer using RWMCMC. It is generally better to use Optim
# Arguments
-  `initial_𝑪`:(vector) initial parameter
- `nll`:(function) negative log-likelihood. The function to minimize
# Keyword Arguments
- `nt`: (int), how many steps of the random walk to take
- `restart`: (int), restart at the optimal value this many times
- `proposal`: (function), proposal function for performing random walk
- `scale`: (real), scale for constructing default proposal
- `filename`: (string), a place to save the JLD2 file for optimization
- `rescale`: (boolean), allows one to rescale the loss function over iterations
- `freq`: how often to output progress, make this larger than nt for no output
- `verbose`: (boolean), outputs optimal values with frequence = freq
# Comments
- This is a "prep step" in mcmc
"""
function optimize(initial_𝑪, nll; nt = 10000, restart = 0, proposal = [], scale = 0.2, filename = [], rescale = true, freq = 10001, verbose = true)
    if proposal == []
        perturbation = closure_proposal(initial_𝑪 * scale)
    else
        perturbation = proposal
    end
    if rescale == true
        scale = nll(initial_𝑪)
    else
        scale = 1.0
    end
    ℒ(𝑪) = nll(𝑪) / scale
    # perform random walk
    tmp_𝑪 = copy(initial_𝑪)
    for i in 1:(restart+1)
        new_𝑪, new_ε = markov_chain(ℒ, tmp_𝑪, perturbation, nt; freq = freq, filename = filename, verbose = verbose)
        # pick out new optimal value
        optimal_index = argmin(new_ε)
        opt_𝑪 = new_𝑪[:, optimal_index]
        tmp_𝑪 = opt_𝑪
        if rescale == true
            ℒ(𝑪) = nll(𝑪) / nll(tmp_𝑪)
        end
    end
    return tmp_𝑪
end


"""
optimize_and_estimate_proposal(initial_𝑪, nll, left_bounds, right_bounds; nt = 1000, restart = 0, proposal = [], scale = 1, filename = [], rescale = true, freq = 1001)
# Description
- A generic optimizer using RWMCMC. It also tries to estimate a new proposal
# Arguments
-  `initial_𝑪`:(vector) initial parameter
- `nll`:(function) negative log-likelihood. The function to minimize
- `left_bounds`: bounds for the proposal
- `right_bounds`: bounds for the proposal
# Keyword Arguments
- `nt`: (int), how many steps of the random walk to take
- `restart`: (int), restart at the optimal value this many times
- `proposal`: (function), proposal function for performing random walk
- `scale`: (real), scale for constructing default proposal
- `filename`: (string), a place to save the JLD2 file for optimization
- `rescale`: (boolean), allows one to rescale the loss function over iterations
- `freq`: how often to output progress, make this larger than nt for no output
- `verbose`: (boolean), outputs optimal values with frequence = freq
# Comments
- This is a "prep step" in mcmc
"""
function optimize_and_estimate_proposal(initial_𝑪, nll, left_bounds, right_bounds; nt = 10000, restart = 0, proposal = [], scale = 0.2, filename = [], rescale = true, freq = 10001, verbose = true)
    if proposal == []
        perturbation = closure_proposal(initial_𝑪 * scale, left_bounds = left_bounds, right_bounds = right_bounds)
    else
        perturbation = proposal
    end
    if rescale == true
        scale = nll(initial_𝑪)
    else
        scale = 1.0
    end
    ℒ(𝑪) = nll(𝑪) / scale
    # perform random walk
    tmp_𝑪 = copy(initial_𝑪)
    Σ = randn(length(initial_𝑪),length(initial_𝑪))
    for i in 1:(restart+1)
        new_𝑪, new_ε = markov_chain(ℒ, tmp_𝑪, perturbation, nt; freq = freq, filename = filename, verbose = verbose)
        # pick out new optimal value
        optimal_index = argmin(new_ε)
        opt_𝑪 = new_𝑪[:, optimal_index]
        tmp_𝑪 = opt_𝑪
        tmp_Σ = cov(new_𝑪')
        println(Σ)
        @. Σ = tmp_Σ
        perturbation = closure_proposal(Σ, left_bounds = left_bounds, right_bounds = right_bounds)
        if rescale == true
            ℒ(𝑪) = nll(𝑪) / nll(tmp_𝑪)
        end
    end
    return tmp_𝑪, Σ
end

function optimize_and_estimate_proposal(initial_𝑪, nll; nt = 10000, restart = 0, proposal = [], scale = 0.2, filename = [], rescale = true, freq = 10001, verbose = true)
    if proposal == []
        perturbation = closure_proposal(initial_𝑪 * scale)
    else
        perturbation = proposal
    end
    if rescale == true
        scale = nll(initial_𝑪)
    else
        scale = 1.0
    end
    ℒ(𝑪) = nll(𝑪) / scale
    # perform random walk
    tmp_𝑪 = copy(initial_𝑪)
    Σ = randn(length(initial_𝑪),length(initial_𝑪))
    for i in 1:(restart+1)
        new_𝑪, new_ε = markov_chain(ℒ, tmp_𝑪, perturbation, nt; freq = freq, verbose = verbose)
        # pick out new optimal value
        optimal_index = argmin(new_ε)
        opt_𝑪 = new_𝑪[:, optimal_index]
        tmp_𝑪 = opt_𝑪
        tmp_Σ = cov(new_𝑪')
        @. Σ = tmp_Σ
        Σ = Σ + sqrt(eps(maximum(Σ))) * I
        perturbation = closure_proposal(Σ)
        if rescale == true
            ℒ(𝑪) = nll(𝑪) / nll(tmp_𝑪)
        end
    end
    return tmp_𝑪, Σ
end

# Defines several functions useful for performing a random walk


"""
accept_reject(Δℒ)
# Description
- Determines the accept or reject criteria for the Monte Carlo method.
# Input: Δℒ
- `Δℒ`: (scalar) Difference of negative log likehood functions
# Output
- Boolean Value: True or False
"""
accept_reject(Δℒ) = log(rand(Uniform(0, 1))) <= Δℒ

"""
markov_link(nll, 𝑪, ε, proposal)
# Description
- Takes a single step in the random walk markov chain monte carlo algorithm and outputs proposal parameters, new parameters, and the evaluate of the loss function
# Arguments
- `nll`: The negative log-likelihood function. In the absence of priors this becomes a loss function
- `𝑪`: (array), current parameter
- `ε`: (scalar), ε = nll(𝑪). The value of negative log-likelihood of the current parameter
- `proposal`: (function), determines the proposal step
# Return
- `new_𝑪`: The value of the accepted 𝑪
- `new_ε`: value of nll(new_𝑪)
- `proposal_𝑪`: The 𝑪 from the "proposal step". Was either rejected or accepted.
- `proposal_ε`: value of nll(test_𝑪)
"""
function markov_link(nll, 𝑪, ε, proposal)
    proposal_𝑪 = proposal(𝑪)
    proposal_ε = nll(proposal_𝑪)
    Δε = (ε - proposal_ε)
    if accept_reject(Δε)
        new_ε = proposal_ε
        new_𝑪 = proposal_𝑪
    else
        new_ε = ε
        new_𝑪 = 𝑪
    end
    return new_𝑪, new_ε, proposal_𝑪, proposal_ε
end



"""
markov_chain_with_save(nll, init_𝑪, proposal, nt, filename, freq)
# Description
- A random walk that computes the posterior distribution
# Arguments
- `nll`: The negative log-likelihood function. In the absence of priors this becomes a loss function
- `init_𝑪`: (Array), initial parameter values
- `proposal`: (function), proposal function for MCMC
- `nt`: (Int) number of markov chain monte carlo steps
- `perturb`: a function that performs a perturbation of 𝑪
# Keyword Arguments
- `filename`: name for output file in JLD2 format
- `freq`: how often to save output (in terms of iterations)
- `verbose`: (bool), if true then print current optimal parameters
# Return
- `param`: The matrix of accepted parameters in the random walk
- `ε`: The array of errors associated with each step in param chain
"""
function markov_chain(nll, initial_𝑪, proposal, nt;
                      filename = [], freq = 1, verbose = false)
    𝑪 = ones(length(initial_𝑪),nt+1)
    @. 𝑪[:,1] = initial_𝑪
    proposal_𝑪 = copy(𝑪)
    ε = ones(nt+1)
    proposal_ε = copy(ε)
    ε[1] = nll(initial_𝑪)
    for i in 1:nt
        new_𝑪, new_ε, proposed_𝑪, proposed_ε = markov_link(nll, 𝑪[:,i], ε[i], proposal)
        @. 𝑪[:,i+1] = new_𝑪
        ε[i+1] = new_ε
        @. proposal_𝑪[:,i+1] = proposed_𝑪
        proposal_ε[i+1] = proposed_ε
        if i%freq==0
            println("saving index " * string(i))
            if !isempty(filename)
                @save filename ε 𝑪 proposal_ε proposal_𝑪
            end
            if verbose==true
                indmin = argmin(ε[1:i])
                println("The current optimal parameters are")
                println(𝑪[:,indmin])
                println("The loss function is " * string(ε[indmin]))
                tmpstrng = string(ε[1] / ε[indmin] )
                println("This is an improvement of " * tmpstrng)
                acceptance_rate = sum(ε[1:i] .== proposal_ε[1:i]) / length(ε[1:i])
                println("The current acceptance rate is $acceptance_rate")
            end
        end
    end
    return 𝑪, ε
end

"""
torus(x, a, b)
# Description
- Takes x ∈ ℝ and outputs torus(x) ∈ [a, b] in a periodic way.
- If a particle is moving to the right then it will pop from b to the point a
# Arguments: x, a, b
- `x`: (scalar). Current location of particle
- `a`: (scalar). left endpoint of interval
- `b`: (scalar). right endpoint of interval
# Output
-  `y`: (scalar). a value in the interval [a,b]
"""
torus(x::Number, a::Number, b::Number) = (((x-a)/(b-a))%1 - 0.5 * (sign((x-a)/(b-a)) - 1) )*(b-a) + a

"""
torus(x, a, b)
# Description
- Takes x ∈ ℝⁿ and outputs torus(x) ∈ ∏[aⁿ, bⁿ] in a periodic way.
- If a particle is moving to the right then it will pop from one part of the box to the oher
# Arguments: x, a, b
- `x`: (array). Current location of particle
- `a`: (array). left endpoint of tensor product interval
- `b`: (array). right endpoint of tensor product interval
# Output
-  `y`: (array). a value in the interval ∏[aⁿ, bⁿ]
"""
function torus(x::AbstractArray, a::AbstractArray, b::AbstractArray)
    N = length(x)
    y = zeros(N)
    for i in 1:N
        y[i] = torus(x[i], a[i], b[i])
    end
    return y
end


"""
closure_proprosal(covariance = Σ; left_bounds = [], right_bounds = []))
# Description
- Constructs a proposal for the Monte Carlo method.
# Arguments
- `covariance`: (vector) proposal parameter
# Keyword Arguments
- `left_bounds`: (array), left bounds for parameters
- `right_bounds`: (array), right bounds for parameters
# Output:
- `proposal`: (function), a function that outputs the proposal parameter
"""
function closure_proposal(Σ; left_bounds = [], right_bounds = [])
    perturbation = MvNormal(Σ)
    function proposal(𝑪)
        proposal_𝑪 = copy(𝑪)
        proposal_𝑪 .+= rand(perturbation)
        # limit ranges for the parameters
        if isempty(left_bounds)
            return proposal_𝑪
        else
            return torus(proposal_𝑪, left_bounds, right_bounds)
        end
        return nothing
    end
    return proposal
end
