# ================================================== #
#                        Core                        #
# ================================================== #

mutable struct Round{N}
    id              ::Int64                     # ID of the round
    assignment      ::Vector{Int64}             # each batches places on the outputs
    batches         ::NTuple{N, Int64}          # list all batches in the Round
end

mutable struct Session
    C               ::Int64                     # Capacity of each outputs
    rounds          ::Vector{Round}             # Round containe in the solution (letter well placed)
    loads           ::Vector{Int64}             # curent state of the outputs
end

mutable struct Solution
    permutation     ::Vector{Int64}             # used permutation to find the solution
    sessions        ::Vector{Session}           # List of solution
end

mutable struct Instance
    C               ::Int64                     # capacity of the bin (Lmax)
    nbOut           ::Int64                     # number of outputs of the machine |O|
    nbRound         ::Int64                     # number of route in the instance |R|
    rounds          ::Vector{Round}             # set of route
end

# ================================================== #
#                         GA                         #
# ================================================== #

mutable struct GenVal
    elite           ::Vector{Float64}           # best 20% of previous population or created by FFD in G0 
    ffr             ::Vector{Float64}           # individuals created using FFD heuristic 
    cross_over      ::Vector{Float64}           # individuals created using Cross-Over
end

mutable struct Generation
    id              ::Int64
    pop             ::Vector{Solution}
end
