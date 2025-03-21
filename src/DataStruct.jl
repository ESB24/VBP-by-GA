# ================================================== #
#                        Core                        #
# ================================================== #

mutable struct Route{N}
    id              ::Int64                     # ID of the route
    assignment      ::Vector{Int64}             # mails positionned on the outputs
    mail            ::NTuple{N, Int64}          # list all mails in the Route
end

mutable struct Session
    Lmax            ::Int64                     # Maximum output Load
    route           ::Vector{Route}             # Route currently in the solution
    load            ::Vector{Int64}             # Current output loads
end

mutable struct Solution
    permutation     ::Vector{Int64}             # used permutation to find the solution
    sessions        ::Vector{Session}           # List of solution
end

mutable struct Instance
    Lmax            ::Int64                     # capacity of the bin (Lmax)
    nbOut           ::Int64                     # number of outputs of the machine |O|
    nbRoute         ::Int64                     # number of route in the instance |R|
    route           ::Vector{Route}             # set of route
end

# ================================================== #
#                         GA                         #
# ================================================== #

call = 0
repairedBuild = 0
improvedOverAll = 0
locked = 0
delta1 = 0.3
delta2 = 0.9

mutable struct GenVal
    elite           ::Vector{Float64}           # best 20% of previous population or created by FFD in G0 
    ffr             ::Vector{Float64}           # individuals created using FFD heuristic 
    cross_over      ::Vector{Float64}           # individuals created using Cross-Over
end

mutable struct Generation
    id              ::Int64
    pop             ::Vector{Solution}
end
