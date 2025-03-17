begin
    using Random
    using Statistics
    using Gurobi
    using JuMP
end

begin
    include("DataStruct.jl")
    include("MyDistribution.jl")
end

## ============================================================================================================== ##
##   ######    ######   ##    ##   #######  ########  #######   ##    ##   ######   ########   ######   #######   ##
##  ##    ##  ##    ##  ###   ##  ##           ##     ##    ##  ##    ##  ##    ##     ##     ##    ##  ##    ##  ##
##  ##        ##    ##  ## ## ##   ######      ##     #######   ##    ##  ##           ##     ##    ##  #######   ##
##  ##    ##  ##    ##  ##   ###        ##     ##     ##  ##    ##    ##  ##    ##     ##     ##    ##  ##  ##    ##
##   ######    ######   ##    ##  #######      ##     ##   ##    ######    ######      ##      ######   ##   ##   ##
## ============================================================================================================== ##

"""
```Julia
    Route(id::Int64, O::Int64, mails::Vector{Int64})
    Route(id::Int64, O::Int64, mails::NTuple{N, Int64} = ())
```
Two constructor of the `Route` struct.
`id` is the route id. `O` is the length of the assignment vector. `mails` contain each mail in the Route (faster data struct than a vector)
"""

function Route(id::Int64, O::Int64, mails::Vector{Int64})
    return Route(id, O, (mails...,))
end

function Route(id::Int64, O::Int64, mails::NTuple{N, Int64} = ()) where N
    newAssignment::Vector{Int64} = zeros(Int64, O)
    isempty(mails) || (newAssignment[1:length(mails)] .= mails)
    return Route(id, newAssignment, mails)
end

## ============================================================================================================== ##
##                      ########  ########  ########  ##    ##  ########   #######   #######                      ##
##                      ##           ##        ##     ###   ##  ##        ##        ##                            ##
##                      #####        ##        ##     ## ## ##  #####      ######    ######                       ##
##                      ##           ##        ##     ##   ###  ##              ##        ##                      ##
##                      ##        ########     ##     ##    ##  ########  #######   #######                       ##
## ============================================================================================================== ##

abstract type RouteFitness end

struct Identity         <: RouteFitness end
struct MaxMin           <: RouteFitness end
struct Max              <: RouteFitness end

struct MailNb           <: RouteFitness end
struct MailVolume       <: RouteFitness end
struct MeanVolume       <: RouteFitness end
struct MeanAssignment   <: RouteFitness end

struct MailStd          <: RouteFitness end
struct AssignmentStd    <: RouteFitness end

@inline specialized(r::Route{N},::Type{Identity}        ) where N   = float(-r.id)
@inline specialized(r::Route{N},::Type{MaxMin}          ) where N   = float(maximum(r.mail) - minimum(r.mail))
@inline specialized(r::Route{N},::Type{Max}             ) where N   = float(maximum(r.mail))

@inline specialized(r::Route{N},::Type{MailNb}          ) where N   = float(length(r.mail))
@inline specialized(r::Route{N},::Type{MailVolume}      ) where N   = sum(r.mail)
@inline specialized(r::Route{N},::Type{MeanVolume}      ) where N   = sum(r.mail) / length(r.mail)
@inline specialized(r::Route{N},::Type{MeanAssignment}  ) where N   = sum(r.mail) / length(r.assignment)

@inline specialized(r::Route{N},::Type{MailStd}         ) where N   = std(r.mail)
@inline specialized(r::Route{N},::Type{AssignmentStd}   ) where N   = std(r.assignment)



"""
```Julia
    fitness(r::Route{N}, TAG::Type{<:RouteFitness}) where N         # one criterion
    fitness(r::Route{N}, TAGs::Tuple{Vararg{DataType, N}}) where N  # multiple criterion
```
Evaluate a given route `r` on one or more criteria.

There are 9 tags each related to criteria to assign a value to a route:
 - `Identity`: return the route id
 - `MaxMin`: difference between the biggest and smallest mail
 - `Max`: biggest mail in the route
 - `StdBatches`: standard deviation of the mail sizes.
 - `StdAssignment`: standard deviation of the assignment vector (including zeros which is not the case with `StdBatches`).
 - `NbBatches`: number of mails in the route
 - `MeanBatches`: average mail size.
 - `MeanAssignment`: average mail size (include zeros of the assignment vector).
 - `SumBatches`: total mail volume in the route
"""

@inline fitness(r::Route{N}, TAG::Type{<:RouteFitness}) where N                     = specialized(r, TAG)
@inline fitness(r::Route{N}, TAGs::Tuple{DataType}) where N                         = (specialized(r, TAGs[1]))::Tuple{Float64}
@inline fitness(r::Route{N}, TAGs::Tuple{DataType, DataType}) where N               = (specialized(r, TAGs[1]), specialized(r, TAGs[2]))::Tuple{Float64, Float64}
@inline fitness(r::Route{N}, TAGs::Tuple{DataType, DataType, DataType}) where N     = (specialized(r, TAGs[1]), specialized(r, TAGs[2]), specialized(r, TAGs[3]))::Tuple{Float64, Float64, Float64}
@inline fitness(r::Route{N}, TAGs::Tuple{Vararg{DataType, N}}) where N              = ntuple(i -> specialized(r, TAGs[i]), N)::NTuple{N, Float64}

## ============================================================================================================== ##
##                      ######    ########   #######  #######   ##         ######   ##    ##                      ##
##                      ##    ##     ##     ##        ##    ##  ##        ##    ##   ##  ##                       ##
##                      ##     #     ##      ######   #######   ##        ########     ##                         ##
##                      ##    ##     ##           ##  ##        ##        ##    ##     ##                         ##
##                      ######    ########  #######   ##        ########  ##    ##     ##                         ##
## ============================================================================================================== ##

function Base.show(io::IO, r::Route)
    print(io, "< id: $(r.id)")
    
    if length(r.assignment) <= 40
        print(io, ", assignment: ")
        for e in r.assignment
            print(io, "$(if e != 0 string(e) else "_" end), ")
        end

        print(io, "mails: $(r.mail)>")
    else
        print(io, " >")
    end
end

## ============================================================================================================== ##
##                                     ########  ########   #######  ########                                     ##
##                                        ##     ##        ##           ##                                        ##
##                                        ##     #####      ######      ##                                        ##
##                                        ##     ##              ##     ##                                        ##
##                                        ##     ########  #######      ##                                        ##
## ============================================================================================================== ##

"""
```Julia
    isRouteValid(r::Route{N}) where N
```
Test the validity of a route: 
 - precedence constraint
 - each mail appears exactly once
 - no extra mail
"""

function isRouteValid(r::Route{N}) where N
    i::Int64 = 1
    for e in r.assignment
        if e != 0
            (i > N) && (return false)
            (e == r.mail[i]) ? (i += 1) : (return false) 
        end
    end
    return (i == N+1)
end

"""
```Julia
    isRouteRaw(r::Route{N}) where N
```
Test if all the mail of the route `r` are allign to the left.
"""
function isRouteRaw(r::Route{N}) where N
    i::Int64 = 1
    left::Bool = true

    isRouteValid(r) || (return false)

    for e in r.assignment
        if e == 0
            left = false
        else # e != 0
            (!left) && (return false)
        end
    end
    return true
end

## ============================================================================================================== ##
##                                     ##    ##  ########   #######   ######                                      ##
##                                     ###  ###     ##     ##        ##    ##                                     ##
##                                     ## ## ##     ##      ######   ##                                           ##
##                                     ##    ##     ##           ##  ##    ##                                     ##
##                                     ##    ##  ########  #######    ######                                      ##
## ============================================================================================================== ##

"""
```Julia
    Random.shuffle!(r::Route)
    Random.shuffle(r::Route)
```
Move mail across the assignment vector randomly. (still following the precedence constraint)
"""

function Random.shuffle!(r::Route)
    n::Int64 = length(r.assignment)
    b::Int64 = length(r.mail)

    # new position for each batch in the route (ordered)
    newPos::Vector{Int64} = sort!(randperm(n)[1:b])
    
    # new assignment vector
    newAssignment::Vector{Int64} = zeros(Int64, n)
    for i=1:b
        newAssignment[newPos[i]] = r.mail[i]
    end

    r.assignment = newAssignment

    return r
end

function Random.shuffle(r::Route)
    return shuffle!(deepcopy(r))
end

"""
```Julia
    randomRoute(id::Int64 = 0, O::Int64 = 20, Br::Union{Int64, Nothing} = nothing, C::Int64 = 20, distrib_fct::Function = distrib_1onN)
```
Create a random route.

`id` is the route id. `O` is the assignment vector length. `Br` is the number of mail wanted in the route. `C` is the maximum size of a mail. `distrib_fct` is the distribution used to generate the mail sizes.

"""

function randomRoute(
        id          ::Int64                 = 0             , 
        O           ::Int64                 = 20            , 
        Br          ::Union{Int64, Nothing} = nothing       ,
        C           ::Int64                 = 20            , 
        distrib_fct ::Function              = distrib_1onN  ,
    )

    (Br === nothing) && (Br = rand(Route(Int64, O/4):ceil(Int64, O/2)))

    distrib = distrib_fct(C)

    mails::Tuple{Vararg{Int64, Br}} = ntuple(i -> rand(distrib), Br)

    assignment::Vector{Int64} = zeros(Int64, O)
    assignment[1:length(mails)] .+= mails

    return Route(id, assignment, mails) 
end

function partitioning_01_KP(routes::Vector{Route}, Lmax::Int64, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())

    @assert !isempty(routes)

    items::Vector{Int64} = map(x::Route -> sum(x.mail), routes)

    I::Int64 = length(items)
    O::Int64 = length(routes[1].assignment)

    model = Model(() -> Gurobi.Optimizer(env))
    set_silent(model)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", tl)

    @variable(model, x[i=1:I], Bin)

    @objective(model, Max, sum([x[i] * items[i] for i=1:I]))  

    @constraint(model, sum([x[i] * items[i] for i=1:I]) <= Lmax * O)

    optimize!(model)

    res = [k for (k,e) in enumerate(value.(x)) if e == 1]

    return model, res, sum([items[e] for e in res])
end

function full_partitioning_01_KP(routes::Vector{Route}, Lmax::Int64, lb::Int64 = 20, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())

    @assert !isempty(routes)

    items::Vector{Int64} = map(x::Route -> sum(x.mail), routes)

    I::Int64 = length(items)
    O::Int64 = length(routes[1].assignment)

    model = Model(() -> Gurobi.Optimizer(env))
    set_silent(model)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", tl)

    @variable(model, x[i=1:I, s=1:lb], Bin)
    @variable(model, y[s=1:lb], Bin)

    @objective(model, Min, sum([y[s] for s=1:lb]))  

    for s=1:lb
        @constraint(model, sum([x[i, s] * items[i] for i=1:I]) <= Lmax * O * y[s])
    end

    for i=1:I
        @constraint(model, sum(x[i, :]) == 1)
    end

    if lb >= 2
        for s=1:lb-1
            @constraint(model, y[s] >= y[s+1])
        end
    end

    optimize!(model)

    if termination_status(model) == MOI.OPTIMAL
        return model, [sum([value(x[i, s]) * items[i] for i=1:I]) for s=1:lb if value(y[s]) == 1], sum(value.(y)), [k for s=1:lb if value(y[s]) == 1 for (k, i) in enumerate(1:I) if value(x[i, s]) == 1]
    else
        print("<x>")
        return nothing, nothing, nothing, randperm(length(items))
    end
end