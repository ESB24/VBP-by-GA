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

# ================================================== #
#                    Constructor                     #
# ================================================== #

"""
```Julia
    Round(id::Int64, O::Int64, batches::Vector{Int64})
    Round(id::Int64, O::Int64, batches::NTuple{N, Int64} = ())
```
Two constructor of the `Round` struct.
`id` is the route id. `O` is the length of the assignment vector. `batches` contain each mail in the round (faster data struct than a vector)
"""

function Round(id::Int64, O::Int64, batches::Vector{Int64})
    println("nop")
    return Round(id, O, (batches...,))
end

function Round(id::Int64, O::Int64, batches::NTuple{N, Int64} = ()) where N
    newAssignment::Vector{Int64} = zeros(Int64, O)
    isempty(batches) || (newAssignment[1:length(batches)] .= batches)
    return Round(id, newAssignment, batches)
end

# ================================================== #
#                       Fitness                      #
# ================================================== #
abstract type FitnessRound end

struct Identity         <: FitnessRound end
struct MaxMin           <: FitnessRound end
struct Max              <: FitnessRound end
struct StdBatches       <: FitnessRound end
struct StdAssignment    <: FitnessRound end
struct NbBatches        <: FitnessRound end
struct MeanBatches      <: FitnessRound end
struct MeanAssignment   <: FitnessRound end
struct SumBatches       <: FitnessRound end

@inline specialized(r::Round{N},::Type{Identity}) where N           = float(-r.id)
@inline specialized(r::Round{N},::Type{MaxMin}) where N             = float(maximum(r.batches) - minimum(r.batches))
@inline specialized(r::Round{N},::Type{Max}) where N                = float(maximum(r.batches))
@inline specialized(r::Round{N},::Type{StdBatches}) where N         = std(r.batches)
@inline specialized(r::Round{N},::Type{StdAssignment}) where N      = std(r.assignment)
@inline specialized(r::Round{N},::Type{NbBatches}) where N          = float(length(r.batches))
@inline specialized(r::Round{N},::Type{MeanBatches}) where N        = sum(r.batches) / length(r.batches)
@inline specialized(r::Round{N},::Type{MeanAssignment}) where N     = sum(r.batches) / length(r.assignment)
@inline specialized(r::Round{N},::Type{SumBatches}) where N         = sum(r.batches)


"""
```Julia
    fitness(r::Round{N}, TAG::Type{<:FitnessRound}) where N         # one criterion
    fitness(r::Round{N}, TAGs::Tuple{Vararg{DataType, N}}) where N  # multiple criterion
```
Evaluate a given route `r` on one or more criteria.

There are 9 tags each related to criteria to assign a value to a route:
 - `Identity`: return the route id
 - `MaxMin`: difference between the biggest and smallest mail
 - `Max`: biggest mail in the route
 - `StdBatches`: standard deviation of the mail sizes.
 - `StdAssignment`: standard deviation of the assignment vector (including zeros which is not the case with `StdBatches`).
 - `NbBatches`: number of batches (mail) in the route
 - `MeanBatches`: average mail size.
 - `MeanAssignment`: average mail size (include zeros of the assignment vector).
 - `SumBatches`: total mail volume in the route
"""

@inline fitness(r::Round{N}, TAG::Type{<:FitnessRound}) where N                     = specialized(r, TAG)
@inline fitness(r::Round{N}, TAGs::Tuple{DataType}) where N                         = (specialized(r, TAGs[1]))::Tuple{Float64}
@inline fitness(r::Round{N}, TAGs::Tuple{DataType, DataType}) where N               = (specialized(r, TAGs[1]), specialized(r, TAGs[2]))::Tuple{Float64, Float64}
@inline fitness(r::Round{N}, TAGs::Tuple{DataType, DataType, DataType}) where N     = (specialized(r, TAGs[1]), specialized(r, TAGs[2]), specialized(r, TAGs[3]))::Tuple{Float64, Float64, Float64}
@inline fitness(r::Round{N}, TAGs::Tuple{Vararg{DataType, N}}) where N              = ntuple(i -> specialized(r, TAGs[i]), N)::NTuple{N, Float64}

# ================================================== #
#                       Display                      #
# ================================================== #

function Base.show(io::IO, r::Round)
    print(io, "< id: $(r.id)")
    
    if length(r.assignment) <= 40
        print(io, ", assignment: ")
        for e in r.assignment
            print(io, "$(if e != 0 string(e) else "_" end), ")
        end

        print(io, "batches: $(r.batches)>")
    else
        print(io, " >")
    end
end

# ================================================== #
#                     Miscelaneous                   #
# ================================================== #

"""
```Julia
    Random.shuffle!(r::Round)
    Random.shuffle(r::Round)
```
Move mail across the assignment vector randomly. (still following the precedence constraint)
"""

function Random.shuffle!(r::Round)
    n::Int64 = length(r.assignment)
    b::Int64 = length(r.batches)

    # new position for each batch in the round (ordered)
    newPos::Vector{Int64} = sort!(randperm(n)[1:b])
    
    # new assignment vector
    newAssignment::Vector{Int64} = zeros(Int64, n)
    for i=1:b
        newAssignment[newPos[i]] = r.batches[i]
    end

    r.assignment = newAssignment

    return r
end

function Random.shuffle(r::Round)
    return shuffle!(deepcopy(r))
end

"""
```Julia
    validRound(r::Round{N}) where N
```
Test the validity of a route: 
 - precenece constraint
 - each mail appears exactly once
 - no extra mail
"""

function validRound(r::Round{N}) where N
    i::Int64 = 1
    for e in r.assignment
        if e != 0
            (i > N) && (return false)
            (e == r.batches[i]) ? (i += 1) : (return false) 
        end
    end
    return (i == N+1)
end

"""
```Julia
    rawRound(r::Round{N}) where N
```
Test if all the mail of the route `r` are allign to the left.
"""
function rawRound(r::Round{N}) where N
    i::Int64 = 1
    left::Bool = true
    for e in r.assignment
        if e == 0
            left = false
        else # e != 0
            (!left) && (return false)
        end
    end
    return true
end

"""
```Julia
    randRound(id::Int64 = 0, O::Int64 = 20, Br::Union{Int64, Nothing} = nothing, C::Int64 = 20, distrib_fct::Function = distrib_1onN)
```
Create a random route.

`id` is the route id. `O` is the assignment vector length. `Br` is the number of mail wanted in the route. `C` is the maximum size of a mail. `distrib_fct` is the distribution used to generate the mail sizes.

"""

function randRound(
        id          ::Int64                 = 0             , 
        O           ::Int64                 = 20            , 
        Br          ::Union{Int64, Nothing} = nothing       ,
        C           ::Int64                 = 20            , 
        distrib_fct ::Function              = distrib_1onN  ,
    )

    (Br === nothing) && (Br = rand(round(Int64, O/4):ceil(Int64, O/2)))

    distrib = distrib_fct(C)

    batches::Tuple{Vararg{Int64, Br}} = ntuple(i -> rand(distrib), Br)

    assignment::Vector{Int64} = zeros(Int64, O)
    assignment[1:length(batches)] .+= batches

    return Round(id, assignment, batches) 
end

function partitioning_01_KP(routes::Vector{Round}, Lmax::Int64, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())

    @assert !isempty(routes)

    items::Vector{Int64} = map(x::Round -> sum(x.batches), routes)

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

function full_partitioning_01_KP(routes::Vector{Round}, Lmax::Int64, lb::Int64 = 20, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())

    @assert !isempty(routes)

    items::Vector{Int64} = map(x::Round -> sum(x.batches), routes)

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