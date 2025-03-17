begin
    using Combinatorics
    using Gurobi
    using GLPK
    using JuMP
end

begin 
    include("Route.jl")
    include("DataStruct.jl")
    include("SessionRebuild.jl")
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
    Session(Lmax::Int64, O::Int64)::Sesssion
    Session(Lmax::Int64, r::Route{N}) where N
    Session(Lmax::Int64, r::Vector{Route{N}}) where N
```
Constructor for the `Session` data structure.
`Lmax` is the maximum capacity per output (Lmax). 
`O` is the length of the assignment vector (the session will start empty). 
Otherwise and replacing `O` to create non empty session you could provide a unique or a vector of route `r`. 

"""

function Session(Lmax::Int64, O::Int64)::Session
    return Session(Lmax, [], zeros(Int64, O))
end

function Session(Lmax::Int64, r::Route{N}) where N
    return Session(Lmax, [r], r.assignment)
end

function Session(Lmax::Int64, r::Vector{Union{Route{N}, Route}}) where N
    return Session(Lmax, r, compute_output(r))
end

"""
```Julia
    compute_output(routes::Vector{Union{Route{N}, Route}})
    compute_output(s::Session)::Vector{Int64}
    compute_output!(s::Session)::Session
```

Compute the output of a given session or list of route. (according to there current mail assignment)
If `compute_output` is called the load vector is computed and returned.
If `compute_output!` then the ouput load is also updated in the given session `s`.
"""

function compute_output(routes::Vector{Union{Route{N}, Route}}) where N
    O::Int64 = length(routes[1].assignment)
    R::Int64 = length(routes)

    return [sum([routes[i].assignment[k] for i=1:R]) for k=1:O]
end

function compute_output(s::Session)::Vector{Int64}
    return isempty(s.route) ? zeros(Int64, length(s.load)) : compute_output(s.route)
end

function compute_output!(s::Session)::Session
    s.load = compute_output(s)
    return s
end

## ============================================================================================================== ##
##                                     ########  ########   #######  ########                                     ##
##                                        ##     ##        ##           ##                                        ##
##                                        ##     #####      ######      ##                                        ##
##                                        ##     ##              ##     ##                                        ##
##                                        ##     ########  #######      ##                                        ##
## ============================================================================================================== ##

isValidLoad(loads::Vector{Int64}, Lmax::Int64)::Bool = sum(loads .> Lmax) == 0
isValidLoad(s::Session)::Bool = isValidLoad(s.load, s.Lmax)

function isSessionValid(s::Session, display::Bool = false)
    valid::Bool = true

    for r in s.route
        isRouteValid(r) || (valid = false; (display && (println(" - Route $(r.id) -> invalid"))))
    end

    compute_output!(s)

    isValidLoad(s.load, s.Lmax) || (valid = false; (display && (println(" - Loads invalid"))))

    return valid
end

"""
```Julia
    canAddRouteToSession(r::Route{N}, s::Session) where N
```
Test if a route `r` could be added to the current session `s` without editing the route nor the session.
"""
canAddRouteToSession(r::Route{N}, s::Session) where N = sum(deepcopy(s.load) + r.assignment .> s.Lmax) == 0

## ============================================================================================================== ##
##                      ########  ########  ########  ##    ##  ########   #######   #######                      ##
##                      ##           ##        ##     ###   ##  ##        ##        ##                            ##
##                      #####        ##        ##     ## ## ##  #####      ######    ######                       ##
##                      ##           ##        ##     ##   ###  ##              ##        ##                      ##
##                      ##        ########     ##     ##    ##  ########  #######   #######                       ##
## ============================================================================================================== ##

abstract type FitnessSession end

struct HollowedPercentage   <: FitnessSession end
struct NonFilledOutputs     <: FitnessSession end
struct LoadSTD              <: FitnessSession end
struct NonNulLoadSTD        <: FitnessSession end

# Minimalistic Session
@inline specialized(loads::Vector{Int64}, Lmax::Int64,::Type{HollowedPercentage})      = (100 * ((Lmax * length(loads)) - sum(loads))) / (Lmax * length(loads))
@inline specialized(loads::Vector{Int64}, Lmax::Int64,::Type{NonFilledOutputs})        = float(count(x -> x >= Lmax, loads))
@inline specialized(loads::Vector{Int64}, Lmax::Int64,::Type{LoadSTD})                 = std(loads)
@inline specialized(loads::Vector{Int64}, Lmax::Int64,::Type{NonNulLoadSTD})           = std([e for e in loads if e != 0])

# Full Session
@inline specialized(s::Session,::Type{HollowedPercentage})  = specialized(s.load, s.Lmax, HollowedPercentage)
@inline specialized(s::Session,::Type{NonFilledOutputs})    = specialized(s.load, s.Lmax, NonFilledOutputs)
@inline specialized(s::Session,::Type{LoadSTD})             = specialized(s.load, s.Lmax, LoadSTD)
@inline specialized(s::Session,::Type{NonNulLoadSTD})       = specialized(s.load, s.Lmax, LoadSTD)


"""
```Julia
    fitness(loads::Vector{Int64}, Lmax::Int64, TAG::Type{<:FitnessSession})
    fitness(loads::Vector{Int64}, Lmax::Int64, TAG::Tuple{Vararg{DataType, N}})
    fitness(s::Session, TAG::Type{<:FitnessSession})
    fitness(s::Session, TAG::Tuple{Vararg{DataType, N}}) where N
```

Apply one or many fitness criteria to the given session `s` *resp.* load vector `loads` along with output capacity `Lmax`.
criterion includes:
 - `HollowedPercentage`: the empty space in the session (in term of mail volume).
 - `NonFilledOutputs`: the number of non full output (very optimistic).
 - `LoadSTD`: standard deviation of the output loads.
 - `NonNulLoadSTD`: standard deviation of the output loads. (don't take into account empty output)
These function will return a `Float64` or a `Tuple{Vararg{Float64, N}}` (N is the number of criterion used), structure which could easily be used along with the `isless` or `sortperm` functions (many Julia soting function rely on `isless()` if the `by` attribute is not overwritten).
"""
# apply multiple criteria on Minimalistic Session
@inline fitness(loads::Vector{Int64}, Lmax::Int64, TAG::Type{<:FitnessSession})::Float64                                       = specialized(loads, Lmax, TAG)
@inline fitness(loads::Vector{Int64}, Lmax::Int64, TAG::Tuple{DataType, DataType})::Tuple{Float64, Float64}                    = (specialized(loads, Lmax, TAG[1]), specialized(loads, Lmax, TAG[2]))
@inline fitness(loads::Vector{Int64}, Lmax::Int64, TAG::Tuple{DataType, DataType, DataType})::Tuple{Float64, Float64, Float64} = (specialized(loads, Lmax, TAG[1]), specialized(loads, Lmax, TAG[2]), specialized(loads, Lmax, TAG[3]))
@inline fitness(loads::Vector{Int64}, Lmax::Int64, TAG::Tuple{Vararg{DataType, N}}) where N                                    = ntuple((i -> specialized(loads, Lmax, TAG[i])), N)::NTuple{N, Float64}

# apply multiple criteria on full Session
@inline fitness(s::Session, TAG::Type{<:FitnessSession})::Float64                                       = specialized(s.load, s.Lmax, TAG)
@inline fitness(s::Session, TAG::Tuple{DataType, DataType})::Tuple{Float64, Float64}                    = (specialized(s.load, s.Lmax, TAG[1]), specialized(s.load, s.Lmax, TAG[2]))
@inline fitness(s::Session, TAG::Tuple{DataType, DataType, DataType})::Tuple{Float64, Float64, Float64} = (specialized(s.load, s.Lmax, TAG[1]), specialized(s.load, s.Lmax, TAG[2]), specialized(s.load, s.Lmax, TAG[3]))
@inline fitness(s::Session, TAG::Tuple{Vararg{DataType, N}}) where N                                    = ntuple((i -> specialized(s.load, s.Lmax, TAG[i])), N)::NTuple{N, Float64}

## ============================================================================================================== ##
##   ######   ########  #######   ########  ########  ########  ########   ######    ######   ########  ########  ##
##  ##    ##  ##        ##    ##     ##        ##     ##           ##     ##    ##  ##    ##     ##     ##        ##
##  ##        #####     #######      ##        ##     #####        ##     ##        ########     ##     #####     ##
##  ##    ##  ##        ##  ##       ##        ##     ##           ##     ##    ##  ##    ##     ##     ##        ##
##   ######   ########  ##   ##      ##     ########  ##        ########   ######   ##    ##     ##     ########  ##
## ============================================================================================================== ##

"""
```Julia
    certificat_CapacityVolume(s::Session, r::Route{N}) where N
```
Volume of mail certificate. Return `true` if the total volume in the session `s` and the route `r` do not exceed the sesson's capacity (computed as the number of output multiplyed by each output capacity). Otherwise return `false`.
"""
function certificat_CapacityVolume(s::Session, r::Route{N}) where N
    return (sum(s.load) + sum(r.mail) <= (s.Lmax * length(s.load)))
end

"""
```Julia
    certificat_MaxBatcheMinLoad(s::Session, r::Route{N}) where N
```
Compare the largest mail in route `r` to the least loaded output of session `s`.
Return `true` if the largest mail of the route is smaller than the remaining space in the least loaded output of the session. Otherwise return `false`.
"""
function certificat_MaxRouteNumber(s::Session, r::Route{N}) where N
    return (length(s.route) + 1 <= length(s.load))
end

"""
```Julia
    certificat_MaxBatcheMinLoad(s::Session, r::Route{N}) where N
```
Compare the largest mail in route `r` to the least loaded output of session `s`.
Return `true` if the largest mail of the route is smaller than the remaining space in the least loaded output of the session. Otherwise return `false`.
"""
function certificat_MaxBatcheMinLoad(s::Session, r::Route{N}) where N
    return (maximum(r.mail) <= (s.Lmax - minimum(s.load)))
end

"""
```Julia
    certificat_LeftAlligned(s::Session, r::Route{N}) where N
```
Left allign assignment under  the certificate format. Will only return a boolean not the actual assignment
"""
function certificat_LeftAlligned(s::Session, r::Route{N}) where N
    certificat_MaxBatcheMinLoad(s, r) || return false

    valid ::Bool = false
    O::Int64 = length(s.load)
    Br::Int64 = length(r.mail)
    currentBatch ::Int64 = 1
    for (out, load) in enumerate(s.load)
        if r.mail[currentBatch] + load <= s.Lmax
            currentBatch += 1

            ((Br - currentBatch) > (O - out)) && (valid = false; break) # not enougth outputs left to fit each mails
            (currentBatch > Br) && (valid = true; break) # no more mails to fit
        end
    end

    return valid
end

function certificate_constrainedLostSpace(s::Session, r::Route{N}) where N
    O                   ::Int64         = length(s.load)
    R                   ::Vector{Route} = [s.route; r]

    total_mail_volume   ::Int64         = sum([sum(sr.mail) for sr in R])
    total_mail_capacity ::Int64         = O * s.Lmax
    max_mail_nb         ::Int64         = maximum([length(sr.mail) for sr in R])

    println("total_mail_volume: $total_mail_volume")
    println("total_mail_capacity: $total_mail_capacity")

    lost_space = 0

    for k=1:max_mail_nb
        # max(0, k * s.Lmax - sum([sum([sr.mail[j] for j=1:length(sr.mail)]) for sr in R]))
        tmp = []
        for sr in R
            tmp2 = [sr.mail[j] for j=1:min(k, length(sr.mail))]
            push!(tmp, tmp2)
            println(tmp2)
        end
        println()

        
    end

    lost_space += maximum([max(0, k * s.Lmax - sum([sum([sr.mail[j] for j=1:min(k, length(sr.mail))]) for sr in R])) for k=1:max_mail_nb])

    return lost_space
end

## ============================================================================================================== ##
##             ######   ######    ######              #######    ######   ##    ##  ########  ########            ##
##            ##    ##  ##    ##  ##    ##            ##    ##  ##    ##  ##    ##     ##     ##                  ##
##            ########  ##     #  ##     #            #######   ##    ##  ##    ##     ##     #####               ##
##            ##    ##  ##    ##  ##    ##            ##  ##    ##    ##  ##    ##     ##     ##                  ##
##            ##    ##  ######    ######              ##   ##    ######    ######      ##     ########            ##
## ============================================================================================================== ##

"""

---

# Inserting a route into a session

Each following function will take at least a session `s` and a route `r`.
The result of each function will be a session along with a boolean.
The boolean will be set to `true` if the input session has been modified and to `false` otherwise.
The returned session is just the session given in the input and will be modified only if the insertion of the new route result in a valid session.
"""

"""
```Julia
    addRoute_RAW!(s::Session, r::Route{N}) where N
```
Add the given route `r` to the session `s` without modifying the route nor the session assignments.
> Note that `canAddRouteToSession` function test if the resulting session of such affectation will be valid (before inserting the route).
"""
function addRoute_RAW!(s::Session, r::Route{N}) where N
    newLoads::Vector{Int64} = s.load + r.assignment

    if isValidLoad(newLoads, s.Lmax)
        s.load = newLoads
        push!(s.route, r)
        return (s, true)::Tuple{Session, Bool}
    else
        return (s, false)::Tuple{Session, Bool}
    end
end

"""
```Julia
    addRoute_SmoothAssigned!(s::Session, r::Route{N}) where N
```
Smooth assigned affectation.
"""
function addRoute_SmoothAssigned!(s::Session, r::Route{N}) where N
    if isempty(s.route)                    # > Empty Session (special case)
        push!(s.route, r)                  #       add route
        s.load += r.assignment             #       init loads
        return s, true
    end

    # ==========< Init >==========
    O::Int64 = length(s.load)
    Br::Int64 = length(r.mail)

    (Br == O) && (return addRoute_RAW!(s, r)) # route is full

    valid::Bool = true
    newPos::Vector{Int64} = sort!(sortperm(s.load)[1:Br]) # O(n log(n))

    # ==========< Update >==========
    newAssignment::Vector{Int64} = zeros(Int64, O)
    i::Int64 = 1
    while i <= Br && valid
        if s.load[newPos[i]] + r.mail[i] <= s.Lmax # updated loads overflow
            newAssignment[newPos[i]] = r.mail[i]
        else
            valid = false
        end
        i += 1
    end

    # ==========< Results >==========
    if valid
        s.load += newAssignment
        push!(s.route, Route(r.id, newAssignment, r.mail))
        return (s, true)::Tuple{Session, Bool}
    else
        return (s, false)::Tuple{Session, Bool}
    end
end

"""
```Julia
    addRoute_SAVANT!(s::Session, r::Route{N}, Δ::Int64 = 2)
```
Smooth assigned affectation variant. 
Rather than considering only |Br| (number of mails in the route) least loaded outputs to assign mails, this method consider |Br| + `Δ` least loaded outputs and return the first valid assignemnent found (if found).
"""
function addRoute_SAVANT!(s::Session, r::Route{N}, Δ::Int64 = 2) where N
    if isempty(s.route)                    # > Empty Session (special case)
        push!(s.route, r)                  #       add route
        s.load += r.assignment             #       init loads
        return (s, true)
    end

    O::Int64 = length(s.load)
    Br::Int64 = length(r.mail)
    valid::Bool = true
    
    certificat_LeftAlligned(s, r) || (return (s, false))  # SAVA-NT impossible proved by Left Alligned certificate
    (Br == O) && (return addRoute_RAW!(s, r))   # route is full
    (Br + Δ > O) && (Δ = O - Br)                # valid Δ 
    
    # least loaded outputs
    ll_out::Vector{Int64} = sort!(sortperm(s.load)[1:Br + Δ]) 

    # all br sized ordered position index in ll_out
    CLC = CoolLexCombinations(Br + Δ, Br)

    # all br sized ordered position subset in ll_out
    all_newPos = [[ll_out[i] for i in e] for e in CLC]

    for newPos in all_newPos
        valid = true
        newLoads::Vector{Int64} = deepcopy(s.load)
        newAssignment::Vector{Int64} = zeros(Int64, O)

        # position the mails to their new location
        for i=1:Br
            newAssignment[newPos[i]] = r.mail[i]
            newLoads[newPos[i]] += r.mail[i]
    
            (newLoads[newPos[i]] <= s.Lmax) || (valid = false; break) # updated loads overflow
        end

        if valid # route fit in the session
            s.load = newLoads
            push!(s.route, Route(r.id, newAssignment, r.mail))
            return (s, true)::Tuple{Session, Bool}
        end
    end
    return (s, false)::Tuple{Session, Bool}
end

"""
```Julia
    addRoute_SAVANT_MINSTD!(s::Session, r::Route{N}, Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)
```
Another smooth assigned affectation variant.
Also considering the |Br| + `Δ` least loaded output of the session, this method will return the best value among all valid assignemnent regarding `LoadSTD`.
"""
function addRoute_SAVANT_MINSTD!(s::Session, r::Route{N}, Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD) where N
    if isempty(s.route)                    # > Empty Session (special case)
        push!(s.route, r)                  #       add route
        s.load += r.assignment             #       init loads
        return (s, true)
    end

    O::Int64 = length(s.load)
    Br::Int64 = length(r.mail)
    valid::Bool = true
    
    certificat_LeftAlligned(s, r) || (return s, false)  # SAVA-NT impossible proved by Left Alligned certificate
    (Br == O) && (return addRoute_RAW!(s, r))   # route is full
    (Br + Δ > O) && (Δ = O - Br)                # valid Δ 
    
    # least loaded outputs
    ll_out::Vector{Int64} = sort!(sortperm(s.load)[1:Br + Δ]) 

    # all br sized ordered position index in ll_out
    CLC = CoolLexCombinations(Br + Δ, Br)

    # all br sized ordered position subset in ll_out
    all_newPos = [[ll_out[i] for i in e] for e in CLC]

    bestLoads::Union{Nothing, Vector{Int64}} = nothing
    bestAssignment::Union{Nothing, Vector{Int64}} = nothing

    for newPos in all_newPos
        valid = true
        newLoads::Vector{Int64} = deepcopy(s.load)
        newAssignment::Vector{Int64} = zeros(Int64, O)

        # position the mails to their new location
        for i=1:Br
            newAssignment[newPos[i]] = r.mail[i]
            newLoads[newPos[i]] += r.mail[i]
    
            (newLoads[newPos[i]] <= s.Lmax) || (valid = false; break) # updated loads overflow
        end

        if valid && (bestLoads === nothing || fitness(bestLoads, s.Lmax, TAG_FitSes) > fitness(newLoads, s.Lmax, TAG_FitSes))
            bestLoads = newLoads
            bestAssignment = newAssignment
        end
    end

    if bestLoads === nothing
        return (s, false)::Tuple{Session, Bool}
    else
        s.load = bestLoads
        push!(s.route, Route(r.id, bestAssignment, r.mail))
        return (s, true)::Tuple{Session, Bool}
    end
end

"""
```Julia
    addRoute_LeftAlligned(s::Session, r::Route{N})  # wont edit s, return a boolean along with the new assignment
    addRoute_LeftAlligned!(s::Session, r::Route{N}) # may edit s
```
Left-alligned algorithm.
"""
function addRoute_LeftAlligned(s::Session, r::Route{N}) where N
    j               ::Int64         = 1
    k               ::Int64         = 1
    O               ::Int64         = length(s.load)
    Br              ::Int64         = length(r.mail)
    valid           ::Bool          = true
    newAssignment   ::Vector{Int64} = zeros(Int64, O)

    while valid
        if r.mail[j] + s.load[k] <= s.Lmax
            newAssignment[k] = r.mail[j]
            j += 1
        end
        k += 1
        (j > Br) && (break)
        (Br - j > O - k) && (valid = false)
    end
    return (valid, newAssignment)::Tuple{Bool, Vector{Int64}}
end

function addRoute_LeftAlligned!(s::Session, r::Route{N}) where N
    j               ::Int64         = 1
    k               ::Int64         = 1
    O               ::Int64         = length(s.load)
    Br              ::Int64         = length(r.mail)
    valid           ::Bool          = true
    newAssignment   ::Vector{Int64} = zeros(Int64, O)

    while valid
        if r.mail[j] + s.load[k] <= s.Lmax
            newAssignment[k] = r.mail[j]
            j += 1
        end
        k += 1
        (j > Br) && (break)
        (Br - j > O - k) && (valid = false)
    end

    if valid
        push!(s.route, r)
        s.load += newAssignment
        r.assignment = newAssignment
    end

    return (valid, newAssignment)::Tuple{Bool, Vector{Int64}}
end

"""
```Julia
    addRoute_SUB_01LP!(s::Session, r::Route{N}, tl::Int64 = 100, env::Gurobi.Env = Gurobi.Env())
```
Matheuristic reliying on a model to found an intresting mail assignment in the added route.

**Binary variables:**
 - ```Latex x \\in \\{0, 1\\} ```
"""
function addRoute_SUB_01LP!(s::Session, r::Route{N}, tl::Int64 = 100, env::Gurobi.Env = Gurobi.Env()) where N
# ===========< Data:                    OK >===========
    O::Int64 = length(r.assignment)

    # compute certificat
    flag::Bool, newAssignment::Vector{Int64} = addRoute_LeftAlligned(s, r)
    flag || (return (s, false)::Tuple{Session, Bool}) # print("<!>");   
    newRoute::Route{N} = Route(r.id, newAssignment, r.mail)

    # Br::Vector{Int64} = collect(1:N)

    # vrj = r.mail 

    Orj::Vector{UnitRange{Int64}} = [(j: O - N + j) for j=1:N]

    # Lmax = s.Lmax

# ===========< Model:                   OK >===========

    model = Model(() -> Gurobi.Optimizer(env)) # () -> GLPK.Optimizer()
    set_silent(model)
    # set_time_limit_sec(model, tl)
    set_optimizer_attribute(model, "TimeLimit", tl)
    set_optimizer_attribute(model, "OutputFlag", 0)

# ===========< Variables:               OK >===========

    @variable(model, x[j=1:N, k=Orj[j]], Bin)

# ===========< Objective:               OK >===========

    @objective(model, Min, sum([sum([((r.mail[j] + s.load[k]) * x[j, k]) for k=Orj[j]]) for j=1:N]))
    
# ===========< Constraint: (4, 7, 10)   OK >===========

    # (4) -> require that each batch j of route r be assigned to exactly one output
    for j=1:N
        @constraint(model, sum([x[j, k] for k=Orj[j]]) == 1)
    end

    # (7) -> precedence mail constraints for each route r
    for j=1:N
        if j != N
            @constraint(model, 1 + sum([k * x[j, k] for k=Orj[j]]) <= sum([k * x[j+1, k] for k=Orj[j + 1]]))
        end
    end

    # (10) -> capacity constraint
    for j=1:N
        @constraint(model, sum([((r.mail[j] + s.load[k]) * x[j, k]) for k=Orj[j]]) <= s.Lmax)
    end

# ===========< Warmup:                  OK >===========

    mail::Int64 = 1
    for k=1:O
        if newRoute.assignment[k] != 0
            set_start_value(x[mail, k], 1)
            mail += 1
            
            (mail > N) && break
        end
    end

# ===========< Results:                 OK >===========
    optimize!(model)

    if termination_status(model) == OPTIMAL || MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0
        newRoute.assignment = zeros(Int64, O)
        for j=1:N
            newPos::Int64 = 0
            for k=Orj[j] 
                (value(x[j,k]) != 0.0) && (newPos = k; break) 
            end
            newRoute.assignment[newPos] = newRoute.mail[j]
        end

        push!(s.route, newRoute)

        s.load += newRoute.assignment
        # termination_status(model) == OPTIMAL ? print("<o>") : print("<0>")
        return s, true
    else
        # print("<x>")
        push!(s.route, newRoute) 
        return s, true # should never trigger but better safe tha sorry (memory overflow)
    end
end

"""
```Julia
    addRoute_PART_SUB_01LP!(s::Session, r::Route{N}, tl::Int64 = 100, env::Gurobi.Env = Gurobi.Env())
```
Matheuristic reliying on a model to found an intresting mail assignment in the added route. 
This model won't enforce the output capacity constraint in the first place hopping that thanks to the objective function no problem will come up.
If the resulting solution isn't valid (due to a capacity overflow), the capacity constraint of each problematic output will be added to the model.
The updated model will be solved again and this process will be repeted again and again until a valid solution is found. 
"""
function addRoute_PART_SUB_01LP!(s::Session, r::Route{N}, tl::Int64 = 100, env::Gurobi.Env = Gurobi.Env()) where N
# ===========< Data:                    OK >===========
    O::Int64 = length(r.assignment)

    # compute certificat
    flag::Bool, newAssignment::Vector{Int64} = addRoute_LeftAlligned(s, r)
    flag || (print("<!>"), return (s, false)::Tuple{Session, Bool}) # print("<c>"), 
    newRoute::Route{N} = Route(r.id, newAssignment, r.mail)

    # Br::Vector{Int64} = collect(1:N)

    # vrj = r.mail 

    Orj::Vector{UnitRange{Int64}} = [(j: O - N + j) for j=1:N]

    # Lmax = s.Lmax

# ===========< Model:                   OK >===========

    model = Model(() -> Gurobi.Optimizer(env))
    set_silent(model)
    set_optimizer_attribute(model, "TimeLimit", tl)
    set_optimizer_attribute(model, "OutputFlag", 0)

# ===========< Variables:               OK >===========

    @variable(model, x[j=1:N, k=Orj[j]], Bin)

# ===========< Objective:               OK >===========

    # minimize the output loads in which a mail is affected
    @objective(model, Min, sum([sum([((r.mail[j] + s.load[k]) * x[j, k]) for k=Orj[j]]) for j=1:N]))
    
# ===========< Constraint: (ctr 4, 7)   OK >===========

    # (4) -> require that each batch j of route r be assigned to exactly one output
    for j=1:N
        @constraint(model, sum([x[j, k] for k=Orj[j]]) == 1)
    end

    # (7) -> precedence mail constraints for each route r
    for j=1:N
        if j != N
            @constraint(model, 1 + sum([k * x[j, k] for k=Orj[j]]) <= sum([k * x[j+1, k] for k=Orj[j + 1]]))
        end
    end

# ===========< Warmup:                  OK >===========

    mail::Int64 = 1
    for k=1:O
        if newRoute.assignment[k] != 0
            set_start_value(x[mail, k], 1)
            mail += 1
            
            (mail > N) && break
        end
    end

# ===========< Results:                 OK >===========
    optimize!(model)

    if termination_status(model) == OPTIMAL || MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0
        newRoute.assignment = zeros(Int64, O)
        for j=1:N
            newPos::Int64 = 0
            for k=Orj[j] 
                (value(x[j,k]) != 0.0) && (newPos = k; break) 
            end
            newRoute.assignment[newPos] = newRoute.mail[j]
        end

        newLoad::Vector{Int64} = s.load + newRoute.assignment

        if isValidLoad(newLoad, s.Lmax)
            termination_status(model) == OPTIMAL ? print("<o>") : print("<0>")
            push!(s.route, newRoute)
            s.load += newRoute.assignment
            return s, true
        else
# ===========< Re-Optimize: (+ ctr 10)  OK >===========
            print("<x")
            for j=1:N
                @constraint(model, sum([((r.mail[j] + s.load[k]) * x[j, k]) for k=Orj[j]]) <= s.Lmax)
            end
            optimize!(model)

            if termination_status(model) == OPTIMAL || MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0
                newRoute.assignment = zeros(Int64, O)
                for j=1:N
                    newPos::Int64 = 0
                    for k=Orj[j] 
                        (value(x[j,k]) != 0.0) && (newPos = k; break) 
                    end
                    newRoute.assignment[newPos] = newRoute.mail[j]
                end

                newLoad = s.load + newRoute.assignment

                if isValidLoad(newLoad, s.Lmax)
                    termination_status(model) == OPTIMAL ? print("o>") : print("0>")
                    push!(s.route, newRoute)
                    s.load += newRoute.assignment
                    return s, true
                else
                    println("x>")
                    return s, false
                end
            end
        end
    else
        return s, true # should never trigger but better safe tha sorry (memory overflow)
    end
end

function addRoute_Rebuild!(s::Session, r::Route{N})::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    (global call += 1)

    ns::Session = Session(s.Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
    push!(ns.route, Route(r.id, deepcopy(r.assignment), r.mail))

    ns, flag = rebuildSession(ns)
    
    (flag) && (global repairedBuild += 1; s.route = ns.route; s.load = ns.load)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRoute_Rebuild_Knapsack!(s::Session, r::Route{N})::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    (global call += 1)

    ns::Session = Session(s.Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
    push!(ns.route, Route(r.id, deepcopy(r.assignment), r.mail))

    ns, flag = rebuildSession_knapSack(ns)
    
    (flag) && (global repairedBuild += 1; s.route = ns.route; s.load = ns.load)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRoute_Rebuild_Knapsack_model!(s::Session, r::Route{N}, env::Gurobi.Env = Gurobi.Env())::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    (global call += 1)

    ns::Session = Session(s.Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
    push!(ns.route, Route(r.id, deepcopy(r.assignment), r.mail))

    ns, flag = rebuildSession_knapSack_model(ns, env)
    
    (flag) && (global repairedBuild += 1; s.route = ns.route; s.load = ns.load)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRoute_Rebuild_Knapsack_model_V2!(s::Session, r::Route{N}, tl::Int64=10, env::Gurobi.Env = Gurobi.Env())::Tuple{Session, Bool} where N
    (global call += 1)

    ns::Session = Session(s.Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
    push!(ns.route, Route(r.id, deepcopy(r.assignment), r.mail))

    ns, flag = rebuildSession_knapSack_model_V2(ns, tl, env)
    
    (flag) && (global repairedBuild += 1; s.route = ns.route; s.load = ns.load)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRoute_Rebuild_Knapsack_model_V3!(s::Session, r::Route{N}, tl::Int64=10, env::Gurobi.Env = Gurobi.Env())::Tuple{Session, Bool} where N
    (global call += 1)

    ns::Session = Session(s.Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
    push!(ns.route, Route(r.id, deepcopy(r.assignment), r.mail))

    ns, flag = rebuildSession_knapSack_model_V3!(ns, tl, env)
    
    (flag) && (global repairedBuild += 1; s.route = ns.route; s.load = ns.load)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRoute_OPTIMOVE_stage1!(s::Session, r::Route{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)

    ns::Session = Session(s.Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
    push!(ns.route, Route(r.id, deepcopy(r.assignment), r.mail))

    ns, nsFit, flag = improvedOptiMove_S1_V2!(ns, TAG_FitSes)
    
    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.route = ns.route; s.load = ns.load)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRoute_OPTIMOVE_Stage1Valid!(s::Session, r::Route{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)

    ns::Session = Session(s.Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
    push!(ns.route, Route(r.id, deepcopy(r.assignment), r.mail))

    ns, nsFit, flag = improvedOptiMove_S1_V3!(ns, TAG_FitSes)
    
    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.route = ns.route; s.load = ns.load)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRoute_OPTIMOVE_Stage1ValidInf!(s::Session, r::Route{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)
    
    ns::Session = Session(s.Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
    push!(ns.route, Route(r.id, deepcopy(r.assignment), r.mail))

    ns, nsFit, flag = improvedOptiMove_S1_V4!(ns, TAG_FitSes)
    
    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.route = ns.route; s.load = ns.load)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRoute_OPTIMOVE_oneMove!(s::Session, r::Route{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)

    ns::Session = Session(s.Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
    push!(ns.route, Route(r.id, deepcopy(r.assignment), r.mail))

    ns, nsFit, flag = improvedOptiMove_S1_V5!(ns, TAG_FitSes)
    
    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.route = ns.route; s.load = ns.load)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRoute_OPTIMOVE_3Stages!(s::Session, r::Route{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)

    ns::Session = Session(s.Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
    push!(ns.route, Route(r.id, deepcopy(r.assignment), r.mail))

    ns, nsFit, flag = improvedOptiMove_V1(ns, TAG_FitSes)
    
    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.route = ns.route; s.load = ns.load)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRoute_OPTIMOVE_Stage1NegSTD!(s::Session, r::Route{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    # println("\n\n$s\n")

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)

    ns::Session = Session(s.Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
    push!(ns.route, Route(r.id, deepcopy(r.assignment), r.mail))

    ns, nsFit, flag = improvedOptiMove_S1_V6!(ns, TAG_FitSes)
    
    # println("\n\n$ns\n")

    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.route = ns.route; s.load = ns.load)
    return (ns, flag)::Tuple{Session, Bool}
end

abstract type SimpleAddRoute end

# MOVE 1 ROUND
struct RAW                  <: SimpleAddRoute end
struct SMOOTH_ASSIGNED      <: SimpleAddRoute end
struct SAVANT               <: SimpleAddRoute end
struct SAVANT_MIN_STD       <: SimpleAddRoute end
struct LEFT_ALLIGNED        <: SimpleAddRoute end
struct SUB_01LP             <: SimpleAddRoute end
struct PART_SUB_01LP        <: SimpleAddRoute end

# MOVE SESSION
struct OPTIMOVE_STAGE1      <: SimpleAddRoute end
struct OPTIMOVE_VALID       <: SimpleAddRoute end
struct OPTIMOVE_VALID_INF   <: SimpleAddRoute end
struct OPTIMOVE_1BATCH      <: SimpleAddRoute end
struct OPTIMOVE_3S          <: SimpleAddRoute end
struct OPTIMOVE_S1_NEGSTD   <: SimpleAddRoute end
struct REBUILD_KP           <: SimpleAddRoute end
struct REBUILD_KP_01LP      <: SimpleAddRoute end
struct REBUILD_KP_01LP_V2   <: SimpleAddRoute end

# MOVE 1 ROUND
@inline addRoute(s::Session, r::Route{N}, ::Type{RAW}                 , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_RAW!(s, r)
@inline addRoute(s::Session, r::Route{N}, ::Type{SMOOTH_ASSIGNED}                , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_SmoothAssigned!(s, r)
@inline addRoute(s::Session, r::Route{N}, ::Type{LEFT_ALLIGNED}                , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_LeftAlligned!(s, r)
@inline addRoute(s::Session, r::Route{N}, ::Type{SAVANT}              , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_SAVANT!(s, r, Δ)
@inline addRoute(s::Session, r::Route{N}, ::Type{SAVANT_MIN_STD}      , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_SAVANT_MINSTD!(s, r, Δ, TAG_FitSes)
@inline addRoute(s::Session, r::Route{N}, ::Type{SUB_01LP}            , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_SUB_01LP!(s, r, tl, env)
@inline addRoute(s::Session, r::Route{N}, ::Type{PART_SUB_01LP}       , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_PART_SUB_01LP!(s, r, tl, env)
@inline addRoute(s::Session, r::Route{N}, ::Type{REBUILD_KP}          , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_Rebuild_Knapsack!(s, r)
@inline addRoute(s::Session, r::Route{N}, ::Type{REBUILD_KP_01LP}     , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_Rebuild_Knapsack_model!(s, r, env)
@inline addRoute(s::Session, r::Route{N}, ::Type{REBUILD_KP_01LP_V2}  , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_Rebuild_Knapsack_model_V2!(s, r, tl, env)

# MOVE SESSION
@inline addRoute(s::Session, r::Route{N}, ::Type{OPTIMOVE_STAGE1}     , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_OPTIMOVE_stage1!(s, r, TAG_FitSes)
@inline addRoute(s::Session, r::Route{N}, ::Type{OPTIMOVE_VALID}      , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_OPTIMOVE_Stage1Valid!(s, r, TAG_FitSes)
@inline addRoute(s::Session, r::Route{N}, ::Type{OPTIMOVE_VALID_INF}  , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_OPTIMOVE_Stage1ValidInf!(s, r, TAG_FitSes)
@inline addRoute(s::Session, r::Route{N}, ::Type{OPTIMOVE_1BATCH}     , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_OPTIMOVE_oneMove!(s, r, TAG_FitSes)
@inline addRoute(s::Session, r::Route{N}, ::Type{OPTIMOVE_3S}         , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_OPTIMOVE_3Stages!(s, r, TAG_FitSes)
@inline addRoute(s::Session, r::Route{N}, ::Type{OPTIMOVE_S1_NEGSTD}  , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRoute_OPTIMOVE_Stage1NegSTD!(s, r, TAG_FitSes)


## ============================================================================================================== ##
##                      ######    ########   #######  #######   ##         ######   ##    ##                      ##
##                      ##    ##     ##     ##        ##    ##  ##        ##    ##   ##  ##                       ##
##                      ##     #     ##      ######   #######   ##        ########     ##                         ##
##                      ##    ##     ##           ##  ##        ##        ##    ##     ##                         ##
##                      ######    ########  #######   ##        ########  ##    ##     ##                         ##
## ============================================================================================================== ##

function Base.show(io::IO, s::Session)
    println(io, "Session: (Valid loads: $(isValidLoad(s))) (Valid session: $(isSessionValid(s))) (filled at $(round(100 - fitness(s, HollowedPercentage), digits=3))%) (Standard deviation: $(round(fitness(s, LoadSTD), digits=3)))\n    Lmax: $(s.Lmax)\n    Routes: (x$(length(s.route)))")
    if length(s.load) <= 40
        for r in s.route
            println(io, "$(r)")
        end
    else
        for (i, r) in enumerate(s.route)
            i%10 == 0 ? println(io, "$(r)") : print(io, "$(r)")
        end
    end
    println(io, "\n    Loads:\n$(s.load)")
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
    Random.shuffle(s::Session)::Session     # wont modify the session
    Random.shuffle!(s::Session)::Session    # will modify the session
```
shuffle the mail position within the session. (will keep the order among mails)
"""
function Random.shuffle!(s::Session)::Session
    for r in s.route
        shuffle!(r)
    end

    return compute_output!(s)
end

function Random.shuffle(s::Session)::Session
    return shuffle!(deepcopy(s))
end

"""
```Julia
    randSession(Lmax::Int64=25, O::Int64=20, nbRoutes::Int64=10)
```
Create a random session. `Lmax` is the maximum capacity of each output. `O` is the number of output of the session and `nbRoutes` is the amount of random routes that the session will bore.
 > :warning: note that this session may be invalid with the wrong set of parameter as most of it is completely random and no safety procedure are implemented.
"""
function randSession(Lmax::Int64=25, O::Int64=20, nbRoutes::Int64=10)
    routes::Vector{Route} = [randRoute(i, O, nothing, Lmax) for i=1:nbRoutes]

    return Session(Lmax, routes, compute_output(routes))
end
