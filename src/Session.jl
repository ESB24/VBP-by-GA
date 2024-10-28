begin
    using Combinatorics
    using Gurobi
    using GLPK
    using JuMP
end

begin 
    include("Round.jl")
    include("DataStruct.jl")
    include("SessionRebuild.jl")
end

# ================================================== #
#                    Constructor                     #
# ================================================== #

"""
```Julia
    Session(C::Int64, O::Int64)::Sesssion
    Session(C::Int64, r::Round{N}) where N
    Session(C::Int64, r::Vector{Round{N}}) where N
```
Constructor for the Session data structure.
`C` is the maximum capacity per output (Lmax). 
`O` is the length of the assignment vector (the session will start empty). 
Otherwise and replacing `O` to create non empty session you could provide a unique or a vector of route `r`. 

"""

function Session(C::Int64, O::Int64)::Session
    return Session(C, [], zeros(Int64, O))
end

function Session(C::Int64, r::Round{N}) where N
    return Session(C, [r], r.assignment)
end

function Session(C::Int64, r::Vector{Round{N}}) where N
    return Session(C, r, compute_output(r))
end

"""
```Julia
    compute_output(rounds::Vector{Union{Round{N}, Round}})
    compute_output(s::Session)::Vector{Int64}
    compute_output!(s::Session)::Session
```

Compute the output of a given session or list of route. (according to there current mail assignment)
If `compute_output` is called the load vector is computed and returned.
If `compute_output!` then the ouput load is also updated in the given session `s`.
"""

function compute_output(rounds::Vector{Union{Round{N}, Round}}) where N
    O::Int64 = length(rounds[1].assignment)
    R::Int64 = length(rounds)

    return [sum([rounds[i].assignment[k] for i=1:R]) for k=1:O]
end

function compute_output(s::Session)::Vector{Int64}
    return isempty(s.rounds) ? zeros(Int64, length(s.loads)) : compute_output(s.rounds)
end

function compute_output!(s::Session)::Session
    s.loads = compute_output(s)
    return s
end

# ================================================== #
#                         Test                       #
# ================================================== #

"""
```Julia
    validLoad(loads::Vector{Int64}, C::Int64)::Bool
    validLoad(s::Session)::Bool
```
Test the validity of a given load vector.
Either with `loads` the load vector and `C` the maximum capacity of each output or with the given session `s`. 
"""
validLoad(loads::Vector{Int64}, C::Int64)::Bool = sum(loads .> C) == 0
validLoad(s::Session)::Bool = validLoad(s.loads, s.C)

"""
```Julia
    validSession(s::Session, display::Bool = true)
```
Test the validity of a session.
Take into accound each route individual validity, and the capacity constraint on each output.
"""
function validSession(s::Session, display::Bool = true)
    valid::Bool = true

    for r in s.rounds
        validRound(r) || (valid = false; (display || (println("-> Round $(r.id) -> invalid"))))
    end

    compute_output!(s)

    validLoad(s.loads, s.C) || (valid = false; (display || (println("-> Loads invalid"))))

    return valid
end

"""
```Julia
    canAddRoundToSession(r::Round{N}, s::Session) where N
```
Test if a route `r` could be added to the current session `s` without editing the route nor the session.
"""
canAddRoundToSession(r::Round{N}, s::Session) where N = sum(deepcopy(s.loads) + r.assignment .> s.C) == 0

# ================================================== #
#                   Fitness Function                 #
# ================================================== #

abstract type FitnessSession end

struct HollowedPercentage   <: FitnessSession end
struct NonFilledOutputs     <: FitnessSession end
struct LoadSTD              <: FitnessSession end
struct NonNulLoadSTD        <: FitnessSession end

# Minimalistic Session
@inline specialized(loads::Vector{Int64}, C::Int64,::Type{HollowedPercentage})      = (100 * ((C * length(loads)) - sum(loads))) / (C * length(loads))
@inline specialized(loads::Vector{Int64}, C::Int64,::Type{NonFilledOutputs})        = count(x -> x >= C, loads)
@inline specialized(loads::Vector{Int64}, C::Int64,::Type{LoadSTD})                 = std(loads)
@inline specialized(loads::Vector{Int64}, C::Int64,::Type{NonNulLoadSTD})           = std([e for e in loads if e != 0])

# Full Session
@inline specialized(s::Session,::Type{HollowedPercentage})  = specialized(s.loads, s.C, HollowedPercentage)
@inline specialized(s::Session,::Type{NonFilledOutputs})    = specialized(s.loads, s.C, NonFilledOutputs)
@inline specialized(s::Session,::Type{LoadSTD})             = specialized(s.loads, s.C, LoadSTD)
@inline specialized(s::Session,::Type{NonNulLoadSTD})       = specialized(s.loads, s.C, LoadSTD)


"""
```Julia
    fitness(loads::Vector{Int64}, C::Int64, TAG::Type{<:FitnessSession})
    fitness(loads::Vector{Int64}, C::Int64, TAG::Tuple{Vararg{DataType, N}})
    fitness(s::Session, TAG::Type{<:FitnessSession})
    fitness(s::Session, TAG::Tuple{Vararg{DataType, N}}) where N
```

Apply one or many fitness criteria to the given session `s` *resp.* load vector `loads` along with output capacity `C`.
criterion includes:
 - `HollowedPercentage`: the empty space in the session (in term of mail volume).
 - `NonFilledOutputs`: the number of non full output (very optimistic).
 - `LoadSTD`: standard deviation of the output loads.
 - `NonNulLoadSTD`: standard deviation of the output loads. (don't take into account empty output)
These function will return a `Float64` or a `Tuple{Vararg{Float64, N}}` (N is the number of criterion), structure which could easily be used along with the `isless` or `sortperm` functions (many Julia soting function rely on `isless()` if the `by` attribute is not overwritten).
"""
# apply multiple criteria on Minimalistic Session
@inline fitness(loads::Vector{Int64}, C::Int64, TAG::Type{<:FitnessSession})::Float64                                       = specialized(loads, C, TAG)
@inline fitness(loads::Vector{Int64}, C::Int64, TAG::Tuple{DataType, DataType})::Tuple{Float64, Float64}                    = (specialized(loads, C, TAG[1]), specialized(loads, C, TAG[2]))
@inline fitness(loads::Vector{Int64}, C::Int64, TAG::Tuple{DataType, DataType, DataType})::Tuple{Float64, Float64, Float64} = (specialized(loads, C, TAG[1]), specialized(loads, C, TAG[2]), specialized(loads, C, TAG[3]))
@inline fitness(loads::Vector{Int64}, C::Int64, TAG::Tuple{Vararg{DataType, N}}) where N                                    = ntuple((i -> specialized(loads, C, TAG[i])), N)::NTuple{N, Float64}

# apply multiple criteria on full Session
@inline fitness(s::Session, TAG::Type{<:FitnessSession})::Float64                                       = specialized(s.loads, s.C, TAG)
@inline fitness(s::Session, TAG::Tuple{DataType, DataType})::Tuple{Float64, Float64}                    = (specialized(s.loads, s.C, TAG[1]), specialized(s.loads, s.C, TAG[2]))
@inline fitness(s::Session, TAG::Tuple{DataType, DataType, DataType})::Tuple{Float64, Float64, Float64} = (specialized(s.loads, s.C, TAG[1]), specialized(s.loads, s.C, TAG[2]), specialized(s.loads, s.C, TAG[3]))
@inline fitness(s::Session, TAG::Tuple{Vararg{DataType, N}}) where N                                    = ntuple((i -> specialized(s.loads, s.C, TAG[i])), N)::NTuple{N, Float64}

function compare_1Criteria(
        s1              ::Session, 
        s2              ::Session; 
        TAG_FitSes      ::Type{<:FitnessSession}    = LoadSTD,
        τ               ::Float64                   = 0.0,
    )::Bool

    return fitness(s1, TAG_FitSes) < fitness(s2, TAG_FitSes) + τ
end

function compare_NCriteria(
        s1              ::Session, 
        s2              ::Session; 
        TAG_FitSes      ::Vector{Type{<:FitnessSession}}        = Type{<:FitnessSession}[LoadSTD, HollowedPercentage, NonFilledOutputs],
        τ               ::Float64       = 0.0,
    )::Bool

    if isempty(TAG_FitSes) 
        return false 
    else
        val_s1 = fitness(s1, TAG_FitSes[1])
        val_s2 = fitness(s2, TAG_FitSes[1]) + τ

        val_s1 == val_s2 ? compare_NCriteria(s1, s2, TAG_FitSes=TAG_FitSes[2:end], τ=τ) : val_s1 < val_s2
    end
end

# ================================================== #
#                      Certificat                    #
# ================================================== #

"""
```Julia
    certificat_CapacityVolume(s::Session, r::Round{N}) where N
```
Volume of mail certificate. Return `true` if the total volume in the session `s` and the route `r` do not exceed the sesson's capacity (computed as the number of output multiplyed by each output capacity). Otherwise return `false`.
"""
function certificat_CapacityVolume(s::Session, r::Round{N}) where N
    return (sum(s.loads) + sum(r.batches) <= (s.C * length(s.loads)))
end

"""
```Julia
    certificat_MaxBatcheMinLoad(s::Session, r::Round{N}) where N
```
Compare the largest mail in route `r` to the least loaded output of session `s`.
Return `true` if the largest mail of the route is smaller than the remaining space in the least loaded output of the session. Otherwise return `false`.
"""
function certificat_MaxRoundNumber(s::Session, r::Round{N}) where N
    return (length(s.rounds) + 1 <= length(s.loads))
end

"""
```Julia
    certificat_MaxBatcheMinLoad(s::Session, r::Round{N}) where N
```
Compare the largest mail in route `r` to the least loaded output of session `s`.
Return `true` if the largest mail of the route is smaller than the remaining space in the least loaded output of the session. Otherwise return `false`.
"""
function certificat_MaxBatcheMinLoad(s::Session, r::Round{N}) where N
    return (maximum(r.batches) <= (s.C - minimum(s.loads)))
end

"""
```Julia
    certificat_NFBA(s::Session, r::Round{N}) where N
```
Left allign assignment under  the certificate format. Will only return a boolean not the actual assignment
"""
function certificat_NFBA(s::Session, r::Round{N}) where N
    certificat_MaxBatcheMinLoad(s, r) || return false

    valid ::Bool = false
    O::Int64 = length(s.loads)
    Br::Int64 = length(r.batches)
    currentBatch ::Int64 = 1
    for (out, load) in enumerate(s.loads)
        if r.batches[currentBatch] + load <= s.C
            currentBatch += 1

            ((Br - currentBatch) > (O - out)) && (valid = false; break) # not enougth outputs left to fit each batches
            (currentBatch > Br) && (valid = true; break) # no more batches to fit
        end
    end

    return valid
end

# ================================================== #
#                      Add Round                     #
# ================================================== #

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
    addRound_RAW!(s::Session, r::Round{N}) where N
```
Add the given route `r` to the session `s` without modifying the route nor the session assignments.
> Note that `canAddRoundToSession` function test if the resulting session of such affectation will be valid (before inserting the route).
"""
function addRound_RAW!(s::Session, r::Round{N}) where N
    newLoads::Vector{Int64} = s.loads + r.assignment

    if validLoad(newLoads, s.C)
        s.loads = newLoads
        push!(s.rounds, r)
        return (s, true)::Tuple{Session, Bool}
    else
        return (s, false)::Tuple{Session, Bool}
    end
end

"""
```Julia
    addRound_SAVA!(s::Session, r::Round{N}) where N
```
Smooth assigned affectation.
"""
function addRound_SAVA!(s::Session, r::Round{N}) where N
    if isempty(s.rounds)                    # > Empty Session (special case)
        push!(s.rounds, r)                  #       add round
        s.loads += r.assignment             #       init loads
        return s, true
    end

    # ==========< Init >==========
    O::Int64 = length(s.loads)
    Br::Int64 = length(r.batches)

    (Br == O) && (return addRound_RAW!(s, r)) # round is full

    valid::Bool = true
    newPos::Vector{Int64} = sort!(sortperm(s.loads)[1:Br]) # O(n log(n))

    # ==========< Update >==========
    newAssignment::Vector{Int64} = zeros(Int64, O)
    i::Int64 = 1
    while i <= Br && valid
        if s.loads[newPos[i]] + r.batches[i] <= s.C # updated loads overflow
            newAssignment[newPos[i]] = r.batches[i]
        else
            valid = false
        end
        i += 1
    end

    # ==========< Results >==========
    if valid
        s.loads += newAssignment
        push!(s.rounds, Round(r.id, newAssignment, r.batches))
        return (s, true)::Tuple{Session, Bool}
    else
        return (s, false)::Tuple{Session, Bool}
    end
end

"""
```Julia
    addRound_SAVANT!(s::Session, r::Round{N}, Δ::Int64 = 2)
```
Smooth assigned affectation variant. 
Rather than considering only |Br| (number of mails in the route) least loaded outputs to assign mails, this method consider |Br| + `Δ` least loaded outputs and return the first valid assignemnent found (if found).
"""
function addRound_SAVANT!(s::Session, r::Round{N}, Δ::Int64 = 2) where N
    if isempty(s.rounds)                    # > Empty Session (special case)
        push!(s.rounds, r)                  #       add round
        s.loads += r.assignment             #       init loads
        return (s, true)
    end

    O::Int64 = length(s.loads)
    Br::Int64 = length(r.batches)
    valid::Bool = true
    
    certificat_NFBA(s, r) || (return (s, false))  # SAVA-NT impossible proved by NFBA certificate
    (Br == O) && (return addRound_RAW!(s, r))   # round is full
    (Br + Δ > O) && (Δ = O - Br)                # valid Δ 
    
    # least loaded outputs
    ll_out::Vector{Int64} = sort!(sortperm(s.loads)[1:Br + Δ]) 

    # all br sized ordered position index in ll_out
    CLC = CoolLexCombinations(Br + Δ, Br)

    # all br sized ordered position subset in ll_out
    all_newPos = [[ll_out[i] for i in e] for e in CLC]

    for newPos in all_newPos
        valid = true
        newLoads::Vector{Int64} = deepcopy(s.loads)
        newAssignment::Vector{Int64} = zeros(Int64, O)

        # position the batches to their new location
        for i=1:Br
            newAssignment[newPos[i]] = r.batches[i]
            newLoads[newPos[i]] += r.batches[i]
    
            (newLoads[newPos[i]] <= s.C) || (valid = false; break) # updated loads overflow
        end

        if valid # round fit in the session
            s.loads = newLoads
            push!(s.rounds, Round(r.id, newAssignment, r.batches))
            return (s, true)::Tuple{Session, Bool}
        end
    end
    return (s, false)::Tuple{Session, Bool}
end

"""
```Julia
    addRound_SAVANT_MINSTD!(s::Session, r::Round{N}, Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)
```
Another smooth assigned affectation variant.
Also considering the |Br| + `Δ` least loaded output of the session, this method will return the best value among all valid assignemnent regarding `LoadSTD`.
"""
function addRound_SAVANT_MINSTD!(s::Session, r::Round{N}, Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD) where N
    if isempty(s.rounds)                    # > Empty Session (special case)
        push!(s.rounds, r)                  #       add round
        s.loads += r.assignment             #       init loads
        return (s, true)
    end

    O::Int64 = length(s.loads)
    Br::Int64 = length(r.batches)
    valid::Bool = true
    
    certificat_NFBA(s, r) || (return s, false)  # SAVA-NT impossible proved by NFBA certificate
    (Br == O) && (return addRound_RAW!(s, r))   # round is full
    (Br + Δ > O) && (Δ = O - Br)                # valid Δ 
    
    # least loaded outputs
    ll_out::Vector{Int64} = sort!(sortperm(s.loads)[1:Br + Δ]) 

    # all br sized ordered position index in ll_out
    CLC = CoolLexCombinations(Br + Δ, Br)

    # all br sized ordered position subset in ll_out
    all_newPos = [[ll_out[i] for i in e] for e in CLC]

    bestLoads::Union{Nothing, Vector{Int64}} = nothing
    bestAssignment::Union{Nothing, Vector{Int64}} = nothing

    for newPos in all_newPos
        valid = true
        newLoads::Vector{Int64} = deepcopy(s.loads)
        newAssignment::Vector{Int64} = zeros(Int64, O)

        # position the batches to their new location
        for i=1:Br
            newAssignment[newPos[i]] = r.batches[i]
            newLoads[newPos[i]] += r.batches[i]
    
            (newLoads[newPos[i]] <= s.C) || (valid = false; break) # updated loads overflow
        end

        if valid && (bestLoads === nothing || fitness(bestLoads, s.C, TAG_FitSes) > fitness(newLoads, s.C, TAG_FitSes))
            bestLoads = newLoads
            bestAssignment = newAssignment
        end
    end

    if bestLoads === nothing
        return (s, false)::Tuple{Session, Bool}
    else
        s.loads = bestLoads
        push!(s.rounds, Round(r.id, bestAssignment, r.batches))
        return (s, true)::Tuple{Session, Bool}
    end
end

"""
```Julia
    addRound_NFBA(s::Session, r::Round{N})  # wont edit s, return a boolean along with the new assignment
    addRound_NFBA!(s::Session, r::Round{N}) # may edit s
```
Left-alligned algorithm.
"""
function addRound_NFBA(s::Session, r::Round{N}) where N
    j               ::Int64         = 1
    k               ::Int64         = 1
    O               ::Int64         = length(s.loads)
    Br              ::Int64         = length(r.batches)
    valid           ::Bool          = true
    newAssignment   ::Vector{Int64} = zeros(Int64, O)

    while valid
        if r.batches[j] + s.loads[k] <= s.C
            newAssignment[k] = r.batches[j]
            j += 1
        end
        k += 1
        (j > Br) && (break)
        (Br - j > O - k) && (valid = false)
    end
    return (valid, newAssignment)::Tuple{Bool, Vector{Int64}}
end

function addRound_NFBA!(s::Session, r::Round{N}) where N
    j               ::Int64         = 1
    k               ::Int64         = 1
    O               ::Int64         = length(s.loads)
    Br              ::Int64         = length(r.batches)
    valid           ::Bool          = true
    newAssignment   ::Vector{Int64} = zeros(Int64, O)

    while valid
        if r.batches[j] + s.loads[k] <= s.C
            newAssignment[k] = r.batches[j]
            j += 1
        end
        k += 1
        (j > Br) && (break)
        (Br - j > O - k) && (valid = false)
    end

    if valid
        push!(s.rounds, r)
        s.loads += newAssignment
        r.assignment = newAssignment
    end

    return (valid, newAssignment)::Tuple{Bool, Vector{Int64}}
end

"""
```Julia
    addRound_SUB_01LP!(s::Session, r::Round{N}, tl::Int64 = 100, env::Gurobi.Env = Gurobi.Env())
```
Matheuristic reliying on a model to found an intresting mail assignment in the added route.
"""
function addRound_SUB_01LP!(s::Session, r::Round{N}, tl::Int64 = 100, env::Gurobi.Env = Gurobi.Env()) where N
# ===========< Data:                    OK >===========
    O::Int64 = length(r.assignment)

    # compute certificat
    _, flag::Bool, newAssignment::Vector{Int64} = addRound_NFBA(s, r)
    flag || (return (s, false)::Tuple{Session, Bool}) # print("<!>");   
    newRound::Round{N} = Round(r.id, newAssignment, r.batches)

    # Br::Vector{Int64} = collect(1:N)

    # vrj = r.batches 

    Orj::Vector{UnitRange{Int64}} = [(j: O - N + j) for j=1:N]

    # Lmax = s.C

# ===========< Model:                   OK >===========

    model = Model(() -> Gurobi.Optimizer(env)) # () -> GLPK.Optimizer()
    set_silent(model)
    # set_time_limit_sec(model, tl)
    set_optimizer_attribute(model, "TimeLimit", tl)
    set_optimizer_attribute(model, "OutputFlag", 0)

# ===========< Variables:               OK >===========

    @variable(model, x[j=1:N, k=Orj[j]], Bin)

# ===========< Objective:               OK >===========

    @objective(model, Min, sum([sum([((r.batches[j] + s.loads[k]) * x[j, k]) for k=Orj[j]]) for j=1:N]))
    
# ===========< Constraint: (4, 7, 10)   OK >===========

    # (4) -> require that each batch j of round r be assigned to exactly oneoutput
    for j=1:N
        @constraint(model, sum([x[j, k] for k=Orj[j]]) == 1)
    end

    # (7) -> precedence mail constraints for each round r
    for j=1:N
        if j != N
            @constraint(model, 1 + sum([k * x[j, k] for k=Orj[j]]) <= sum([k * x[j+1, k] for k=Orj[j + 1]]))
        end
    end

    # (10) -> capacity constraint
    for j=1:N
        @constraint(model, sum([((r.batches[j] + s.loads[k]) * x[j, k]) for k=Orj[j]]) <= s.C)
    end

# ===========< Warmup:                  OK >===========

    batche::Int64 = 1
    for k=1:O
        if newRound.assignment[k] != 0
            set_start_value(x[batche, k], 1)
            batche += 1
            
            (batche > N) && break
        end
    end

# ===========< Results:                 OK >===========
    optimize!(model)

    if termination_status(model) == OPTIMAL || MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0
        newRound.assignment = zeros(Int64, O)
        for j=1:N
            newPos::Int64 = 0
            for k=Orj[j] 
                (value(x[j,k]) != 0.0) && (newPos = k; break) 
            end
            newRound.assignment[newPos] = newRound.batches[j]
        end

        push!(s.rounds, newRound)

        s.loads += newRound.assignment
        # termination_status(model) == OPTIMAL ? print("<o>") : print("<0>")
        return s, true
    else
        # print("<x>")
        push!(s.rounds, newRound) 
        return s, true # should never trigger but better safe tha sorry (memory overflow)
    end
end

"""
```Julia
    addRound_PART_SUB_01LP!(s::Session, r::Round{N}, tl::Int64 = 100, env::Gurobi.Env = Gurobi.Env())
```
Matheuristic reliying on a model to found an intresting mail assignment in the added route. 
This model won't enforce the output capacity constraint in the first place hopping that thanks to the objective function no problem will come up.
If the resulting solution isn't valid (due to a capacity overflow), the capacity constraint of each problematic output will be added to the model.
The updated model will be solved again and this process will be repeted again and again until a valid solution is found. 
"""
function addRound_PART_SUB_01LP!(s::Session, r::Round{N}, tl::Int64 = 100, env::Gurobi.Env = Gurobi.Env()) where N
# ===========< Data:                    OK >===========
    O::Int64 = length(r.assignment)

    # compute certificat
    _, flag::Bool, newAssignment::Vector{Int64} = addRound_NFBA(s, r)
    flag || (print("<!>"), return (s, false)::Tuple{Session, Bool}) # print("<c>"), 
    newRound::Round{N} = Round(r.id, newAssignment, r.batches)

    # Br::Vector{Int64} = collect(1:N)

    # vrj = r.batches 

    Orj::Vector{UnitRange{Int64}} = [(j: O - N + j) for j=1:N]

    # Lmax = s.C

# ===========< Model:                   OK >===========

    model = Model(() -> Gurobi.Optimizer(env))
    set_silent(model)
    set_optimizer_attribute(model, "TimeLimit", tl)
    set_optimizer_attribute(model, "OutputFlag", 0)

# ===========< Variables:               OK >===========

    @variable(model, x[j=1:N, k=Orj[j]], Bin)

# ===========< Objective:               OK >===========

    # minimize the output loads in which a mail is affected
    @objective(model, Min, sum([sum([((r.batches[j] + s.loads[k]) * x[j, k]) for k=Orj[j]]) for j=1:N]))
    
# ===========< Constraint: (ctr 4, 7)   OK >===========

    # (4) -> require that each batch j of round r be assigned to exactly one output
    for j=1:N
        @constraint(model, sum([x[j, k] for k=Orj[j]]) == 1)
    end

    # (7) -> precedence mail constraints for each round r
    for j=1:N
        if j != N
            @constraint(model, 1 + sum([k * x[j, k] for k=Orj[j]]) <= sum([k * x[j+1, k] for k=Orj[j + 1]]))
        end
    end

# ===========< Warmup:                  OK >===========

    batche::Int64 = 1
    for k=1:O
        if newRound.assignment[k] != 0
            set_start_value(x[batche, k], 1)
            batche += 1
            
            (batche > N) && break
        end
    end

# ===========< Results:                 OK >===========
    optimize!(model)

    if termination_status(model) == OPTIMAL || MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0
        newRound.assignment = zeros(Int64, O)
        for j=1:N
            newPos::Int64 = 0
            for k=Orj[j] 
                (value(x[j,k]) != 0.0) && (newPos = k; break) 
            end
            newRound.assignment[newPos] = newRound.batches[j]
        end

        newLoad::Vector{Int64} = s.loads + newRound.assignment

        if validLoad(newLoad, s.C)
            termination_status(model) == OPTIMAL ? print("<o>") : print("<0>")
            push!(s.rounds, newRound)
            s.loads += newRound.assignment
            return s, true
        else
# ===========< Re-Optimize: (+ ctr 10)  OK >===========
            print("<x")
            for j=1:N
                @constraint(model, sum([((r.batches[j] + s.loads[k]) * x[j, k]) for k=Orj[j]]) <= s.C)
            end
            optimize!(model)

            if termination_status(model) == OPTIMAL || MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0
                newRound.assignment = zeros(Int64, O)
                for j=1:N
                    newPos::Int64 = 0
                    for k=Orj[j] 
                        (value(x[j,k]) != 0.0) && (newPos = k; break) 
                    end
                    newRound.assignment[newPos] = newRound.batches[j]
                end

                newLoad = s.loads + newRound.assignment

                if validLoad(newLoad, s.C)
                    termination_status(model) == OPTIMAL ? print("o>") : print("0>")
                    push!(s.rounds, newRound)
                    s.loads += newRound.assignment
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

function addRound_Rebuild!(s::Session, r::Round{N})::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    (global call += 1)

    ns::Session = Session(s.C, [Round(cr.id, deepcopy(cr.assignment), cr.batches) for cr in s.rounds], s.loads + r.assignment)
    push!(ns.rounds, Round(r.id, deepcopy(r.assignment), r.batches))

    ns, flag = rebuildSession(ns)
    
    (flag) && (global repairedBuild += 1; s.rounds = ns.rounds; s.loads = ns.loads)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRound_Rebuild_Knapsack!(s::Session, r::Round{N})::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    (global call += 1)

    ns::Session = Session(s.C, [Round(cr.id, deepcopy(cr.assignment), cr.batches) for cr in s.rounds], s.loads + r.assignment)
    push!(ns.rounds, Round(r.id, deepcopy(r.assignment), r.batches))

    ns, flag = rebuildSession_knapSack(ns)
    
    (flag) && (global repairedBuild += 1; s.rounds = ns.rounds; s.loads = ns.loads)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRound_Rebuild_Knapsack_model!(s::Session, r::Round{N}, env::Gurobi.Env = Gurobi.Env())::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    (global call += 1)

    ns::Session = Session(s.C, [Round(cr.id, deepcopy(cr.assignment), cr.batches) for cr in s.rounds], s.loads + r.assignment)
    push!(ns.rounds, Round(r.id, deepcopy(r.assignment), r.batches))

    ns, flag = rebuildSession_knapSack_model(ns, env)
    
    (flag) && (global repairedBuild += 1; s.rounds = ns.rounds; s.loads = ns.loads)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRound_Rebuild_Knapsack_model_V2!(s::Session, r::Round{N}, tl::Int64=10, env::Gurobi.Env = Gurobi.Env())::Tuple{Session, Bool} where N
    (global call += 1)

    ns::Session = Session(s.C, [Round(cr.id, deepcopy(cr.assignment), cr.batches) for cr in s.rounds], s.loads + r.assignment)
    push!(ns.rounds, Round(r.id, deepcopy(r.assignment), r.batches))

    ns, flag = rebuildSession_knapSack_model_V2(ns, tl, env)
    
    (flag) && (global repairedBuild += 1; s.rounds = ns.rounds; s.loads = ns.loads)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRound_Rebuild_Knapsack_model_V3!(s::Session, r::Round{N}, tl::Int64=10, env::Gurobi.Env = Gurobi.Env())::Tuple{Session, Bool} where N
    (global call += 1)

    ns::Session = Session(s.C, [Round(cr.id, deepcopy(cr.assignment), cr.batches) for cr in s.rounds], s.loads + r.assignment)
    push!(ns.rounds, Round(r.id, deepcopy(r.assignment), r.batches))

    ns, flag = rebuildSession_knapSack_model_V3!(ns, tl, env)
    
    (flag) && (global repairedBuild += 1; s.rounds = ns.rounds; s.loads = ns.loads)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRound_OPTIMOVE_stage1!(s::Session, r::Round{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)

    ns::Session = Session(s.C, [Round(cr.id, deepcopy(cr.assignment), cr.batches) for cr in s.rounds], s.loads + r.assignment)
    push!(ns.rounds, Round(r.id, deepcopy(r.assignment), r.batches))

    ns, nsFit, flag = improvedOptiMove_S1_V2!(ns, TAG_FitSes)
    
    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.rounds = ns.rounds; s.loads = ns.loads)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRound_OPTIMOVE_Stage1Valid!(s::Session, r::Round{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)

    ns::Session = Session(s.C, [Round(cr.id, deepcopy(cr.assignment), cr.batches) for cr in s.rounds], s.loads + r.assignment)
    push!(ns.rounds, Round(r.id, deepcopy(r.assignment), r.batches))

    ns, nsFit, flag = improvedOptiMove_S1_V3!(ns, TAG_FitSes)
    
    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.rounds = ns.rounds; s.loads = ns.loads)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRound_OPTIMOVE_Stage1ValidInf!(s::Session, r::Round{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)
    
    ns::Session = Session(s.C, [Round(cr.id, deepcopy(cr.assignment), cr.batches) for cr in s.rounds], s.loads + r.assignment)
    push!(ns.rounds, Round(r.id, deepcopy(r.assignment), r.batches))

    ns, nsFit, flag = improvedOptiMove_S1_V4!(ns, TAG_FitSes)
    
    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.rounds = ns.rounds; s.loads = ns.loads)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRound_OPTIMOVE_oneMove!(s::Session, r::Round{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)

    ns::Session = Session(s.C, [Round(cr.id, deepcopy(cr.assignment), cr.batches) for cr in s.rounds], s.loads + r.assignment)
    push!(ns.rounds, Round(r.id, deepcopy(r.assignment), r.batches))

    ns, nsFit, flag = improvedOptiMove_S1_V5!(ns, TAG_FitSes)
    
    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.rounds = ns.rounds; s.loads = ns.loads)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRound_OPTIMOVE_3Stages!(s::Session, r::Round{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)

    ns::Session = Session(s.C, [Round(cr.id, deepcopy(cr.assignment), cr.batches) for cr in s.rounds], s.loads + r.assignment)
    push!(ns.rounds, Round(r.id, deepcopy(r.assignment), r.batches))

    ns, nsFit, flag = improvedOptiMove_V1(ns, TAG_FitSes)
    
    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.rounds = ns.rounds; s.loads = ns.loads)
    return (ns, flag)::Tuple{Session, Bool}
end

function addRound_OPTIMOVE_Stage1NegSTD!(s::Session, r::Round{N}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD)::Tuple{Session, Bool} where N
    # session capacity
    # certificat_CapacityVolume(s, r) || (return (s, false)::Tuple{Session, Bool}) # print("<oc>"),  

    # println("\n\n$s\n")

    sFit::Float64 = fitness(s, TAG_FitSes)
    (global call += 1)

    ns::Session = Session(s.C, [Round(cr.id, deepcopy(cr.assignment), cr.batches) for cr in s.rounds], s.loads + r.assignment)
    push!(ns.rounds, Round(r.id, deepcopy(r.assignment), r.batches))

    ns, nsFit, flag = improvedOptiMove_S1_V6!(ns, TAG_FitSes)
    
    # println("\n\n$ns\n")

    (nsFit < sFit) && (global improvedOverAll += 1)
    (flag) && (global repairedBuild += 1; s.rounds = ns.rounds; s.loads = ns.loads)
    return (ns, flag)::Tuple{Session, Bool}
end

abstract type SimpleAddRound end

# MOVE 1 ROUND
struct RAW                  <: SimpleAddRound end
struct SAVA                 <: SimpleAddRound end
struct SAVANT               <: SimpleAddRound end
struct SAVANT_MIN_STD       <: SimpleAddRound end
struct NFBA                 <: SimpleAddRound end
struct SUB_01LP             <: SimpleAddRound end
struct PART_SUB_01LP        <: SimpleAddRound end

# MOVE SESSION
struct OPTIMOVE_STAGE1      <: SimpleAddRound end
struct OPTIMOVE_VALID       <: SimpleAddRound end
struct OPTIMOVE_VALID_INF   <: SimpleAddRound end
struct OPTIMOVE_1BATCH      <: SimpleAddRound end
struct OPTIMOVE_3S          <: SimpleAddRound end
struct OPTIMOVE_S1_NEGSTD   <: SimpleAddRound end
struct REBUILD_KP           <: SimpleAddRound end
struct REBUILD_KP_01LP      <: SimpleAddRound end
struct REBUILD_KP_01LP_V2   <: SimpleAddRound end

# MOVE 1 ROUND
@inline addRound(s::Session, r::Round{N}, ::Type{RAW}                 , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_RAW!(s, r)
@inline addRound(s::Session, r::Round{N}, ::Type{SAVA}                , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_SAVA!(s, r)
@inline addRound(s::Session, r::Round{N}, ::Type{NFBA}                , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_NFBA!(s, r)
@inline addRound(s::Session, r::Round{N}, ::Type{SAVANT}              , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_SAVANT!(s, r, Δ)
@inline addRound(s::Session, r::Round{N}, ::Type{SAVANT_MIN_STD}      , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_SAVANT_MINSTD!(s, r, Δ, TAG_FitSes)
@inline addRound(s::Session, r::Round{N}, ::Type{SUB_01LP}            , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_SUB_01LP!(s, r, env, tl)
@inline addRound(s::Session, r::Round{N}, ::Type{PART_SUB_01LP}       , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_PART_SUB_01LP!(s, r, env, tl)
@inline addRound(s::Session, r::Round{N}, ::Type{REBUILD_KP}          , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_Rebuild_Knapsack!(s, r)
@inline addRound(s::Session, r::Round{N}, ::Type{REBUILD_KP_01LP}     , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_Rebuild_Knapsack_model!(s, r, env)
@inline addRound(s::Session, r::Round{N}, ::Type{REBUILD_KP_01LP_V2}  , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_Rebuild_Knapsack_model_V2!(s, r, tl, env)

# MOVE SESSION
@inline addRound(s::Session, r::Round{N}, ::Type{OPTIMOVE_STAGE1}     , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_OPTIMOVE_stage1!(s, r, TAG_FitSes)
@inline addRound(s::Session, r::Round{N}, ::Type{OPTIMOVE_VALID}      , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_OPTIMOVE_Stage1Valid!(s, r, TAG_FitSes)
@inline addRound(s::Session, r::Round{N}, ::Type{OPTIMOVE_VALID_INF}  , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_OPTIMOVE_Stage1ValidInf!(s, r, TAG_FitSes)
@inline addRound(s::Session, r::Round{N}, ::Type{OPTIMOVE_1BATCH}     , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_OPTIMOVE_oneMove!(s, r, TAG_FitSes)
@inline addRound(s::Session, r::Round{N}, ::Type{OPTIMOVE_3S}         , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_OPTIMOVE_3Stages!(s, r, TAG_FitSes)
@inline addRound(s::Session, r::Round{N}, ::Type{OPTIMOVE_S1_NEGSTD}  , Δ::Int64 = 2, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, tl::Int64 = 2, env::Gurobi.Env = Gurobi.Env()) where N = addRound_OPTIMOVE_Stage1NegSTD!(s, r, TAG_FitSes)


# ================================================== #
#                       Display                      #
# ================================================== #

function Base.show(io::IO, s::Session)
    println(io, "Session: (Valid loads: $(validLoad(s))) (Valid session: $(validSession(s))) (filled at $(round(100 - fitness(s, HollowedPercentage), digits=3))%) (Standard deviation: $(round(fitness(s, LoadSTD), digits=3)))\n    C: $(s.C)\n    Rounds: (x$(length(s.rounds)))")
    if length(s.loads) <= 40
        for r in s.rounds
            println(io, "$(r)")
        end
    else
        for (i, r) in enumerate(s.rounds)
            i%10 == 0 ? println(io, "$(r)") : print(io, "$(r)")
        end
    end
    println(io, "\n    Loads:\n$(s.loads)")
end

# ================================================== #
#                    Miscelaneous                    #
# ================================================== #

"""
```Julia
    Random.shuffle(s::Session)::Session     # wont modify the session
    Random.shuffle!(s::Session)::Session    # will modify the session
```
shuffle the mail position within the session. (will keep the order among mails)
"""
function Random.shuffle!(s::Session)::Session
    for r in s.rounds
        shuffle!(r)
    end

    return compute_output!(s)
end

function Random.shuffle(s::Session)::Session
    return shuffle!(deepcopy(s))
end

"""
```Julia
    randSession(C::Int64=25, O::Int64=20, nbRounds::Int64=10)
```
Create a random session. `C` is the maximum capacity of each output. `O` is the number of output of the session and `nbRounds` is the amount of random routes that the session will bore.
 > :warning: note that this session may be invalid with the wrong set of parameter as most of it is completely random and no safety procedure are implemented.
"""
function randSession(C::Int64=25, O::Int64=20, nbRounds::Int64=10)
    rounds::Vector{Round} = [randRound(i, O, nothing, C) for i=1:nbRounds]

    return Session(C, rounds, compute_output(rounds))
end