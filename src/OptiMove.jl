begin # Library
    using Statistics
end

begin # Files
    include("Instance.jl")
end

# ================================================== #
#                      Getters                       #
# ================================================== #

"""
```Julia
function Ers(r::Route)::Vector{Int64}
```

Return the set of position corresponding to the empty output in route `r`.

`component of the Opti-Move heristic`

# Example
```jldoctest
julia> r = Route(0, [1, 2, 0, 4, 0, 0, 7, 8, 0, 0])
Route(0, [1, 2, 0, 4, 0, 0, 7, 8, 0, 0])

julia> Ers(r)
5-element Vector{Int64}:
  3
  5
  6
  9
 10
```
"""
function Ers(r::Route)::Vector{Int64}
    return findall(iszero, r.assignment)
end

"""
```Julia
function Frs(r::Route)::Vector{Int64}
```

Return the set of position corresponding to the non-empty output in route `r`.

`component of the Opti-Move heristic`

# Example
```jldoctest
julia> r = Route(0, [1, 2, 0, 4, 0, 0, 7, 8, 0, 0])
Route(0, [1, 2, 0, 4, 0, 0, 7, 8, 0, 0])

julia> Frs(r)
5-element Vector{Int64}:
 1
 2
 4
 7
 8
```
"""
function Frs(r::Route)::Vector{Int64}
    return findall(!iszero, r.assignment)
end

function adapt_EF(r::Route)
    open        ::Bool                                          = false
    groups      ::Vector{Tuple{Vector{Int64}, Vector{Int64}}}   = []
    lastBatch   ::Union{Nothing, Int64}                         = nothing

    for (i, b) in enumerate(r.assignment)
        if b == 0
            if open
                push!(groups[end][2], i)
            else
                open = true
                (lastBatch === nothing) ? push!(groups, ([], [i])) : push!(groups, ([lastBatch], [i]))
            end
        else
            if open
                open = false
                push!(groups[end][1], i)
            end
            lastBatch = i
        end
    end

    return groups
end

function select_ef(s::Session, r::Route)
    groups::Vector{Tuple{Vector{Int64}, Vector{Int64}}} = adapt_EF(r)

    f       ::Union{Nothing, Int64} = nothing
    group_f ::Vector{Int64} = []
    
    for (id, (F, E)) in enumerate(groups)
        for i in F
            ((f === nothing) || (s.load[f] < s.load[i])) && (f = i; group_f = [id])
            (f === i) && (push!(group_f, id))
        end
    end

    e       ::Union{Nothing, Int64} = nothing

    if !isempty(group_f)
        for group in group_f
            for i in groups[group][2]
                ((e === nothing) || (s.load[e] > s.load[i])) && (e = i)
            end
        end
    end

    return (f, e)
end

# ================================================== #
#                EmptyMove / ij-shift                #
# ================================================== #

# ijshift_V1!
function ijshift_V1!(
        s::Session,     # Session to edit
        rId::Int64,     # index of the route
        f::Int64,       # position of the most loaded letter in the route
        e::Int64        # position of the least loaded empty spot in the route
    )::Session
    # println("ijshift(r = ", rId, ", f = ", f, ", e = ", e, ")")

    r::Route = s.route[rId] # get the route
    interval::Vector{Int64} = (e < f) ? collect(e:f) : reverse(collect(f:e)) # interval between "f" and "e" ("e" must be in the first position of the vector)
    pl::Int64 = popfirst!(interval) # Previous letter <=> empty spot e (will be overwritten by the first found letter)

    for i in interval 
        if r.assignment[i] != 0 && i != e # if this position isn't empty
            # println("replace r.assignment[$pl] = $(r.assignment[pl]) by r.assignment[$i] = $(r.assignment[i])")
            r.assignment[pl] = r.assignment[i] # move the content of the current position to the previous letter position
            
            s.load[pl] += r.assignment[i] # add current value to precedent letter output
            s.load[i] -= r.assignment[i] # remove current value to curent letter output
            
            pl = i # set precedent spot to current location
        end
    end
    r.assignment[f] = 0 # delete the last letter (had been moved to last pl)

    return s
end

# ijshift_V1! but faster
function ijshift_V2(
        s       ::Session       ,     # Session to edit
        r       ::Route{N}         ,     # index of the route
        f       ::Int64         ,     # position of the most loaded letter in the route
        e       ::Int64         ,     # position of the least loaded empty spot in the route
    ) where N

    interval::StepRange{Int64, Int64}= e<f ? ((e+1):1:f) : ((e-1):-1:f) # interval between "f" and "e" ("e" must be in the first position of the vector)
    newPos::Int64 = e # Previous letter <=> empty spot e (will be overwritten by the first found letter)
    newAssignment::Vector{Int64} = deepcopy(r.assignment)
    newLoad::Vector{Int64} = deepcopy(s.load)


    for oldPos=interval 
        if newAssignment[oldPos] != 0 # && i != e # if this position isn't empty
            newAssignment[newPos] = newAssignment[oldPos]

            newLoad[newPos] += newAssignment[newPos]
            newLoad[oldPos] -= newAssignment[newPos]
            
            newPos = oldPos
        end
    end
    newAssignment[f] = 0
    return (newLoad, newAssignment)
end

# ijshift_V2 only valid movement
function ijshift_V3(
        s       ::Session       ,     # Session to edit
        r       ::Route{N}      ,     # index of the route
        f       ::Int64         ,     # position of the most loaded letter in the route
        e       ::Int64         ,     # position of the least loaded empty spot in the route
    ) where N

    interval::StepRange{Int64, Int64}= e<f ? ((e+1):1:f) : ((e-1):-1:f) # interval between "f" and "e" ("e" must be in the first position of the vector)
    newPos::Int64 = e # Previous letter <=> empty spot e (will be overwritten by the first found letter)
    newAssignment::Vector{Int64} = deepcopy(r.assignment)
    newLoad::Vector{Int64} = deepcopy(s.load)

    valid::Bool = true

    for oldPos=interval 
        if newAssignment[oldPos] != 0 # && i != e # if this position isn't empty
            newAssignment[newPos] = newAssignment[oldPos]

            (newLoad[oldPos] <= s.Lmax && newLoad[newPos] <= s.Lmax && newLoad[newPos] + newAssignment[oldPos] > s.Lmax) && (global locked += 1; valid = false; break) # print("!"); 

            newLoad[newPos] += newAssignment[newPos]
            newLoad[oldPos] -= newAssignment[newPos]
            
            newPos = oldPos
        end
    end
    newAssignment[f] = 0
    return (newLoad, newAssignment, valid)
end

# ================================================== #
#                  ij-shift Neighbor                 #
# ================================================== #

function Neighbor_V1(
        s::Session, # Session to edit
        τ::Float64 = 0.0, # tolerance
        TAG_FitSes::Type{T} = LoadSTD # fistness function
    )::Session where {T<:FitnessSession}
    # ==========< Step 1 >==========

    for (rId, r) in enumerate(s.route)
        # ==========< Step 2 >==========
        F::Vector{Int64} = Frs(r) # set of non-empty output in route r
        E::Vector{Int64} = Ers(r) # set of empty output in route r
    
        # ==========< Step 3 >==========
        if !isempty(E) && !isempty(F) # if none is empty
            f::Int64 = F[1] # most loaded output of F (k in the Adrien's paper)
            for i=2:length(F)
                if s.load[f] < s.load[F[i]]
                    f = F[i]
                end
            end

            e::Int64 = E[1] # least loaded output of E (q in the Adrien's paper)
            for i=2:length(E)
                if s.load[e] > s.load[E[i]]
                    e = E[i]
                end
            end
            
            # ==========< Step 4 >==========
            st::Session = ijshift_V1!(deepcopy(s), rId, f, e)

            if compare_1Criteria(st, s, TAG_FitSes=TAG_FitSes, τ=τ)
               s = st
            end
        end
    end

    return s
end

function Neighbor_V2(
        s::Session, # Session to edit
        τ::Float64 = 0.0, # tolerance
        TAG_FitSes::Type{<:FitnessSession} = LoadSTD # fistness function
    )
    # ==========< Step 1 >==========

    R               ::Int64                     = length(s.route)          # number of route in the session
    edit            ::Bool                      = false                     # true if we change s (made a copy of it)
    sFit            ::Float64                   = fitness(s, TAG_FitSes)    # fitness of the current session
    newFit          ::Union{Float64, Nothing}   = nothing                   # fitness of a new session

    for rId=1:R
        r           ::Route                     = s.route[rId]             # current batch

        # ==========< Step 2 >==========
        e           ::Union{Nothing, Int64}     = nothing                   # most loaded output of s with a non nul batch in r
        f           ::Union{Nothing, Int64}     = nothing                   # least loaded output of s with a nul batch in r
        
        for (k, b) in enumerate(r.assignment)
            if b == 0
                ((e === nothing) || (s.load[e] > s.load[k])) && (e = k)
            else
                ((f === nothing) || (s.load[f] < s.load[k])) && (f = k)
            end
        end

        # ==========< Step 3 >==========
        if !(e === nothing || f === nothing) # if none is empty
           
            # ==========< Step 4 >==========
            newLoad::Vector{Int64}, newAssignment::Vector{Int64} = ijshift_V2(s, r, f, e)
            newFit = fitness(newLoad, s.Lmax, TAG_FitSes)

            if newFit < (sFit + τ)
                sFit = newFit
                if edit # Update the current Session (outer scop session untouch because a copy of it is made if no edit as been made before)
                    s.load         = newLoad
                    s.route[rId]   = Route(r.id, newAssignment, r.mail)
                else    # Create a new session with the correct updates
                    s               = Session(s.Lmax, [(cr.id == r.id) ? Route(cr.id, newAssignment, cr.mail) : cr for cr in s.route], newLoad)
                    edit            = true
                end
            end
        end
    end

    return s, edit, sFit
end

function Neighbor_V2!(
        s::Session, # Session to edit
        τ::Float64 = 0.0, # tolerance
        TAG_FitSes::Type{<:FitnessSession} = LoadSTD # fistness function
    )
    # ==========< Step 1 >==========

    R               ::Int64                     = length(s.route)          # number of route in the session
    edit            ::Bool                      = false                     # true if we change s (made a copy of it)
    sFit            ::Float64                   = fitness(s, TAG_FitSes)    # fitness of the current session
    newFit          ::Union{Float64, Nothing}   = nothing                   # fitness of a new session

    for rId=1:R
        r           ::Route                     = s.route[rId]             # current batch

        # ==========< Step 2 >==========
        e           ::Union{Nothing, Int64}     = nothing                   # most loaded output of s with a non nul batch in r
        f           ::Union{Nothing, Int64}     = nothing                   # least loaded output of s with a nul batch in r
        
        for (k, b) in enumerate(r.assignment)
            if b == 0
                ((e === nothing) || (s.load[e] > s.load[k])) && (e = k)
            else
                ((f === nothing) || (s.load[f] < s.load[k])) && (f = k)
            end
        end

        # ==========< Step 3 >==========
        if !(e === nothing || f === nothing) # if none is empty
        
            # ==========< Step 4 >==========
            newLoad::Vector{Int64}, newAssignment::Vector{Int64} = ijshift_V2(s, r, f, e)
            newFit = fitness(newLoad, s.Lmax, TAG_FitSes)

            if newFit < (sFit + τ)
                edit = true
                sFit = newFit
                s.load         = newLoad
                r.assignment    = newAssignment
            end
        end
    end

    return s, edit, sFit
end

function Neighbor_V3(
        s::Session, # Session to edit
        τ::Float64 = 0.0, # tolerance
        TAG_FitSes::Type{<:FitnessSession} = LoadSTD # fistness function 
    )
    # ==========< Step 1 >==========

    R               ::Int64                     = length(s.route)          # number of route in the session
    edit            ::Bool                      = false                     # true if we change s (made a copy of it)
    sFit            ::Float64                   = fitness(s, TAG_FitSes)    # fitness of the current session
    newFit          ::Union{Float64, Nothing}   = nothing                   # fitness of a new session

    for rId=1:R
        r           ::Route                     = s.route[rId]             # current batch

        # ==========< Step 2 >==========
        e           ::Union{Nothing, Int64}     = nothing                   # most loaded output of s with a non nul batch in r
        f           ::Union{Nothing, Int64}     = nothing                   # least loaded output of s with a nul batch in r
        
        for (k, b) in enumerate(r.assignment)
            if b == 0
                ((e === nothing) || (s.load[e] > s.load[k])) && (e = k)
            else
                ((f === nothing) || (s.load[f] < s.load[k])) && (f = k)
            end
        end

        # ==========< Step 3 >==========
        if !(e === nothing || f === nothing) # if none is empty
        
            # ==========< Step 4 >==========
            newLoad::Vector{Int64}, newAssignment::Vector{Int64}, valid::Bool = ijshift_V3(s, r, f, e)
            if valid
                newFit = fitness(newLoad, s.Lmax, TAG_FitSes)

                if newFit < (sFit + τ)
                    sFit = newFit
                    if edit # Update the current Session (outer scop session untouch because a copy of it is made if no edit as been made before)
                        s.load         = newLoad
                        s.route[rId]   = Route(r.id, newAssignment, r.mail)
                    else    # Create a new session with the correct updates
                        s               = Session(s.Lmax, [(cr.id == r.id) ? Route(cr.id, newAssignment, cr.mail) : cr for cr in s.route], newLoad)
                        edit            = true
                    end
                end
            end
        end
    end

    return s, edit, sFit
end

function Neighbor_V3!(
        s::Session, # Session to edit
        τ::Float64 = 0.0, # tolerance
        TAG_FitSes::Type{<:FitnessSession} = LoadSTD # fistness function
    )
    # ==========< Step 1 >==========

    R               ::Int64                     = length(s.route)          # number of route in the session
    edit            ::Bool                      = false                     # true if we change s (made a copy of it)
    sFit            ::Float64                   = fitness(s, TAG_FitSes)    # fitness of the current session
    newFit          ::Union{Float64, Nothing}   = nothing                   # fitness of a new session

    for rId=1:R
        r           ::Route                     = s.route[rId]             # current batch

        # ==========< Step 2 >==========
        e           ::Union{Nothing, Int64}     = nothing                   # most loaded output of s with a non nul batch in r
        f           ::Union{Nothing, Int64}     = nothing                   # least loaded output of s with a nul batch in r
        
        for (k, b) in enumerate(r.assignment)
            if b == 0
                ((e === nothing) || (s.load[e] > s.load[k])) && (e = k)
            else
                ((f === nothing) || (s.load[f] < s.load[k])) && (f = k)
            end
        end

        # ==========< Step 3 >==========
        if !(e === nothing || f === nothing) # if none is empty
        
            # ==========< Step 4 >==========
            newLoad::Vector{Int64}, newAssignment::Vector{Int64}, valid::Bool = ijshift_V3(s, r, f, e)
            if valid
                newFit = fitness(newLoad, s.Lmax, TAG_FitSes)

                if newFit < (sFit + τ)
                    edit = true
                    sFit = newFit
                    s.load         = newLoad
                    r.assignment    = newAssignment
                end
            end
        end
    end

    return s, edit, sFit
end

function Neighbor_V4!(
        s::Session, # Session to edit
        τ::Float64 = 0.0, # tolerance
        TAG_FitSes::Type{<:FitnessSession} = LoadSTD # fistness function
    )
    # ==========< Step 1 >==========

    R               ::Int64                     = length(s.route)          # number of route in the session
    edit            ::Bool                      = false                     # true if we change s (made a copy of it)
    sFit            ::Float64                   = fitness(s, TAG_FitSes)    # fitness of the current session
    newFit          ::Union{Float64, Nothing}   = nothing                   # fitness of a new session

    for rId=1:R
        r           ::Route                     = s.route[rId]             # current batch

        # ==========< Step 2 >==========
        e           ::Union{Nothing, Int64}     = nothing                   # most loaded output of s with a non nul batch in r
        f           ::Union{Nothing, Int64}     = nothing                   # least loaded output of s with a nul batch in r
        
        for (k, b) in enumerate(r.assignment)
            if b == 0
                ((e === nothing) || (s.load[e] > s.load[k])) && (e = k)
            else
                ((f === nothing) || (s.load[f] < s.load[k])) && (f = k)
            end
        end

        # ==========< Step 3 >==========
        if !(e === nothing || f === nothing) # if none is empty
        
            # ==========< Step 4 >==========
            newLoad::Vector{Int64}, newAssignment::Vector{Int64}, valid::Bool = ijshift_V3(s, r, f, e)
            if valid
                newFit = fitness(newLoad, s.Lmax, TAG_FitSes)
                edit = true
                sFit = newFit
                s.load         = newLoad
                r.assignment    = newAssignment
            end
        end
    end

    return s, edit, sFit
end

function Neighbor_V5(
        s::Session, # Session to edit
        τ::Float64 = 0.0, # tolerance
        TAG_FitSes::Type{<:FitnessSession} = LoadSTD # fistness function
    )
    # ==========< Step 1 >==========
    N               ::Int64                     = 10
    C               ::Int64                     = N
    R               ::Int64                     = length(s.route)          # number of route in the session
    edit            ::Bool                      = false                     # true if we change s (made a copy of it)
    keepRunning     ::Bool                      = false 
    improvement     ::Bool                      = false 
    sFit            ::Float64                   = fitness(s, TAG_FitSes)    # fitness of the current session
    newFit          ::Union{Float64, Nothing}   = nothing                   # fitness of a new session

    for (rId, r) in enumerate(s.route)
        # println(r)
        newLoad         ::Vector{Int64}             = deepcopy(s.load)
        newAssignment   ::Vector{Int64}             = deepcopy(r.assignment)
        # println(newAssignment)
        
        improvement = false
        keepRunning = true
        while keepRunning
            f::Union{Nothing, Int64} , e::Union{Nothing, Int64} = select_ef(s, r)
            # print("-> e = $e, f = $f -> ")

            # Non empty and no more un valid outputs
            if !(e === nothing || f === nothing) && (newLoad[e] + newAssignment[f] <= s.Lmax || !(newLoad[f] <= s.Lmax) || !(newLoad[e] <= s.Lmax))

                startFit = fitness(newLoad, s.Lmax, TAG_FitSes)

                # Updates
                newAssignment[e] = newAssignment[f]
                newLoad[e] += newAssignment[f]
                newLoad[f] -= newAssignment[f]
                newAssignment[f] = 0

                r.assignment = newAssignment

                newFit = fitness(newLoad, s.Lmax, TAG_FitSes)

                improvement = true
                if newFit < startFit
                    # print("o")
                    keepRunning = true
                else 
                    # print("x")
                    C -= 1
                    keepRunning = (C >= 0)
                end
            else
                keepRunning = false
            end
        end
        # println(newAssignment)

        if improvement
            # print("-")
            sFit = newFit
            if edit # Update the current Session (outer scop session untouch because a copy of it is made if no edit as been made before)
                s.load         = newLoad
                s.route[rId]   = Route(r.id, newAssignment, r.mail)
            else    # Create a new session with the correct updates
                s               = Session(s.Lmax, [(cr.id == r.id) ? Route(cr.id, newAssignment, cr.mail) : cr for cr in s.route], newLoad)
                edit            = true
                # print("deepcopy")
            end
        end
    end

    return s, edit, sFit
end

function Neighbor_V6!(
        s::Session, # Session to edit
        τ::Float64 = 0.0, # tolerance
        TAG_FitSes::Type{<:FitnessSession} = LoadSTD # fistness function
    )
    # ==========< Step 1 >==========

    R               ::Int64                     = length(s.route)          # number of route in the session
    edit            ::Bool                      = false                     # true if we change s (made a copy of it)
    sFit            ::Float64                   = fitness(s, TAG_FitSes)    # fitness of the current session
    newFit          ::Union{Float64, Nothing}   = nothing                   # fitness of a new session

    for rId=1:R
        r           ::Route                     = s.route[rId]             # current batch

        # ==========< Step 2 >==========
        e           ::Union{Nothing, Int64}     = nothing                   # most loaded output of s with a non nul batch in r
        f           ::Union{Nothing, Int64}     = nothing                   # least loaded output of s with a nul batch in r
        
        for (k, b) in enumerate(r.assignment)
            if b == 0
                ((e === nothing) || (s.load[e] > s.load[k])) && (e = k)
            else
                ((f === nothing) || (s.load[f] < s.load[k])) && (f = k)
            end
        end

        # ==========< Step 3 >==========
        if !(e === nothing || f === nothing) # if none is empty
        
            # ==========< Step 4 >==========
            capacityOverflow::Int64 = count(x -> x >= s.Lmax, s.load)
            
            newLoad::Vector{Int64}, newAssignment::Vector{Int64}, valid::Bool = ijshift_V3(s, r, f, e)
            newFit = fitness(newLoad, s.Lmax, TAG_FitSes)

            if valid || (newFit > (sFit + τ) && capacityOverflow >= count(x -> x >= s.Lmax, newLoad))
                edit = true
                sFit = newFit
                s.load         = newLoad
                r.assignment    = newAssignment
            end
        end
    end

    return s, edit, sFit
end

# ================================================== #
#                  Improved OptiMove                 #
# ================================================== #

function improvedOptiMove_S1_V1(
        s::Session, 
        # N::Int64 = 10, # not used in Stage 1
        TAG_FitSes::Type{T} = LoadSTD # fitness function
    )::Session where {T<:FitnessSession}
    # ==========< Stage 1 >==========
    # =====< 1.1 >=====
    τ::Float64 = 0
    # print("|")

    # =====< 1.2 >=====
    keepRunning = true
    while keepRunning
        # print("^")
        sc::Session = Neighbor_V1(s, τ, TAG_FitSes)
        compare_1Criteria(sc, s, TAG_FitSes=TAG_FitSes, τ=τ) ? s = sc : (s = sc; keepRunning = false)
        validSession(s) && (return (s)) # print("!"); 
    end

    return s
end

function improvedOptiMove_S1_V2!(
        s::Session, 
        # N::Int64 = 10, # not used in Stage 1
        TAG_FitSes::Type{T} = LoadSTD # fitness function
    ) where {T<:FitnessSession}
    # ==========< Stage 1 >==========
    # =====< 1.1 >=====
    # print(" ($(s.route[end].id)) -> ")
    τ           ::Float64   = 0.0
    sFit        ::Float64   = fitness(s, TAG_FitSes)
    keepRunning ::Bool      = true                      # true until no modification were made with Neighbor_V2!

    # =====< 1.2 >=====
    while keepRunning
        # print("^")
        s, keepRunning, sFit = Neighbor_V2!(s, τ, TAG_FitSes)
        validSession(s) && (return (s, sFit, true)) # print("!"); 
    end

    return s, sFit, false
end

function improvedOptiMove_S1_V3!(
        s::Session, 
        # N::Int64 = 10, # not used in Stage 1
        TAG_FitSes::Type{T} = LoadSTD # fitness function
    ) where {T<:FitnessSession}
    # ==========< Stage 1 >==========
    # =====< 1.1 >=====
    # print(" ($(s.route[end].id)) -> ")
    τ           ::Float64   = 0.0
    sFit        ::Float64   = fitness(s, TAG_FitSes)
    keepRunning ::Bool      = true                      # true until no modification were made with Neighbor_V2!
    
    # =====< 1.2 >=====
    while keepRunning
        # print("^")
        s, keepRunning, sFit = Neighbor_V3!(s, τ, TAG_FitSes)
        validSession(s) && (return (s, sFit, true)) # print("!"); 
    end

    return s, sFit, false
end

function improvedOptiMove_V1(
        s::Session, 
        TAG_FitSes::Type{T} = LoadSTD, # fitness function
        N::Int64 = 10, # not used in Stage 1
    ) where {T<:FitnessSession}
    # print("|")
    # ==========< Stage 1 >==========
    # =====< 1.1 >=====
    sFit::Float64 = fitness(s, TAG_FitSes)
    τ::Float64 = 0.0

    # =====< 1.2 >=====
    keepRunning::Bool = true
    while keepRunning
        # print("^")
        s, keepRunning, sFit = Neighbor_V2!(s, τ, TAG_FitSes)
        validSession(s) && (return (s, sFit, true)) # print("!"); 
    end

    # print(">")

    # ==========< Stage 2 >==========
    # =====< 2.1 >=====
    # start::Float64 = fitness(sb, TAG_FitSes)
    τ2::Float64 = 0.05 * fitness(s, TAG_FitSes)
    Δ::Float64 = 0.02
    C::Int64 = N

    bsFit::Float64 = sFit
    nsFit::Float64 = sFit
    bs::Session = s
    ns::Session = s
    # start = bsFit

    # =====< 2.2 >=====
    while τ2 > 0
        # print("τ")
        C -= 1
        ns, _, nsFit = Neighbor_V3(ns, τ2, TAG_FitSes)
        keepRunning = nsFit < sFit
        (nsFit < bsFit) && (bsFit = nsFit, bs = ns)
        
        while keepRunning
            # print("-")
            C -= 1
            ns, _, nsFit = Neighbor_V3(ns, τ2, TAG_FitSes)
            validSession(ns) && (return (ns, nsFit, true)) # print("!"); 
            if nsFit < bsFit
                keepRunning = C > -N
                (nsFit < bsFit) && (bsFit = nsFit, bs = ns)
            else
                keepRunning = false
            end
        end

        if C <= 0 # =====< 2.3 >=====
            keepRunning = true
            while keepRunning
                # print("^")
                ns, keepRunning, nsFit = Neighbor_V2!(ns, τ, TAG_FitSes)
                validSession(ns) && (return (ns, nsFit, true)) # print("!"); 
                (nsFit < bsFit) && (bsFit = nsFit; bs = ns)               # new best solution
            end
            C = N
        else # =====< 2.4 >=====
            τ2 -= Δ
            s = ns
            # (nsFit < sFit) ? (ns = s) : (s = ns)
        end
    end
    # print(">")

    # ==========< Stage 3 >==========
    # =====< 3.1 >=====
    # =====< 3.2 >=====

    keepRunning = true
    while keepRunning
        # print("^")
        bs, keepRunning, bsFit = Neighbor_V2!(bs, τ, TAG_FitSes)
        validSession(bs) && (return (bs, bsFit, true)) # print("!"); 
    end

    return bs, bsFit, false
end

function improvedOptiMove_S1_V4!(
        s::Session, 
        # N::Int64 = 10, # not used in Stage 1
        TAG_FitSes::Type{T} = LoadSTD # fitness function
    ) where {T<:FitnessSession}
    # ==========< Stage 1 >==========
    # =====< 1.1 >=====
    # print(" ($(s.route[end].id)) -> ")
    τ           ::Float64   = 0.0
    sFit        ::Float64   = fitness(s, TAG_FitSes)
    keepRunning ::Bool      = true                      # true until no modification were made with Neighbor_V2!
    i = 1

    # =====< 1.2 >=====
    while keepRunning
        # print("^")
        s, keepRunning, sFit = Neighbor_V4!(s, τ, TAG_FitSes)
        validSession(s) && (return (s, sFit, true)) #  print("!");
        i += 1
        (i > 50) && (keepRunning = false)
    end

    return s, sFit, false
end

function improvedOptiMove_S1_V5!(
        s::Session, 
        # N::Int64 = 10, # not used in Stage 1
        TAG_FitSes::Type{T} = LoadSTD # fitness function
    ) where {T<:FitnessSession}
    # ==========< Stage 1 >==========
    # =====< 1.1 >=====
    # print(" ($(s.route[end].id)) -> ")
    τ           ::Float64   = 0.0
    sFit        ::Float64   = fitness(s, TAG_FitSes)
    keepRunning ::Bool      = true                      # true until no modification were made with Neighbor_V2!

    # =====< 1.2 >=====
    i = 1
    while keepRunning
        # print("^")
        s, keepRunning, sFit = Neighbor_V5(s, τ, TAG_FitSes)
        validSession(s) && (return (s, sFit, true)) # print("!"); 
        i += 1
        (i > 50) && (keepRunning = false)
    end

    return s, sFit, false
end

function improvedOptiMove_S1_V6!(
        s::Session, 
        # N::Int64 = 10, # not used in Stage 1
        TAG_FitSes::Type{T} = LoadSTD # fitness function
    ) where {T<:FitnessSession}
    # ==========< Stage 1 >==========
    # =====< 1.1 >=====
    # print(" ($(s.route[end].id)) -> ")
    τ           ::Float64   = 0.0
    sFit        ::Float64   = fitness(s, TAG_FitSes)
    keepRunning ::Bool      = true                      # true until no modification were made with Neighbor_V2!
    i           ::Int64     = 0

    # =====< 1.2 >=====
    while keepRunning
        # print("^")
        s, keepRunning, sFit = Neighbor_V6!(s, τ, TAG_FitSes)
        validSession(s) && (return (s, sFit, true)) # print("!");
        (i >= length(s.load)) && (keepRunning = false)
        i += 1
    end

    return s, sFit, false
end