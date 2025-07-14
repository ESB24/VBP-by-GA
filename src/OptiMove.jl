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
        r       ::Route{N}      ,     # index of the route
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

            (newLoad[oldPos] <= s.Lmax && newLoad[newPos] <= s.Lmax && newLoad[newPos] + newAssignment[oldPos] > s.Lmax) && (valid = false; break) # print("!"); global locked += 1; 

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
        r           ::Route                     = s.route[rId]             # current route

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

    R               ::Int64                     = length(s.route)           # number of route in the session
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
        isSessionValid(s) && (return (s)) # print("!"); 
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
        isSessionValid(s) && (return (s, sFit, true)) # print("!"); 
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
        isSessionValid(s) && (return (s, sFit, true)) # print("!"); 
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
        isSessionValid(s) && (return (s, sFit, true)) # print("!"); 
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
        (nsFit < bsFit) && (bsFit = nsFit; bs = ns)
        
        while keepRunning
            # print("-")
            C -= 1
            ns, _, nsFit = Neighbor_V3(ns, τ2, TAG_FitSes)
            isSessionValid(ns) && (return (ns, nsFit, true)) # print("!"); 
            if nsFit < bsFit
                keepRunning = C > -N
                (nsFit < bsFit) && (bsFit = nsFit; bs = ns)
            else
                keepRunning = false
            end
        end

        if C <= 0 # =====< 2.3 >=====
            keepRunning = true
            while keepRunning
                # print("^")
                ns, keepRunning, nsFit = Neighbor_V2!(ns, τ, TAG_FitSes)
                isSessionValid(ns) && (return (ns, nsFit, true)) # print("!"); 
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
        isSessionValid(bs) && (return (bs, bsFit, true)) # print("!"); 
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
        isSessionValid(s) && (return (s, sFit, true)) #  print("!");
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
        isSessionValid(s) && (return (s, sFit, true)) # print("!"); 
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
        isSessionValid(s) && (return (s, sFit, true)) # print("!");
        (i >= length(s.load)) && (keepRunning = false)
        i += 1
    end

    return s, sFit, false
end

## ============================================================================================================== ##
##       ########  ##    ##  #######   ########  ##    ##            ##    ##   ######   ##    ##  ########       ##
##       ##        ###  ###  ##    ##     ##      ##  ##             ###  ###  ##    ##  ##    ##  ##             ##
##       #####     ## ## ##  #######      ##        ##               ## ## ##  ##    ##   ##  ##   #####          ##
##       ##        ##    ##  ##           ##        ##               ##    ##  ##    ##   ##  ##   ##             ##
##       ########  ##    ##  ##           ##        ##               ##    ##   ######      ##     ########       ##
## ============================================================================================================== ##

function EmptyMove(
        s       ::Session               ,     # Session to edit
        r       ::Route{N}              ;     # index of the Route
        τ       ::Float64       = 0.0   ,
        obj     ::Type{<:FitnessSession}    = OutputVariance    ,   # regarded objective
    ) where N
    
    e::Int64 = findall(iszero, r.assignment)[argmin(map(k -> s.load[k], findall(iszero, r.assignment)))]    # least loaded empty output
    f::Int64 = findall(!iszero, r.assignment)[argmax(map(k -> s.load[k], findall(!iszero, r.assignment)))]    # least loaded empty output
    range::StepRange{Int64, Int64} = (e<f) ? ((e+1):f) : ((e-1):-1:f)           # interval from e to f
    
    newPos          ::Int64         = e                         # ne position for moved mail
    newAssignment   ::Vector{Int64} = deepcopy(r.assignment)    # updated positions
    newLoad         ::Vector{Int64} = deepcopy(s.load)          # updated loads

    valid::Bool = true

    # println("$e, $f, $range", ANSI_cyan)

    for oldPos = range
        if newAssignment[oldPos] != 0
            newAssignment[newPos] = newAssignment[oldPos]

            if (newLoad[oldPos] <= s.Lmax) && (newLoad[newPos] <= s.Lmax) && (newLoad[newPos] + newAssignment[oldPos] > s.Lmax)
                # a constraint violation has been created from valid position
                valid = false
                # break
            end

            newLoad[newPos] += newAssignment[newPos]
            newLoad[oldPos] -= newAssignment[newPos]
            
            newPos = oldPos
        end
    end
    
    newAssignment[f] = 0

    # print_verbose("<EM $(valid ? "V" : "X") - f=$f to e=$e >", valid ? ANSI_green : ANSI_red)
    
    return (newLoad, newAssignment, valid)
end

function EmptyMove_V2!(
        s       ::Session                                       ,   # Session to edit
        r       ::Route{N}                                      ,   # index of the Route
        obj     ::Type{<:FitnessSession}    = OutputVariance    ,   # regarded objective
    ) where N
    
    load::Vector{Int64} = s.load
    assignment::Vector{Int64} = r.assignment
    O::Int64 = length(load)
    e::Int64 = -1 # least loaded empty output
    f::Int64 = -1 # most loaded filled output
    for k=1:O
        if assignment[k] == 0
            ((e == -1) || (load[e] > load[k])) && (e = k)
        else
            ((f == -1) || (load[f] < load[k])) && (f = k)
        end
    end
    range::StepRange{Int64, Int64} = (e<f) ? ((e+1):f) : ((e-1):-1:f)   # interval from e to f
    newPos::Int64 = e                         # ne position for moved mail
    valid::Bool = true

    if ((e != -1) && (f != -1))
        for oldPos=range
            if assignment[oldPos] != 0
                assignment[newPos] = assignment[oldPos]
                load[newPos] += assignment[newPos]
                load[oldPos] -= assignment[newPos]
                
                newPos = oldPos
            end
        end
        assignment[f] = 0
    end
end

function EmptyMove_V3!(
        s       ::Session                                       ,   # Session to edit
        r       ::Route{N}                                      ,   # index of the Route
        obj     ::Type{<:FitnessSession}    = OutputVariance    ,   # regarded objective
    ) where N

    load::Vector{Int64} = s.load
    assignment::Vector{Int64} = r.assignment
    O::Int64 = length(load)
    Lmax::Int64 = s.Lmax
    e::Int64 = -1 # least loaded empty output
    f::Int64 = -1 # most loaded filled output
    for k=1:O
        if assignment[k] == 0
            ((e == -1) || (load[e] > load[k])) && (e = k)
        else
            ((f == -1) || (load[f] < load[k])) && (f = k)
        end
    end
    range::StepRange{Int64, Int64} = (e<f) ? ((e+1):f) : ((e-1):-1:f)   # interval from e to f
    newPos::Int64 = e                         # ne position for moved mail
    edited::Vector{Int64} = zeros(Int64, O)

    # println("e=$e, f=$f, range=$range")

    for oldPos=range
        if assignment[oldPos] != 0
            if ((load[newPos] <= Lmax) && (load[newPos] + assignment[newPos] > Lmax))
                edited[newPos] += 1 
            end
            if ((load[oldPos] > Lmax) && (load[oldPos] - assignment[oldPos] <= Lmax))
                edited[oldPos] -= 1
            end
            assignment[newPos] = assignment[oldPos]
            load[newPos] += assignment[newPos]
            load[oldPos] -= assignment[newPos]
            newPos = oldPos
        end
    end
    assignment[f] = 0

    return edited
end

function EmptyMove_V4(
        s       ::Session                                       ,   # Session to edit
        r       ::Route{N}                                      ,   # index of the Route
        obj     ::Type{<:FitnessSession}    = OutputVariance    ,   # regarded objective
    ) where N
    
    load::Vector{Int64} = deepcopy(s.load)
    assignment::Vector{Int64} = deepcopy(r.assignment)
    O::Int64 = length(load)
    e::Int64 = -1 # least loaded empty output
    f::Int64 = -1 # most loaded filled output
    for k=1:O
        if assignment[k] == 0
            ((e == -1) || (load[e] > load[k])) && (e = k)
        else
            ((f == -1) || (load[f] < load[k])) && (f = k)
        end
    end
    range::StepRange{Int64, Int64} = (e<f) ? ((e+1):f) : ((e-1):-1:f)   # interval from e to f
    newPos::Int64 = e                         # ne position for moved mail
    valid::Bool = true

    if ((e != -1) && (f != -1))
        for oldPos=range
            if assignment[oldPos] != 0
                assignment[newPos] = assignment[oldPos]
                load[newPos] += assignment[newPos]
                load[oldPos] -= assignment[newPos]
                
                newPos = oldPos
            end
        end
        assignment[f] = 0
    end

    return assignment, load
end

## ============================================================================================================== ##
##                 ##    ##  ########  ########   ######   ##    ##  ######     ######   #######                  ##
##                 ###   ##  ##           ##     ##        ##    ##  ##    ##  ##    ##  ##    ##                 ##
##                 ## ## ##  #####        ##     ##  ###   ########  #######   ##    ##  #######                  ##
##                 ##   ###  ##           ##     ##    ##  ##    ##  ##    ##  ##    ##  ##  ##                   ##
##                 ##    ##  ########  ########   ######   ##    ##  #######    ######   ##   ##                  ##
## ============================================================================================================== ##

function EmptyMoveNeighborhood_cascade(
        s       ::Session                                       ;   # Session to edit
        τ       ::Float64                   = 0.0               ,   # Fault tolerance
        obj     ::Type{<:FitnessSession}    = OutputVariance    ,   # regarded objective
    )
    s = Session(s.Lmax, [Route(r.id, deepcopy(r.assignment), r.mail) for r in s.route], deepcopy(s.load))

    # print_verbose("<EMN Cascade { ", ANSI_cyan)

    tmpObj = round(fitness(s, obj), digits=3)

    movement::Int64 = 0

    for r in s.route
        newLoad, newAssignment, valid = EmptyMove(s, r, τ=τ, obj=obj)
        if true # valid    # Valid movement
            s.load          = newLoad       # update session's loads
            r.assignment    = newAssignment # update route's assignment
            movement += 1
        end
    end

    return s, movement
end

function EmptyMoveNeighborhood_cascade_V2(
        s       ::Session                                       ,   # initial session
        obj     ::Type{<:FitnessSession}    = OutputVariance    ,   # regarded objective
    )
    s = Session(s.Lmax, [Route(r.id, deepcopy(r.assignment), r.mail) for r in s.route], deepcopy(s.load))

    for r in s.route
        EmptyMove_V2!(s, r, obj)
    end

    return s
end

function EmptyMoveNeighborhood_randOrder(
        s       ::Session                                       ;   # Session to edit
        τ       ::Float64                   = 0.0               ,   # Fault tolerance
        obj     ::Type{<:FitnessSession}    = OutputVariance    ,   # regarded objective
    )
    s = Session(s.Lmax, [Route(r.id, deepcopy(r.assignment), r.mail) for r in s.route], deepcopy(s.load))

    tmpObj = round(fitness(s, obj), digits=3)

    movement::Int64 = 0

    for r_id in randperm(length(s.route))
        r = s.route[r_id]
        newLoad, newAssignment, valid = EmptyMove(s, r, τ=τ, obj=obj)
        if true # valid    # Valid movement
            s.load          = newLoad       # update session's loads
            r.assignment    = newAssignment # update route's assignment
            movement += 1
        end
    end

    return s, movement
end

function EmptyMoveNeighborhood_randOrder_V2(
        s       ::Session                                       ,   # Session to edit
        obj     ::Type{<:FitnessSession}    = OutputVariance    ,   # regarded objective
    )
    s::Session = Session(s.Lmax, [Route(r.id, deepcopy(r.assignment), r.mail) for r in s.route], deepcopy(s.load))
    perm = randperm(length(s.route))
    for i=1:length(s.route)
        EmptyMove(s, s.route[perm[i]])
    end

    return s
end

function EmptyMoveNeighborhood_randDistrib(
        s       ::Session                                       ,   # Session to edit
        obj     ::Type{<:FitnessSession}    = OutputVariance    ,   # regarded objective
    )
    s = Session(s.Lmax, [Route(r.id, deepcopy(r.assignment), r.mail) for r in s.route], deepcopy(s.load))

    distrib = [Float64(sum(r.mail)) for r in s.route]
    for k=1:length(s.load)
        if s.load[k] >= s.Lmax
            for (r_id, r) in enumerate(s.route)
                if r.assignment[k] != 0
                    distrib[r_id] += s.Lmax
                end
            end
        end
    end
    distrib ./= sum(distrib)
    for _=1:length(s.route)
        valid = true
        iter = 0
        while valid && iter < 10
            r = s.route[rand(Categorical(distrib))]
            newLoad, newAssignment, valid = EmptyMove(s, r, obj=obj)
            if true # valid    # Valid movement
                s.load          = newLoad       # update session's loads
                r.assignment    = newAssignment # update route's assignment
            end
            iter += 1
        end
    end

    return s
end

function EmptyMoveNeighborhood_randDistrib_V2(
        s       ::Session                                       ,   # Session to edit
        obj     ::Type{<:FitnessSession}    = OutputVariance    ,   # regarded objective
    )
    s = Session(s.Lmax, [Route(r.id, deepcopy(r.assignment), r.mail) for r in s.route], deepcopy(s.load))

    distrib = ones(Float64, length(s.route))
    for k=1:length(s.load)
        if s.load[k] >= s.Lmax
            for (r_id, r) in enumerate(s.route)
                if r.assignment[k] != 0
                    distrib[r_id] += 1
                end
            end
        end
    end

    iter = 0
    while iter <= 2 * length(s.route)
        r = s.route[rand(Categorical(distrib ./ sum(distrib)))]
        edited = EmptyMove_V3!(s, r, obj)

        if iter < length(s.route)
            for (k, e) in enumerate(edited)
                if e != 0
                    for (r_id, r) in enumerate(s.route)
                        if r.assignment[k] != 0
                            distrib[r_id] += e
                            if distrib[r_id] < 1
                                distrib[r_id] = 1
                            end
                        end
                    end
                end
            end
        end

        iter += 1
    end

    return s
end

function EmptyMoveNeighborhood_oneMove(
        s       ::Session                                       ,   # initial session
        r       ::Int64                                         ,   # route on which to do the move on 
        obj     ::Type{<:FitnessSession}    = OutputVariance    ,   # regarded objective
    )
    s = Session(s.Lmax, [Route(r.id, deepcopy(r.assignment), r.mail) for r in s.route], deepcopy(s.load))

    EmptyMove_V2!(s, s.route[r], obj)

    return s
end

## ============================================================================================================== ##
##                                                #######   ######                                                ##
##                                               ##        ##    ##                                               ##
##                                                ######   ########                                               ##
##                                                     ##  ##    ##                                               ##
##                                               #######   ##    ##                                               ##
## ============================================================================================================== ##

function SimulatedAnnealing_V2(
        s_init          ::Session                                   ; # initial solution
        obj             ::Type{<:FitnessSession} = OutputVariance   , # regarded objective
        display_plot    ::Bool = false                              , # generate plots
        display_state   ::Bool = false                              , # debug info
    )

    Lmax    ::Int64     = s_init.Lmax

    τ_T     ::Float64   = .35                                   # allowed error %
    T       ::Float64   = fitness(s_init, obj) * .2 * τ_T       # initial temperature
    α       ::Float64   = 0.96                                  # cooling ratio

    Emax    ::Int64     = round(Int64, length(s_init.load))/2   # maximum number of epochs
    j       ::Int64     = 1                                     # current epoch
    Mmax    ::Int64     = round(Int64, length(s_init.load))/4   # maximum number of moves per epoch
    i       ::Int64     = 1                                     # iter in epoch
    
    s_best  ::Session   = s_init                                # best solution
    s_curr  ::Session   = s_init                                # current solution

    nogood  ::Int64     = 0                                     # number of consecutive rejected solution
    nobest  ::Int64     = 0                                     # number of consecutive rejected solution

    # =====< Results tracking >=====
    # [1, 1:Emax]  -> T
    # [2, 1:Emax]  -> TODO
    # [3, 1:Emax]  -> # kept sol
    # [4, 1:Emax]  -> ∑ kept sol obj
    # [5, 1:Emax]  -> # degrading sol
    # [6, 1:Emax]  -> # kept degrading sol
    # [7, 1:Emax]  -> ∑ degrading sol obj kept
    # [8, 1:Emax]  -> obj at 0.6% accept chance regarding average s_curr
    # [9, 1:Emax] -> obj at 50% accept chance regarding average s_curr
    # [10, 1:Emax] -> obj(s_best)
    # [11, 1:Emax] -> worst s* of the eppoch
    # [12, 1:Emax] -> best s* of the eppoch
    # [13, 1:Emax] -> s_current at end of epoch
    resEpoch::Array{Union{Float64, Nothing}} = Array{Union{Float64, Nothing}}(nothing, 13, Emax)
    resEpoch[1:10, :] .= 0.

    # =====< Results tracking >=====
    # [1, 1:Emax]  -> s_star
    # [2, 1:Emax]  -> s_curr
    # [3, 1:Emax]  -> s_best
    resIter::Array{Union{Float64, Nothing}} = Array{Union{Float64, Nothing}}(nothing, 3, Emax * Mmax)

    while true
        if i > Mmax
            resEpoch[1, j] = T                       # store epoche temperature
            resEpoch[8, j] = 5 * T + (resEpoch[4, j] / resEpoch[3, j])
            resEpoch[9, j] = .69 * T + (resEpoch[4, j] / resEpoch[3, j])
            resEpoch[10, j] = fitness(s_best, obj)   # store epoche obj(s_best)
            resEpoch[13, j] = fitness(s_curr, obj)   # store epoche obj(s_best)

            (display_state) && (print("\x1b[3F\x1b[1GEpoch $(j)/$(Emax):\n - temperature: $(round(T, digits=3))          \n - avg kept sol obj: $((resEpoch[3, j] != 0) ? (round(resEpoch[4, j] / resEpoch[3, j], digits=3)) : "∅")          \n - best sol obj: $(round(resEpoch[10, j], digits=3))          "))
            
            if j < Emax
                j += 1                                          # next epoch
                i = 1
                
                T *= α                                          # temperature cooling
            else
                break
            end
        else
            s_star = s_curr
            # s_star, _ = EmptyMoveNeighborhood_cascade(s_star, obj=obj) # movement in neighborhood
            # s_star, _ = EmptyMoveNeighborhood_randOrder(s_star, obj=obj)
            # s_star, _ = EmptyMoveNeighborhood_randDistrib(s_star, obj=obj)
            s_star = EmptyMoveNeighborhood_randOrder_V2(s_star, obj) # movement in neighborhood
            s_star = EmptyMoveNeighborhood_randDistrib_V2(s_star, obj) # movement in neighborhood

            ((resEpoch[11, j] === nothing) || (resEpoch[11, j] > fitness(s_star, obj))) && (resEpoch[11, j] = fitness(s_star, obj))
            ((resEpoch[12, j] === nothing) || (resEpoch[12, j] < fitness(s_star, obj))) && (resEpoch[12, j] = fitness(s_star, obj))

            if fitness(s_star, obj) < fitness(s_best, obj)  # computed solution better than best solution
                s_best = s_star                             # update best solution
                nobest = 0
            else
                nobest += 1
            end
            
            Δ = fitness(s_star, obj) - fitness(s_curr, obj)
            if (Δ ≤ 0) || (rand() < exp(-Δ/T))
                s_curr = s_star
                nogood = 0

                resEpoch[3, j] += 1.
                resEpoch[4, j] += fitness(s_star, obj)

                if (Δ > 0)
                    resEpoch[6, j] += 1.
                    resEpoch[7, j] += fitness(s_star, obj)
                end
            else
                nogood += 1
            end
            if (Δ > 0)
                resEpoch[5, j] += 1.
            end

            resIter[1, (j-1)*Mmax + i] = fitness(s_star, obj)
            resIter[2, (j-1)*Mmax + i] = fitness(s_curr, obj)
            resIter[3, (j-1)*Mmax + i] = fitness(s_best, obj)

            if isSessionValid(s_star) # feasibility check
                break
            end

            if nogood >= Mmax                                   # temperature is too low
                if nobest >= Mmax * Emax * .3
                    τ_T *= 1.15                                 # 125 % of previous temperature
                    T = fitness(s_init, obj) * .2 * τ_T         # warmum temperature
                    #α -= .01                                    # decrease cooling ration

                    nobest = 0
                else
                    τ_T *= 0.85                                 # 85 % of previous temperature
                    T = fitness(s_init, obj) * .2 * τ_T         # warmum temperature
                    α -= .005                                    # increase cooling ration
                end
                nogood = 0
            end

            i += 1
        end
    end

    if display_plot
        replace!(resEpoch, nothing=>0.)
        replace!(resIter, nothing=>0.)

        x = 1:Emax

        p1 = plot(x, resEpoch[1, :], lc="blue", lw=1.5, title="Temperature", label="T")
        (isSessionValid(s_best)) && (plot!([j, j], [minimum(resEpoch[1, :]), maximum(resEpoch[1, :])], lc="black", lw=1.5, label="opt reached"))

        p2 = plot(x, resEpoch[6, :] ./ resEpoch[5, :] .* 100, title = "Acceptance ratio", label="% of degrading sol used", lw=2)
        plot!(x, resEpoch[3, :] ./ Emax .* 100, label="% of sol used", lw=2)
        # plot!(x, resEpoch[5, :] ./ Emax .* 100, label="% of degrading sol", lw=2)
        (isSessionValid(s_best)) && (plot!([j, j], [minimum(resEpoch[6, :] ./ resEpoch[5, :] .* 100), maximum(resEpoch[6, :] ./ resEpoch[5, :] .* 100)], lc="black", lw=1.5, label="opt reached"))

        x = 2:Emax
        p3 = plot(title="Avg epoch solutions", label="valid sol threshold")
        plot!(x, [resEpoch[11, 2:end], resEpoch[12, 2:end]], fillrange = Vector{Float64}(resEpoch[12, 2:end]), fa = .1, fc="blue", lc="blue", lw=1.5, label=["range of s* values" ""])
        plot!(x, [resEpoch[9, 2:end], (resEpoch[4, 2:end] ./ resEpoch[3, 2:end])], fillrange = Vector{Float64}((resEpoch[4, 2:end] ./ resEpoch[3, 2:end])), fa = .3, fc="orange", lc="orange", lw=1.5, label=["obj with [50, 100]% acceptation chance" ""])
        plot!(x, resEpoch[4, 2:end] ./ resEpoch[3, 2:end], lc="green", lw=1.5, label="avg s* kept")
        # plot!(x, resEpoch[13, 2:end], lc="brown", lw=1.5, label="s_curr at end of epoch")
        plot!(x, resEpoch[10, 2:end], lc="purple", lw=1.5, label="s_best")
        (isSessionValid(s_best)) && (plot!([j, j], [minimum(resEpoch[10, 2:end]), maximum(resEpoch[11, 2:end])], lc="black", lw=1.5, label="opt reached"))

        x = Mmax:(Emax*Mmax)
        p4 = plot(x, [resIter[1, Mmax:end] resIter[2, Mmax:end] resIter[3, Mmax:end]], lw=1.5, title="Iter solutions", label=["s*" "s_curr" "s_best"])
        (isSessionValid(s_best)) && (plot!([j, j], [minimum(resIter[3, Mmax:end]), maximum(resIter[1, Mmax:end])], lc="black", lw=1.5, label="opt reached"))


        display(plot(p1, p2, p3, p4, layout=(4,1), size = (1200, 1600)))
    end

    return s_best, isSessionValid(s_best), (resEpoch, resIter)
end

function SimulatedAnnealing_V3(
        s_init          ::Session                                   ; # initial solution
        obj             ::Type{<:FitnessSession} = LoadSTD          , # regarded objective
        display_plot    ::Bool = false                              , # generate plots
        display_state   ::Bool = false                              , # debug info
        τ               ::Float64 = .5                              , # initial allowed error %
        α               ::Float64 = .98                             , # # cooling ratio
        Emax            ::Int64   = round(Int64, length(s_init.load))*2 , # maximum number of epochs
        Mmax            ::Int64   = round(Int64, length(s_init.load))*2 , # maximum number of moves per epoch
    )

    Lmax    ::Int64     = s_init.Lmax
    R       ::Int64     = length(s_init.route)

    τ_T     ::Float64   = τ                                  # allowed error %
    T       ::Float64   = fitness(s_init, obj) * .2 * τ_T         # initial temperature
    # α       ::Float64   = 0.98                                # cooling ratio

    j       ::Int64     = 1                                     # current epoch
    i       ::Int64     = 1                                     # iter in epoch

    s_best  ::Session   = s_init                                # best solution
    s_curr  ::Session   = s_init                                # current solution

    nogood  ::Int64     = 0                                     # number of consecutive rejected solution

    route_id::Int64     = R+1
    order::Vector{Int64} = []

    # println("NEW RUN -> Emax = $Emax, Mmax=$Mmax, T=$T, i = $i, j = $j, obj = $(round(fitness(s_init, obj), digits=3))")

    # =====< Results tracking >=====
    # [1, 1:Emax]  -> T
    # [2, 1:Emax]  -> TODO
    # [3, 1:Emax]  -> # kept sol
    # [4, 1:Emax]  -> ∑ kept sol obj
    # [5, 1:Emax]  -> # degrading sol
    # [6, 1:Emax]  -> # kept degrading sol
    # [7, 1:Emax]  -> ∑ degrading sol obj kept
    # [8, 1:Emax]  -> obj at 0.6% accept chance regarding average s_curr
    # [9, 1:Emax] -> obj at 50% accept chance regarding average s_curr
    # [10, 1:Emax] -> obj(s_best)
    # [11, 1:Emax] -> worst s* of the eppoch
    # [12, 1:Emax] -> best s* of the eppoch
    # [13, 1:Emax] -> s_current at end of epoch
    resEpoch::Array{Union{Float64, Nothing}} = Array{Union{Float64, Nothing}}(nothing, 13, Emax)
    resEpoch[1:10, :] .= 0.

    # =====< Results tracking >=====
    # [1, 1:Emax]  -> s_star
    # [2, 1:Emax]  -> s_curr
    # [3, 1:Emax]  -> s_best
    resIter::Array{Union{Float64, Nothing}} = Array{Union{Float64, Nothing}}(nothing, 3, Emax * Mmax)

    (display_state) && (println("\n\n"))

    while true
        if i > Mmax
            resEpoch[1, j] = T                       # store epoche temperature
            resEpoch[8, j] = 5 * T + (resEpoch[4, j] / resEpoch[3, j])
            resEpoch[9, j] = .69 * T + (resEpoch[4, j] / resEpoch[3, j])
            resEpoch[10, j] = fitness(s_best, obj)   # store epoche obj(s_best)
            resEpoch[13, j] = fitness(s_curr, obj)   # store epoche obj(s_best)

            (display_state) && (print("\x1b[3F\x1b[1GEpoch $(j)/$(Emax):\n - temperature: $(round(T, digits=3))          \n - avg kept sol obj: $((resEpoch[3, j] != 0) ? (round(resEpoch[4, j] / resEpoch[3, j], digits=3)) : "∅")          \n - best sol obj: $(round(resEpoch[10, j], digits=3))          "))
            
            if j < Emax
                j += 1                                          # next epoch
                i = 1
                
                T *= α                                          # temperature cooling
            else
                # println("\n\nEND ITER\n")
                break
            end
        else
            # route_id = (route_id % R) +1

            if route_id > R-1
                order = randperm(R)
                route_id = 1
            else
                route_id += 1
            end
            s_star = EmptyMoveNeighborhood_oneMove(s_curr, order[route_id], obj)
            
            ((resEpoch[11, j] === nothing) || (resEpoch[11, j] > fitness(s_star, obj))) && (resEpoch[11, j] = fitness(s_star, obj))
            ((resEpoch[12, j] === nothing) || (resEpoch[12, j] < fitness(s_star, obj))) && (resEpoch[12, j] = fitness(s_star, obj))

            if fitness(s_star, obj) < fitness(s_best, obj)  # computed solution better than best solution
                s_best = s_star                             # update best solution
            end
            
            Δ = fitness(s_star, obj) - fitness(s_curr, obj)
            if (Δ ≤ 0) || (rand() < exp(-Δ/T))
                s_curr = s_star
                nogood = 0

                resEpoch[3, j] += 1.
                resEpoch[4, j] += fitness(s_star, obj)

                if (Δ > 0)
                    resEpoch[6, j] += 1.
                    resEpoch[7, j] += fitness(s_star, obj)
                end
            else
                nogood += 1
            end
            if (Δ > 0)
                resEpoch[5, j] += 1.
            end

            resIter[1, (j-1)*Mmax + i] = fitness(s_star, obj)
            resIter[2, (j-1)*Mmax + i] = fitness(s_curr, obj)
            resIter[3, (j-1)*Mmax + i] = fitness(s_best, obj)

            # if isSessionValid(s_star) # feasibility check
            #     s_best = s_star
            #     # println("\n\nEND VALID\n")
            #     break
            # end

            if nogood >= Emax# 2*R # || T <= 0.001                                   # temperature is too low
                τ_T *= 0.85                                 # 85 % of previous temperature
                T = fitness(s_init, obj) * .2 * τ_T         # warmum temperature
                α -= .005                                    # increase cooling ration

                nogood = 0
            end

            i += 1
        end
    end

    if display_plot
        replace!(resEpoch, nothing=>0.)
        replace!(resIter, nothing=>0.)

        x = 1:Emax

        p1 = plot(x, resEpoch[1, :], lc="blue", lw=1.5, title="Temperature", label="T")
        (isSessionValid(s_best)) && (plot!([j, j], [minimum(resEpoch[1, :]), maximum(resEpoch[1, :])], lc="black", lw=1.5, label="opt reached"))

        p2 = plot(x, resEpoch[6, :] ./ resEpoch[5, :] .* 100, title = "Acceptance ratio", label="% of degrading sol used", lw=2)
        plot!(x, resEpoch[3, :] ./ Emax .* 100, label="% of sol used", lw=2)
        # plot!(x, resEpoch[5, :] ./ Emax .* 100, label="% of degrading sol", lw=2)
        (isSessionValid(s_best)) && (plot!([j, j], [minimum(resEpoch[6, :] ./ resEpoch[5, :] .* 100), maximum(resEpoch[6, :] ./ resEpoch[5, :] .* 100)], lc="black", lw=1.5, label="opt reached"))

        x = 2:Emax
        p3 = plot(title="Avg epoch solutions", label="valid sol threshold")
        plot!(x, [resEpoch[11, 2:end], resEpoch[12, 2:end]], fillrange = Vector{Float64}(resEpoch[12, 2:end]), fa = .1, fc="blue", lc="blue", lw=1.5, label=["range of s* values" ""])
        plot!(x, [resEpoch[9, 2:end], (resEpoch[4, 2:end] ./ resEpoch[3, 2:end])], fillrange = Vector{Float64}((resEpoch[4, 2:end] ./ resEpoch[3, 2:end])), fa = .3, fc="orange", lc="orange", lw=1.5, label=["obj with [50, 100]% acceptation chance" ""])
        plot!(x, resEpoch[4, 2:end] ./ resEpoch[3, 2:end], lc="green", lw=1.5, label="avg s* kept")
        # plot!(x, resEpoch[13, 2:end], lc="brown", lw=1.5, label="s_curr at end of epoch")
        plot!(x, resEpoch[10, 2:end], lc="purple", lw=1.5, label="s_best")
        (isSessionValid(s_best)) && (plot!([j, j], [minimum(resEpoch[10, 2:end]), maximum(resEpoch[11, 2:end])], lc="black", lw=1.5, label="opt reached"))

        x = Mmax:(Emax*Mmax)
        p4 = plot(x, [resIter[1, Mmax:end] resIter[2, Mmax:end] resIter[3, Mmax:end]], lw=1.5, title="Iter solutions", label=["s*" "s_curr" "s_best"])
        (isSessionValid(s_best)) && (plot!([j, j], [minimum(resIter[3, Mmax:end]), maximum(resIter[1, Mmax:end])], lc="black", lw=1.5, label="opt reached"))


        display(plot(p1, p2, p3, p4, layout=(4,1), size = (1200, 1600)))
    end

    return s_best, isSessionValid(s_best), (resEpoch, resIter)
end

function SA_V4(
        s_init          ::Session                                   ; # initial solution
        obj1            ::Type{<:FitnessSession} = LoadSTD          , # regarded objective
        obj2            ::Type{<:FitnessSession} = MostLoadedOut    , # regarded objective
        display_plot    ::Bool = false                              , # generate plots
        display_state   ::Bool = false                              , # debug info
        τ               ::Float64 = .5                              , # initial allowed error %
        α               ::Float64 = .96                             , # # cooling ratio
        Emax            ::Int64   = round(Int64, length(s_init.load))*2 , # maximum number of epochs
        Mmax            ::Int64   = round(Int64, length(s_init.load))*2 , # maximum number of moves per epoch
    )

    Lmax    ::Int64     = s_init.Lmax
    R       ::Int64     = length(s_init.route)

    τ_T     ::Float64   = τ                                     # allowed error %
    T       ::Float64   = fitness(s_init, obj1) * .2 * τ_T      # initial temperature
    # α       ::Float64   = 0.98                                # cooling ratio

    j       ::Int64     = 1                                     # current epoch
    i       ::Int64     = 1                                     # iter in epoch

    s_best1 ::Session   = s_init                                # best solution
    s_best2 ::Session   = s_init                                # best solution
    s_curr  ::Session   = s_init                                # current solution

    nogood  ::Int64     = 0                                     # number of consecutive rejected solution

    route_id::Int64     = R+1
    order::Vector{Int64} = []

    # println("NEW RUN -> Emax = $Emax, Mmax=$Mmax, T=$T, i = $i, j = $j, obj = $(round(fitness(s_init, obj), digits=3))")

    # =====< Results tracking >=====
    # [1, 1:Emax]  -> T
    # [2, 1:Emax]  -> sc
    # [3, 1:Emax]  -> # kept sol
    # [4, 1:Emax]  -> ∑ kept sol obj
    # [5, 1:Emax]  -> # degrading sol
    # [6, 1:Emax]  -> # kept degrading sol
    # [7, 1:Emax]  -> ∑ degrading sol obj kept
    # [8, 1:Emax]  -> obj at 0.6% accept chance regarding average s_curr
    # [9, 1:Emax]  -> obj at 50% accept chance regarding average s_curr
    # [10, 1:Emax] -> obj(s_best)
    # [11, 1:Emax] -> worst s* of the eppoch
    # [12, 1:Emax] -> best s* of the eppoch
    # [13, 1:Emax] -> s_current at end of epoch
    resEpoch::Array{Union{Float64, Nothing}} = Array{Union{Float64, Nothing}}(nothing, 13, Emax)
    resEpoch[1:10, :] .= 0.

    # =====< Results tracking >=====
    # [1, 1:Emax]  -> s_star
    # [2, 1:Emax]  -> s_curr
    # [3, 1:Emax]  -> s_best
    resIter::Array{Union{Float64, Nothing}} = Array{Union{Float64, Nothing}}(nothing, 3, Emax * Mmax)

    (display_state) && (println("\n\n"))

    while true
        if i > Mmax
            resEpoch[1, j] = T                       # store epoche temperature
            resEpoch[8, j] = 5 * T + (resEpoch[4, j] / resEpoch[3, j])
            resEpoch[9, j] = .69 * T + (resEpoch[4, j] / resEpoch[3, j])
            resEpoch[10, j] = fitness(s_best1, obj1)   # store epoche obj(s_best)
            resEpoch[13, j] = fitness(s_curr, obj1)   # store epoche obj(s_best)

            (display_state) && (print("\x1b[3F\x1b[1GEpoch $(j)/$(Emax):\n - temperature: $(round(T, digits=3))          \n - avg kept sol obj: $((resEpoch[3, j] != 0) ? (round(resEpoch[4, j] / resEpoch[3, j], digits=3)) : "∅")          \n - best sol obj: $(round(resEpoch[10, j], digits=3))          "))
            
            if j < Emax
                j += 1                                          # next epoch
                i = 1
                
                T *= α                                          # temperature cooling
            else
                # println("\n\nEND ITER\n")
                break
            end
        else
            # route_id = (route_id % R) +1

            if route_id > R-1
                order = randperm(R)
                route_id = 1
            else
                route_id += 1
            end
            s_star = EmptyMoveNeighborhood_oneMove(s_curr, order[route_id], obj1)
            
            ((resEpoch[11, j] === nothing) || (resEpoch[11, j] > fitness(s_star, obj1))) && (resEpoch[11, j] = fitness(s_star, obj1))
            ((resEpoch[12, j] === nothing) || (resEpoch[12, j] < fitness(s_star, obj1))) && (resEpoch[12, j] = fitness(s_star, obj1))

            if fitness(s_star, obj1) < fitness(s_best1, obj1)  # computed solution better than best solution
                s_best1 = s_star                             # update best solution
            end

            if fitness(s_star, obj2) < fitness(s_best2, obj2)  # computed solution better than best solution
                s_best2 = s_star                             # update best solution
            end
            
            Δ = fitness(s_star, obj1) - fitness(s_curr, obj1)
            if (Δ ≤ 0) || (rand() < exp(-Δ/T))
                s_curr = s_star
                nogood = 0

                resEpoch[3, j] += 1.
                resEpoch[4, j] += fitness(s_star, obj1)

                if (Δ > 0)
                    resEpoch[6, j] += 1.
                    resEpoch[7, j] += fitness(s_star, obj1)
                end
            else
                nogood += 1
            end
            if (Δ > 0)
                resEpoch[5, j] += 1.
            end

            resIter[1, (j-1)*Mmax + i] = fitness(s_star, obj1)
            resIter[2, (j-1)*Mmax + i] = fitness(s_curr, obj1)
            resIter[3, (j-1)*Mmax + i] = fitness(s_best1, obj1)

            if isSessionValid(s_star) # feasibility check
                s_best1 = s_star
                s_best2 = s_star
                # println("\n\nEND VALID\n")
                break
            end

            if nogood >= Mmax # 2*R # || T <= 0.001                                   # temperature is too low
                τ_T *= 0.85                                 # 85 % of previous temperature
                T = fitness(s_init, obj1) * .2 * τ_T         # warmum temperature
                α -= .005                                    # increase cooling ration

                nogood = 0
            end

            i += 1
        end
    end

    if display_plot
        replace!(resEpoch, nothing=>0.)
        replace!(resIter, nothing=>0.)

        x = 1:Emax

        p1 = plot(x, resEpoch[1, :], lc="blue", lw=1.5, title="Temperature", label="T")
        (isSessionValid(s_best1)) && (plot!([j, j], [minimum(resEpoch[1, :]), maximum(resEpoch[1, :])], lc="black", lw=1.5, label="opt reached"))

        p2 = plot(x, resEpoch[6, :] ./ resEpoch[5, :] .* 100, title = "Acceptance ratio", label="% of degrading sol used", lw=2)
        plot!(x, resEpoch[3, :] ./ Emax .* 100, label="% of sol used", lw=2)
        # plot!(x, resEpoch[5, :] ./ Emax .* 100, label="% of degrading sol", lw=2)
        (isSessionValid(s_best1)) && (plot!([j, j], [minimum(resEpoch[6, :] ./ resEpoch[5, :] .* 100), maximum(resEpoch[6, :] ./ resEpoch[5, :] .* 100)], lc="black", lw=1.5, label="opt reached"))

        x = 2:Emax
        p3 = plot(title="Avg epoch solutions", label="valid sol threshold")
        plot!(x, [resEpoch[11, 2:end], resEpoch[12, 2:end]], fillrange = Vector{Float64}(resEpoch[12, 2:end]), fa = .1, fc="blue", lc="blue", lw=1.5, label=["range of s* values" ""])
        plot!(x, [resEpoch[9, 2:end], (resEpoch[4, 2:end] ./ resEpoch[3, 2:end])], fillrange = Vector{Float64}((resEpoch[4, 2:end] ./ resEpoch[3, 2:end])), fa = .3, fc="orange", lc="orange", lw=1.5, label=["obj with [50, 100]% acceptation chance" ""])
        plot!(x, resEpoch[4, 2:end] ./ resEpoch[3, 2:end], lc="green", lw=1.5, label="avg s* kept")
        # plot!(x, resEpoch[13, 2:end], lc="brown", lw=1.5, label="s_curr at end of epoch")
        plot!(x, resEpoch[10, 2:end], lc="purple", lw=1.5, label="s_best")
        (isSessionValid(s_best1)) && (plot!([j, j], [minimum(resEpoch[10, 2:end]), maximum(resEpoch[11, 2:end])], lc="black", lw=1.5, label="opt reached"))

        x = Mmax:(Emax*Mmax)
        p4 = plot(x, [resIter[1, Mmax:end] resIter[2, Mmax:end] resIter[3, Mmax:end]], lw=1.5, title="Iter solutions", label=["s*" "s_curr" "s_best"])
        (isSessionValid(s_best1)) && (plot!([j, j], [minimum(resIter[3, Mmax:end]), maximum(resIter[1, Mmax:end])], lc="black", lw=1.5, label="opt reached"))


        display(plot(p1, p2, p3, p4, layout=(4,1), size = (1200, 1600)))
    end

    return s_best1, s_best2, isSessionValid(s_best1), isSessionValid(s_best2), (resEpoch, resIter)
end

function EM_VE(
        s               ::Session                                   ,   # initial session
        obj1            ::Type{<:FitnessSession} = LoadSTD          , # regarded objective
        obj2            ::Type{<:FitnessSession} = MostLoadedOut    , # regarded objective
    )

    i       ::Int64     = 0
    Lmax    ::Int64     = s.Lmax
    F       ::Bool      = true
    sb      ::Session   = s

    while F
        F = false
        for r in s.route
            assignment, load = EmptyMove_V4(s, r, obj1)

            if fitness(load, Lmax, obj1) < fitness(s.load, Lmax, obj1)
                i += 1
                s.load = load
                r.assignment = assignment

                F = true

                if fitness(load, Lmax, obj2) < fitness(s.load, Lmax, obj2)
                    sb = Session(s.Lmax, [Route(r.id, deepcopy(r.assignment), r.mail) for r in s.route], deepcopy(s.load))
                end
            end
        end
    end

    return s, sb, i
end

