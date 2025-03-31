begin
    import Base:copy
    import Random
end

begin # Files 
    include("Instance.jl")
    include("DataStruct.jl")  
    include("OptiMove.jl")  
end

# ================================================== #
#                    Constructor                     #
# ================================================== #

function buildSolution_FF(
        instance            ::Instance                                                                          ; 
        TAG_AddRoute        ::Vector{Type{<:SimpleAddRoute}}    = Type{<:SimpleAddRoute}[SMOOTH_ASSIGNED, LEFT_ALLIGNED, REBUILD_KP_01LP_V2]  , 
        Δ                   ::Int64                             = 2                                             ,
        TAG_FitSes          ::Type{<:FitnessSession}            = LoadSTD                                       ,
        env                 ::Gurobi.Env                        = Gurobi.Env()                                  ,
        tl                  ::Int64                             = 100                                           , 
    )
    return buildSolution_FFD(instance, randperm(instance.nbRoute), Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRoute=TAG_AddRoute, env=env, tl=tl)
end

function buildSolution_FFD(
        instance            ::Instance                                                                          , 
        perm                ::Vector{Int64}                                                                     ;
        TAG_AddRoute        ::Vector{Type{<:SimpleAddRoute}}    = Type{<:SimpleAddRoute}[SAVA, NFBA, OPTIMOVE_STAGE1]  ,  
        Δ                   ::Int64                             = 2                                             ,
        TAG_FitSes          ::Type{<:FitnessSession}            = LoadSTD                                       ,
        env                 ::Gurobi.Env                        = Gurobi.Env()                                  ,
        tl                  ::Int64                             = 100                                           , 
    )

    sol::Solution = Solution(perm, [Session(instance.Lmax, instance.nbOut) for _=1:minSession(instance)])
    # upgrades::Vector{Int64} = zeros(Int64, length(TAG_AddRoute))

    for routeId in perm
        println(" - route: $routeId")
        sId         ::Int64     = 1
        added       ::Bool      = false
        r           ::Route     = instance.route[routeId]

        while !added
            if sId <= length(sol.sessions)
                println("   - session: $sId")
                # print("(sava_nt), ")
                # st          ::Union{Session, Nothing} = nothing
                # flag        ::Bool      = false
                flagId      ::Int64     = 1

                while !added && flagId <= length(TAG_AddRoute)
                    # st, added = addRoute(sol.sessions[sId], r, TAG_AddRoute[flagId], Δ=Δ, TAG_FitSes=TAG_FitSes, tl=tl, env=env)
                    println("     - add TAG : $(TAG_AddRoute[flagId])")
                    res, _ = addRoute(sol.sessions[sId], r, TAG_AddRoute[flagId], Δ, TAG_FitSes, tl, env)

                    # flag && (upgrades[flagId] += 1)
                    flagId += 1
                end
                sId += 1
            else # create new Session
                push!(sol.sessions, Session(instance.Lmax, [r], deepcopy(r.assignment)))

                added = true
            end
        end
    end

    return sol # , upgrades
end

function buildSolution_FF_final(instance::Instance, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    return buildSolution_FFD_final(instance, randperm(instance.nbRoute), tl, env)
end

function buildSolution_FFD_final(instance::Instance, perm::Vector{Int64}, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    sol::Solution = Solution(perm, [Session(instance.Lmax, instance.nbOut) for _=1:minSession(instance)])
    for routeId in perm
        r::Route = instance.route[routeId]
        sol = InsertRoute_Procedure!(sol, r, tl, env)
    end
    return sol
end

function buildSolution_BF_final(instance::Instance, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    return buildSolution_BFD_final(instance, randperm(instance.nbRoute), tl, env)
end

function buildSolution_BFD_final(instance::Instance, perm::Vector{Int64}, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    sol::Solution = Solution(perm, [Session(instance.Lmax, instance.nbOut)])
    for routeId in perm
        r::Route = instance.route[routeId]
        # println("new order -> $(sortperm(sol.sessions, by=x -> (100 - fitness(x, HollowedPercentage))))")
        sort!(sol.sessions, by=x -> (fitness(x, HollowedPercentage)))
        sol = InsertRoute_Procedure!(sol, r, tl, env)
    end
    return sol
end

function buildSolution_WF_final(instance::Instance, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    return buildSolution_WFD_final(instance, randperm(instance.nbRoute), tl, env)
end

function buildSolution_WFD_final(instance::Instance, perm::Vector{Int64}, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    sol::Solution = Solution(perm, [Session(instance.Lmax, instance.nbOut)])
    for routeId in perm
        r::Route = instance.route[routeId]
        # println("new order -> $(sortperm(sol.sessions, by=x -> (100 - fitness(x, HollowedPercentage))))")
        sort!(sol.sessions, by=x -> (-fitness(x, HollowedPercentage)))
        sol = InsertRoute_Procedure!(sol, r, tl, env)
    end
    return sol
end

"""
function InsertRoute_Procedure!(sol::Solution, r::Route, env::Gurobi.Env = Gurobi.Env(), tl::Int64 = 10)

Full procedure to insert a route r into a solution s, following a First-fit like behavior.
"""
@inline function InsertRoute_Procedure!(sol::Solution, r::Route, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    sId         ::Int64     = 1
    added       ::Bool      = false

    while !added                                                                                                # While r hasn't been added to any sorting session
        if sId <= length(sol.sessions)                                                                          #  |  if any more session
            s::Session = sol.sessions[sId]                                                                      #  |   |
            if length(s.route) < length(s.load) && certificat_CapacityVolume(s, r)                              #  |   | if Volume and number of route certificates
                newAssignment::Vector{Int64}, added = addRoute_LeftAlligned(s, r)                               #  |   |   |  Try adding r using left alligned assignment
                if added                                                                                        #  |   |   |  if Left alligned assignment valid
                    s, added = addRoute_SmoothAssigned!(s, r)                                                   #  |   |   |   |  Try adding r using a smoother assignment
                    if !added                                                                                   #  |   |   |   |  if smoother assignment not valid
                        r = Route(r.id, newAssignment, r.mail)                                                  #  |   |   |   |   |  use Left alligned assignment
                        s.load += newAssignment                                                                 #  |   |   |   |   |
                        push!(s.route, r)                                                                       #  |   |   |   |   |
                        added = true                                                                            #  |   |   |   |   |
                    end                                                                                         #  |   |   |   |  end
                else                                                                                            #  |   |   |  else
                    s, added = addRoute_Rebuild_Knapsack_model_V3!(s, r, tl, env)                               #  |   |   |   |  rebuild session
                end                                                                                             #  |   |   |  endif
            end                                                                                                 #  |   |  endif
            sId += 1                                                                                            #  |   |  Pass to next Session (relevent only r has not been added yet)
        else                                                                                                    #  |  else
            push!(sol.sessions, Session(sol.sessions[1].Lmax, [r], deepcopy(r.assignment)))                     #  |   |  Create a new Session with r inside (will end the while loop)
            added = true                                                                                        #  |   |
        end                                                                                                     #  |  end
    end                                                                                                         # end
    return sol                                                                                                  # return solution
end

@inline function insert_route!(s::Session, r::Route)::Bool
    added::Bool = false
    
    if length(s.route) < length(s.load) && certificat_CapacityVolume(s, r)                             #  if Volume and number of route certificates
        added, newAssignment::Vector{Int64} = addRoute_NFBA(s, r)                                       #   |  Try adding r using left alligned assignment
        if added                                                                                        #   |  if Left alligned assignment valid
            s, added = addRoute_SAVA!(s, r)                                                             #   |   |  Try adding r using a smoother assignment
            if !added                                                                                   #   |   |  if smoother assignment not valid
                r = Route(r.id, newAssignment, r.mail)                                                  #   |   |   |  use Left alligned assignment
                s.load += newAssignment                                                                #   |   |   |
                push!(s.route, r)                                                                       #   |   |   |
                added = true                                                                            #   |   |   |
            end                                                                                         #   |   |  end
        else                                                                                            #   |  else
            s, added = addRoute_Rebuild_Knapsack_model_V3!(s, r, tl, env)                               #   |   |  rebuild session
        end                                                                                             #   |  endif
    end                                                                                                 #  endif

    return added
end

## ============================================================================================================== ##
##                           ######    ########  ######              ########  ##    ##                           ##
##                           ##    ##  ##        ##    ##            ##        ###  ###                           ##
##                           #######   #####     ##     #    ####    #####     ## ## ##                           ##
##                           ##    ##  ##        ##    ##            ##        ##    ##                           ##
##                           #######   ##        ######              ########  ##    ##                           ##
## ============================================================================================================== ##

function BFD_EmptyMove(
        Routes  ::Vector{Route}     , # Routes of the instance
        Lmax    ::Int64             , # Maximum capacity of an output
        O       ::Int64 = nothing   , # Number of output of each session
        R       ::Int64 = nothing   , # Number of route to sort
    )::Vector{Session}
    
    # Solution (Initialisation with an empty session) 
    S::Vector{Session} = [Session(Lmax, O)]

    # permutation = order to consider routes
    perm::Vector{Int64} = sortperm([(-fitness(r, MaxMin), -fitness(r, MailNb)) for r in Routes])

    for route_id in perm

        # Sort sessions from most to least loaded (mail volume)
        sort!(S, by=x -> (fitness(x, HollowedPercentage)))
        
        added               ::Bool  = false
        current_session_id  ::Int64 = 1
        r                   ::Route = Routes[route_id]

        while !added
            if current_session_id <= length(S)
                
                # try inserting the route into the current session
                s::Session = S[current_session_id]

                if length(s.route) < length(s.load) && certificat_CapacityVolume(s, r)

                    newAssignment::Vector{Int64}, added = addRoute_LeftAlligned(s, r)

                    if added

                        s, added = addRoute_SmoothAssigned!(s, r)

                        if !added
                            # Smooth assigned failed use the left alligned valid assignemnent
                            r = Route(r.id, newAssignment, r.mail)
                            s.load += newAssignment
                            push!(s.route, r)
                            added = true
                        end
                    else
                        # heavy procedure EMPTY-MOVE
                        new_session::Session = Session(Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
                        push!(new_session.route, Route(r.id, deepcopy(r.assignment), r.mail))

                        # println("\n\x1b[31m=====< Before >=====\n\x1b[0m$new_session")

                        new_session, _, added = improvedOptiMove_V1(new_session)

                        isSessionValid(new_session) ? print("\x1b[31mo\x1b[0m") : print("\x1b[31m-\x1b[0m")

                        # if Empty move managed to insert the route update s
                        (added) && (S[current_session_id] = new_session)
                    end
                end
            else
                # open a new session with the current route inside
                push!(S, Session(Lmax, [r], deepcopy(r.assignment)))
                added = true
            end

            current_session_id += 1
        end
    end

    return S
end

## ============================================================================================================== ##
##                           ######    ########  ######               ######   #######                            ##
##                           ##    ##  ##        ##    ##            ##        ##    ##                           ##
##                           #######   #####     ##     #    ####    ##  ###   #######                            ##
##                           ##    ##  ##        ##    ##            ##    ##  ##  ##                             ##
##                           #######   ##        ######               ######   ##   ##                            ##
## ============================================================================================================== ##

function BFD_GreadyRebuild(
        Routes  ::Vector{Route}             , # Routes of the instance
        Lmax    ::Int64                     , # Maximum capacity of an output
        O       ::Int64 = nothing           , # Number of output of each session
        R       ::Int64 = nothing           , # Number of route to sort
        tl      ::Int64 = 10                , # Model time limite
        env     ::Gurobi.Env = Gurobi.Env() , # Gurobi environement (for display purpose)
    )::Vector{Session}

    # Solution (Initialisation with an empty session) 
    S::Vector{Session} = [Session(Lmax, O)]

    # permutation = order to consider routes
    perm::Vector{Int64} = sortperm([(-fitness(r, MaxMin), -fitness(r, MailNb)) for r in Routes])

    for route_id in perm

        # Sort sessions from most to least loaded (mail volume)
        sort!(S, by=x -> (fitness(x, HollowedPercentage)))
        
        added               ::Bool  = false
        current_session_id  ::Int64 = 1
        r                   ::Route = Routes[route_id]

        while !added
            if current_session_id <= length(S)
                
                # try inserting the route into the current session
                s::Session = S[current_session_id]

                if length(s.route) < length(s.load) && certificat_CapacityVolume(s, r)

                    newAssignment::Vector{Int64}, added = addRoute_LeftAlligned(s, r)

                    if added

                        s, added = addRoute_SmoothAssigned!(s, r)

                        if !added
                            # Smooth assigned failed use the left alligned valid assignemnent
                            r = Route(r.id, newAssignment, r.mail)
                            s.load += newAssignment
                            push!(s.route, r)
                            added = true
                        end
                    else
                        # heavy procedure EMPTY-MOVE
                        new_session::Session = Session(Lmax, [Route(cr.id, deepcopy(cr.assignment), cr.mail) for cr in s.route], s.load + r.assignment)
                        push!(new_session.route, Route(r.id, deepcopy(r.assignment), r.mail))

                        new_session, added = rebuildSession_knapSack_model_V3!(new_session, tl, env)

                        isSessionValid(new_session) ? print("\x1b[32mo\x1b[0m") : print("\x1b[32m-\x1b[0m")

                        # if Empty move managed to insert the route update s
                        (added) && (S[current_session_id] = new_session)
                    end
                end
            else
                # open a new session with the current route inside
                push!(S, Session(Lmax, [r], deepcopy(r.assignment)))
                added = true
            end

            current_session_id += 1
        end
    end

    return S
end

## ============================================================================================================== ##
##             /###     ######              ######    ########  ######               ######   #######             ##
##             # ##     ##    ##            ##    ##  ##        ##    ##    ##      ##        ##    ##            ##
##               ##     ##     #    ####    #######   #####     ##     #      ##    ##  ###   #######             ##
##               ##     ##    ##            ##    ##  ##        ##    ##    ##      ##    ##  ##  ##              ##
##             ######   ######              #######   ##        ######               ######   ##   ##             ##
## ============================================================================================================== ##

function BF_1D(
        inst    ::tInstance                                 ; # Instance 
        lb      ::Union{Int64, Nothing}         = nothing   , # Lower bound
        order   ::Union{Vector{Int64}, Nothing} = nothing   , # Order to consider items
    )::Vector{tBin}

    I       ::Vector{tItem} = inst.I                                                    # Set of item to partition
    C       ::Int64         = inst.C                                                    # Maximum bin capacity
    order   ::Vector{Int64} = randperm(length(inst.I))                                  # Order to consider item
    (lb == nothing) && (lb = ceil(Int64, sum([i.size for i::tItem in I]) / C))          # LB on the bin number

    # Solution (set of bin ≡ item partition)
    sol     ::Vector{tBin}  = tBin[tBin(Vector{tItem}(undef, 0), 0, C) for _=1:lb]

    for k in order
        i       ::tItem     = inst.I[k] # Current Item
        added   ::Bool      = false     # Is the current item inserted into a bin ?
        b_id    ::Int64     = 1         # Current session id
        
        sort!(sol, by=(x -> -x.L))

        while !added
            if b_id ≤ length(sol)
                b::tBin = sol[b_id]     # Current bin

                if i.size + b.L ≤ C
                    # Insert item in the current bin
                    push!(b.I, i)
                    b.L += i.size

                    added = true
                else
                    # Pass to the next bin
                    b_id += 1
                end
            else
                # Open a new bin
                push!(sol, tBin([i], i.size, C))

                added = true
            end
        end
    end

    return sol
end

function BF_1D(I::Vector{Int64}, C::Int64)::Vector{tBin}
    return BF_1D(tInstance(enumerate(I), C))
end

function BFD_1D(inst::tInstance)::Vector{tBin}
    return BF_1D(inst, order=(sortperm(inst.I, by=(x -> x.size))))
end

function BFD_1D(I::Vector{Int64}, C::Int64)::Vector{tBin}
    return BFD_1D(tInstance([tItem(id, size) for (id, size) in enumerate(I)], C))
end

function BFD_1D_into_GreadyRebuild(
        Routes  ::Vector{Route}             , # Routes of the instance
        Lmax    ::Int64                     , # Maximum capacity of an output
        O       ::Int64 = nothing           , # Number of output of each session
        R       ::Int64 = nothing           , # Number of route to sort
        tl      ::Int64 = 10                , # Time limit
        env     ::Gurobi.Env = Gurobi.Env() , # Gurobi environement (for display purpose)
    )::Tuple{Vector{Session}, Bool}


    (R == nothing) && (R = length(Routes))
    (O == nothing) && (O = length(Routes[1].assignment))

    sol_relax::Vector{tBin} = BFD_1D([sum(r.mail) for r in Routes], Lmax * O)

    sol::Vector{Session} = Vector{Session}(undef, length(sol_relax))

    valid::Bool = true

    # println(sol_relax)

    for (s_id::Int64, b::tBin) in enumerate(sol_relax)

        # println("b $s_id -> $b")
        ses = Session(Lmax, Route[Routes[i.id] for i in b.I])

        ses, valid = rebuildSession_knapSack_model_V3!(ses, tl, env)

        if valid
            sol[s_id] = ses
        else
            break
        end
    end

    return (sol, valid)
end

## ============================================================================================================== ##
##                  /###     ######              ######    #######              ######   #######                  ##
##                  # ##     ##    ##            ##    ##  ##    ##    ##      ##        ##    ##                 ##
##                    ##     ##     #    ####    #######   #######       ##    ##  ###   #######                  ##
##                    ##     ##    ##            ##    ##  ##          ##      ##    ##  ##  ##                   ##
##                  ######   ######              #######   ##                   ######   ##   ##                  ##
## ============================================================================================================== ##

function setup_1DBP(
        V   ::Vector{Int64}                                 , 
        C   ::Int64                                         ; 
        env ::Gurobi.Env                    = Gurobi.Env()  , 
        tl  ::Int64                         = 10            ,
        sol ::Union{Vector{tBin}, Nothing}  = nothing       ,
    )::Union{Vector{tBin}, Nothing}

    len_I::Int64 = length(V)
    UB  ::Int64 = (sol == nothing) ? (len_I) : (length(sol))

    I = collect(1:len_I)
    B = collect(1:UB)

    md = Model(() -> Gurobi.Optimizer(env))
    set_silent(md)
    set_optimizer_attribute(md, "OutputFlag", 0)
    set_optimizer_attribute(md, "TimeLimit", tl)

    # x[i, b] ∈ {0, 1}, equals 1 iff item i is in bin b, 0 otherwise.
    @variable(md, x[i in I, b in B], Bin)
    
    # y[b] ∈ {0, 1}, equals 1 iff bin b is open, 0 otherwise.
    @variable(md, y[b in B], Bin)

    # OBJ: minimise the number of opened bin
    @objective(md, Min, sum(y)) 

    # CST 1: Each item mus be assigned to a bin
    for i in I
        @constraint(md, sum(x[i, :]) == 1)
    end

    # CST 2: Capacity constraint
    for b in B
        @constraint(md, sum([x[i, b] * V[i] for i in I]) ≤ C * y[b])
    end

    # CST 3: open bin b before b+1 (to break symetrie)
    for b in B[1:end-1]    
        @constraint(md, y[b] >= y[b+1])
    end 

    if sol != nothing
        for b in B
            for i in sol[b].I
                set_start_value(x[i.id, b], 1)
            end
            set_start_value(y[b], 1)
        end
    end

    optimize!(md)

    if termination_status(md) == MOI.OPTIMAL
        x = round.(Int64, value.(md[:x]))
        y = round.(Int64, value.(md[:y]))

        obj = sum(y)
        sol = Vector{tBin}(undef, round(Int64, obj))

        println(x)
        println(y)
        println("obj -> $obj, $(sum(y)), $(y[1])")

        for b=1:obj
            sol[b] = tBin([tItem(i, V[i]) for i in I if (x[i, b] == 1)], sum([x[i, b] * V[i] for i in I]), C)
        end

        println(sol)

        return sol
    else
        return nothing
    end
end

function BP_1D_into_GreadyRebuild(
        Routes  ::Vector{Route}             , # Routes of the instance
        Lmax    ::Int64                     , # Maximum capacity of an output
        O       ::Int64 = nothing           , # Number of output of each session
        R       ::Int64 = nothing           , # Number of route to sort
        tl      ::Int64 = 100               , # Time limit
        env     ::Gurobi.Env = Gurobi.Env() , # Gurobi environement (for display purpose)
    )::Tuple{Vector{Session}, Bool}


    (R == nothing) && (R = length(Routes))
    (O == nothing) && (O = length(Routes[1].assignment))

    relaxed_instance    ::Vector{Int64} = [sum(r.mail) for r in Routes]
    sol_heuristic_relax ::Vector{tBin}  = BFD_1D(relaxed_instance, Lmax * O)

    sol_relax::Vector{tBin} = setup_1DBP(relaxed_instance, Lmax * O, sol=sol_heuristic_relax, tl=tl, env=env)

    # println(sol_heuristic_relax)
    # println(sol_relax)

    (sol_relax == nothing) && (sol_relax = sol_heuristic_relax)

    sol::Vector{Session} = Vector{Session}(undef, length(sol_relax))

    valid::Bool = true

    println(sol_relax)

    for (s_id::Int64, b::tBin) in enumerate(sol_relax)

        ses = Session(Lmax, Route[Routes[i.id] for i in b.I])

        ses, valid = rebuildSession_knapSack_model_V3!(ses, tl, env)
        # println("b $s_id -> $b , rebuild success ? -> $valid")

        if valid
            sol[s_id] = ses
        else
            break
        end
    end

    return (sol, valid)
end


# ================================================== #
#                       Display                      #
# ================================================== #

function Base.show(io::IO, sol::Solution)
    println(io, "Solution:\n    Permutation: $(sol.permutation)\n    \nSessions: (x$(length(sol.sessions)))\n")

    for e in sol.sessions
        println(e)
    end

    println("Average load STD: $(round(fitness(sol, SessionLoadSTD)/length(sol.sessions), digits=2))%")
end

# ================================================== #
#                   Fitness Function                 #
# ================================================== #

abstract type FitnessSolution end

struct WeightedHollowedPercentage   <: FitnessSolution end
struct SessionNonFilledOut          <: FitnessSolution end
struct SessionLoadSTD               <: FitnessSolution end
struct NbSession                    <: FitnessSolution end
struct ObjGA                        <: FitnessSolution end
struct ObjGA_2                      <: FitnessSolution end

# Fitness on Full Solution
@inline specialized(sol::Solution,::Type{WeightedHollowedPercentage}) = sum([1/(i * (100 - fitness(s, HollowedPercentage)))  for (i, s) in enumerate(sol.sessions)])
@inline specialized(sol::Solution,::Type{SessionNonFilledOut})        = sum([fitness(s, NonFilledOutputs) for s in sol.sessions])
@inline specialized(sol::Solution,::Type{SessionLoadSTD})             = sum([fitness(s, LoadSTD) for s in sol.sessions])
@inline specialized(sol::Solution,::Type{NbSession})                  = length(sol.sessions)
@inline specialized(sol::Solution,::Type{ObjGA})                      = fitness(sol, NbSession) + (fitness(sol, SessionLoadSTD) / (fitness(sol, NbSession) * sol.sessions[1].Lmax))
@inline specialized(sol::Solution,::Type{ObjGA_2})                    = fitness(sol, NbSession) + fitness(sol, WeightedHollowedPercentage)

# Multiple criteria fitness on Solution
@inline fitness(sol::Solution, TAG::Type{<:FitnessSolution}) = specialized(sol, TAG)
@inline fitness(sol::Solution, TAG::Tuple{DataType, DataType}) = (specialized(sol, TAG[1]), specialized(sol, TAG[2]))
@inline fitness(sol::Solution, TAG::Tuple{DataType, DataType, DataType}) = (specialized(sol, TAG[1]), specialized(sol, TAG[2]), specialized(sol, TAG[3]))
@inline fitness(sol::Solution, TAG::Tuple{Vararg{DataType, N}}) where N = ntuple(i -> specialized(sol, TAG[i]), N)::NTuple{N, Float64}

@inline fitness(sol::Solution, TAGs::Vector{Type{<:FitnessSolution}}) = [specialized(sol, TAG) for TAG in TAGs]

# ================================================== #
#                     Miscelaneous                   #
# ================================================== #

function validSolution(instance::Instance, sol::Solution, display::Bool = true)
    usedRoute_vect::Vector{Int64} = zeros(Int64, instance.nbRoute)
    valid::Bool = true

    for s::Session in sol.sessions
        isSessionValid(s) || (valid = false)
        for r in s.route
            usedRoute_vect[r.id] += 1
        end
    end

    for (k, v) in enumerate(usedRoute_vect)
        (v == 1) || (valid = false; (display && println(" -> Route $k: appears $v times.")))
    end

    return valid
end

function isSolutionValid(instance::Instance, sol::Vector{Session}, display::Bool = true)
    usedRoute_vect::Vector{Int64} = zeros(Int64, instance.nbRoute)
    valid::Bool = true

    for s::Session in sol
        isSessionValid(s) || (valid = false)
        for r in s.route
            isRouteValid(r) || (valid = false)
            (instance.route[r.id].mail == r.mail) || (valid = false)
            usedRoute_vect[r.id] += 1
        end
    end

    for (k, v) in enumerate(usedRoute_vect)
        (v == 1) || (valid = false; (display && println(" -> Route $k: appears $v times.")))
    end

    return valid
end

function isSolutionValid(instance::Instance, sol::Solution, display::Bool = true)
    return isSolutionValid(instance, sol.route, display)
end
