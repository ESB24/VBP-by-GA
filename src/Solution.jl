begin
    import Base:copy
    import Random
end

begin # Files 
    include("Instance.jl")
    include("DataStruct.jl")    
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
        sId         ::Int64     = 1
        added       ::Bool      = false
        r           ::Route     = instance.route[routeId]

        while !added
            if sId <= length(sol.sessions)
                # print("(sava_nt), ")
                # st          ::Union{Session, Nothing} = nothing
                # flag        ::Bool      = false
                flagId      ::Int64     = 1

                while !added && flagId <= length(TAG_AddRoute)
                    # st, added = addRoute(sol.sessions[sId], r, TAG_AddRoute[flagId], Δ=Δ, TAG_FitSes=TAG_FitSes, tl=tl, env=env)
                    _, added = addRoute(sol.sessions[sId], r, TAG_AddRoute[flagId], Δ, TAG_FitSes, tl, env)
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
            if length(s.route) < length(s.loads) && certificat_CapacityVolume(s, r)                             #  |   | if Volume and number of route certificates
                added, newAssignment::Vector{Int64} = addRoute_NFBA(s, r)                                       #  |   |   |  Try adding r using left alligned assignment
                if added                                                                                        #  |   |   |  if Left alligned assignment valid
                    s, added = addRoute_SAVA!(s, r)                                                             #  |   |   |   |  Try adding r using a smoother assignment
                    if !added                                                                                   #  |   |   |   |  if smoother assignment not valid
                        r = Route(r.id, newAssignment, r.mail)                                                  #  |   |   |   |   |  use Left alligned assignment
                        s.loads += newAssignment                                                                #  |   |   |   |   |
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
    
    if length(s.route) < length(s.loads) && certificat_CapacityVolume(s, r)                             #  if Volume and number of route certificates
        added, newAssignment::Vector{Int64} = addRoute_NFBA(s, r)                                       #   |  Try adding r using left alligned assignment
        if added                                                                                        #   |  if Left alligned assignment valid
            s, added = addRoute_SAVA!(s, r)                                                             #   |   |  Try adding r using a smoother assignment
            if !added                                                                                   #   |   |  if smoother assignment not valid
                r = Route(r.id, newAssignment, r.mail)                                                  #   |   |   |  use Left alligned assignment
                s.loads += newAssignment                                                                #   |   |   |
                push!(s.route, r)                                                                       #   |   |   |
                added = true                                                                            #   |   |   |
            end                                                                                         #   |   |  end
        else                                                                                            #   |  else
            s, added = addRoute_Rebuild_Knapsack_model_V3!(s, r, tl, env)                               #   |   |  rebuild session
        end                                                                                             #   |  endif
    end                                                                                                 #  endif

    return added
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
