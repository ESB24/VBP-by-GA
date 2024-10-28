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
        TAG_AddRound        ::Vector{Type{<:SimpleAddRound}}    = Type{<:SimpleAddRound}[SAVA, SPCF, OPTIMOVE]  , 
        Δ                   ::Int64                             = 2                                             ,
        TAG_FitSes          ::Type{<:FitnessSession}            = LoadSTD                                       ,
        env                 ::Gurobi.Env                        = Gurobi.Env()                                  ,
        tl                  ::Int64                             = 100                                           , 
    )
    return buildSolution_FFD(instance, randperm(instance.nbRound), Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRound=TAG_AddRound, env=env, tl=tl)
end

function buildSolution_FFD(
        instance            ::Instance                                                                          , 
        perm                ::Vector{Int64}                                                                     ;
        TAG_AddRound        ::Vector{Type{<:SimpleAddRound}}    = Type{<:SimpleAddRound}[SAVA, NFBA, OPTIMOVE_STAGE1]  ,  
        Δ                   ::Int64                             = 2                                             ,
        TAG_FitSes          ::Type{<:FitnessSession}            = LoadSTD                                       ,
        env                 ::Gurobi.Env                        = Gurobi.Env()                                  ,
        tl                  ::Int64                             = 100                                           , 
    )

    sol::Solution = Solution(perm, [Session(instance.C, instance.nbOut) for _=1:minSession(instance)])
    # upgrades::Vector{Int64} = zeros(Int64, length(TAG_AddRound))

    for roundId in perm
        sId         ::Int64     = 1
        added       ::Bool      = false
        r           ::Round     = instance.rounds[roundId]

        while !added
            if sId <= length(sol.sessions)
                # print("(sava_nt), ")
                # st          ::Union{Session, Nothing} = nothing
                # flag        ::Bool      = false
                flagId      ::Int64     = 1

                while !added && flagId <= length(TAG_AddRound)
                    # st, added = addRound(sol.sessions[sId], r, TAG_AddRound[flagId], Δ=Δ, TAG_FitSes=TAG_FitSes, tl=tl, env=env)
                    _, added = addRound(sol.sessions[sId], r, TAG_AddRound[flagId], Δ, TAG_FitSes, tl, env)
                    # flag && (upgrades[flagId] += 1)
                    flagId += 1
                end
                sId += 1
            else # create new Session
                push!(sol.sessions, Session(instance.C, [r], deepcopy(r.assignment)))

                added = true
            end
        end
    end

    return sol # , upgrades
end

function buildSolution_FF_final(instance::Instance, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    return buildSolution_FFD_final(instance, randperm(instance.nbRound), tl, env)
end

function buildSolution_FFD_final(instance::Instance, perm::Vector{Int64}, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    sol::Solution = Solution(perm, [Session(instance.C, instance.nbOut) for _=1:minSession(instance)])
    for roundId in perm
        r::Round = instance.rounds[roundId]
        sol = InsertRound_Procedure!(sol, r, tl, env)
    end
    return sol
end

function buildSolution_BF_final(instance::Instance, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    return buildSolution_BFD_final(instance, randperm(instance.nbRound), tl, env)
end

function buildSolution_BFD_final(instance::Instance, perm::Vector{Int64}, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    sol::Solution = Solution(perm, [Session(instance.C, instance.nbOut)])
    for roundId in perm
        r::Round = instance.rounds[roundId]
        # println("new order -> $(sortperm(sol.sessions, by=x -> (100 - fitness(x, HollowedPercentage))))")
        sort!(sol.sessions, by=x -> (fitness(x, HollowedPercentage)))
        sol = InsertRound_Procedure!(sol, r, tl, env)
    end
    return sol
end

function buildSolution_WF_final(instance::Instance, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    return buildSolution_WFD_final(instance, randperm(instance.nbRound), tl, env)
end

function buildSolution_WFD_final(instance::Instance, perm::Vector{Int64}, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    sol::Solution = Solution(perm, [Session(instance.C, instance.nbOut)])
    for roundId in perm
        r::Round = instance.rounds[roundId]
        # println("new order -> $(sortperm(sol.sessions, by=x -> (100 - fitness(x, HollowedPercentage))))")
        sort!(sol.sessions, by=x -> (-fitness(x, HollowedPercentage)))
        sol = InsertRound_Procedure!(sol, r, tl, env)
    end
    return sol
end

"""
function InsertRound_Procedure!(sol::Solution, r::Round, env::Gurobi.Env = Gurobi.Env(), tl::Int64 = 10)

Full procedure to insert a round r into a solution s, following a First-fit like behavior.
"""
@inline function InsertRound_Procedure!(sol::Solution, r::Round, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    sId         ::Int64     = 1
    added       ::Bool      = false

    while !added                                                                                                # While r hasn't been added to any sorting session
        if sId <= length(sol.sessions)                                                                          #  |  if any more session
            s::Session = sol.sessions[sId]                                                                      #  |   |
            if length(s.rounds) < length(s.loads) && certificat_CapacityVolume(s, r)                            #  |   | if Volume and number of round certificates
                added, newAssignment::Vector{Int64} = addRound_NFBA(s, r)                                       #  |   |   |  Try adding r using left alligned assignment
                if added                                                                                        #  |   |   |  if Left alligned assignment valid
                    s, added = addRound_SAVA!(s, r)                                                             #  |   |   |   |  Try adding r using a smoother assignment
                    if !added                                                                                   #  |   |   |   |  if smoother assignment not valid
                        r = Round(r.id, newAssignment, r.batches)                                               #  |   |   |   |   |  use Left alligned assignment
                        s.loads += newAssignment                                                                #  |   |   |   |   |
                        push!(s.rounds, r)                                                                      #  |   |   |   |   |
                        added = true                                                                            #  |   |   |   |   |
                    end                                                                                         #  |   |   |   |  end
                else                                                                                            #  |   |   |  else
                    s, added = addRound_Rebuild_Knapsack_model_V3!(s, r, tl, env)                               #  |   |   |   |  rebuild session
                end                                                                                             #  |   |   |  endif
            end                                                                                                 #  |   |  endif
            sId += 1                                                                                            #  |   |  Pass to next Session (relevent only r has not been added yet)
        else                                                                                                    #  |  else
            push!(sol.sessions, Session(sol.sessions[1].C, [r], deepcopy(r.assignment)))                        #  |   |  Create a new Session with r inside (will end the while loop)
            added = true                                                                                        #  |   |
        end                                                                                                     #  |  end
    end                                                                                                         # end
    return sol                                                                                                  # return solution
end

@inline function insert_route!(s::Session, r::Round)::Bool
    added::Bool = false
    
    if length(s.rounds) < length(s.loads) && certificat_CapacityVolume(s, r)                            #  if Volume and number of round certificates
        added, newAssignment::Vector{Int64} = addRound_NFBA(s, r)                                       #   |  Try adding r using left alligned assignment
        if added                                                                                        #   |  if Left alligned assignment valid
            s, added = addRound_SAVA!(s, r)                                                             #   |   |  Try adding r using a smoother assignment
            if !added                                                                                   #   |   |  if smoother assignment not valid
                r = Round(r.id, newAssignment, r.batches)                                               #   |   |   |  use Left alligned assignment
                s.loads += newAssignment                                                                #   |   |   |
                push!(s.rounds, r)                                                                      #   |   |   |
                added = true                                                                            #   |   |   |
            end                                                                                         #   |   |  end
        else                                                                                            #   |  else
            s, added = addRound_Rebuild_Knapsack_model_V3!(s, r, tl, env)                               #   |   |  rebuild session
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
@inline specialized(sol::Solution,::Type{ObjGA})                      = fitness(sol, NbSession) + (fitness(sol, SessionLoadSTD) / (fitness(sol, NbSession) * sol.sessions[1].C))
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
    usedRound_vect::Vector{Int64} = zeros(Int64, instance.nbRound)
    valid::Bool = true

    for s::Session in sol.sessions
        validSession(s) || (valid = false)
        for r in s.rounds
            usedRound_vect[r.id] += 1
        end
    end

    for (k, v) in enumerate(usedRound_vect)
        (v == 1) || (valid = false; (display && println(" -> Round $k: appears $v times.")))
    end

    return valid
end
