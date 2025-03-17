begin
    # data structures
    include("DataStruct.jl")
    include("Route.jl")
    include("Session.jl")
    include("Solution.jl")
    include("Instance.jl")
end

begin
    # using PyPlot
    using Random
end

# ================================================== #
#                   Son Acceptance                   #
# ================================================== #

abstract type AcceptSon end

struct LandingPercent   <: AcceptSon end
struct BestOneParent    <: AcceptSon end
struct Always           <: AcceptSon end

"""
```Julia
    function acceptSon(p1::Solution, p2::Solution, offspring::Solution, TAG_FitSol::Type{<:FitnessSolution},::Type{<:AcceptSon})::Bool
```
AcceptSon function aim to determine if a offspring is good enough to be keept in the current generation regarding it's parents objective value.

`p1` and `p2` are the parents of the `offspring`.
`TAG_FitSol` is the tag linked to the objective function you want tu use to compare p1, p2 and offspring.
The last attribute is a tag to choose which criteria you want to use to determine if an offspring should be keept:
 - `LandingPercent`: the offspring will be keept if he is as good as one of its parent with a `landing`% tolerance. (`landing` is an optional attribute set to 0.2% as default)
 - `BestOneParent`:  the offspring will be keept if he is as good as one of its parent.
 - `Always`: the offspring will always be keept.
""" 
@inline function acceptSon(p1::Solution, p2::Solution, offspring::Solution, TAG_FitSol::Type{<:FitnessSolution},::Type{LandingPercent}; landing::Float64 = 0.2, args...)::Bool

    fit_offspring   ::Float64   = fitness(offspring, TAG_FitSol) # compute all objective value
    fit_p1          ::Float64   = fitness(p1, TAG_FitSol)
    fit_p2          ::Float64   = fitness(p2, TAG_FitSol)

    return ((100. * fit_offspring) / fit_p1 <= (100 + landing)) && ((100. * fit_offspring) / fit_p2 <= (100 + landing))
end

@inline function acceptSon(p1::Solution, p2::Solution, offspring::Solution, TAG_FitSol::Type{<:FitnessSolution},::Type{BestOneParent})::Bool

    fit_offspring ::Float64   = fitness(offspring, TAG_FitSol)
    fit_p1  ::Float64   = fitness(p1, TAG_FitSol)
    fit_p2  ::Float64   = fitness(p2, TAG_FitSol)

    return fit_offspring <= fit_p1 || fit_offspring <= fit_p2
end

@inline acceptSon(p1::Solution, p2::Solution, offspring::Solution, TAG_FitSol::Type{<:FitnessSolution},::Type{Always})::Bool = true

# ================================================== #
#                   Parent Selection                 #
# ================================================== #

"""
```Julia
    function selectParent(previousGen::Vector{Solution}, ::Type{<:SelectParent})::Tuple{Solution, Solution}
```
Select 2 parents in order to create an offspring.
`previousGen` is the set of individuals among which the parents could be selected.
The last tag is again a way to determine a way to select the parents:
 - `SelectRandom`: will choose 2 different parent completely randomly.
 - `SelectElite`: will choose 2 different parent randomly but will value parents with a lower amount of sorting session.
"""

abstract type SelectParent end

struct SelectRandom           <: SelectParent end
struct SelectElite            <: SelectParent end

@inline function selectParent(previousGen::Vector{Solution}, ::Type{SelectRandom})::Tuple{Solution, Solution}
    id1, id2 = randperm(length(previousGen))[1:2]

    return previousGen[id1], previousGen[id2]
end

@inline function selectParent(previousGen::Vector{Solution}, ::Type{SelectElite})::Tuple{Solution, Solution}
    function biased_random(n::Int64)
        return ceil(Int64, n * (1 - sqrt(rand())))
    end

    id1 = biased_random(length(previousGen))
    id2 = biased_random(length(previousGen))

    # Assurer que les deux nombres sont différents
    while id1 == id2
        id2 = biased_random(length(previousGen))
    end

    return previousGen[id1], previousGen[id2]
end

# ================================================== #
#                      CrossOver                     #
# ================================================== #

"""
```Julia
    crossOver(instance::Instance, p1::Solution, p2::Solution, ::Type{Point1}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, env::Gurobi.Env=Gurobi.Env())
```
Apply the selected Cross-Over on two given parents.
`instance` is the treated instance.
`p1` and `p2` are the parents of the future offspring.
A tag to choose among tree types of Cross-Over:
 - `Point1`: Cross-Over one point. (on permutation)
 - `PointN`: Cross-Over n point. (on permutation)
 - `Edit`: will modify the p1 solution with some p2 features. (modify the solution directly not the permutation which is faster)

`TAG_FitSes` and `env` are two optional attribute to choose the used objective function and provide a Gurobi environement (display purpose) respectively.
"""

abstract type CrossOver end

struct Point1           <: CrossOver end
struct PointN           <: CrossOver end
struct Edit             <: CrossOver end

function crossOver(instance::Instance, p1::Solution, p2::Solution, ::Type{Point1}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, env::Gurobi.Env=Gurobi.Env())# ::Solution
    O::Int64 = length(p1.permutation)
    c::Int64 = rand(1:O-1)

    newPerm::Vector{Int64} = zeros(Int64, O)

    # copy first segment from p1
    newPerm[1:c] = p1.permutation[1:c]

    # deduce second segment from p2
    i::Int64 = c+1
    for v in p2.permutation
        if !(v in newPerm)
            newPerm[i] = v
            i += 1
        end
    end

    count(x -> x == 0, newPerm) > 0 && println("\n\np1 -> $(p1.permutation)\np2 -> $(p2.permutation)\nson -> $(newPerm)")

    return buildSolution_FFD_final(instance, newPerm, TAG_FitSes, env)
end

function crossOver(instance::Instance, p1::Solution, p2::Solution, ::Type{PointN}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, env::Gurobi.Env=Gurobi.Env())::Solution
    m::Int64 = length(p1.permutation)
    n::Int64 = rand(1:ceil(Int64, m/4))

    p::Vector{Int64} = sort!(push!(randperm(m)[1:n], 0, m)) # n points

    interval::Vector{Tuple{Int64, UnitRange{Int64}}} = [(if (k%2 == 0) 1 else 2 end, p[k]+1:p[k+1]) for k=1:n+2 if k < n+2] # interval between the n points

    newPerm::Vector{Int64} = []
    for (s, i) in interval
        if s == 1 
            for k in i
                push!(newPerm, p1.permutation[k]) 
            end
        else 
            for k in i
                push!(newPerm, 0) 
            end
        end 
    end

    for (s, i) in interval # TODO reverse the adding in the new Permutation (new perm)
        if s == 2
            for pos=i
                for (k, v) in enumerate(p2.permutation)
                    if !(v in newPerm)
                        newPerm[pos] = v
                    end
                end
            end
        end
    end

    return buildSolution_FFD_final(instance, newPerm, TAG_FitSes, env)
end

function crossOver(instance::Instance, p1::Solution, p2::Solution, ::Type{Edit}, TAG_FitSes::Type{<:FitnessSession} = LoadSTD, env::Gurobi.Env=Gurobi.Env())
    
    son::Solution = Solution(deepcopy(p1.permutation), [Session(s.Lmax, [Route(r.id, deepcopy(r.assignment), r.mail) for r in s.route], deepcopy(s.loads)) for s in p1.sessions])

    removedRoute::Vector{Int64} = [e.id for e in p1.sessions[end].route]

    nonFullSession::Vector{Int64} = [i for (i, s) in enumerate(son.sessions[1:end-1]) if (fitness(s, HollowedPercentage) >= 1)]

    emptySession::Vector{Int64} = [length(son.sessions)]

    nbRemovedRoute::Int64 = 0

    while nbRemovedRoute <= round(Int64, instance.nbRoute/4) && length(emptySession) < length(nonFullSession)
        sId::Int64 = rand(nonFullSession)
        s::Session = son.sessions[sId]

        if !isempty(s.route)
            r = rand(s.route)

            filter!(e -> e.id != r.id, son.sessions[sId].route)
            filter!(e -> e != r.id, son.permutation)

            push!(removedRoute, r.id)

            # println("r = $(r.id), ")

            if isempty(s.route)
                push!(emptySession, sId) 
            end

            nbRemovedRoute += 1
        end
    end

    for s in son.sessions
        isempty(s.route) ? s.loads == zeros(Int64, length(s.loads)) : compute_output!(s)
    end
    son.sessions = son.sessions[1:end-1]
    filter!(e -> sum(e.loads) != 0, son.sessions)
    
    # round id to add 
    addOrder::Vector{Int64} = filter(e -> e in removedRoute, p2.permutation)
    filter!(e-> !(e in removedRoute), son.permutation)
    son.permutation = [son.permutation; addOrder]

    # for r in addOrder
    #     son = InsertRoute_Procedure!(instance, son, r, TAG_FitSes, env)
    # end

    for roundId in addOrder
        r::Route = instance.route[roundId]
        # println("new order -> $(sortperm(sol.sessions, by=x -> (100 - fitness(x, HollowedPercentage))))")
        sort!(son.sessions, by=x -> (fitness(x, HollowedPercentage)))
        son = InsertRoute_Procedure!(son, r, 10, env)
    end

    # println("removedRoute -> $(removedRoute)")
    # println("addOrder -> $(addOrder)")

    return son
end

# ================================================== #
#                       Mutation                     #
# ================================================== #

"""
```Julia
    mutation_swap!(perm::Vector{Int64})
```
A swap on a given mutation

"""

function mutation_swap!(perm::Vector{Int64})
    p1::Int64, p2::Int64 = randperm(length(perm))[1:2]
    perm[p1], perm[p2] = perm[p2], perm[p1]
    return perm
end

# write a raw version of the OnePoint Cross-Over to be faster. (ex: reduce type inference due to optional variable)
function crossOver_onePoint(p1_perm::Vector{Int64}, p2_perm::Vector{Int64})# ::Solution
    O::Int64 = length(p1_perm)
    c::Int64 = rand(1:O-1)

    son_perm::Vector{Int64} = zeros(Int64, O)

    # copy first segment from p1
    son_perm[1:c] = p1_perm[1:c]

    # deduce second segment from p2
    i::Int64 = c+1
    for v in p2_perm
        !(v in son_perm) && (son_perm[i] = v; i += 1)
    end

    return son_perm
end

# same as OnePoint Cross-Over but with Edit Cross-Over
function crossOver_edit(instance::Instance, p1::Solution, p2::Solution, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env())
    son::Solution = Solution(deepcopy(p1.permutation), [Session(s.Lmax, [Route(r.id, deepcopy(r.assignment), r.mail) for r in s.route], deepcopy(s.loads)) for s in p1.sessions])

    removedRoute::Vector{Int64} = [e.id for e in p1.sessions[end].route]

    nonFullSession::Vector{Int64} = [i for (i, s) in enumerate(son.sessions[1:end-1]) if (fitness(s, HollowedPercentage) >= 1)]

    emptySession::Vector{Int64} = [length(son.sessions)]

    nbRemovedRoute::Int64 = 0

    while nbRemovedRoute <= round(Int64, instance.nbRoute/4) && length(emptySession) < length(nonFullSession)
        sId::Int64 = rand(nonFullSession)
        s::Session = son.sessions[sId]

        if !isempty(s.route)
            r = rand(s.route)

            filter!(e -> e.id != r.id, son.sessions[sId].route)
            filter!(e -> e != r.id, son.permutation)

            push!(removedRoute, r.id)

            # println("r = $(r.id), ")

            if isempty(s.route)
                push!(emptySession, sId) 
            end

            nbRemovedRoute += 1
        end
    end

    for s in son.sessions
        isempty(s.route) ? s.loads == zeros(Int64, length(s.loads)) : compute_output!(s)
    end
    son.sessions = son.sessions[1:end-1]
    filter!(e -> sum(e.loads) != 0, son.sessions)
    
    # round id to add 
    addOrder::Vector{Int64} = filter(e -> e in removedRoute, p2.permutation)
    filter!(e-> !(e in removedRoute), son.permutation)
    son.permutation = [son.permutation; addOrder]

    # for r in addOrder
    #     son = InsertRoute_Procedure!(instance, son, r, TAG_FitSes, env)
    # end

    for roundId in addOrder
        r::Route = instance.route[roundId]
        # println("new order -> $(sortperm(sol.sessions, by=x -> (100 - fitness(x, HollowedPercentage))))")
        sort!(son.sessions, by=x -> (fitness(x, HollowedPercentage)))
        son = InsertRoute_Procedure!(son, r, 10, env)
    end

    # println("removedRoute -> $(removedRoute)")
    # println("addOrder -> $(addOrder)")

    return son
end

"""
```Julia
    GA(instance::Instance, max_pop::Int64, max_gen::Int64, max_time::Float64, fit::Function, crossover::Function, objective::Function)
```
> GA: Genetic Algoritm
"""
function GA(
            instance            ::Instance                                                  ;
            max_time            ::Float64                       = 120.                      ,   # time limit
            α_elite             ::Float64                       = 0.4                       ,   # pourcentage of elite in each core
            TAG_FitGA           ::Type{<:FitnessSolution}       = ObjGA                     ,
        )
# ====================< INIT >====================

    println("================================================================================")
    println("                            > Genetic Algorithm ")
    println("================================================================================")

    R::Int64 = instance.nbRoute
    O::Int64 = instance.nbOut

    max_pop::Int64 = max(R, 50)
    max_gen::Int64 = R * 10

    elite_size::Int64 = max(round(Int64, max_pop * 2α_elite), 50)       # space occupied by the elite of the previous generation in the current core

    # Exit Parameters
    opti_sol    ::Int64 = minSession(instance)  # Best known lower bound
    flag_opti   ::Bool  = false                 # Best bound had been reached

    start_time ::Float64 = time()               # Store start time (to allow a time limited forced exit)
    flag_tl    ::Bool    = false                # Time limit had been reached

    lastElite     ::Union{Float64, Nothing} = nothing   # Store the score of the last best generation. Allow an early exit if no improvement were made since last generation
    no_improv     ::Int64                   = 0         # number of generation without any improvement
    max_no_improv ::Int64                   = 4         # maximum number of generation without any improvement (for early exit)
    flag_improv   ::Bool                    = false     # no improvement for more than `max_no_imporv` generation

    mod            ::Int64 = 2       # Current building method 1 -> FF, 2 -> BF
    mod_switch     ::Int64 = 0       # Counter of mod switching
    max_mod_switch ::Int64 = 20      # Maximum number of mod switching

    # results & geneneration management
    gen_val ::Vector{GenVal} = [GenVal([], [], [])]  # Result of the GA (objective value of each individuals. Intended for plotting)
    gen     ::Generation     = Generation(0, [])     # Current generation

    # Building procedure Parameters
    env ::Gurobi.Env = Gurobi.Env()     # in order to remove gurobi licence display when a model is called
    tl  ::Int64      = 10               # models time limit
    tl_perm ::Int64  = 600

    # Sorting criteria for FFD (used to create the first generation)
    all_TAG = (MaxMin, NbBatches, Max, StdBatches, MeanBatches, SumBatches)
    all_sort_crit::Vector{Tuple} = []
    for TAG1 in all_TAG
        for TAG2 in all_TAG
            (TAG1 != TAG2) && push!(all_sort_crit, (TAG1, TAG2))
        end
    end

# ====================< GEN LOOP >====================

    while (gen.id <= max_gen) && !flag_opti && !flag_tl && !flag_improv

# ====================< Initial Population/ Elite selection >====================

        if gen.id == 0 
            print("\n\n <> Generation $(gen.id):\n     - FFD: (x$(length(all_sort_crit)+1))\n        ")

            Tags         = all_sort_crit[length(gen.pop)+1]                                                     # get sorting criteria tags (for FFD)
            perm         = sortperm([(-fitness(r, Tags[1]), -fitness(r, Tags[2])) for r in instance.route])    # compute used permutation (according to sorting criteria)
            
            sol  = buildSolution_BFD_final(instance, perm, tl, env)   # compute solution usinf FFD heuristic
            val  = fitness(sol, TAG_FitGA)                            # compute objective value

            push!(gen.pop, sol)                                                 # Store Solution
            push!(gen_val[gen.id+1].elite, val)                                 # Store objective value
            
            flag_opti = length(sol.sessions) <= opti_sol                        # test optimality
            flag_tl   = time() - start_time > max_time                          # test time limit

            print("($(all_sort_crit[length(gen.pop)]): $(round(val, digits=3)))$(if length(gen_val[gen.id+1].elite)%5 == 0 && length(gen_val[gen.id+1].elite) != 0 "\n        " else ", " end)")

            if !flag_opti && !flag_tl
                Tags         = all_sort_crit[length(gen.pop)+1]                                                     # get sorting criteria tags (for FFD)
                _, _, tmp_bound, perm = full_partitioning_01_KP(instance.route, instance.Lmax, length(sol.sessions), tl_perm, env)   # compute used permutation (according to sorting criteria)

                (tmp_bound != nothing) && (opti_sol = max(opti_sol, tmp_bound))
                (length(sol.sessions) <= opti_sol) && (return sol, true, gen_val, time() - start_time)
                
                sol = buildSolution_BFD_final(instance, perm, tl, env)              # compute solution usinf FFD heuristic
                val = fitness(sol, TAG_FitGA)                                       # compute objective value

                push!(gen.pop, sol)                                                 # Store Solution
                push!(gen_val[gen.id+1].elite, val)                                 # Store objective value
                
                flag_opti = length(sol.sessions) <= opti_sol                        # test optimality
                flag_tl   = time() - start_time > max_time                          # test time limit

                print("(BP perm: $(round(val, digits=3)))$(if length(gen_val[gen.id+1].elite)%5 == 0 && length(gen_val[gen.id+1].elite) != 0 "\n        " else ", " end)")
            end

            while length(gen.pop)+1 < length(all_sort_crit) && !flag_opti && !flag_tl
                Tags    = all_sort_crit[length(gen.pop)]                                                     # get sorting criteria tags (for FFD)
                perm    = sortperm([(-fitness(r, Tags[1]), -fitness(r, Tags[2])) for r in instance.route])    # compute used permutation (according to sorting criteria)
                
                sol     =  buildSolution_BF_final(instance, tl, env)                # buildSolution_BFD_final(instance, perm, tl, env)   # compute solution usinf FFD heuristic
                val     = fitness(sol, TAG_FitGA)                                   # compute objective value

                push!(gen.pop, sol)                                                 # Store Solution
                push!(gen_val[gen.id+1].elite, val)                                 # Store objective value
                
                flag_opti = length(sol.sessions) <= opti_sol                        # test optimality
                flag_tl   = time() - start_time > max_time                          # test time limit

                print("($(all_sort_crit[length(gen.pop)-1]): $(round(val, digits=3)))$(if length(gen_val[gen.id+1].elite)%5 == 0 && length(gen_val[gen.id+1].elite) != 0 "\n        " else ", " end)")
            end
        else
            print("\n\n <> Generation $(gen.id):")
            gen.pop = sort!(gen.pop, by = x -> fitness(x, TAG_FitGA))[1:elite_size]     # next gen initialisation
            push!(gen_val, GenVal([fitness(e, TAG_FitGA) for e in gen.pop], [], []))    # Store objective value

            if lastElite === nothing                                        # If gen n° == 1
                lastElite = sum(gen_val[end].elite)                         #  |  store generation score
            else                                                            # else
                currentElite = sum(gen_val[end].elite)                      #  |  compute current generation score
                if lastElite > currentElite                                 #  |  If current generation better than last generation
                    lastElite = currentElite                                #  |   |  update current generation score
                else                                                        #  |  else
                    no_improv += 1                                          #  |   |  
                    flag_improv = no_improv > max_no_improv

                    # (mod == 1) ? (mod = 2) : (mod = 1)
                    mod_switch += 1                 

                    print(" (no improvement, mod switched to $(mod==1 ? "FF" : "BF"))")
                end                                                         #  |  endif
            end                                                             # endif
        
            print("\n     - Elite: (x$(elite_size))\n        ")
            for (k, e) in enumerate(gen_val[gen.id+1].elite)
                print("($(round(e, digits=3)))$(if k%10 == 0 && k != 0 "\n        " else ", " end)")
            end
        end

# ====================< Random individuals >====================

        print("\n     - $(mod==1 ? "FF" : "BF"): (x$(if gen.id == 0 max_pop - length(all_sort_crit) else max_pop - elite_size end))\n        ")
        while length(gen.pop) < max_pop && !flag_opti && !flag_tl && !flag_improv
            sol::Solution = (mod == 1 ? buildSolution_FF_final(instance, tl, env) : buildSolution_BF_final(instance, tl, env))            # solution computed usinf FF heuristic
            val::Float64  = fitness(sol, TAG_FitGA)                           # compute objective value

            push!(gen.pop, sol)                          # Store Solution
            push!(gen_val[gen.id+1].ffr, val)            # Store objective value

            flag_opti = length(sol.sessions) <= opti_sol # test optimality
            flag_tl   = time() - start_time > max_time      # test time limit

            print("($(round(val, digits=3)))$(if length(gen_val[gen.id+1].ffr)%10 == 0 && length(gen_val[gen.id+1].ffr) != 0 "\n        " else ", " end)")
        end

# ====================< Extended Population -> Cross-Over, Mutation >====================

        print("\n     - Cross-Over: (x$(max_pop))\n        ")
        parent_order::Vector{Int64} = randperm(length(gen.pop))
        while length(gen.pop) < 2max_pop && !flag_opti && !flag_tl && !flag_improv
            p1::Solution = gen.pop[popfirst!(parent_order)] # Get the next two parents
            p2::Solution = gen.pop[popfirst!(parent_order)]

            perm_son1::Vector{Int64} = mutation_swap!(crossOver_onePoint(p1.permutation, p1.permutation))   # compute the first son permutation CrossOver 1 point + swap
            perm_son2::Vector{Int64} = mutation_swap!(crossOver_onePoint(p2.permutation, p1.permutation))   # compute the second son permutation CrossOver 1 point + swap

            sol_son1::Solution = (mod == 1 ? buildSolution_FFD_final(instance, perm_son1, tl, env) : buildSolution_BFD_final(instance, perm_son1, tl, env))
            val_son1::Float64 = fitness(sol_son1, TAG_FitGA)                                                # compute first son objective value

            push!(gen.pop, sol_son1)                          # add first son to the population
            push!(gen_val[gen.id+1].cross_over, val_son1)     # add the first son objective value to the generation value

            flag_opti = length(sol_son1.sessions) <= opti_sol
            flag_tl   = time() - start_time > max_time
            if flag_opti || flag_tl
                print("($(round(fitness(p1, TAG_FitGA), digits=3)) & $(round(fitness(p2, TAG_FitGA), digits=3)) | CrossOver 1Point + Swap ->  $(round(val_son1, digits=3)))$(if length(gen_val[gen.id+1].cross_over)%4 == 0 && length(gen_val[gen.id+1].cross_over) != 0 "\n        " else ", " end)")
            else
                sol_son2::Solution = (mod == 1 ? buildSolution_FFD_final(instance, perm_son2, tl, env) : buildSolution_BFD_final(instance, perm_son2, tl, env))                      # compute second son solution
                val_son2::Float64 = fitness(sol_son2, TAG_FitGA)                                                # compute second son objective value

                push!(gen.pop, sol_son2)                          # add second son to the population
                push!(gen_val[gen.id+1].cross_over, val_son2)     # add the second son objective value to the generation value

                # update other component
                flag_opti = length(sol_son2.sessions) <= opti_sol
                flag_tl   = time() - start_time > max_time

                print("($(round(fitness(p1, TAG_FitGA), digits=3)) & $(round(fitness(p2, TAG_FitGA), digits=3)) | CrossOver 1Point + Swap + $(mod==1 ? "FF" : "BF") ->  $(round(val_son1, digits=3)) and $(round(val_son2, digits=3)))$(if length(gen_val[gen.id+1].cross_over)%4 == 0 && length(gen_val[gen.id+1].cross_over) != 0 "\n        " else ", " end)")
            end
        end

        println("\n <> Elapsed time: $(round(time() - start_time, digits=3))s")
        gen.id += 1
    end

# ====================< End results >====================
    # sort the last generation computed 
    sort!(gen.pop, by = x -> fitness(x, TAG_FitGA))
    # handle flags (special stop condition)
    if flag_opti
        println("\n <:D> Optimum reached: $(round(fitness(gen.pop[1], TAG_FitGA), digits=3)) <:D>\n")
    elseif flag_tl
        println("\n <!> Time limit reached: $(time() - start_time) <!>\n")
    elseif flag_improv
        println("\n <!> No improvements in $(max_no_improv) generations <!>\n")
    end

    return gen.pop[1], flag_opti, gen_val, time() - start_time
end

# copy of GA() to test random manipulation without breaking everything.
function GA_V2(
        instance            ::Instance                                                  ;
        max_time            ::Float64                       = 120.                      ,   # time limit
        α_elite             ::Float64                       = 0.4                       ,   # pourcentage of elite in each core
        TAG_FitGA           ::Type{<:FitnessSolution}       = ObjGA                     ,
    )
# ====================< INIT >====================

    println("================================================================================")
    println("                            > Genetic Algirithm ")
    println("================================================================================")

    R::Int64 = instance.nbRoute
    O::Int64 = instance.nbOut

    max_pop::Int64 = R
    max_gen::Int64 = R * 10

    elite_size::Int64 = round(Int64, max_pop * 2α_elite)       # space occupied by the elite of the previous generation in the current core

    # Exit Parameters
    opti_sol    ::Int64 = minSession(instance)  # Best known lower bound
    flag_opti   ::Bool  = false                 # Best bound had been reached

    start_time ::Float64 = time()               # Store start time (to allow a time limited forced exit)
    flag_tl    ::Bool    = false                # Time limit had been reached

    lastElite     ::Union{Float64, Nothing} = nothing   # Store the score of the last best generation. Allow an early exit if no improvement were made since last generation
    no_improv     ::Int64                   = 0         # number of generation without any improvement
    max_no_improv ::Int64                   = 4         # maximum number of generation without any improvement (for early exit)
    flag_improv   ::Bool                    = false     # no improvement for more than `max_no_imporv` generation

    mod            ::Int64 = 1       # Current building method 1 -> FF, 2 -> BF

    # results & geneneration management
    gen_val ::Vector{GenVal} = [GenVal([], [], [])]  # Result of the GA (objective value of each individuals. Intended for plotting)
    gen     ::Generation     = Generation(0, [])     # Current generation

    # Building procedure Parameters
    env ::Gurobi.Env = Gurobi.Env()     # in order to remove gurobi licence display when a model is called
    tl  ::Int64      = 10               # models time limit

    # Sorting criteria for FFD (used to create the first generation)
    all_TAG = (MaxMin, Max, StdBatches, StdAssignment, NbBatches, MeanBatches, MeanAssignment)
    all_sort_crit::Vector{Tuple} = []
    for TAG1 in all_TAG
        for TAG2 in all_TAG
            (TAG1 != TAG2) && push!(all_sort_crit, (TAG1, TAG2))
        end
    end

# ====================< GEN LOOP >====================

    while (gen.id <= max_gen) && !flag_opti && !flag_tl && !flag_improv

# ====================< Initial Population/ Elite selection >====================

        if gen.id == 0 
            print("\n\n <> Generation $(gen.id):\n     - FFD: (x$(length(all_sort_crit)))\n        ")
            while length(gen.pop) < length(all_sort_crit) && !flag_opti && !flag_tl
                Tags ::Tuple         = all_sort_crit[length(gen.pop)+1]                                                     # get sorting criteria tags (for FFD)
                perm ::Vector{Int64} = sortperm([(-fitness(r, Tags[1]), -fitness(r, Tags[2])) for r in instance.route])    # compute used permutation (according to sorting criteria)
                
                sol ::Solution = buildSolution_FFD_final(instance, perm, tl, env)   # compute solution usinf FFD heuristic
                val ::Float64  = fitness(sol, TAG_FitGA)                            # compute objective value

                push!(gen.pop, sol)                                                 # Store Solution
                push!(gen_val[gen.id+1].elite, val)                                 # Store objective value
                
                flag_opti = length(sol.sessions) <= opti_sol                        # test optimality
                flag_tl   = time() - start_time > max_time                          # test time limit

                print("($(all_sort_crit[length(gen.pop)]): $(round(val, digits=3)))$(if length(gen_val[gen.id+1].elite)%5 == 0 && length(gen_val[gen.id+1].elite) != 0 "\n        " else ", " end)")
            end
        else
            print("\n\n <> Generation $(gen.id):")
            gen.pop = sort!(gen.pop, by = x -> fitness(x, TAG_FitGA))[1:elite_size]     # next gen initialisation
            push!(gen_val, GenVal([fitness(e, TAG_FitGA) for e in gen.pop], [], []))    # Store objective value

            if lastElite === nothing                                        # If gen n° == 1
                lastElite = sum(gen_val[end].elite)                         #  |  store generation score
            else                                                            # else
                currentElite = sum(gen_val[end].elite)                      #  |  compute current generation score
                if lastElite > currentElite                                 #  |  If current generation better than last generation
                    lastElite = currentElite                                #  |   |  update current generation score
                else                                                        #  |  else
                    no_improv += 1                                          #  |   |  
                    flag_improv = no_improv > max_no_improv

                    mod = (mod + 1)%3 + 1             

                    print(" (no improvement, mod switched to $(mod==1 ? "FF" : mod==3 ? "BF" : "WF"))")
                end                                                         #  |  endif
            end                                                             # endif

            print("\n     - Elite: (x$(elite_size))\n        ")
            for (k, e) in enumerate(gen_val[gen.id+1].elite)
                print("($(round(e, digits=3)))$(if k%10 == 0 && k != 0 "\n        " else ", " end)")
            end
        end

# ====================< Random individuals >====================

    print("\n     - $(mod==1 ? "FF" : mod==3 ? "BF" : "WF"): (x$(if gen.id == 0 max_pop - length(all_sort_crit) else max_pop - elite_size end))\n        ")
    while length(gen.pop) < max_pop && !flag_opti && !flag_tl && !flag_improv
        sol::Solution = (mod == 1 ? buildSolution_FF_final(instance, tl, env) : mod == 3 ? buildSolution_BF_final(instance, tl, env) : buildSolution_WF_final(instance, tl, env))            # solution computed usinf FF heuristic
        val::Float64  = fitness(sol, TAG_FitGA)                           # compute objective value

        push!(gen.pop, sol)                          # Store Solution
        push!(gen_val[gen.id+1].ffr, val)            # Store objective value

        flag_opti = length(sol.sessions) <= opti_sol # test optimality
        flag_tl   = time() - start_time > max_time      # test time limit

        print("($(round(val, digits=3)))$(if length(gen_val[gen.id+1].ffr)%10 == 0 && length(gen_val[gen.id+1].ffr) != 0 "\n        " else ", " end)")
    end

# ====================< Extended Population -> Cross-Over, Mutation >====================

    print("\n     - Cross-Over: (x$(max_pop))\n        ")
    parent_order::Vector{Int64} = randperm(length(gen.pop))
    while length(gen.pop) < 2max_pop && !flag_opti && !flag_tl && !flag_improv
        p1::Solution = gen.pop[popfirst!(parent_order)] # Get the next two parents
        p2::Solution = gen.pop[popfirst!(parent_order)]

        perm_son1::Vector{Int64} = mutation_swap!(crossOver_onePoint(p1.permutation, p1.permutation))   # compute the first son permutation CrossOver 1 point + swap
        perm_son2::Vector{Int64} = mutation_swap!(crossOver_onePoint(p2.permutation, p1.permutation))   # compute the second son permutation CrossOver 1 point + swap

        sol_son1::Solution = (mod == 1 ? buildSolution_FF_final(instance, tl, env) : mod == 3 ? buildSolution_BF_final(instance, tl, env) : buildSolution_WF_final(instance, tl, env))
        val_son1::Float64 = fitness(sol_son1, TAG_FitGA)                                                # compute first son objective value

        push!(gen.pop, sol_son1)                          # add first son to the population
        push!(gen_val[gen.id+1].cross_over, val_son1)     # add the first son objective value to the generation value

        flag_opti = length(sol_son1.sessions) <= opti_sol
        flag_tl   = time() - start_time > max_time
        if flag_opti || flag_tl
            print("($(round(fitness(p1, TAG_FitGA), digits=3)) & $(round(fitness(p2, TAG_FitGA), digits=3)) | CrossOver 1Point + Swap ->  $(round(val_son1, digits=3)))$(if length(gen_val[gen.id+1].cross_over)%4 == 0 && length(gen_val[gen.id+1].cross_over) != 0 "\n        " else ", " end)")
        else
            sol_son2::Solution = (mod == 1 ? buildSolution_FF_final(instance, tl, env) : mod == 3 ? buildSolution_BF_final(instance, tl, env) : buildSolution_WF_final(instance, tl, env))                      # compute second son solution
            val_son2::Float64 = fitness(sol_son2, TAG_FitGA)                                                # compute second son objective value

            push!(gen.pop, sol_son2)                          # add second son to the population
            push!(gen_val[gen.id+1].cross_over, val_son2)     # add the second son objective value to the generation value

            # update other component
            flag_opti = length(sol_son2.sessions) <= opti_sol
            flag_tl   = time() - start_time > max_time

            print("($(round(fitness(p1, TAG_FitGA), digits=3)) & $(round(fitness(p2, TAG_FitGA), digits=3)) | CrossOver 1Point + Swap + $(mod==1 ? "FF" : mod==3 ? "BF" : "WF") ->  $(round(val_son1, digits=3)) and $(round(val_son2, digits=3)))$(if length(gen_val[gen.id+1].cross_over)%4 == 0 && length(gen_val[gen.id+1].cross_over) != 0 "\n        " else ", " end)")
        end
    end

    println("\n <> Elapsed time: $(round(time() - start_time, digits=3))s")
    gen.id += 1
    end

# ====================< End results >====================
    # sort the last generation computed 
    sort!(gen.pop, by = x -> fitness(x, TAG_FitGA))
    # handle flags (special stop condition)
    if flag_opti
    println("\n <:D> Optimum reached: $(round(fitness(gen.pop[1], TAG_FitGA), digits=3)) <:D>\n")
    elseif flag_tl
    println("\n <!> Time limit reached: $(time() - start_time) <!>\n")
    elseif flag_improv
    println("\n <!> No improvements in $(max_no_improv) generations <!>\n")
    end

    return gen.pop[1], flag_opti, gen_val, time() - start_time
end
