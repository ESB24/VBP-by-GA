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
    GA(instance::Instance, Pmax::Int64, Gmax::Int64, tl::Float64, fit::Function, crossover::Function, objective::Function)
```
> GA: Genetic Algoritm
"""


function GA_V3(
        instance            ::Instance                                                  ;
        tl                  ::Float64                       = 600.                      ,   # time limit
        TAG_FitGA           ::Type{<:FitnessSolution}       = ObjGA_2                     ,
        display_plot    ::Bool = false                              , # generate plots
        display_state   ::Bool = true                              , # debug info
    )
# ====================< INIT >====================

    R::Int64 = instance.nbRoute
    O::Int64 = instance.nbOut
    Lmax::Int64 = instance.Lmax
    α::Float64 = 0.4

    Pmax::Int64 = R
    Gmax::Int64 = R * 10

    elite_size::Int64 = round(Int64, Pmax * 2α)       # space occupied by the elite of the previous generation in the current core

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

    strat::Int64 = 1 # 1=BF, 2=FF, 3=NF, 4=WF
    strat_name::String = "?"
    f::Union{Function, Nothing} = nothing
    if strat == 1
        f = BFD_EM
        strat_name = "BFD_EM"
    elseif strat == 2
        f = FFD_EM
        strat_name = "FFD_EM"
    elseif strat == 3
        f = NFD_EM
        strat_name = "NFD_EM"
    elseif strat == 4
        f = WFD_EM
        strat_name = "WFD_EM"
    elseif strat == 5
        f = BFD_SA
        strat_name = "BFD_SA"
    elseif strat == 6
        f = FFD_SA
        strat_name = "FFD_SA"
    elseif strat == 7
        f = NFD_SA
        strat_name = "NFD_SA"
    elseif strat == 8
        f = WFD_SA
        strat_name = "WFD_SA"
    end

    # results & geneneration management
    gen_val ::Vector{GenVal} = [GenVal([], [], [])]  # Result of the GA (objective value of each individuals. Intended for plotting)
    gen     ::Generation     = Generation(0, [])     # Current generation

    # Sorting criteria for FFD (used to create the first generation)
    all_TAG = (MaxMin, Max, MailNb, MailVolume, MeanVolume, MeanAssignment, MailStd, AssignmentStd)
    all_sort_crit::Vector{Tuple} = []
    for TAG1 in all_TAG
        for TAG2 in all_TAG
            (TAG1 != TAG2) && push!(all_sort_crit, (TAG1, TAG2))
        end
    end

# ====================< GEN LOOP >====================

    while (gen.id <= Gmax) && !flag_opti && !flag_tl && !flag_improv
        display_state && (print("\n\n \x1b[36m<-> Generation $(gen.id):\x1b[38;5;250m"))

# ====================< Initial Population/ Elite selection >====================

        if gen.id == 0 
            display_state && (print("\n    \x1b[36m<◿> Sorted perm $(strat_name): (x$(length(all_sort_crit)))\x1b[38;5;250m\n        "))
            while length(gen.pop) < length(all_sort_crit) && !flag_opti && !flag_tl
                Tags ::Tuple         = all_sort_crit[length(gen.pop)+1]                                                     # get sorting criteria tags (for FFD)
                perm ::Vector{Int64} = sortperm([(-fitness(r, Tags[1]), -fitness(r, Tags[2])) for r in instance.route])    # compute used permutation (according to sorting criteria)
                
                tmp, _ = f(instance.route, Lmax, O, R, perm, display_state=false)
                sol::Solution = Solution(perm, tmp)
                val ::Float64  = fitness(sol, TAG_FitGA)

                push!(gen.pop, sol)                                                 # Store Solution
                push!(gen_val[gen.id+1].elite, val)                                 # Store objective value
                
                flag_opti = length(sol.sessions) <= opti_sol                        # test optimality
                flag_tl   = time() - start_time > tl                          # test time limit

                display_state && (print("($(all_sort_crit[length(gen.pop)]): $(round(val, digits=3)))$(if length(gen_val[gen.id+1].elite)%5 == 0 && length(gen_val[gen.id+1].elite) != 0 "\n        " else ", " end)"))
            end
        else
            gen.pop = sort!(gen.pop, by = x -> fitness(x, TAG_FitGA))[1:elite_size]     # next gen initialisation
            push!(gen_val, GenVal([fitness(e, TAG_FitGA) for e in gen.pop], [], []))    # Store objective value

            if lastElite === nothing                                        # If gen n° == 1
                lastElite = sum(gen_val[end].elite)                         #  |  store generation score
            else                                                            # else
                currentElite = sum(gen_val[end].elite)                      #  |  compute current generation score
                if lastElite > currentElite                                 #  |  If current generation better than last generation
                    lastElite = currentElite                                #  |   |  update current generation score
                else                                                        #  |  else
                    # no_improv += 1                                          #  |   |  
                    flag_improv = true # no_improv > max_no_improv

                    # mod = (mod + 1)%3 + 1             

                    # print(" (no improvement, mod switched to $(mod==1 ? "FF" : mod==3 ? "BF" : "WF"))")
                end                                                         #  |  endif
            end                                                             # endif

            display_state && (print("\n    \x1b[36m<♔> Elite selection: (x$(elite_size))\x1b[38;5;250m\n        "))
            for (k, e) in enumerate(gen_val[gen.id+1].elite)
                display_state && (print("($(round(e, digits=3)))$(if k%10 == 0 && k != 0 "\n        " else ", " end)"))
            end
        end

# ====================< Random individuals >====================

    print("\n    \x1b[36m<🎲> $(strat_name): (x$(if gen.id == 0 Pmax - length(all_sort_crit) else Pmax - elite_size end))\x1b[38;5;250m\n        ")
    while length(gen.pop) < Pmax && !flag_opti && !flag_tl && !flag_improv
        perm = randperm(R)
        tmp, _ = f(instance.route, Lmax, O, R, perm, display_state=false)
        sol::Solution = Solution(perm, tmp)
        val ::Float64  = fitness(sol, TAG_FitGA)

        push!(gen.pop, sol)                          # Store Solution
        push!(gen_val[gen.id+1].ffr, val)            # Store objective value

        flag_opti = length(sol.sessions) <= opti_sol # test optimality
        flag_tl   = time() - start_time > tl      # test time limit

        print("($(round(val, digits=3)))$(if length(gen_val[gen.id+1].ffr)%10 == 0 && length(gen_val[gen.id+1].ffr) != 0 "\n        " else ", " end)")
    end

# ====================< Extended Population -> Cross-Over, Mutation >====================

        print("\n    \x1b[36m<🧬> Cross-Over: (x$(Pmax))\x1b[38;5;250m\n        ")
        parent_order::Vector{Int64} = randperm(length(gen.pop))
        while length(gen.pop) < 2Pmax && !flag_opti && !flag_tl && !flag_improv
            p1::Solution = gen.pop[popfirst!(parent_order)] # Get the next two parents
            p2::Solution = gen.pop[popfirst!(parent_order)]

            perm_son1::Vector{Int64} = mutation_swap!(crossOver_onePoint(p1.permutation, p1.permutation))   # compute the first son permutation CrossOver 1 point + swap
            perm_son2::Vector{Int64} = mutation_swap!(crossOver_onePoint(p2.permutation, p1.permutation))   # compute the second son permutation CrossOver 1 point + swap

            tmp, _ = f(instance.route, Lmax, O, R, perm_son1, display_state=false)
            sol_son1::Solution = Solution(perm_son1, tmp)
            val_son1 ::Float64  = fitness(sol_son1, TAG_FitGA)

            push!(gen.pop, sol_son1)                          # add first son to the population
            push!(gen_val[gen.id+1].cross_over, val_son1)     # add the first son objective value to the generation value

            flag_opti = length(sol_son1.sessions) <= opti_sol
            flag_tl   = time() - start_time > tl
            if flag_opti || flag_tl
                print("($(round(fitness(p1, TAG_FitGA), digits=3)) & $(round(fitness(p2, TAG_FitGA), digits=3)) | CrossOver 1Point + Swap ->  $(round(val_son1, digits=3)))$(if length(gen_val[gen.id+1].cross_over)%4 == 0 && length(gen_val[gen.id+1].cross_over) != 0 "\n        " else ", " end)")
            else
                tmp, _ = f(instance.route, Lmax, O, R, perm_son2, display_state=false)
                sol_son2::Solution = Solution(perm_son2, tmp)
                val_son2 ::Float64  = fitness(sol_son2, TAG_FitGA)

                push!(gen.pop, sol_son2)                          # add second son to the population
                push!(gen_val[gen.id+1].cross_over, val_son2)     # add the second son objective value to the generation value

                # update other component
                flag_opti = length(sol_son2.sessions) <= opti_sol
                flag_tl   = time() - start_time > tl

                print("($(round(fitness(p1, TAG_FitGA), digits=3)) & $(round(fitness(p2, TAG_FitGA), digits=3)) | CrossOver 1Point + Swap + $(strat_name) ->  $(round(val_son1, digits=3)) and $(round(val_son2, digits=3)))$(if length(gen_val[gen.id+1].cross_over)%4 == 0 && length(gen_val[gen.id+1].cross_over) != 0 "\n        " else ", " end)")
            end
        end

        println("\n    \x1b[36m<◷> Elapsed time: $(round(time() - start_time, digits=3))s\x1b[0m")
        gen.id += 1
    end

# ====================< End results >====================
    # sort the last generation computed 
    sort!(gen.pop, by = x -> fitness(x, TAG_FitGA))
    # handle flags (special stop condition)
    if flag_opti
        display_state && (println("\n \x1b[36m<:D> Optimum reached: $(round(fitness(gen.pop[1], TAG_FitGA), digits=3)) <:D>\x1b[0m\n"))
    elseif flag_tl
        display_state && (println("\n \x1b[36m<!> Time limit reached: $(time() - start_time) <!>\x1b[0m\n"))
    elseif flag_improv
        display_state && (println("\n \x1b[36m<!> No improvements in $(max_no_improv) generations <!>\x1b[0m\n"))
    end

# ====================< Plot >====================
    if display_plot
        plot1 = plot(title="Genetic algorithm", size = (1000, 1000))

        scatter!([1 + 0.1 for _=1:length(gen_val[1].cross_over)], gen_val[1].cross_over, label="cross over", markercolor=:blue, lw=1.5)
        scatter!([1 for _=1:length(gen_val[1].ffr)], gen_val[1].ffr, label="random perm", markercolor=:orange, lw=1.5)
        scatter!([1 - 0.1 for _=1:length(gen_val[1].elite)], gen_val[1].elite, label="best α% and sorted perm", markercolor=:green, lw=1.5)

        for i=2:length(gen_val)
            scatter!([i + 0.1 for _=1:length(gen_val[i].cross_over)], gen_val[i].cross_over, label="", markercolor=:blue, lw=1.5)
            scatter!([i for _=1:length(gen_val[i].ffr)], gen_val[i].ffr, label="", markercolor=:orange, lw=1.5)
            scatter!([i - 0.1 for _=1:length(gen_val[i].elite)], gen_val[i].elite, label="", markercolor=:green, lw=1.5)
        end

        plot!(
            collect(Iterators.flatten([[i-0.1, Float64(i), i+0.1] for i=1:length(gen_val)])), 
            collect(Iterators.flatten([[minimum(gen_val[i].elite), minimum([gen_val[i].elite; gen_val[i].ffr]), minimum([gen_val[i].elite; gen_val[i].ffr; gen_val[i].cross_over])] for i=1:length(gen_val)])), 
            lc="green", label="best solution", lw=1.5)

        display(plot1)
    end

    return gen.pop[1], flag_opti, gen_val, time() - start_time
end
