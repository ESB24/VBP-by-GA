
begin # Files
    include("Route.jl")
    include("Session.jl")
    include("DataStruct.jl")
    include("MyDistribution.jl")
end

# ================================================== #
#                        Test                        #
# ================================================== #

function minSession(instance::Instance)::Int64
    return max(ceil(Int64, instance.nbRoute/instance.nbOut), ceil(Int64, sum([sum(r.mail) for r in instance.route]) / (instance.Lmax * instance.nbOut)))
end

function rawInstance(instance::Instance)
    for r in instance.route
        rawRoute(r) ? print("-") : print("x")
    end
end


# ================================================== #
#                    Initialisation                  #
# ================================================== #

function Instance(mat::Matrix{Int64}, λ::Float64 = 0.25)
    R::Int64, O::Int64 = size(mat)

    # maximum weighted sum of each round
    Kmax = maximum(mat * [λ for _=1:O])

    # Biggest batch in the instance
    vmax = maximum(mat)

    # machine outputs capacity
    C::Int64 = ceil(max(vmax , Kmax))

    return Instance(C, O, R, [Route(i, O, (filter(x -> x != 0, mat[i, :])..., )) for i=1:R])
end

# ================================================== #
#                       Display                      #
# ================================================== #

function Base.show(io::IO, i::Instance)
    print(io, "Instance:\n     C : $(i.Lmax)\n    |O|: $(i.nbOut)\n    |R|: $(i.nbRoute)\n     R :\n")

    for r in i.route
        println(io, r)
    end
end

function printASCII(i::Instance)
    print("Instance:\n     C : $(i.Lmax)\n    |O|: $(i.nbOut)\n    |R|: $(i.nbRoute)\n     R :\n")
    for r::Route in i.route
        print("$(r.id): ")
        for b::Int64 in r.assignment
            if b == 0
                print(".")
            else
                if b <= i.Lmax/5
                    print("_")
                else
                    print("|")
                end
            end
        end
        println()
    end
end

# ================================================================================= #
#                           ######   #######  ######   ##   ##                      #
#                           #     #  #        #     #  # # # #                      #
#                           ######   ####     ######   #  #  #                      #
#                           #        #        #    #   #     #                      #
#                           #        #######  #     #  #     #                      #
# ================================================================================= #

function compute_perm(instance::Instance, N::Int64, nb_sol::Int64 = 10, tl::Int64 = 60)
    R = collect(1:instance.nbRoute)
    V = [sum(instance.route[r].mail) for r in R]
    O = instance.nbOut
    Lmax = instance.Lmax
    S = collect(1:N)

    # println("========\nR = $R\nB = $B\nV = $V\nJ = $J\nC = $C\nO = $O\n")
    env=Gurobi.Env()
    model = Model(() -> (Gurobi.Optimizer(env)))
    # set_silent(model)
    # set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", tl)

    @variable(model, y[s in S, r in R], Bin)
    @variable(model, z[s in S], Bin)

    # warmup
    bound = minSession(instance)
    for s=1:bound
        set_start_value(z[s], 1)
    end

    @objective(model, Min, sum(z))

    @constraint(model, sum(z) >= bound)

    # activation y
    for r in R
        @constraint(model, sum([y[s, r] for s in S]) == 1)
    end

    # Capacity + activation of z
    for s in S
        @constraint(model, sum([y[s, r] * V[r] for r in R]) <= O * Lmax)
    end

    # number of route per session
    for s in S
        @constraint(model, sum(y[s, :]) <= O * z[s])
    end

    # use z[s] before z[s+1]
    # for s in S[1:end-1]
    #     @constraint(model, z[s] >= z[s+1])
    # end

    optimize!(model)

    if termination_status(model) == OPTIMAL || MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0
        res = [[] for _ in S]
        loads = zeros(Int64, length(S))

        for s in S
            for r in R
                (value(y[s, r]) != 0) && (push!(res[s], r); loads[s] += V[r])
            end
        end
        println("==========\n    O * Lmax = $(O * Lmax)\n    status ? $(termination_status(model))\n    loads -> $(loads)\n==========")

        return res, model
    else
        println("==========\n    O * Lmax = $(O * Lmax)\n    no solution. \n==========")

        return nothing, model
    end
end

# ================================================================================= #
#  ######   #######   #####            ######    #####   #######   #####   #     #  #
#  #     #     #     #                 #     #  #     #     #     #     #  #     #  #
#  ######      #     #  ###            ######   #######     #     #        #######  #
#  #     #     #     #     #           #     #  #     #     #     #     #  #     #  #
#  #######  #######   #####            #######  #     #     #      #####   #     #  #
# ================================================================================= #

function createInstanceV1(
        nbRoute             ::Int64                         = 200                   ,   # number of round in the instance
        nbOutputs           ::Int64                         = 200                   ,   # number of outputs per round
        λ                   ::Float64                       = 0.25                  # parameter for the capacity of the outputs
    )::Instance
    mat                     ::Matrix{Int64}                 = Matrix{Int64}(undef, nbRoute, nbOutputs)

    minBatchesNb            ::Int64                         = round(Int64, 2*nbOutputs/3)    # minimum number of batches per round
    maxBatchesNb            ::Int64                         = nbOutputs             # maximum number of batches per round
    BatchesSizeProba        ::Vector{Float64}               = [1. / convert(Float64, nbOutputs*1.5), 1. - (1. / convert(Float64, nbOutputs*1.5))]            # distribution of the batches size
    BatchesSizeInter        ::Vector{UnitRange{Int64}}      = [nbOutputs*3:round(Int64, nbOutputs*3.2), 1:round(Int64, nbOutputs/10)]          # size of the batches according to the distribution

    for i=1:nbRoute
        nbBatches::Int64 = rand(minBatchesNb:maxBatchesNb)
        for j=1:nbOutputs
            if j <= nbBatches
                mat[i, j] = rand(BatchesSizeInter[rand(Categorical(BatchesSizeProba))])
            else
                mat[i, j] = 0
            end
        end 
    end

    return Instance(mat, λ)
end

function createInstanceV2(
        nbRoute             ::Int64                         = 200                   ,   # number of round in the instance
        nbOutputs           ::Int64                         = 200                   ,
        λ                   ::Float64                       = 0.25                     # number of outputs per round
    )::Instance
    mat                     ::Matrix{Int64}                 = Matrix{Int64}(undef, nbRoute, nbOutputs)

    minBatchesNb            ::Int64                         = round(Int64, 2*nbOutputs/3)    # minimum number of batches per round
    maxBatchesNb            ::Int64                         = nbOutputs             # maximum number of batches per round
    BatchesSizeProba        ::Vector{Float64}               = [1, 0]            # distribution of the batches size
    BatchesSizeInter        ::Vector{UnitRange{Int64}}      = [round(Int64, nbOutputs*0.95):nbOutputs, 1:1]          # size of the batches according to the distribution

    nbBigBatches::Int64 = 3

    for i=1:nbRoute
        nbBatches::Int64 = rand(minBatchesNb:maxBatchesNb)

        bigPos::Vector{Int64} = randperm(nbBatches)[1:nbBigBatches]

        for j=1:nbOutputs
            if j <= nbBatches
                if j in bigPos
                    mat[i, j] = rand(BatchesSizeInter[rand(Categorical(BatchesSizeProba))])
                else
                    mat[i, j] = 1
                end
            else
                mat[i, j] = 0
            end
        end 
    end

    return Instance(mat, λ)
end

function generateRoutes_BigBatche(
        nbRoute         ::Int64             ,
        smallBatches    ::UnitRange{Int64}  ,
        bigBatches      ::UnitRange{Int64}  ,
        amontBigBatches ::UnitRange{Int64}  ,
        minBatches      ::Int64             ,
        maxBatches      ::Int64             ,
        Lmax            ::Int64             , 
        O               ::Int64             ,
        currentId       ::Int64             ,
    )
# ==========< Generate Batches >==========
    loads::Vector{Int64} = [Lmax for _=1:O]

    batches::Vector{Vector{Int64}} = [Int64[] for _=1:O]
    nbBatches::Int64 = 0

    while sum(loads) > 0
        isBig::Bool = (rand(amontBigBatches) == 1)
        maxLoad = maximum(loads)

        if isBig
            if maxLoad >= bigBatches[1] # can create a big batch
                
                size = bigBatches[1]:(maxLoad >= bigBatches[end] ? bigBatches[end] : maxLoad)
                pos = sortperm(loads, rev=true)[1]

                b = rand(size)

                loads[pos] - b < 0 && (println("error1"); global error += 1)

                push!(batches[pos], b)
                loads[pos] -= b
                nbBatches += 1
            else
                isBig = false
            end
        end

        if !isBig
            size = smallBatches[1]:(maxLoad >= smallBatches[end] ? smallBatches[end] : maxLoad)
            pos = sortperm(loads, rev=true)[1]

            b = rand(size)

            loads[pos] - b < 0 && (println("error2"); global error += 1)

            push!(batches[pos], b)
            loads[pos] -= b
            nbBatches += 1
        end
    end
# ==========< Assign Batches >==========
    
    roundAssignment::Vector{Vector{Int64}} = [zeros(Int64, O) for _=1:nbRoute]
    roundBatchesNumber::Vector{Int64} = zeros(Int64, nbRoute)

    while nbBatches > 0
        pos = sortperm(batches, by= x -> length(x), rev=true)[1]

        subRoute::Vector{Int64} = [] # sub set of batches that can receave the current batch
        priority::Bool = false
        for (k, v) in enumerate(roundAssignment)
            if v[pos] == 0
                if roundBatchesNumber[k] < minBatches
                    priority || (priority = true; subRoute = [])
                    push!(subRoute, k)    
                elseif !priority
                    push!(subRoute, k)
                end
            end
        end

        if !isempty(subRoute)
            shuffle!(subRoute)

            r = subRoute[1]

            roundAssignment[r][pos] += popfirst!(batches[pos])
            roundBatchesNumber[r] += 1
            nbBatches -= 1
        else
            assigned = false
            while !assigned
                r = rand(1:nbRoute)
                if roundAssignment[r][pos] + batches[pos][1] <= Lmax
                    roundAssignment[r][pos] += popfirst!(batches[pos])
                    nbBatches -= 1
                    assigned = true
                end
            end
        end
    end

# ==========< Build Routes >==========

    rounds::Vector{Route} = []
    for rid=1:nbRoute
        roundBatches::Vector{Int64} = filter(x -> x != 0, roundAssignment[rid])
        Br::Int64 = length(roundBatches)
        r::Route{Br} = Route(currentId, roundAssignment[rid], ntuple(i -> roundBatches[i], Br))
        push!(rounds, r)
        currentId += 1
    end

    return (rounds, currentId)
end

function generateSolution_BigBatche(
        nbSession           ::Int64             ,
        nbRoute             ::Int64             ,
        smallBatches        ::UnitRange{Int64}  ,
        bigBatches          ::UnitRange{Int64}  ,
        amontBigBatches     ::UnitRange{Int64}  ,
        minBatches          ::Int64             ,
        maxBatches          ::Int64             ,
        Lmax                ::Int64             , 
        O                   ::Int64             ,
    )

    totalRoutes             ::Int64 = 0
    currentId               ::Int64 = 1
    instanceRoutes          ::Vector{Vector{Int64}} = []
    sol                     ::Solution = Solution([], [])

    for _=1:nbSession
        rounds, currentId = generateRoutes_BigBatche(nbRoute, smallBatches, bigBatches, amontBigBatches, minBatches, maxBatches, Lmax, O, currentId)
        totalRoutes += length(rounds)
        
        push!(sol.sessions, Session(Lmax, rounds, compute_output(rounds)))
    end

    # shuffle indexes of each round
    perm                    ::Vector{Int64} = randperm(totalRoutes)
    sol.permutation = perm
    i = 1
    for s in sol.sessions
        for r in s.route
            r.id = perm[i]
            i += 1
        end
    end

    return sol
end

function writeSolution_BigBatche(path::String, sol::Solution, id::Int64 = 0)
    O::Int64 = length(sol.sessions[1].loads)
    nbRoutes::Int64 = sum([length(s.route) for s in sol.sessions])
    Lmax::Int64 = sol.sessions[1].Lmax
    nbSession::Int64 = length(sol.sessions)
    

    filename = "instanceBigBatche_$(id)_O$(O)_R$(nbRoutes)_C$(Lmax)_opt_$(nbSession).txt"
        
    # Create file
    println("instance generated : "*filename)
    fd = open("$(path)/$(filename)", "w+")

    write(fd, "\n$O\n$nbRoutes\n$Lmax\n$nbSession\n\n")

    # for e in sol.permutation
    #     write(fd, "$(e), ")
    # end
    # write(fd, "\n\n")

    for s in sol.sessions
        for r in s.route
            write(fd, "$(r.id): ")

            for a in r.assignment
                write(fd, "$(a), ")
            end

            write(fd, ": ")

            for b in r.mail
                write(fd, "$(b), ")
            end
            write(fd, "\n")
        end
        write(fd, "\n\n")
    end
    
    close(fd)
end

function writeInstanceBatterie_BigBatche(
        path                ::String            ,
        nbInstance          ::Int64             ,
        nbSession           ::Int64             ,
        nbRoute             ::Int64             ,
        O                   ::Int64             ,
        Lmax                ::Int64             , 
        minBatches          ::Int64             ,
        maxBatches          ::Int64             ,
        smallBatches        ::UnitRange{Int64}  ,
        bigBatches          ::UnitRange{Int64}  ,
        amontBigBatches     ::UnitRange{Int64}  ,
    )

    for i=1:nbInstance
        sol = generateSolution_BigBatche(nbSession, nbRoute, smallBatches, bigBatches, amontBigBatches, minBatches, maxBatches, Lmax, O)
        writeSolution_BigBatche(path, sol, i)
    end
end

# ================================================================================= #
#           ######   #######   ######  #######  ######   #######  ######            #
#           #     #     #     #           #     #     #     #     #     #           #
#           #     #     #      #####      #     ######      #     ######            #
#           #     #     #           #     #     #   #       #     #     #           #
#           ######   #######  ######      #     #     #  #######  #######           #
# ================================================================================= #

# begin
#     mat = parseMailXLSX("TER/data/Trafic/trafic_200_200_35_1.xlsx")
#     R, O = size(mat)
    

#     total = zeros(Float64, maximum(mat))

#     for i=1:maximum(mat)
#         total[i] = count(x -> x == i, mat)
#     end
# end

# begin 
#     fd = open("TER/data/Trafic/distrib_trafic.txt", "a")

#     write(fd, "O$(O): ")

#     for e in total
#         write(fd, "$(e), ")
#     end
#     write(fd, "\n")

#     close(fd)
# end

# begin
#     normalize([float(e) for e in total])
# end

# ================================================================================= #
#                #######  ######    #####     #######  #######   ######             #
#                   #     #     #  #     #    #           #     #                   #
#                   #     ######   #######    ####        #     #                   #
#                   #     #    #   #     #    #           #     #                   #
#                   #     #     #  #     #    #        #######   ######             #
# ================================================================================= #

function generateRoutes_Distrib_easy(R::Int64, O::Int64)

    distrib_sizeBatch   = distrib_trafic_batchSize(O)
    distrib_nbBatch     = distrib_trafic_batchCount(O)

    rounds::Vector{Route} = Vector{Route}(undef, R)

    for i=1:R
        nbBatch = rand(distrib_nbBatch)

        rounds[i] = Route(i, O, ntuple(i -> rand(distrib_sizeBatch), nbBatch))
    end

    return rounds
end

function writeRoutes_Distrib_easy(rounds::Vector{Route}, id::Int64 = 0)
    O::Int64 = length(rounds[1].assignment)
    R::Int64 = length(rounds)


    filename = "myTrafic_$(id)_O$(O)_R$(R).txt"
    path::String = "TER/data/MyTrafic/" 
        
    # Create file
    println("instance generated : "*filename)
    fd = open("$(path)/$(filename)", "w+")
    
    # Head
    write(fd, "$O\n$R\n")

    # write rounds
    for r in rounds
        write(fd, "$(r.id): ")

        for a in r.assignment
            write(fd, "$(a), ")
        end

        write(fd, ": ")

        for b in r.mail
            write(fd, "$(b), ")
        end
        write(fd, "\n")
    end
    write(fd, "\n\n")

    close(fd)
end

function writeInstanceBatterie_Distrib_easy(nbInstance::Int64, R::Int64, O::Int64)
    for i=1:nbInstance
        rounds = generateRoutes_Distrib_easy(R, O)
        writeRoutes_Distrib_easy(rounds, i)
    end
end

function parseMyInstance_easy(path::String)
    fd  = open(path, "r")

    O::Int64 = parse(Int64, readline(fd))
    R::Int64 = parse(Int64, readline(fd))

    mat::Matrix{Int64} = zeros(Int64, R, O)

    for i=1:R
        line = readline(fd)
        while line == ""
            line = readline(fd)
        end

        id, _, batches = split(line, ": ")

        id = parse(Int64, id)
        batches = parse.(Int64, split(batches, ", ")[1:end-1])

        mat[id, 1:length(batches)] = batches
    end

    close(fd)
    return mat
end



# ================================================================================= #
#                           #     #   #####   ######   ######                       #
#                           #     #  #     #  #     #  #     #                      #
#                           #######  #######  ######   #     #                      #
#                           #     #  #     #  #    #   #     #                      #
#                           #     #  #     #  #     #  ######                       #
# ================================================================================= #

function generateRoutes_Distrib(
        nbRoute         ::Int64             ,
        smallerBatches  ::Int64             ,
        biggestBatches  ::Int64             ,
        minBatches      ::Int64             ,
        maxBatches      ::Int64             ,
        Lmax            ::Int64             , 
        O               ::Int64             ,
        currentId       ::Int64             ,
        fct_distrib     ::Function          ,
    )

# ==========< Generate Batches >==========
    
    distrib = fct_distrib((biggestBatches+1) - smallerBatches)
    batcheValues = [smallerBatches + i for i=0:(biggestBatches - smallerBatches)]

    loads::Vector{Int64} = [Lmax for _=1:O]

    batches::Vector{Vector{Int64}} = [Int64[] for _=1:O]
    nbBatches::Int64 = 0

    while sum(loads) > 0
        # batche size = minimum value between smallest load left and randomly generated batch (according to the distribution)
        size = min(batcheValues[rand(distrib)], maximum(loads))

        # remove the batches from the most loaded output
        pos = sortperm(loads, rev=true)[1]

        # test validity (should not trigger)
        loads[pos] - size < 0 && (println("error1"); global error += 1)

        push!(batches[pos], size)
        loads[pos] -= size
        nbBatches += 1
    end

# ==========< Assign Batches >==========

    roundAssignment::Vector{Vector{Int64}} = [zeros(Int64, O) for _=1:nbRoute]
    roundBatchesNumber::Vector{Int64} = zeros(Int64, nbRoute)

    while nbBatches > 0
        pos = sortperm(batches, by= x -> length(x), rev=true)[1]

        subRoute::Vector{Int64} = [] # sub set of batches that can receave the current batch
        priority::Bool = false
        for (k, v) in enumerate(roundAssignment)
            if v[pos] == 0
                if roundBatchesNumber[k] < minBatches
                    priority || (priority = true; subRoute = [])
                    push!(subRoute, k)    
                elseif !priority
                    push!(subRoute, k)
                end
            end
        end

        if !isempty(subRoute)
            shuffle!(subRoute)

            r = subRoute[1]

            roundAssignment[r][pos] += popfirst!(batches[pos])
            roundBatchesNumber[r] += 1
            nbBatches -= 1
        else
            assigned = false
            while !assigned
                r = rand(1:nbRoute)
                if roundAssignment[r][pos] + batches[pos][1] <= Lmax
                    roundAssignment[r][pos] += popfirst!(batches[pos])
                    nbBatches -= 1
                    assigned = true
                end
            end
        end
    end

# ==========< Build Routes >==========

    rounds::Vector{Route} = []
    for rid=1:nbRoute
        roundBatches::Vector{Int64} = filter(x -> x != 0, roundAssignment[rid])
        Br::Int64 = length(roundBatches)
        r::Route{Br} = Route(currentId, roundAssignment[rid], ntuple(i -> roundBatches[i], Br))
        push!(rounds, r)
        currentId += 1
    end

    return (rounds, currentId)
end

function generateSolution_Distrib(
        nbSession       ::Int64             ,
        nbRoute         ::Int64             ,
        smallerBatches  ::Int64             ,
        biggestBatches  ::Int64             ,
        minBatches      ::Int64             ,
        maxBatches      ::Int64             ,
        Lmax            ::Int64             , 
        O               ::Int64             ,
        fct_distrib     ::Function          ,
    )

    totalRoutes             ::Int64 = 0
    currentId               ::Int64 = 1
    instanceRoutes          ::Vector{Vector{Int64}} = []
    sol                     ::Solution = Solution([], [])

    for _=1:nbSession
        rounds, currentId = generateRoutes_Distrib(nbRoute, smallerBatches, biggestBatches, minBatches, maxBatches, Lmax, O, currentId, fct_distrib)
        totalRoutes += length(rounds)
        
        push!(sol.sessions, Session(Lmax, rounds, compute_output(rounds)))
    end

    # shuffle indexes of each round
    perm                    ::Vector{Int64} = randperm(totalRoutes)
    sol.permutation = perm
    i = 1
    for s in sol.sessions
        for r in s.route
            r.id = perm[i]
            i += 1
        end
    end

    return sol
end

function writeSolution_Distrib(path::String, sol::Solution, id::Int64 = 0)
    O::Int64 = length(sol.sessions[1].loads)
    nbRoutes::Int64 = sum([length(s.route) for s in sol.sessions])
    Lmax::Int64 = sol.sessions[1].Lmax
    nbSession::Int64 = length(sol.sessions)


    filename = "instanceDistribHard_$(id)_O$(O)_R$(nbRoutes)_C$(Lmax)_opt_$(nbSession).txt"
        
    # Create file
    println("instance generated : "*filename)
    fd = open("$(path)/$(filename)", "w+")

    write(fd, "\n$O\n$nbRoutes\n$Lmax\n$nbSession\n\n")

    # for e in sol.permutation
    #     write(fd, "$(e), ")
    # end
    # write(fd, "\n\n")

    for s in sol.sessions
        for r in s.route
            write(fd, "$(r.id): ")

            for a in r.assignment
                write(fd, "$(a), ")
            end

            write(fd, ": ")

            for b in r.mail
                write(fd, "$(b), ")
            end
            write(fd, "\n")
        end
        write(fd, "\n\n")
    end

    close(fd)
end

function writeInstanceBatterie_Distrib(
        path            ::String            ,
        nbInstance      ::Int64             ,
        nbSession       ::Int64             ,
        nbRoute         ::Int64             ,
        O               ::Int64             ,
        Lmax            ::Int64             , 
        minBatches      ::Int64             ,
        maxBatches      ::Int64             ,
        smallerBatches  ::Int64             ,
        biggestBatches  ::Int64             ,
        fct_distrib     ::Function          ,
    )

    for i=1:nbInstance
        sol = generateSolution_Distrib(nbSession, nbRoute, smallerBatches, biggestBatches, minBatches, maxBatches, Lmax, O, fct_distrib)
        writeSolution_Distrib(path, sol, i)
    end
end

# ================================================================================= #
#                           #     #   #####   ######   ######                       #
#                           #     #  #     #  #     #  #     #                      #
#                           #######  #######  ######   #     #                      #
#                           #     #  #     #  #    #   #     #                      #
#                           #     #  #     #  #     #  ######                       #
# ================================================================================= #

function generateRoutes_Distrib_V2(
        nbRoute         ::Int64             ,
        smallerBatches  ::Int64             ,
        biggestBatches  ::Int64             ,
        minBatches      ::Int64             ,
        maxBatches      ::Int64             ,
        Lmax            ::Int64             , 
        O               ::Int64             ,
        currentId       ::Int64             ,
        fct_distribVolume::Function          ,
        fct_distribNumber::Function          ,
    )

# ==========< Generate Batches >==========

    distribVolume = fct_distribVolume((biggestBatches+1) - smallerBatches)
    batcheVolume = [smallerBatches + i for i=0:(biggestBatches - smallerBatches)]

    loads::Vector{Int64} = [Lmax for _=1:O]

    batches::Vector{Vector{Int64}} = [Int64[] for _=1:O]
    nbBatches::Int64 = 0

    while sum(loads) > 0
        # batche size = minimum value between smallest load left and randomly generated batch (according to the distribution)
        size = min(batcheVolume[rand(distribVolume)], maximum(loads))

        # remove the batches from the most loaded output
        pos = sortperm(loads, rev=true)[1]

        # test validity (should not trigger)
        loads[pos] - size < 0 && (println("error1"); global error += 1)

        push!(batches[pos], size)
        loads[pos] -= size
        nbBatches += 1
    end

# ==========< Assign Batches >==========
    
    distribNumber = fct_distribNumber((maxBatches+1) - minBatches)
    batcheNumber = [minBatches + i for i=0:(maxBatches - minBatches)]

    roundAssignment::Vector{Vector{Int64}} = [zeros(Int64, O) for _=1:nbRoute]
    roundBatchesNumber::Vector{Int64} = [batcheNumber[rand(distribNumber)] for _=1:nbRoute]

    while sum(roundBatchesNumber) > nbBatches
        pos = rand(1:nbRoute)

        if (roundBatchesNumber[pos] > minBatches) || (maximum(roundBatchesNumber) <= minBatches)
            roundBatchesNumber[pos] -= 1
        end
    end

    while sum(roundBatchesNumber) < nbBatches
        pos = rand(1:nbRoute)

        if (roundBatchesNumber[pos] < maxBatches) || (minimum(roundBatchesNumber) >= maxBatches)
            roundBatchesNumber[pos] += 1
        end
    end

    for r=1:nbRoute
        outputs = randperm(O)
        assigned = 0

        while assigned < roundBatchesNumber[r] && !isempty(outputs)
            k = popfirst!(outputs)

            if !isempty(batches[k])
                b = popfirst!(batches[k])

                roundAssignment[r][k] += b 
                assigned += 1
            end
        end
    end

    for k=1:O
        while !isempty(batches[k])
            roundAssignment[rand(1:nbRoute)][k] += popfirst!(batches[k])
        end
    end

# ==========< Build Routes >==========

    rounds::Vector{Route} = []
    for rid=1:nbRoute
        roundBatches::Vector{Int64} = filter(x -> x != 0, roundAssignment[rid])
        Br::Int64 = length(roundBatches)
        r::Route{Br} = Route(currentId, roundAssignment[rid], ntuple(i -> roundBatches[i], Br))
        push!(rounds, r)
        currentId += 1
    end

    return (rounds, currentId)
end

function generateSolution_Distrib_V2(
        nbSession       ::Int64             ,
        nbRoute         ::Int64             ,
        smallerBatches  ::Int64             ,
        biggestBatches  ::Int64             ,
        minBatches      ::Int64             ,
        maxBatches      ::Int64             ,
        Lmax            ::Int64             , 
        O               ::Int64             ,
        fct_distrib     ::Function          ,
        fct_distribNumber::Function         ,
    )

    totalRoutes             ::Int64 = 0
    currentId               ::Int64 = 1
    instanceRoutes          ::Vector{Vector{Int64}} = []
    sol                     ::Solution = Solution([], [])

    for _=1:nbSession
        rounds, currentId = generateRoutes_Distrib_V2(nbRoute, smallerBatches, biggestBatches, minBatches, maxBatches, Lmax, O, currentId, fct_distrib, fct_distribNumber)
        totalRoutes += length(rounds)
        
        push!(sol.sessions, Session(Lmax, rounds, compute_output(rounds)))
    end

    # shuffle indexes of each round
    perm                    ::Vector{Int64} = randperm(totalRoutes)
    sol.permutation = perm
    i = 1
    for s in sol.sessions
        for r in s.route
            r.id = perm[i]
            i += 1
        end
    end

    return sol
end

function writeSolution_Distrib_V2(path::String, sol::Solution, id::Int64 = 0)
    O::Int64 = length(sol.sessions[1].loads)
    nbRoutes::Int64 = sum([length(s.route) for s in sol.sessions])
    Lmax::Int64 = sol.sessions[1].Lmax
    nbSession::Int64 = length(sol.sessions)


    filename = "instanceDistribV2_$(id)_O$(O)_R$(nbRoutes)_C$(Lmax)_opt_$(nbSession).txt"
        
    # Create file
    println("instance generated : "*filename)
    fd = open("$(path)/$(filename)", "w+")

    write(fd, "\n$O\n$nbRoutes\n$Lmax\n$nbSession\n\n")

    # for e in sol.permutation
    #     write(fd, "$(e), ")
    # end
    # write(fd, "\n\n")

    for s in sol.sessions
        for r in s.route
            write(fd, "$(r.id): ")

            for a in r.assignment
                write(fd, "$(a), ")
            end

            write(fd, ": ")

            for b in r.mail
                write(fd, "$(b), ")
            end
            write(fd, "\n")
        end
        write(fd, "\n\n")
    end

    close(fd)
end

function writeInstanceBatterie_Distrib_V2(
        path            ::String            ,
        nbInstance      ::Int64             ,
        nbSession       ::Int64             ,
        nbRoute         ::Int64             ,
        O               ::Int64             ,
        Lmax            ::Int64             , 
        minBatches      ::Int64             ,
        maxBatches      ::Int64             ,
        fct_distribNumber::Function         ,
        smallerBatches  ::Int64             ,
        biggestBatches  ::Int64             ,
        fct_distrib     ::Function          ,
    )

    for i=1:nbInstance
        sol = generateSolution_Distrib_V2(nbSession, nbRoute, smallerBatches, biggestBatches, minBatches, maxBatches, Lmax, O, fct_distrib, fct_distribNumber)
        writeSolution_Distrib(path, sol, i)
    end
end

# ================================================================================= #
#        #####    #####   #     #   ######   ######  #######   #####   ##    #      #
#       #        #     #  #     #  #        #           #     #     #  # #   #      #
#       #  ###   #######  #     #   #####    #####      #     #######  #  #  #      #
#       #     #  #     #  #     #        #        #     #     #     #  #   # #      #
#        #####   #     #   #####   ######   ######   #######  #     #  #    ##      #
# ================================================================================= #

using PyPlot
using Random

# Define the Gaussian function
function gaussian(x::Float64; μ::Float64 = 0.0, σ::Float64 = 1.0)
    return (1 / (σ * sqrt(2 * π))) * exp(-((x - μ)^2) / (2 * σ^2))
end

"""
    return the image of x according to the plot which is a sum of gaussian with the all_σ and all_μ parameters 
"""
function gaussianSum(x::Float64, all_σ, all_μ)
    return sum([(i%2 == 0 ? 1 : -1 ) * gaussian(x, μ=μ, σ=σ) for (i ,(σ, μ)) in enumerate(zip(all_σ, all_μ))])
end

"""
    Create n plot wich are the sum of rangeGaussian randomly generated gaussian with O samples.
    return a matrix with the n*O points corresponding of the O points of each n plots normalized in [0, 1].
"""
function generateMultipleGaussianSum(O::Int64, n::Int64, rangeNbGaussian::UnitRange{Int64} = 6:8)
    matrix  ::Matrix{Float64} = zeros(Float64, n, O)
    x_values = 0.0:(10/O):(10-10/O)                                                     # range from 0 to near 10 with O values

    for i=1:n
        nbGaussian = rand(rangeNbGaussian)                                              # random number of gaussian

        all_σ = [rand(0.2:0.1:0.6) for _=1:nbGaussian]                                  # σ parameters for each gaussian
        all_μ = (randperm(100)[1:nbGaussian]./10)                                       # µ parameters for each gaussian

        matrix[i, :] = normalize([gaussianSum(x, all_σ, all_μ) for x in x_values])      # Calculate y values for Gaussian with μ=0, σ=1
    end

    return matrix
end

"""
Create a matrix containing the coordonate multiple plots (randomly generated combinaison of gaussian)
#### Parameters
 - O::Int64 : number of column of the matrix (j coordonates) <=> number of outputs per mmachine
 - nbPlots::Int64 : number of plot contained in the matrix
    → j →
 ↓  |######################################################################| Lmax
 i  |***                                         *****                     |
 ↓  |   **             ********               ***     ***       ***        |
    |     ****     ****        *****      ****           *    **   **   ***|
    |         *****                 ******                *  *       ***   |
    |                                                      **              |
    |######################################################################|
    |        *****                                            ***          |
    |     ***     **                **********              **   **        |
    |*****          *           ****          ****       ***       **      |
    |                ***     ***                  *******            ***   |
    |                   *****                                           ***|
    |######################################################################| 0
"""
function plotGaussianMatrix(O::Int64 = 40, nbPlots::Int64 = 5, rangeNbGaussian::UnitRange{Int64} = 16:20, Lmax::Int64 = 100, minVolume::Int64 = 2)
    tmp = generateMultipleGaussianSum(O, nbPlots, rangeNbGaussian)
    plotRange::Float64 = ((Lmax - minVolume) / nbPlots) - (minVolume) # Full range for each plot

    matrix::Matrix{Int64} = zeros(Int64, nbPlots+2, O)

    matrix[2:(nbPlots+1), :] = floor.(plotRange .* tmp) # scale and round each plot
    matrix[nbPlots+2, :] .= Lmax

    for i=2:nbPlots+1
        matrix[i, :] .+= round((i - 1) * minVolume + (i - 2) * plotRange)
    end

    batches::Vector{Int64} = []

    for k=1:O
        for i=1:nbPlots+1
            push!(batches, matrix[i+1, k] - matrix[i, k])
        end
    end

    for i=1:nbPlots+2
        plot(collect(1:O), matrix[i, :], label="μ=0, σ=1")
    end
end

function generateRoutes_Gaussian(
        nbRoute         ::Int64             ,
        minVolume       ::Int64             ,
        minBatches      ::Int64             ,
        maxBatches      ::Int64             ,
        nbPlots         ::Int64             ,
        nbGaussian      ::UnitRange{Int64}  ,
        Lmax            ::Int64             , 
        O               ::Int64             ,
        currentId       ::Int64             ,
    )
# ==========< Generate Batches >==========

    tmp = generateMultipleGaussianSum(O, nbPlots, nbGaussian)
    plotRange::Float64 = ((Lmax - minVolume) / nbPlots) - (minVolume) # Full range for each plot

    matrix::Matrix{Int64} = zeros(Int64, nbPlots+2, O)

    matrix[2:(nbPlots+1), :] = floor.(plotRange .* tmp) # scale and round each plot
    matrix[nbPlots+2, :] .= Lmax

    for i=2:nbPlots+1
        matrix[i, :] .+= round((i - 1) * minVolume + (i - 2) * plotRange)
    end

    batches::Vector{Int64} = []

    for k=1:O
        for i=1:nbPlots+1
            matrix[i+1, k] != matrix[i, k] && push!(batches, matrix[i+1, k] - matrix[i, k])
        end
    end

# ==========< Assign Batches >==========
    roundBatches::Vector{Vector{Int64}} = [[] for _=1:nbRoute]
    roundAssignments::Vector{Vector{Int64}} = [zeros(Int64, O) for _=1:nbRoute]

    batchesId::Vector{Int64} = randperm(length(batches))

    # assign minBatches batches to each rounds
    for i=1:nbRoute
        while length(roundBatches[i]) < minBatches
            id = popfirst!(batchesId)

            pos = round(Int64, 1+floor((id-1)/(nbPlots+1)))
            if roundAssignments[i][pos] != 0
                push!(batchesId, id)
            else
                roundAssignments[i][pos] = batches[id]
                push!(roundBatches[i], id)
            end
        end
    end

    # assign up to maxBatches batches to each rounds
    nonFullRoutes = collect(1:nbRoute)
    while !isempty(batchesId)
        bid = popfirst!(batchesId) # batch id
        rid = rand(nonFullRoutes)  # round id

        pos = round(Int64, 1+floor((bid-1)/(nbPlots+1)))
        if roundAssignments[rid][pos] != 0
            push!(batchesId, bid)
        else
            roundAssignments[rid][pos] = batches[bid]
            push!(roundBatches[rid], bid)

            (length(roundBatches[rid]) >= maxBatches) && (filter!(x -> x != rid, nonFullRoutes))
        end
    end

# ==========< Build Routes >==========
    rounds::Vector{Route} = []
    for rid=1:nbRoute
        Br::Int64 = length(roundBatches[rid])
        sort!(roundBatches[rid])
        r::Route{Br} = Route(currentId, roundAssignments[rid], ntuple(i -> batches[roundBatches[rid][i]], Br))
        push!(rounds, r)
        currentId += 1
    end

    return (rounds, currentId)
end

# tmp, _ = generateRoutes_Gaussian(30, 2, 8, 16:20, 50, 20, 1)

# for e in tmp
#     println(e)                                                                                                                                                                                   
# end

function generateSolution_Gaussian(
        nbSession           ::Int64 = 5,
        O                   ::Int64 = 40,
        Lmax                ::Int64 = 50,
        roundPerSession     ::Int64 = 20,
        minBatches          ::Int64 = 2,
        maxBatches          ::Int64 = 20,
        minVolume           ::Int64 = 2,
        nbPlots             ::Int64 = 5,
        nbGaussian          ::UnitRange{Int64} = 16:20,
    )

    totalRoutes             ::Int64 = 0
    currentId               ::Int64 = 1
    instanceRoutes          ::Vector{Vector{Int64}} = []
    sol                     ::Solution = Solution([], [])

    for _=1:nbSession
        rounds, currentId = generateRoutes_Gaussian(roundPerSession, minVolume, minBatches, maxBatches, nbPlots, nbGaussian, Lmax, O, currentId)
        totalRoutes += length(rounds)
        
        push!(sol.sessions, Session(Lmax, rounds, compute_output(rounds)))
    end

    perm                    ::Vector{Int64} = randperm(totalRoutes)

    sol.permutation = perm

    i = 1
    for s in sol.sessions
        for r in s.route
            r.id = perm[i]
            i += 1
        end
    end

    return sol
end

function writeSolution_Gaussian(path::String, sol::Solution, id::Int64 = 0)
    O::Int64 = length(sol.sessions[1].loads)
    nbRoutes::Int64 = sum([length(s.route) for s in sol.sessions])
    Lmax::Int64 = sol.sessions[1].Lmax
    nbSession::Int64 = length(sol.sessions)
    

    filename = "instanceGaussian_$(id)_O$(O)_R$(nbRoutes)_C$(Lmax)_opt_$(nbSession).txt"
        
    # Create file
    println("instance generated : "*filename)
    fd = open("$(path)/$(filename)", "w+")

    write(fd, "\n$O\n$nbRoutes\n$Lmax\n$nbSession\n\n")

    # for e in sol.permutation
    #     write(fd, "$(e), ")
    # end
    # write(fd, "\n\n")

    for s in sol.sessions
        for r in s.route
            write(fd, "$(r.id): ")

            for a in r.assignment
                write(fd, "$(a), ")
            end

            write(fd, ": ")

            for b in r.mail
                write(fd, "$(b), ")
            end
            write(fd, "\n")
        end
        write(fd, "\n\n")
    end
    
    close(fd)
end

function writeInstanceBatterie_Gaussian(
        path                ::String,
        nbInstance          ::Int64,
        nbSession           ::Int64,
        nbRoute             ::Int64,
        O                   ::Int64,
        Lmax                ::Int64,
        minBatches          ::Int64,
        maxBatches          ::Int64,
        minVolume           ::Int64,
        nbPlots             ::Int64,
        nbGaussian          ::UnitRange{Int64},
    )

    for i=1:nbInstance
        sol = generateSolution_Gaussian(nbSession, O, Lmax, nbRoute, minBatches, maxBatches, minVolume, nbPlots, nbGaussian)
        writeSolution_Gaussian(path, sol, i)
    end
end

# TODO if needed
# function parseSolution_Gaussian()
# end

function parseMyInstance(path::String)
    fd  = open(path, "r")

    _ = readline(fd) # ""

    O::Int64 = parse(Int64, readline(fd))
    nbRoutes::Int64 = parse(Int64, readline(fd))
    Lmax::Int64 = parse(Int64, readline(fd))
    nbSession::Int64 = parse(Int64, readline(fd))

    _ = readline(fd) # ""

    mat::Matrix{Int64} = zeros(Int64, nbRoutes, O)

    for i=1:nbRoutes
        # eat each empty line
        line = readline(fd)
        while line == ""
            line = readline(fd)
        end
        # println(line)

        id, _, batches = split(line, ": ")

        id = parse(Int64, id)
        batches = parse.(Int64, split(batches, ", ")[1:end-1])

        mat[id, 1:length(batches)] = batches
    end

    close(fd)
    return Lmax, mat, nbSession
end

function parseMyInstance_completed(path::String)
    fd  = open(path, "r")

    _ = readline(fd) # ""

    O::Int64 = parse(Int64, readline(fd))
    nbRoutes::Int64 = parse(Int64, readline(fd))
    Lmax::Int64 = parse(Int64, readline(fd))
    nbSession::Int64 = parse(Int64, readline(fd))

    _ = readline(fd) # ""

    mat::Matrix{Int64} = zeros(Int64, nbRoutes, O)

    for i=1:nbRoutes
        # eat each empty line
        line = readline(fd)
        while line == ""
            line = readline(fd)
        end
        # println(line)

        id, assignment, batches = split(line, ": ")

        id = parse(Int64, id)
        assignment = parse.(Int64, split(assignment, ", ")[1:end-1])

        mat[id, 1:length(assignment)] = assignment
    end

    close(fd)
    return Lmax, mat, nbSession
end


# ================================================================================= #
#           ######    #####   #######   ######   ######   #####   ##    #           #
#           #     #  #     #     #     #        #        #     #  # #   #           #
#           ######   #     #     #      #####    #####   #     #  #  #  #           #
#           #        #     #     #           #        #  #     #  #   # #           #
#           #         #####   #######  ######   ######    #####   #    ##           #
# ================================================================================= #

function parseOptiInstance(path::String)
    
    file = open(path,"r")

    _ = readline(file) # ""
    numberOfOutputs = parse(Int64, readline(file))
    numberOfRoutes = parse(Int64, readline(file))
    numberOfSessions = parse(Int64, readline(file))
    Lmax = parse(Int64, readline(file))
    LmaxLastBin = parse(Int64, readline(file))

    sessions = []

    readline(file) # Skips the first empty line
    # Iterate over the different bins of the 
    for sessionIndex in 1:numberOfSessions
        session = Vector{Vector{Int64}}(undef, 0)
        line = readline(file)
        while !isempty(line)
            println(line)

            valuesAsString = split(line, ",")
            a = [parse(Int64, valuesAsString[i]) for i in 1:length(valuesAsString)-1]
            push!(session, a)
            line = readline(file)
        end
        push!(sessions, session)
        println()
    end

    return numberOfOutputs, numberOfRoutes, numberOfSessions, sessions, Lmax
end

function normalizeOptiInstance(sessions::Vector)
    roundsLeftAligned = []
    rounds = reduce(vcat, sessions) # Puts all rounds of the sessions into a single unique session
    zeroes = [0 for i in 1:length(first(rounds))]

    for round in rounds
        z = deepcopy(zeroes)
        for (indexValue, value) in enumerate([e for e in round if e != 0])
            z[indexValue] = value
        end
        push!(roundsLeftAligned, z)
    end

    return roundsLeftAligned
end

function readjustSupportVector(vp::Vector, index::Int64)
    v2::Vector{Float64} = Vector{Int64}(undef, 0)
    pk = vp[index]
    for (i,p) in enumerate(vp)
        if i != index
            p2 = p + (p/(1-pk))*pk
            push!(v2, p2)
        else
            push!(v2, 0.0)
        end
    end
    vp = v2

    return vp
end

function readjustSupportVector(vx::Vector, indexes::Vector{Int64})
    indexes2 = sort(indexes, rev=true)
    vp = deepcopy(vx)

    v2::Vector{Float64} = Vector{Int64}(undef, 0)
    for index in reverse!(indexes2)
        v2 = Vector{Int64}(undef, 0)
        pk = vp[index]
        for (i,p) in enumerate(vp)
            if i != index
                p2 = p + (p/(1-pk))*pk
                push!(v2, p2)
            else
                push!(v2, 0.0)
            end
        end
        vp = v2
    end

    return vp
end

function countCat(c, n::Int64)
    res = [0 for i in 1:length(c.p)]

    for i in 1:n
        k = rand(c)
        res[k] += 1
    end

    return res
end

function generateHardInstance(Lmax::Int64, n_outputs::Int64, d_batchNumbers::Categorical, d_batchPositions::Categorical, d_batchValue::Categorical)
    sessionRoutes = Vector{Vector{Int64}}(undef,0)
    a = [Lmax for i in 1:n_outputs]

    roundNumber = 1
    while (a'*[1 for i in eachindex(a)]) != 0
        
        nonActiveIndexes = [i for (i,e) in enumerate(a) if e == 0] # The indexes of the vector a that have a value of 0 
        # START CREATION OF THE NEXT ROUND
        round = [0 for i in 1:n_outputs]

        # Number of batches in the round
        newSupport = readjustSupportVector(d_batchNumbers.p, [i for i in (n_outputs-length(nonActiveIndexes)+1):n_outputs])
        n_batchNumber = rand(Categorical(newSupport))

        # Positions of the batches in the round
        positions = Vector{Int64}(undef,0)
        for j in 1:n_batchNumber
            newSupport = readjustSupportVector(d_batchPositions.p, nonActiveIndexes)
            position = rand(Categorical(newSupport))
            push!(nonActiveIndexes, position)
            push!(positions, position)
        end

        # Value of each batch at each position
        for position in positions
            nonActiveIndexesForPosition = [i for i in a[position]+1:Lmax] # Set the list of impossible values for the round
            value = rand(Categorical(readjustSupportVector(d_batchValue.p, nonActiveIndexesForPosition))) # Pick a random value
            round[position] = value # Update the round
            a[position] -= value # Update the aggregate vector
        end

        push!(sessionRoutes, round)

        roundNumber += 1
    end

    return sessionRoutes
end

function vect2D_to_matrix(v)
    mat::Matrix{Int64} = Matrix{Int64}(undef, length(v), length(v[1]))

    for i=1:length(v)
        for j=1:length(v[1])
            mat[i, j] = v[i][j]
        end
    end

    return mat
end

function test_gen_instance_0()
    Lmax = 10
    n_outputs = 10
    d_batchNumbers = Categorical([.1 for i in 1:n_outputs])
    d_batchPositions = Categorical([.1 for i in 1:n_outputs])
    d_batchValue = Categorical([0.4, 0.15, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])

    return generateHardInstance(Lmax, n_outputs, d_batchNumbers, d_batchPositions, d_batchValue)
end

function restrict_distribution(f::Function, n::Int64)
    d = [f(k) for k in 1:n]
    push!(d, 1-(d'*[1 for i in 1:length(d)]))

    return readjustSupportVector(d, length(d))
end

function poissonToCategorical(λ::Float64, n::Int64)
    support = Vector{Float64}(undef,0)

    for k in 1:n
        push!(support, (λ^k * exp(-λ))/factorial(big(k)))
    end

    sum = support'*[1 for i in 1:length(support)]
    push!(support, 1-sum)

    support = readjustSupportVector(support, length(support))
    deleteat!(support, length(support))

    return Categorical(support)
end

function uniformToCategorical(n::Int64)
    support = [1/n for i in 1:n]

    return Categorical(support)
end

function binomialToCategorical(n::Int64, p::Float64)
    n = n-1
    f = factorial
    a = f(big(n))
    support = [(a / (f(big(k))*f(big(n-k)))) * (p^k)*((1-p)^(n-k)) for k in 0:n]
    return Categorical(support)
end

function exponentialToCategorical(λ::Float64, n::Int64)
    n -= 1
    support = Vector{Float64}(undef,0)

    for k in 0:n
        push!(support, exp(-λ*k) - exp(-λ*(k+1)))
    end

    sum = support'*[1 for i in 1:length(support)]
    push!(support, 1-sum)

    support = readjustSupportVector(support, length(support))
    deleteat!(support, length(support))

    return Categorical(support)
end

function gen_instance_uniform_poissonBatchValue(n::Int64, Lmax)
    n_outputs = n
    d_batchNumbers = uniformToCategorical(n)
    d_batchPositions = uniformToCategorical(n)
    d_batchValue = poissonToCategorical(1., Lmax)

    return generateHardInstance(Lmax, n_outputs, d_batchNumbers, d_batchPositions, d_batchValue)
end

function gen_instance_uniform_exponBatchValue(n::Int64, Lmax)
    n_outputs = n
    d_batchNumbers = uniformToCategorical(n)
    d_batchPositions = uniformToCategorical(n)
    d_batchValue = exponentialToCategorical(.5, Lmax)

    return generateHardInstance(Lmax, n_outputs, d_batchNumbers, d_batchPositions, d_batchValue)
end

function generateBins(numberOfBins::Int64, numberOfOutputs, Lmax::Int64, LmaxLastBin::Int64)
    bins = []

    # Bins
    for i in 1:numberOfBins-1
        push!(bins, gen_instance_uniform_exponBatchValue(numberOfOutputs, Lmax))
    end

    # Last bin with LmaxLastBin
    push!(bins, gen_instance_uniform_exponBatchValue(numberOfOutputs, LmaxLastBin))

    return bins
end


"""
    generateInstances(path, numberOfInstances, numberOfBins, numberOfOutputs, Lmax, LmaxLastBin)

**Creates a set of instances randomly generated to have a non trivial optimal solution.**

# Parameters :
- path::String : The directory on which to create the text files of the instances
- numberOfInstances::Int64 : The number of instances to create with the given parameters
- numberOfOutputs::Int64 : The number of outputs for each instance
- Lmax::Int64 : The highest value each output can tolerate
- LmaxLastBin::Int64 : The highest value each output of the last bin can tolerate
"""
function generateInstances(path::String, numberOfInstances::Int64, numberOfBins::Int64, numberOfOutputs::Int64, Lmax::Int64, LmaxLastBin::Int64)
    for instanceNumber in 1:numberOfInstances
        bins = generateBins(numberOfBins, numberOfOutputs, Lmax, LmaxLastBin)
        numberOfRoutes = sum([length(bin) for bin in bins])
        filename = "instance_" * string(instanceNumber) * "_" * string(numberOfOutputs) * "_" * string(numberOfRoutes) * "_opt_" * string(numberOfBins) * "_" * string(Lmax) * "_" * string(LmaxLastBin) * ".txt"
        
        # Create file
        println("instance generated : "*filename)
        file = open(path*"/"*filename, "w")
        
        # Instance infos
        write(file, string(numberOfOutputs)*"\n")
        write(file, string(numberOfRoutes)*"\n")
        write(file, string(numberOfBins)*"\n")
        write(file, string(Lmax)*"\n")
        write(file, string(LmaxLastBin)*"\n")
        write(file, "\n")
        
        for bin in bins
            for round in bin
                for value in round
                    write(file, string(value)*",")
                end
                write(file, "\n")
            end
            write(file, "\n")
        end

        close(file)
    end
end

# generateInstances("TER/data/hard/test", 10, 10, 200, 50, 50)
