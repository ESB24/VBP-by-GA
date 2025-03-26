function generateRounds(
        nbRound         ::Int64             ,
        minVolume       ::Int64             ,
        minBatches      ::Int64             ,
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

    matrix[2:(nbPlots+1), :] = round.(plotRange .* tmp) # scale and round each plot
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

# ==========< Assign Batches >==========

    roundBatches::Vector{Vector{Int64}} = [[] for _=1:nbRound]

    for k=1:O
        order::Vector{Int64} = randperm(nbRound+1)
        for i=1:nbPlot+1
            push!(roundBatches[order[1]], k)
        end
    end

# ==========< Build Rounds >==========
    rounds::Vector{Round} = []    

    for batchesId::Vector{Int64} in roundBatches
        Br::Int64 = length(batchesId)
        r::Round{Br} = Round(currentId, zeros(Int64, O), ntuple(i -> batches[batchesId[i]], Br))
        
        for i=1:Br
            r.assignment[round(Int64, 1+floor((batchesId[i]-1) / (nbPlots+1)))] = batches[batchesId[i]]
        end

        push!(rounds, r)
        currentId += 1
    end

    return (rounds, currentId)
end

begin # Library
    # Test
	import BenchmarkTools
	import ProfileCanvas

    # Miscelaneous
    using Random
end
    
begin # Files
    include("parsing.jl")

    include("Route.jl")
    include("Session.jl")
    include("Solution.jl")

    include("Instance.jl")
    include("MILP.jl")
    include("Miscelaneous.jl")
    
    # include("OptiMove.jl")
    
    include("GeneticAlgorithm.jl")

    include("SessionRebuild.jl")
end

# ================================================================================= #
#  #######  ##    #   ######  #######   #####   ##    #   ######  #######   ######  #
#     #     # #   #  #           #     #     #  # #   #  #        #        #        #
#     #     #  #  #   #####      #     #######  #  #  #  #        ####      #####   #
#     #     #   # #        #     #     #     #  #   # #  #        #              #  #
#  #######  #    ##  ######      #     #     #  #    ##   ######  #######  ######   #
# ================================================================================= #

function test()
    instance, nbSession = parseAnyInstance("instanceSkewed_1_O200_R20_C60_opt_1.txt")

    R = instance.route
    B = [length(r.mail) for r in R]
    J = [rand(1:B[r]) for r=1:length(R)]
    V = [r.mail[J[k]] for (k, r) in enumerate(R)]

    m = maximum(V)
    b = maximum([B[r] - J[r] for r=1:length(R)])

    w1 = zeros(Float64, length(R))
    w2 = zeros(Float64, length(R))

    for r=1:length(R)
        obj1_1 = (.3/(2 * length(R))) * ((B[r]-J[r]) / (b+1))
        obj1_2 = (.9/(2 * length(R))) * (V[r]/(m+1))^2

        w1[r] += V[r] + obj1_1 + obj1_2
        
        obj2_1 = (.61) * ((B[r]-J[r]) / (1+sum(B-J)))
        obj2_2 = (.39) * ((V[r]^2) / (1+sum(V.^2)))

        w2[r] += V[r] + obj2_1 + obj2_2
    end

    (w1 == w2) ? print("o") : print("-")
end




function rebuild1Session(inst_name::String = "instanceSkewed_1_O200_R20_C60_opt_1.txt")
    tmp = VERBOSE
    global VERBOSE = true

    instance, nbSession = parseAnyInstance(inst_name)

    glob_s = Session(instance.Lmax, instance.route)
    glob_tl = 10
    glob_env = Gurobi.Env()

    s, tag = rebuildSession_knapSack_model_V4!(glob_s, glob_tl, glob_env)
    println_verbose("$s")
    println_verbose("VALID = $tag", ANSI_cyan)
    global VERBOSE = tmp
    return s, tag
end

# begin
#     cpt = 0
#     nb_test = 10
#     total_time = 0.
#     for i=1:nb_test
#         start = time()  
#         _, tag = rebuild1Session("instanceSkewed_$(i)_O200_R20_C60_opt_1.txt")
#         total_time +=  time() - start 
#         sleep(5) 
#         tag && (cpt += 1)
#     end
#     global VERBOSE = true
#     println_verbose("optimal = $cpt/$nb_test took in average $(round(total_time/nb_test, digits=3))", ANSI_cyan)
#     global VERBOSE = false
# end

function EM_GR(instance_name::String = "instanceSkewed_1_O200_R20_C60_opt_1.txt")
    instance, nbSession = parseAnyInstance(instance_name) # trafic_200_200_35_1.xlsx
    Lmax = instance.Lmax
    O = instance.nbOut
    R = instance.nbRoute

    start = time()
    sol_BFD_EmptyMove       = BFD_EmptyMove(instance.route, Lmax, O, R)
    time_BFD_EmptyMove      = round(time() - start, digits=3)

    start = time()
    sol_BFD_GreadyRebuild   = BFD_GreadyRebuild(instance.route, Lmax, O, R)
    time_BFD_GreadyRebuild  = round(time() - start, digits=3)

    println("$(ANSI_cyan) ==========< Solution BFD EMPTY-MOVE     ($(length(sol_BFD_EmptyMove)) Session(s) in $(time_BFD_EmptyMove)s, valid?= $(isSolutionValid(instance, sol_BFD_EmptyMove, false))) >========== $(ANSI_reset)")
    println("$sol_BFD_EmptyMove\n\n")

    println("$(ANSI_cyan) ==========< Solution BFD GREADY REBUILD ($(length(sol_BFD_GreadyRebuild)) Session(s) in $(time_BFD_GreadyRebuild)s, valid?= $(isSolutionValid(instance, sol_BFD_GreadyRebuild, false))) >========== $(ANSI_reset)")
    println("$sol_BFD_GreadyRebuild\n\n")
end

function getStdInstance(inst_name::String = "instanceChunk_1_O200_R20_C100_opt_1.txt")
    i, _ = parseAnyInstance(inst_name)

    all_mail = []
    for r in i.route
        all_mail = [all_mail; collect(r.mail)]
    end

    return std(all_mail)
end
