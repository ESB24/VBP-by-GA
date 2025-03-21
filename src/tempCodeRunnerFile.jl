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

# tmp_O::Int64 = 10
# tmp_Lmax::Int64 = 10
# routes = [Route(1, tmp_O, (1, tmp_Lmax, 1)) ;[Route(i, tmp_O, (tmp_Lmax, 1)) for i = 2:(tmp_O-1)]]

# begin
#     glob_tl::Int64 = 10
#     glob_env::Gurobi.Env = Gurobi.Env()

#     all_δ = zeros(Int64, 8, 5)

#     for i=1:100
#         print("<i = $i> ")
#         for (p1, δ1) in enumerate([.19])
#             for (p2, δ2) in enumerate([.15])
#                 Lmax, mat, nbSession = parseMyInstance("../data/Chunk/instanceChunk_$(i)_O20_R20_C100_opt_1.txt")
#                 instance::Instance = Instance(mat, 0.25)
#                 instance.Lmax = Lmax

#                 glob_s = Session(Lmax, instance.route)

#                 _, glob_res = rebuildSession_knapSack_model_V4!(glob_s, glob_tl, glob_env, [δ1, δ2])
#                 print("$(glob_res ? "-" : "x")")
#                 glob_res && (all_δ[p1, p2] += 1)
#             end
#         end
#         println()
#     end

#     print("δ -> $all_δ, ∑ -> $(maximum(all_δ))")
# end

begin
    Lmax, mat, nbSession = parseMyInstance("../data/Chunk/instanceChunk_1_O200_R20_C100_opt_1.txt")
    instance::Instance = Instance(mat, 0.25)
    instance.Lmax = Lmax

    glob_s = Session(Lmax, instance.route)

    rebuildSession_knapSack_model_V4!(glob_s, glob_tl, glob_env, [10., .1, 1.])
end