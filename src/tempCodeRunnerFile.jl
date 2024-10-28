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