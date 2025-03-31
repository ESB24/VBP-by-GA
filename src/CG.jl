begin
    using Random
end

begin
    include("DataStruct.jl")
    include("Solution.jl")
end

## ============================================================================================================== ##
##        ######    ######   ##        ##    ##  ##    ##  ##    ##             ######   ########  ##    ##       ##
##       ##    ##  ##    ##  ##        ##    ##  ###  ###  ###   ##            ##        ##        ###   ##       ##
##       ##        ##    ##  ##        ##    ##  ## ## ##  ## ## ##            ##  ###   #####     ## ## ##       ##
##       ##    ##  ##    ##  ##        ##    ##  ##    ##  ##   ###            ##    ##  ##        ##   ###       ##
##        ######    ######   ########   ######   ##    ##  ##    ##             ######   ########  ##    ##       ##
## ============================================================================================================== ##

function BF_1D(
        inst    ::tInstance                                 ; # Instance 
        lb      ::Union{Int64, Nothing}         = nothing   , # Lower bound
        order   ::Union{Vector{Int64}, Nothing} = nothing   , # Order to consider items
    )

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

function BF_1D(I::Vector{Int64}, C::Int64)
    return BF_1D(tInstance(enumerate(I), C))
end

function BFD_1D(inst::tInstance)
    return BF_1D(inst, order=(sortperm(inst.I, by=(x -> x.size))))
end

function BFD_1D(I::Vector{Int64}, C::Int64)
    return BFD_1D(tInstance(enumerate(I), C))
end

function FF_1D(
        inst    ::tInstance                                 ; # Instance 
        lb      ::Union{Int64, Nothing}         = nothing   , # Lower bound
        order   ::Union{Vector{Int64}, Nothing} = nothing   , # Order to consider items
    )

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

function FFD_1D(inst::tInstance)
    return FF_1D(inst, order=(sortperm(inst.I, by=(x -> x.size))))
end

function randBF()
    I::Vector{tItem} = [tItem(i, rand(1:10)) for i=1:200]
    C::Int64 = 100
    inst::tInstance = tInstance(I, C)

    BF_1D(inst)
end

function randFF()
    I::Vector{tItem} = [tItem(i, rand(1:10)) for i=1:200]
    C::Int64 = 100
    inst::tInstance = tInstance(I, C)

    FF_1D(inst)
end

function heuristic_fill(
        b       ::tBin          , 
        inst    ::tInstance     ,
    )::tBin

    # Update instance with item not in b
    inst = tInstance([i for i in inst.I if i ∉ b.I], inst.C)

    # random order to consider items
    order = randperm(length(inst.I))

    # smallest item not in b
    smallest_item_volume = minimum(x -> x.size, inst.I)

    # current item considered (order id)
    current = 1

    println(inst)
    println(order)

    # while there is still items smaller than the remaining capacity 
    while smallest_item_volume < b.C - b.L && current <= length(order)

        # get current item
        i = inst.I[order[current]]

        # if size fit insert item
        if i.size + b.L ≤ b.C
            push!(b.I, i)
            b.L += i.size
        end

        # move to next item
        current += 1
    end

    return b
end

# function test_heuristic_fill()
#     I = [tItem(i, rand(1:10)) for i=1:10]
#     b = tBin([I[1], I[2]], 18, 30)

#     println("b -> $b")
#     println("I -> $I")

#     b = heuristic_fill(b, tInstance(I, 30))

#     println("updated b -> $b")
# end


## ============================================================================================================== ##
##        ######    ######   ##        ##    ##  ##    ##  ##    ##             ######   #######   ########       ##
##       ##    ##  ##    ##  ##        ##    ##  ###  ###  ###   ##            ##    ##  ##    ##     ##          ##
##       ##        ##    ##  ##        ##    ##  ## ## ##  ## ## ##            ##    ##  #######      ##          ##
##       ##    ##  ##    ##  ##        ##    ##  ##    ##  ##   ###            ##    ##  ##           ##          ##
##        ######    ######   ########   ######   ##    ##  ##    ##             ######   ##           ##          ##
## ============================================================================================================== ##

function H(
        b       ::tBin              ,
        inst    ::tInstance         ,
        a       ::Vector{Int64}     ,
        p       ::Int64  , 
        α       ::Int64             , 
        δ       ::Int64             ,
    )

    return sum([a[p+j * δ] * 2^(α - j) + sum([a[k] * ceil(Int64, (k/2)+1)^2 for k=p+j*δ+1:p+(j+1)*δ-1]) for j=0:α-1])
end

function v(
        b       ::tBin                                      ,
        inst    ::tInstance                                 ;
        α       ::Union{Int64, Nothing}            = nothing, 
        δ       ::Union{Int64, Nothing}            = nothing,
    )

    n = length(inst.I)

    a = [(i ∈ b.I) ? (1) : (0) for i in inst.I]

    (α == nothing) && (α = 20)
    (δ == nothing) && (δ = floor(Int64, n/α))

    return sum([H(b, inst, a, p, α, δ) for p=1:4]) / 3
end

function hashing(b::tBin, columns::Vector{Tuple{Float64, Vector{tBin}}})

    v_b = v(b)

    𝓛

    
    
end

function create_column(𝓛v)

function CG(
        instance
        tl
    )

    𝓛 = []

    LB
    UB

    start_time = time()

    while  
    
end

