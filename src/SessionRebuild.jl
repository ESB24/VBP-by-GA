begin
    using Gurobi
    using JuMP
end

begin
    using BenchmarkTools
    using ProfileCanvas
end

begin
    include("DataStruct.jl")
    include("MyDistribution.jl")
end

function rebuildSession(s::Session)
    O::Int64 = length(s.loads)
    R::Int64 = length(s.rounds)
    C::Int64 = s.C

    s = Session(s.C, [Round(r.id, zeros(Int64, O), r.batches) for r in s.rounds], zeros(Int64, O))

    flag    ::Bool          = true
    rounds  ::Vector{Tuple} = [(Ref(1), r) for r in s.rounds]

    for k=1:O
        sort!(rounds, by= x -> (x[1][] > length(x[2].batches)) ? (1, 1, 1) : (
                                ((length(x[2].batches) - x[1][] >= O - k) ? 0 : 1),     
                                -(x[2].batches[x[1][]]), 
                                -(length(x[2].batches) - x[1][])
                                ))

        for (j, r) in rounds
            if (j[] <= length(r.batches)) && (r.batches[j[]] + s.loads[k] <= C)
                s.loads[k] += r.batches[j[]]        # update loads
                r.assignment[k] = r.batches[j[]]    # Affect batch
                j[] += 1                            # Pass to next batch
            end
        end
    end

    for (j, r) in rounds
        if j[] <= length(r.batches)
            # println("r -> $(r.id)")
            flag = false
            # batchesLeft = r.batches[j[]:end] # batches left to assign

            # toRemove = Set(findall(x -> x==0, r.assignment)[end-length(batchesLeft)+1:end]) # Last n zeros of the assignment vector with n the number of batches left to assign

            # r.assignment = [b for (k, b) in enumerate(r.assignment) if !(k in toRemove)] # assignment vector without the last n zeros
            # append!(r.assignment, batchesLeft) # add the left out batches back
        end
    end
    
    flag || compute_output!(s)

    return s, flag
end

function rebuildSession_knapSack_model(s::Session, env::Gurobi.Env = Gurobi.Env())
    O::Int64 = length(s.loads)
    R::Int64 = length(s.rounds)
    C::Int64 = s.C

    s = Session(s.C, [Round(r.id, zeros(Int64, O), r.batches) for r in s.rounds], zeros(Int64, O))

    flag    ::Bool          = true
    rounds  ::Vector{Tuple} = [(Ref(1), r) for r in s.rounds]

    for k=1:O
        b   ::Vector{Int64} = [(j[] <= length(r.batches)) ? r.batches[j[]] : 0 for (j, r) in rounds]
        p   ::Vector{Int64} = [length(r.batches) - j[] for (j, r) in rounds]
        tl  ::Int64         = 10

        affected, _ = model_O1KP(b, p, C, O, tl, env)
        (affected === nothing) && (flag = false; break)

        for (i, (v, (j, r))) in enumerate(zip(affected, rounds)) 
            if v == 1 && j[] <= length(r.batches)
                batch = r.batches[j[]]
                rounds[i][2].assignment[k] = batch
                rounds[i][1][] += 1
                s.loads[k] += batch
            end
        end

        # println("$k, $(s.loads[k]), $(affected)")
    end

    i=1
    while flag && (i <= length(rounds))
        (rounds[i][1][] <= length(rounds[i][2].batches)) && (flag = false)
        i+=1
    end
    
    # flag || compute_output!(s)
    return s, flag
end

function rebuildSession_knapSack(s::Session)
    O::Int64 = length(s.loads)
    R::Int64 = length(s.rounds)
    C::Int64 = s.C

    s = Session(s.C, [Round(r.id, zeros(Int64, O), r.batches) for r in s.rounds], zeros(Int64, O))

    flag    ::Bool          = true
    rounds  ::Vector{Tuple} = [(Ref(1), r) for r in s.rounds]

    for k=1:O
        println("$k")
        b   ::Vector{Int64} = [(j[] <= length(r.batches)) ? r.batches[j[]] : 0 for (j, r) in rounds]
        p   ::Vector{Int64} = [length(r.batches) - j[] for (j, r) in rounds]
        tl  ::Int64         = 10

        _, affected = knapSack(C, b)

        for v in affected 
            j = rounds[v][1]
            r = rounds[v][2]

            if j[] <= length(r.batches)
                batch = r.batches[j[]]
                r.assignment[k] = batch
                j[] += 1
                s.loads[k] += batch
            end
        end
    end

    i=1
    while flag && (i <= length(rounds))
        (rounds[i][1][] <= length(rounds[i][2].batches)) && (flag = false)
        i+=1
    end
    
    # flag || compute_output!(s)
    return s, flag
end

function model_O1KP(b::Vector{Int64}, p::Vector{Int64}, C::Int64, O::Int64, tl::Int64=20, env::Gurobi.Env = Gurobi.Env())
    R::Int64 = length(b)
    # print("*")

    model = Model(() -> Gurobi.Optimizer(env))
    set_silent(model)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", tl)

    # 1 if the j-th batch of round r is assigned to the output k ∈ O(r)j, 0 otherwise
    @variable(model, x[1:R], Bin)

    @objective(model, Max, sum([x[r] * (b[r] + (p[r] / O)) for r=1:R]))

    @constraint(model, sum([x[r] * b[r] for r=1:R]) <= C)

    set_start_value(x[argmax(b)], 1)

    optimize!(model)

    if MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0
        return value.(x), model
    else
        return nothing, model
    end
end

function rebuildSession_knapSack_model_V2(s::Session, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env(), δ::Vector{Float64} = [1., 1., 1., 1., 1.])
# ==========< Init >==========
    O       ::Int64         = length(s.loads)
    Rs      ::Int64         = length(s.rounds)
    Lmax    ::Int64         = s.C

    s       ::Session       = Session(s.C, [Round(r.id, zeros(Int64, O), r.batches) for r in s.rounds], zeros(Int64, O))

    minCapa ::Int64         = sum([sum(r.batches) for r in s.rounds])
    maxCapa ::Int64         = O * Lmax
    flag    ::Bool          = true
    
    J       ::Vector{Int64} = ones(Int64, Rs)
    B       ::Vector{Int64} = [length(r.batches) for r in s.rounds]

    for k=1:O
# ==========< Forced round >==========
        V   ::Vector{Int64} = [(j <= length(r.batches)) ? r.batches[j] : 0 for (j, r) in zip(J, s.rounds)]
        F   ::Vector{Int64} = [r for r=1:Rs if ((J[r] <= B[r]) && (B[r] - J[r] >= O - k))]
        R   ::Vector{Int64} = [r for r=1:Rs if ((J[r] <= B[r]) && (B[r] - J[r] < O - k))]

        for r in F
            s.rounds[r].assignment[k] = V[r]
            s.loads[k] += V[r]
            J[r] += 1
        end
        
        C   ::Int64         = Lmax - s.loads[k]
        (C < 0) && (flag = false; break)
        ((isempty(R)) || (C == 0)) && (continue)

# ==========< Lauch Model >==========

        affected, _ = model_O1KP_V2(R, B, V, J, C, O, Rs, tl, env, δ)   
        (affected === nothing) && (flag = false; break)

# ==========< Update Session >==========
        for (r, v) in enumerate(affected) 
            if v == 1 && J[r] <= B[r]
                s.rounds[r].assignment[k] = V[r]
                s.loads[k] += V[r]
                J[r] += 1
            end
        end
        
        if s.loads[k] < C
            maxCapa -= C - s.loads[k]
            (maxCapa < minCapa) && (flag = false; break)
        end
    end

    r=1
    while flag && (r <= Rs)
        (J[r] <= B[r]) && (flag = false)
        r+=1
    end
    
    return s, flag
end

function model_O1KP_V2(R, B, V, J, C, O, Rs, tl, env::Gurobi.Env = Gurobi.Env(), δ::Vector{Float64} = [1., 1., 1., 1., 1.])
    m = maximum(V)
    b = maximum([B[r] - J[r] > 0 ? B[r] - J[r] : 1 for r=1:length(R)])

    # println("========\nR = $R\nB = $B\nV = $V\nJ = $J\nC = $C\nO = $O\n")

    model = Model(() -> Gurobi.Optimizer(env))
    set_silent(model)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", tl)

    @variable(model, x[r in R], Bin)

    @objective(model, Max, sum([x[r] * (V[r] + (δ[1]/(length(δ) * length(R))) * ((B[r]-J[r]) / (b+1)) + (δ[2]/(length(δ) * length(R))) * (V[r]/(m+1))^2) for r in R])) # + (δ[1]/length(δ)) * ((B[r]-J[r]) / (b+1)) #  + (δ[2]/length(δ)) * (V[r]/(m+1))^2 #  + (δ[3]/length(δ)) * (((m - V[r])/m)^2) # + (δ[4]/length(δ)) * (sum(x)/length(R))

    @constraint(model, sum([x[r] * V[r] for r in R]) <= C)

    optimize!(model)

    if MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0
        res = zeros(Int64, Rs)
        for r in R
            res[r] = value(x[r])
        end
        return res, model
    else
        return nothing, model
    end
end

function knapSack(C::Int64, b::Vector{Int64}, n::Int64=length(b))::Tuple{Int64, Vector{Int64}}
    (n==0 || C==0 || isempty(b)) && (return (0, []))

    if b[n] > C || b[n] == 0
        print("-")
        return knapSack(C, b, n-1)
    else
        print("<")
        S1 = knapSack(C - b[n], b, n-1)
        V1 = S1[1] + b[n]

        S2 = knapSack(C, b, n-1)
        V2 = S2[1]

        (V1 < V2) && (return S2)
        (V1 > V2) && (return (S1[1] + b[n], push!(S1[2], n)))

        maxS1 = isempty(S1[2]) ? 0 : maximum([b[e] for e in S1[2]])
        maxS2 = isempty(S2[2]) ? 0 : maximum([b[e] for e in S2[2]])
        
        (maxS1 < maxS2) && (return S2)
        (maxS1 > maxS2) && (return (S1[1] + b[n], push!(S1[2], n)))

        (length(S1[2]) < length(S2[2])) && (return S2)
        # (length(S1[2]) > length(S2[2])) && (return (S1[1] + b[n], push!(S1[2], n)))
        return (S1[1] + b[n], push!(S1[2], n))
    end
end

function knapSack_all(C::Int64, b::Vector{Int64}, n::Int64=length(b))::Vector{Tuple{Int64, Vector{Int64}}}
    (n==0 || C==0) && (return [(0, [])])

    if b[n] > C
        return knapSack_all(C, b, n-1)
    else
        S1 = knapSack_all(C - b[n], b, n-1)
        V1 = S1[1][1] + b[n]

        S2 = knapSack_all(C, b, n-1)
        V2 = S2[1][1]

        if V1 < V2
            return S2
        elseif V1 > V2
            return [(o + b[n], push!(p, n)) for (o, p) in S1]
        else
            return [[(o + b[n], push!(p, n)) for (o, p) in S1];S2]
        end
    end
end

function knapSack_all_τ(C::Int64, b::Vector{Int64}, τ::Int64=0, n::Int64=length(b))::Vector{Tuple{Int64, Vector{Int64}}}
    (n==0 || C==0) && (return [(0, [])])

    if b[n] > C
        return knapSack_all_τ(C, b, τ, n-1)
    else
        S1 = knapSack_all_τ(C - b[n], b, τ, n-1)
        V1 = maximum([e[1] for e in S1]) + b[n]

        S2 = knapSack_all_τ(C, b, τ, n-1)
        V2 = S2[1][1]

        if V1+τ < V2
            return S2
        elseif V1 > V2+τ
            return [(o + b[n], push!(p, n)) for (o, p) in S1 if (o+b[n]+τ >= V1)]
        else
            return [[(o + b[n], push!(p, n)) for (o, p) in S1 if (o+b[n]+τ >= V1)];S2]
        end
    end
end

# env = Gurobi.Env()

# BenchmarkTools.@benchmark begin
#     R = 200
#     p = [rand(1:150) for _=1:R]
#     b = [rand(1:50) for _=1:R]

#     C = 300
#     O = 200

#     tl = 20

#     # res, md = model_O1KP(b, p, C, O, tl, env)
#     obj, res = knapSack(C, b)

#     # ((termination_status(md) != OPTIMAL) || (res === nothing)) && (println("not opti!"))

#     # print("$(sum(res .* b)), ")
# end

function rebuild_full_model(s::Session, tl::Int64 = 100, env::Gurobi.Env = Gurobi.Env())::Tuple{Model, Union{Nothing, Solution}}
# ===========< Parameters: >===========

    # is the set of all considered rounds
    R::Vector{Int64} = collect(1:length(s.rounds))

    # is the set of batches for the round r ∈ R
    Br::Vector{Vector{Int64}} = [collect(1:length(r.batches)) for r in s.rounds]

    # is the volume of the j-th batch in the round r ∈ R. Here, j ∈ B(r)
    vrj::Vector{Vector{Int64}} = [[b for b in r.batches] for r in s.rounds]

    # is the set of available machine outputs
    O::Vector{Int64} = collect(1:length(s.loads)) 

    # is the interval of potential outputs for the j-th batch in the round r ∈ R
    Orj::Vector{Vector{Vector{Int64}}} = [[collect(j: length(O) - length(Br[r]) + j) for j in Br[r]] for r in R]
    
    # is the set of mail bathes of round r, which can be potentially assigned to the output k ∈ O
    Urk::Vector{Vector{Vector{Int64}}} = [[[k for (k, v) in enumerate(Orj[r]) if o in v] for o in O] for r in R]
 
    # is the maximal load per output for any sorting session
    Lmax = s.C

# ===========< Model: >===========

    model = Model(() -> Gurobi.Optimizer(env))
    # set_silent(model)
    # set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", tl)
    set_optimizer_attribute(model, "MIPFocus", 1) # Focus on finding feasible solutions
    set_optimizer_attribute(model, "SolutionLimit", 1) # Stop after finding the first feasible solution

# ===========< Variable: >===========

    # 1 if the j-th batch of round r is assigned to the output k ∈ O(r)j, 0 otherwise
    @variable(model, x[r in R, j in Br[r], k in Orj[r][j]], Bin)

# ===========< Objective: TODO >===========

    # (1) -> inimizing the number of used session
    @objective(model, Min, maximum([sum([sum([vrj[r][j] * x[r, j, k] for j in Urk[r][k]]) for k in O]) for r in R]))

# ===========< Constraint: TODO >===========

    # (4) -> require that each batch j of round r be assigned to exactly oneoutput
    for r in R
        for j in Br[r]
            @constraint(model, sum([x[r, j, k] for k in Orj[r][j]]) == 1)
        end
    end

    # (7) -> precedence mail constraints for each round r
    for r in R
        for j in Br[r]
            if j != length(Br[r])
                @constraint(model, sum([k * x[r, j, k] for k in Orj[r][j]]) <= sum([k * x[r, j+1, k] for k in Orj[r][j + 1]]) - 1.0)
            end
        end
    end

    # (8) -> the total load for each output of each session, which can not be greater than Lmax if the session is considered as not empty
    for k in O
        @constraint(model, sum([sum([(vrj[r][j] * x[r, j, k]) for j in Urk[r][k]]) for r in R]) <= Lmax)
    end

    optimize!(model)

# ==========< Results: TODO >==========

    if MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0 || termination_status(model) == OPTIMAL
        res::Session = Session(perm, [Round(r.id, zeros(Int64, length(r.assignment)), r.batches) for r in s.rounds])
        for r in R
            for j in Br[r]
                for k in Orj[r][j]
                    (value(x[r, j, k]) == 1) && (res.rounds[r].StdAssignment[k] = res.rounds[r].batches[j])
                end
            end
        end
        return model, res
    else
        return model, nothing
    end
end



















function setup_O1KP_V3(R, B, V, J, C, O, tl, env::Gurobi.Env = Gurobi.Env(), δ::Vector{Float64} = [1., 1.])
    m = maximum(V)
    b = maximum([B[r] - J[r] > 0 ? B[r] - J[r] : 1 for r=1:length(R)])

    # println("========\nR = $R\nB = $B\nV = $V\nJ = $J\nC = $C\nO = $O\n")

    model = Model(() -> Gurobi.Optimizer(env))
    set_silent(model)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", tl)

    @variable(model, x[r in R], Bin)

    @objective(model, Max, sum([x[r] * (V[r] + (delta1/(length(δ) * length(R))) * ((B[r]-J[r]) / (b+1)) + (delta2/(length(δ) * length(R))) * (V[r]/(m+1))^2) for r in R])) # + (δ[3]/length(δ)) * (((m - V[r])/m)^2)  + (δ[4]/length(δ)) * (sum(x)/length(R))) #  

    @constraint(model, sum([x[r] * V[r] for r in R]) <= C)

    return model
end

function solve_O1KP_V3(model)
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        solution = value.(model[:x])
        return solution
    else
        return nothing
    end
end

function exclude_solution(model, solution, M, V)
    obj = sum([solution[r] * V[r] for r in M])

    @constraint(model, sum([(1 - solution[r]) * model[:x][r] for r in M]) + sum([solution[r] * (1 - model[:x][r]) for r in M]) >= 1)

    @constraint(model, sum([model[:x][r] * V[r] for r in M]) >= obj)

    # println("new cst -> obj >= $(obj)")
end

function Forced_mail_affectation!(s, k, F, V, J)
    Lmax = s.C
    
    ((isempty(F) ? 0 : sum([V[r] for r in F])) > Lmax) && (return false)

    for r in F
        s.rounds[r].assignment[k] = V[r]
        s.loads[k] += V[r]
        J[r] += 1
    end

    return true
end

function Model_mail_afffectation!(s, k, M, V, J, B, Lneed, model = nothing, tl = 10, env = Gurobi.Env(), δ = [0.3, 0.9, 0., 0., 0.])
    Lmax = s.C
    O = length(s.loads)
    L = s.loads[k]
    C = Lmax - L

    (C == 0) && (return (nothing, true)) # forced mail complitely filled the output print(" - C == 0(first) - "); 

    ((isempty(M)) || minimum([V[r] for r in M]) > C || (Lneed - sum(V) > Lmax * (O - k))) && (return (nothing, Lneed - L <= Lmax * (O - k))) # either no mail or to big to do a model print(" - empty(M)/min(M)>C - "); 

    (model === nothing) && (model = setup_O1KP_V3(M, B, V, J, C, O, tl, env, δ))

    solution = solve_O1KP_V3(model)

    (solution === nothing) && (return (nothing, Lneed - L <= Lmax * (O - k))) # print(" - no more sol - "); 

    for (r, a) in zip(M, solution)
        if a == 1
            s.rounds[r].assignment[k] = V[r]
            s.loads[k] += V[r]
            J[r] += 1
        end
    end

    L = s.loads[k]
    C = Lmax - L

    (C == 0) && (return (model, true)) # useless print(" - C == 0 (second) - "); 

    # print(" - C < 0, Lneed=$Lneed, L=$L, Lleft=$(Lmax * (O - k)) - ")
    return (model, Lneed - L <= Lmax * (O - k)) # C > 0
end

function rebuildSession_knapSack_model_V3!(s::Session, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env(), δ::Vector{Float64} = [0.3, 0.9, 0., 0., 0.])
    # print("<")

    O       ::Int64             = length(s.loads)
    Rs      ::Int64             = length(s.rounds)
    Lmax    ::Int64             = s.C

    s       ::Session           = Session(s.C, [Round(r.id, zeros(Int64, O), r.batches) for r in s.rounds], zeros(Int64, O))

    Lneed   ::Int64             = sum([sum(r.batches) for r in s.rounds])
    J       ::Vector{Int64}     = ones(Int64, Rs)
    B       ::Vector{Int64}     = [length(r.batches) for r in s.rounds]
    stack   ::Vector{Tuple}     = []

    reopt   ::Int64 = 0
    max_reopt::Int64 = ceil(Int64, 2O) # 2O

    last_k  ::Int64 = 0
    repeated_k ::Int64 = 0

    k::Int64 = 1
    while k <= O && reopt <= max_reopt
        if length(stack) < k
            V   ::Vector{Int64} = [(0 < j <= length(r.batches)) ? r.batches[j] : 0 for (j, r) in zip(J, s.rounds)]
            F   ::Vector{Int64} = [r for r=1:Rs if ((J[r] <= B[r]) && (B[r] - J[r] >= O - k))]
            M   ::Vector{Int64} = [r for r=1:Rs if ((J[r] <= B[r]) && (B[r] - J[r] < O - k))]
            model = nothing

            if Forced_mail_affectation!(s, k, F, V, J)
                model, valid = Model_mail_afffectation!(s, k, M, V, J, B, Lneed, nothing, tl, env, δ)
                if valid
                    # println("Success Model: (k=$k)")
                    push!(stack, (deepcopy(F), deepcopy(M), deepcopy(V), model))
                    Lneed -= s.loads[k]
                    (Lneed == 0) ? (k=O+1) : (k += 1)
                else
                    # println("Failed Model: (k=$k)")
                    Lneed -= s.loads[k]
                    k -= 1
                end
            else
                # println("Failed Forced: (k=$k)")
                k -= 1
            end
        else
            isempty(stack) && (return s, false)

            reopt += 1

            # print("<k=$(k), repeat=$(repeated_k)>")

            (last_k != k) ? (last_k = k; repeated_k = 0) : (repeated_k += 1) # || last_k != k+1 || last_k != k-1


            rollback = ((repeated_k >= 5 && k >= 5) ? rand(1:min(ceil(Int64, O/20), length(stack)-1)) : (0)) # 

            # println("rollback = $rollback")

            # println("    -> rollback=$rollback - ")
            for r=1:Rs
                if s.rounds[r].assignment[k+1] != 0
                    s.rounds[r].assignment[k+1] = 0
                    J[r] -= 1
                end
            end
            Lneed += s.loads[k+1]
            s.loads[k+1] = 0
            F, M, V, model = ([], [], [], nothing)

            for i=0:rollback
                (i != 0) && (k -= 1)

                F, M, V, model = pop!(stack)

                for r=1:Rs
                    if s.rounds[r].assignment[k] != 0
                        s.rounds[r].assignment[k] = 0
                        J[r] -= 1
                    end
                end
                Lneed += s.loads[k]
                s.loads[k] = 0

                # print(" -> k=$k -> ")
                # for r in s.rounds
                #     print("$(r.assignment[k]), ")
                # end
                # print(":$(s.loads[k])\n - J -> ")
                # for r=1:Rs
                #     print("$(J[r]), ")
                # end
                # println()
            end

            if model === nothing
                # println("    ->  No Model: (k=$k)")
                k -= 1
            else
                # print("    ->  Model: (k=$k)")
                solution = value.(model[:x])

                Forced_mail_affectation!(s, k, F, V, J)

                exclude_solution(model, solution, M, V)

                model, valid = Model_mail_afffectation!(s, k, M, V, J, B, Lneed, model, tl, env, δ)

                if valid
                    Lneed -= s.loads[k]

                    (Lneed == 0) ? (k=O+1) : (k += 1)
                end
                # println("Re-opt $(valid ? "Success" : "Failed")$(model === nothing ? " ! lost model !" : "")")

                push!(stack, (deepcopy(F), deepcopy(M), deepcopy(V), model))
            end
        end
    end

    # println("==========")
    # println("    re-opt:$(reopt)")
    # println("==========")

    # print(">")
    
    return s, k == O+1
end


