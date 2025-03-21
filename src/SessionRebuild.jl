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
    include("GlobalDisplay.jl")
end

function rebuildSession(s::Session)
    O::Int64 = length(s.load)
    R::Int64 = length(s.route)
    C::Int64 = s.Lmax

    s = Session(s.Lmax, [Route(r.id, zeros(Int64, O), r.mail) for r in s.route], zeros(Int64, O))

    flag    ::Bool          = true
    routes  ::Vector{Tuple} = [(Ref(1), r) for r in s.route]

    for k=1:O
        sort!(routes, by= x -> (x[1][] > length(x[2].mail)) ? (1, 1, 1) : (
                                ((length(x[2].mail) - x[1][] >= O - k) ? 0 : 1),     
                                -(x[2].mail[x[1][]]), 
                                -(length(x[2].mail) - x[1][])
                                ))

        for (j, r) in routes
            if (j[] <= length(r.mail)) && (r.mail[j[]] + s.load[k] <= C)
                s.load[k] += r.mail[j[]]        # update loads
                r.assignment[k] = r.mail[j[]]    # Affect batch
                j[] += 1                            # Pass to next batch
            end
        end
    end

    for (j, r) in routes
        if j[] <= length(r.mail)
            # println("r -> $(r.id)")
            flag = false
            # batchesLeft = r.mail[j[]:end] # batches left to assign

            # toRemove = Set(findall(x -> x==0, r.assignment)[end-length(batchesLeft)+1:end]) # Last n zeros of the assignment vector with n the number of batches left to assign

            # r.assignment = [b for (k, b) in enumerate(r.assignment) if !(k in toRemove)] # assignment vector without the last n zeros
            # append!(r.assignment, batchesLeft) # add the left out batches back
        end
    end
    
    flag || compute_output!(s)

    return s, flag
end

function rebuildSession_knapSack_model(s::Session, env::Gurobi.Env = Gurobi.Env())
    O::Int64 = length(s.load)
    R::Int64 = length(s.route)
    C::Int64 = s.Lmax

    s = Session(s.Lmax, [Route(r.id, zeros(Int64, O), r.mail) for r in s.route], zeros(Int64, O))

    flag    ::Bool          = true
    routes  ::Vector{Tuple} = [(Ref(1), r) for r in s.route]

    for k=1:O
        b   ::Vector{Int64} = [(j[] <= length(r.mail)) ? r.mail[j[]] : 0 for (j, r) in routes]
        p   ::Vector{Int64} = [length(r.mail) - j[] for (j, r) in routes]
        tl  ::Int64         = 10

        affected, _ = model_O1KP(b, p, C, O, tl, env)
        (affected === nothing) && (flag = false; break)

        for (i, (v, (j, r))) in enumerate(zip(affected, routes)) 
            if v == 1 && j[] <= length(r.mail)
                batch = r.mail[j[]]
                routes[i][2].assignment[k] = batch
                routes[i][1][] += 1
                s.load[k] += batch
            end
        end

        # println("$k, $(s.load[k]), $(affected)")
    end

    i=1
    while flag && (i <= length(routes))
        (routes[i][1][] <= length(routes[i][2].mail)) && (flag = false)
        i+=1
    end
    
    # flag || compute_output!(s)
    return s, flag
end

function rebuildSession_knapSack(s::Session)
    O::Int64 = length(s.load)
    R::Int64 = length(s.route)
    C::Int64 = s.Lmax

    s = Session(s.Lmax, [Route(r.id, zeros(Int64, O), r.mail) for r in s.route], zeros(Int64, O))

    flag    ::Bool          = true
    routes  ::Vector{Tuple} = [(Ref(1), r) for r in s.route]

    for k=1:O
        # println("$k")
        b   ::Vector{Int64} = [(j[] <= length(r.mail)) ? r.mail[j[]] : 0 for (j, r) in routes]
        p   ::Vector{Int64} = [length(r.mail) - j[] for (j, r) in routes]
        tl  ::Int64         = 10

        _, affected = knapSack(C, b)

        for v in affected 
            j = routes[v][1]
            r = routes[v][2]

            if j[] <= length(r.mail)
                batch = r.mail[j[]]
                r.assignment[k] = batch
                j[] += 1
                s.load[k] += batch
            end
        end
    end

    i=1
    while flag && (i <= length(routes))
        (routes[i][1][] <= length(routes[i][2].mail)) && (flag = false)
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

    # 1 if the j-th batch of route r is assigned to the output k ∈ O(r)j, 0 otherwise
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
    O       ::Int64         = length(s.load)
    Rs      ::Int64         = length(s.route)
    Lmax    ::Int64         = s.Lmax

    s       ::Session       = Session(s.Lmax, [Route(r.id, zeros(Int64, O), r.mail) for r in s.route], zeros(Int64, O))

    minCapa ::Int64         = sum([sum(r.mail) for r in s.route])
    maxCapa ::Int64         = O * Lmax
    flag    ::Bool          = true
    
    J       ::Vector{Int64} = ones(Int64, Rs)
    B       ::Vector{Int64} = [length(r.mail) for r in s.route]

    for k=1:O
# ==========< Forced route >==========
        V   ::Vector{Int64} = [(j <= length(r.mail)) ? r.mail[j] : 0 for (j, r) in zip(J, s.route)]
        F   ::Vector{Int64} = [r for r=1:Rs if ((J[r] <= B[r]) && (B[r] - J[r] >= O - k))]
        R   ::Vector{Int64} = [r for r=1:Rs if ((J[r] <= B[r]) && (B[r] - J[r] < O - k))]

        for r in F
            s.route[r].assignment[k] = V[r]
            s.load[k] += V[r]
            J[r] += 1
        end
        
        C   ::Int64         = Lmax - s.load[k]
        (C < 0) && (flag = false; break)
        ((isempty(R)) || (C == 0)) && (continue)

# ==========< Lauch Model >==========

        affected, _ = model_O1KP_V2(R, B, V, J, C, O, Rs, tl, env, δ)   
        (affected === nothing) && (flag = false; break)

# ==========< Update Session >==========
        for (r, v) in enumerate(affected) 
            if v == 1 && J[r] <= B[r]
                s.route[r].assignment[k] = V[r]
                s.load[k] += V[r]
                J[r] += 1
            end
        end
        
        if s.load[k] < C
            maxCapa -= C - s.load[k]
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
        # print("-")
        return knapSack(C, b, n-1)
    else
        # print("<")
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

    # is the set of all considered routes
    R::Vector{Int64} = collect(1:length(s.route))

    # is the set of batches for the route r ∈ R
    Br::Vector{Vector{Int64}} = [collect(1:length(r.mail)) for r in s.route]

    # is the volume of the j-th batch in the route r ∈ R. Here, j ∈ B(r)
    vrj::Vector{Vector{Int64}} = [[b for b in r.mail] for r in s.route]

    # is the set of available machine outputs
    O::Vector{Int64} = collect(1:length(s.load)) 

    # is the interval of potential outputs for the j-th batch in the route r ∈ R
    Orj::Vector{Vector{Vector{Int64}}} = [[collect(j: length(O) - length(Br[r]) + j) for j in Br[r]] for r in R]
    
    # is the set of mail bathes of route r, which can be potentially assigned to the output k ∈ O
    Urk::Vector{Vector{Vector{Int64}}} = [[[k for (k, v) in enumerate(Orj[r]) if o in v] for o in O] for r in R]
 
    # is the maximal load per output for any sorting session
    Lmax = s.Lmax

# ===========< Model: >===========

    model = Model(() -> Gurobi.Optimizer(env))
    # set_silent(model)
    # set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", tl)
    set_optimizer_attribute(model, "MIPFocus", 1) # Focus on finding feasible solutions
    set_optimizer_attribute(model, "SolutionLimit", 1) # Stop after finding the first feasible solution

# ===========< Variable: >===========

    # 1 if the j-th batch of route r is assigned to the output k ∈ O(r)j, 0 otherwise
    @variable(model, x[r in R, j in Br[r], k in Orj[r][j]], Bin)

# ===========< Objective: TODO >===========

    # (1) -> inimizing the number of used session
    @objective(model, Min, maximum([sum([sum([vrj[r][j] * x[r, j, k] for j in Urk[r][k]]) for k in O]) for r in R]))

# ===========< Constraint: TODO >===========

    # (4) -> require that each batch j of route r be assigned to exactly oneoutput
    for r in R
        for j in Br[r]
            @constraint(model, sum([x[r, j, k] for k in Orj[r][j]]) == 1)
        end
    end

    # (7) -> precedence mail constraints for each route r
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
        res::Session = Session(perm, [Route(r.id, zeros(Int64, length(r.assignment)), r.mail) for r in s.route])
        for r in R
            for j in Br[r]
                for k in Orj[r][j]
                    (value(x[r, j, k]) == 1) && (res.route[r].StdAssignment[k] = res.route[r].mail[j])
                end
            end
        end
        return model, res
    else
        return model, nothing
    end
end


## ============================================================================================================== ##
##                                               ##    ##   ######                                                ##
##                                               ##    ##        ##                                               ##
##                                       ####     ##  ##     #####     ####                                       ##
##                                                ##  ##         ##                                               ##
##                                                  ##      ######                                                ##
## ============================================================================================================== ##


function setup_O1KP_V3(R, B, V, J, C, O, tl, env::Gurobi.Env = Gurobi.Env(), δ::Vector{Float64} = [1., 1.])
    m = maximum(V)
    b = maximum([B[r] - J[r] > 0 ? B[r] - J[r] : 1 for r=1:length(R)])

    # println("========\nR = $R\nB = $B\nV = $V\nJ = $J\nC = $C\nO = $O\n")

    model = Model(() -> Gurobi.Optimizer(env))
    set_silent(model)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", tl)

    @variable(model, x[r in R], Bin)

    @objective(model, Max, sum([x[r] * (V[r] + (δ[1]/(length(δ) * length(R))) * ((B[r]-J[r]) / (b+1)) + (δ[2]/(length(δ) * length(R))) * (V[r]/(m+1))^2) for r in R])) 

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
    Lmax = s.Lmax
    
    if VERBOSE
        if isempty(F)
            print_verbose("<FA-PASS: F is empty> ", ANSI_yellow)
        else
            if sum([V[r] for r in F]) > Lmax
                print_verbose("<FA-FAIL: output overflow> ", ANSI_red)
            end
        end
    end


    ((isempty(F) ? 0 : sum([V[r] for r in F])) > Lmax) && (return false)

    for r in F
        s.route[r].assignment[k] = V[r]
        s.load[k] += V[r]
        J[r] += 1
    end

    isempty(F) || print_verbose("<FA-SUCC: F=$(join(F, ","))> ", ANSI_green)

    return true
end

function Model_mail_afffectation!(
        s::Session, 
        k::Int64, 
        M, 
        V::Vector{Int64}, 
        J::Vector{Int64}, 
        B::Vector{Int64}, 
        Lneed, 
        model = nothing, 
        tl = 10, 
        env = Gurobi.Env(), 
        δ = [0.3, 0.9, 0., 0., 0.]
    )
    Lmax = s.Lmax
    O = length(s.load)
    L = s.load[k]
    C = Lmax - L

    (C == 0) && (return (nothing, true)) # forced mail complitely filled the output print(" - C == 0(first) - "); 

    if VERBOSE
        if isempty(M) || isempty([r for r in M if V[r] <= C])
            if (Lneed - L <= Lmax * (O - k))
                print_verbose("<MA-SUCC: M is empty> ", ANSI_green)
            else
                print_verbose("<MA-FAIL: M is empty> ", ANSI_red)
            end
        else
            if (Lneed - L - sum([V[r] for r in M]) > Lmax * (O - k))
                print_verbose("<MA-FAIL: not enought volume> ", ANSI_red)
            end
        end
    end

    ((isempty(M)) || minimum([V[r] for r in M]) > C || (Lneed - sum(V) > Lmax * (O - k))) && (return (nothing, Lneed - L <= Lmax * (O - k))) # either no mail or to big to do a model print(" - empty(M)/min(M)>C - "); 

    (model === nothing) && (model = setup_O1KP_V3(M, B, V, J, C, O, tl, env, δ))

    solution = solve_O1KP_V3(model)

    (solution === nothing) && (return (nothing, Lneed - L <= Lmax * (O - k))) # print(" - no more sol - "); 

    for (r, x) in zip(M, solution)
        if x == 1
            s.route[r].assignment[k] = V[r]
            s.load[k] += V[r]
            J[r] += 1
        end
    end

    L = s.load[k]
    C = Lmax - L

    # (C == 0) && (return (model, true)) # useless print(" - C == 0 (second) - "); 

    if VERBOSE
        if (Lneed - L <= Lmax * (O - k))
            print_verbose("<MA-SUCC: x=$(join([r for (r, x) in zip(M, solution) if x == 1], ","))> ", ANSI_green)
        else
            print_verbose("<MA-FAIL: not enought inserted>", ANSI_red)
        end
    end

    # print(" - C < 0, Lneed=$Lneed, L=$L, Lleft=$(Lmax * (O - k)) - ")
    return (model, Lneed - L <= Lmax * (O - k)) # C > 0
end

function rebuildSession_knapSack_model_V3!(s::Session, tl::Int64 = 10, env::Gurobi.Env = Gurobi.Env(), δ::Vector{Float64} = [0.3, 0.9, 0., 0., 0.])
    # print("<")

    O       ::Int64             = length(s.load)
    Rs      ::Int64             = length(s.route)
    Lmax    ::Int64             = s.Lmax

    s       ::Session           = Session(s.Lmax, [Route(r.id, zeros(Int64, O), r.mail) for r in s.route], zeros(Int64, O))

    Lneed   ::Int64             = sum([sum(r.mail) for r in s.route])
    J       ::Vector{Int64}     = ones(Int64, Rs)
    B       ::Vector{Int64}     = [length(r.mail) for r in s.route]
    stack   ::Vector{Tuple}     = []

    reopt   ::Int64 = 0
    max_reopt::Int64 = ceil(Int64, 2O) # 2O

    last_k  ::Int64 = 0
    repeated_k ::Int64 = 0

    k::Int64 = 1
    while k <= O && reopt <= max_reopt

        print_verbose("<k = $k> ")

        if length(stack) < k
            V   ::Vector{Int64} = [(0 < j <= length(r.mail)) ? r.mail[j] : 0 for (j, r) in zip(J, s.route)]
            F   ::Vector{Int64} = [r for r=1:Rs if ((J[r] <= B[r]) && (B[r] - J[r] >= O - k))]
            tmp_C = Lmax - sum([V[r] for r in F])
            M   ::Vector{Int64} = [r for r=1:Rs if ((J[r] <= B[r]) && (B[r] - J[r] < O - k) && (V[r] <= tmp_C))]
            model = nothing

            # println_verbose("F = $(join(F, ",")), M = $(join(M, ","))")

            if Forced_mail_affectation!(s, k, F, V, J)
                model, valid = Model_mail_afffectation!(s, k, M, V, J, B, Lneed, nothing, tl, env, δ)
                if valid
                    # println("Success Model: (k=$k)")
                    push!(stack, (deepcopy(F), deepcopy(M), deepcopy(V), model))
                    Lneed -= s.load[k]
                    (Lneed == 0) ? (k=O+1) : (k += 1)
                else
                    # println("Failed Model: (k=$k)")
                    Lneed -= s.load[k]
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


            rollback = ((repeated_k >= 5 && k >= 5) ? rand(1:min(ceil(Int64, O/20), length(stack)-1)) : (0))

            # println("rollback = $rollback")

            # println("    -> rollback=$rollback - ")
            for r=1:Rs
                if s.route[r].assignment[k+1] != 0
                    s.route[r].assignment[k+1] = 0
                    J[r] -= 1
                end
            end
            Lneed += s.load[k+1]
            s.load[k+1] = 0
            F, M, V, model = ([], [], [], nothing)

            for i=0:rollback
                (i != 0) && (k -= 1)

                F, M, V, model = pop!(stack)

                for r=1:Rs
                    if s.route[r].assignment[k] != 0
                        s.route[r].assignment[k] = 0
                        J[r] -= 1
                    end
                end
                Lneed += s.load[k]
                s.load[k] = 0

                # print(" -> k=$k -> ")
                # for r in s.route
                #     print("$(r.assignment[k]), ")
                # end
                # print(":$(s.load[k])\n - J -> ")
                # for r=1:Rs
                #     print("$(J[r]), ")
                # end
                # println()
            end

            if model === nothing
                print_verbose("<No model> ", ANSI_magenta)
                k -= 1
            else
                print_verbose("<Re-optimisation> ", ANSI_magenta)
                solution = value.(model[:x])

                Forced_mail_affectation!(s, k, F, V, J)

                exclude_solution(model, solution, M, V)

                model, valid = Model_mail_afffectation!(s, k, M, V, J, B, Lneed, model, tl, env, δ)

                if valid
                    Lneed -= s.load[k]

                    (Lneed == 0) ? (k=O+1) : (k += 1)
                end
    
                push!(stack, (deepcopy(F), deepcopy(M), deepcopy(V), model))
            end
        end
        println_verbose()
        # println_verbose("$(s.load)")
        # println_verbose("$s")
    end

    return s, k == O+1
end



## ============================================================================================================== ##
##                                               ##    ##     ####                                                ##
##                                               ##    ##    ## ##                                                ##
##                                       ####     ##  ##    ##  ##     ####                                       ##
##                                                ##  ##   ########                                               ##
##                                                  ##          ##                                                ##
## ============================================================================================================== ##


global REBUILD_WEIGHT       = 2 # 1-2
global REBUILD_MODEL        = 1 # 1-2
global REBUILD_BACKTRACK    = 2 # 0-1-2

function setup_O1KP_V4(
        M       ::Vector{Int64}                     , 
        B       ::Vector{Int64}                     , 
        V       ::Vector{Int64}                     , 
        J       ::Vector{Int64}                     , 
        C       ::Int64                             ,
        Lneed   ::Int64                             ,
        Lmax    ::Int64                             ,
        O       ::Int64                             ,
        k       ::Int64                             ,
        τ       ::Int64                             , 
        δ       ::Vector{Float64}                   ,
        tl      ::Int64             = 10            , 
        env     ::Gurobi.Env        = Gurobi.Env()  , 
    )

    function compute_weight(
            M::Vector{Int64},  
            B::Vector{Int64}, 
            V::Vector{Int64}, 
            J::Vector{Int64},
            δ::Vector{Float64},
        )
        
        weights::Vector{Float64} = zeros(Float64, length(M))


        if REBUILD_WEIGHT == 1
            max_mail_left   ::Int64 = maximum([B[r] - J[r] for r in M])
            max_mail_volume ::Int64 = maximum([V[r] for r in M])

            # Main objective: Maximize mail volume
            for (i, r) in enumerate(M)
                weights[i] += V[r]
            end

            # Secondary objective (1/2): priority to route with a lot of mail left to assign
            for (i, r) in enumerate(M)
                weights[i] += (δ[1] / (2 * length(M))) * ((B[r] - J[r]) / (max_mail_left + 1))
            end

            # Secondary objective (2/2): priority to larger mail
            for (i, r) in enumerate(M)
                weights[i] += (δ[2] / (2 * length(M))) * (V[r] / (max_mail_volume+1))^2
            end

            sum_obj1 = round(sum([(δ[1] / (2 * length(M))) * ((B[r] - J[r]) / (max_mail_left + 1))  for (i, r) in enumerate(M)]), digits=10)
            sum_obj2 = round(sum([(δ[2] / (2 * length(M))) * (V[r] / (max_mail_volume+1))^2         for (i, r) in enumerate(M)]), digits=10)

            print_verbose("<obj1>", ANSI_green) #  = $(round(sum_obj1 - sum_obj2, digits=3))
        elseif REBUILD_WEIGHT == 2

            # Main objective: Maximize mail volume
            for (i, r) in enumerate(M)
                weights[i] += V[r]
            end

            # Secondary objective (1/2): priority to route with a lot of mail left to assign
            tmp_sec_obj1 = 1 + sum([B[r] - J[r] for r in M])
            for (i, r) in enumerate(M)
                weights[i] += δ[1] * (B[r] - J[r]) / tmp_sec_obj1
            end

            # Secondary objective (2/2): priority to larger mail
            tmp_sec_obj2 = 1 + sum([(V[r])^2 for r in M])
            for (i, r) in enumerate(M)
                weights[i] += δ[2] * (V[r])^2 / tmp_sec_obj2
            end

            print_verbose("<obj2>", ANSI_green)
        end

        return weights
    end
    
    # Compute variable weight for objective fonction
    w           ::Vector{Float64}   = compute_weight(M, B, V, J, δ)
    Required    ::Int64             = Lneed - ((O - k) * Lmax)

    # println("w = $w")

    model = Model(() -> Gurobi.Optimizer(env))
    set_silent(model)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", tl)

    if REBUILD_MODEL == 1
        # VAR (1/1): x ∈ {0, 1} equals 1 if first non assigned mail of route r is placed in output, 0 otherwise.
        @variable(model, x[r in M], Bin)

        # OBJ (1/1): Maximize ~ mail volume (∈ ℕ) mail selection (∈ [0, 1[) 
        @objective(model, Max, sum([x[r] * w[i] for (i, r) in enumerate(M)]))

        # CST (1/2): Capacity cst
        @constraint(model, sum([x[r] * V[r] for r in M]) ≤ C)

        # CST (2/2): Required mail volume cst
        if (0 < Required)
            @constraint(model, Required ≤ sum([x[r] * V[r] for r in M]))
        end

        print_verbose("<cst1>", ANSI_green)
    elseif REBUILD_MODEL == 2
        # VAR (1/2): x ∈ {0, 1} equals 1 if first non assigned mail of route r is placed in output, 0 otherwise.
        @variable(model, x[r in M], Bin)

        # VAR (2/2): P ∈ ℕ Penality for not reaching the average required volume
        @variable(model, 0 ≤ P ≤ C, Int)

        # OBJ (1/1): Maximize ~ mail volume (∈ [0, 1]), use big mail (∈ [0, 1]), mail with fullest route (∈ [0, 1]) -> weighted using δ
        @objective(model, Max, sum([x[r] * w[i] for (i, r) in enumerate(M)]) - P)

        # CST (1/3): Capacity cst
        @constraint(model, sum([x[r] * V[r] for r in M]) ≤ C)

        # CST (2/3): Required mail volume cst
        if (0 < Required)
            @constraint(model, Required ≤ sum([x[r] * V[r] for r in M]))
        end

        # CST (3/3): Penality cst
        if (Required < (Lmax - τ) < C)
            @constraint(model, (Lmax - (τ + sum([x[r] * V[r] for r in M]))) ≤ P)
        end

        print_verbose("<cst2>", ANSI_green)
    end

    return model
end

function solve_O1KP_V4(model::Model)::Union{Vector{Int64}, Nothing}
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        solution = value.(model[:x])
        return solution
    else
        return nothing
    end
end

function update_01KP_V4!(
        M       ::Vector{Int64},    # Id of the route considered by the model 
        V       ::Vector{Int64},    # Volume of the first non assigned mail of each route in session s
        model   ::Model        ,    # Model to update  
    )
    solution = value.(model[:x])

    @constraint(model, sum([(1 - solution[r]) * model[:x][r] for r in M]) + sum([solution[r] * (1 - model[:x][r]) for r in M]) >= 1)

    # TODO -> update model obj
    # TODO -> update min volume cst ?

    obj = sum([solution[r] * V[r] for r in M])
    @constraint(model, sum([model[:x][r] * V[r] for r in M]) >= obj)
end

function backtrack(
        k                   ::Int64         , # current output
        total_count_reopt   ::Int64         , # Total number of re-optimisation
        total_max_reopt     ::Int64         , # Maximum number of re-optimisation
        count_reopt         ::Vector{Int64} , # Count for each output how many time the model has been re-optimized
        max_reopt           ::Int64         , # Maximum reopt before further backtracking
        threshold_reopt     ::Int64         , # min output index before big backtrack step
        last_k              ::Int64         , # last output that required reoptimisation
        repeated_k          ::Int64         , # how many times last_k has been reoptimisation
        O                   ::Int64         , # Number of output
    )

    total_count_reopt += 1

    if REBUILD_BACKTRACK == 0 # Only backtrack exactly 1 output

        return (k, total_count_reopt, last_k, repeated_k, 0)

    elseif REBUILD_BACKTRACK == 1 # Backtrack more than 1 iff threshold is past and same output is optimized multiple time in a row

        (last_k != k) ? (last_k = k; repeated_k = 0) : (repeated_k += 1)

        rollback = ((repeated_k >= 5 && k >= 5) ? rand(1:min(ceil(Int64, O/20), k)) : (0))

        return (k, total_count_reopt, last_k, repeated_k, rollback)

    elseif REBUILD_BACKTRACK == 2 # Backtrack more than 1 iff threshold is past and same output has already been reoptimized multiple times

        count_reopt[k] += 1
        rollback = 0

        if (threshold_reopt < k) && (max_reopt < count_reopt[k])
            rollback +=  rand(1:min(ceil(Int64, O/20), k))
        end

        return (k, total_count_reopt, last_k, repeated_k, rollback)
    end
end

function Forced_mail_affectation_V4!(
        s::Session          , # Rebuild session 
        k::Int64            , # Current output
        F::Vector{Int64}    , # Route which will force a mail into the current output k
        V::Vector{Int64}    , # Volumes of the mails
        J::Vector{Int64}    , # Indexes of the first non assigned mail
        Lneed::Int64        ,
    )::Tuple{Bool, Int64}

    # Test: no forced affectations
    if isempty(F)
        print_verbose("<FA-PASS: F is empty> ", ANSI_yellow)
        return (true, Lneed)
    end

    # Test: forced affectation make the output overflow
    if sum([V[r] for r in F]) > s.Lmax
        print_verbose("<FA-FAIL: output overflow> ", ANSI_red)
        return (false, Lneed)
    end

    # Force first non assigned mail of each route in F to output k
    for r in F
        s.route[r].assignment[k] = V[r]     # Update mail assignemnent
        s.load[k] += V[r]                   # Update output load
        J[r] += 1                           # Update first non assigned mail
    end

    print_verbose("<FA-SUCC: F=$(join(F, ","))> ", ANSI_green)

    return (true, Lneed - s.load[k])
end

function Model_mail_affectation_V4!(
        s       ::Session                               , # Rebuild session 
        k       ::Int64                                 , # Current output
        O       ::Int64                                 , # Number of output 
        M       ::Vector{Int64}                         , # Set of route to be considered by the model
        V       ::Vector{Int64}                         , # Volumes of the mails
        J       ::Vector{Int64}                         , # Indexes of the first non assigned mail
        B       ::Vector{Int64}                         , # Number of mail in each route of session s
        Lneed   ::Int64                                 , # Required load for output k (to avoid failure)
        τ       ::Int64                                 , # Average tolerable hollowness for an output 
        δ       ::Vector{Float64}                       , # weight objective for the model
        model   ::Union{Model, Nothing} = nothing       , # Already updated model to use
        tl      ::Int64                 = 10            , # Time limite for the model
        env     ::Gurobi.Env            = Gurobi.Env()  , # Gurobi environement (display purpose)
    )::Tuple{Union{Model, Nothing}, Bool, Int64}
    
    Lmax    ::Int64 = s.Lmax            # Maximum capacity of output k
    C       ::Int64 = Lmax - s.load[k]      # Residual capacity left after forced affectation

    # Test: Model has no mail to insert
    if isempty(M)
        if VERBOSE
            if (Lneed <= Lmax * (O - k))
                print_verbose("<MA-SUCC: M is empty> ", ANSI_yellow)
            else
                print_verbose("<MA-FAIL: M is empty> ", ANSI_red)
            end
        end
        return (nothing, Lneed <= Lmax * (O - k), Lneed)
    end
    
    # Test: Model cannot insert enought mail to satisfy the volume constraints
    if (Lneed - sum([V[r] for r in M]) > Lmax * (O - k))
        print_verbose("<MA-FAIL: not enought volume> ", ANSI_red)
        return (nothing, false, Lneed)
    end

    (model == nothing) && (model = setup_O1KP_V4(M, B, V, J, C, Lneed, Lmax, O, k, τ, δ, tl, env))

    solution = solve_O1KP_V4(model)

    (solution === nothing) && (return (nothing, Lneed <= Lmax * (O - k), Lneed))

    L::Int64 = 0 # mail volume added by the model

    for (r, x) in zip(M, solution)
        if x == 1
            s.route[r].assignment[k] = V[r]
            s.load[k] += V[r]
            J[r] += 1

            L += V[r]
        end
    end
    
    if VERBOSE
        if (Lneed - L <= Lmax * (O - k))
            print_verbose("<MA-SUCC: x=$(join([r for (r, x) in zip(M, solution) if x == 1], ","))> ", ANSI_green)
        else
            print_verbose("<MA-FAIL: not enought inserted>", ANSI_red)
        end
    end

    return (model, Lneed - L <= Lmax * (O - k), Lneed - L)
end

function rebuildSession_knapSack_model_V4!(
        s   ::Session                           , # Rebuilded session
        tl  ::Int64             = 10            , # Time limit to execute the model
        env ::Gurobi.Env        = Gurobi.Env()  , # Gurobi model environement (display purpose)
        δ   ::Vector{Float64}   = [[.3, .9], [.2, .15]][REBUILD_MODEL] # Weights for model objective
    )

    # ====================< Miscelaneous >====================

    O       ::Int64         = length(s.load)                        # Number of output of the session
    R       ::Int64         = length(s.route)                       # Number of route to fit in the session
    Lmax    ::Int64         = s.Lmax                                # Maximum capacity of each output
    Lneed   ::Int64         = sum([sum(r.mail) for r in s.route])   # Volume of mail to fit in the session
    τ       ::Int64         = ceil(Int64, ((O * Lmax) - Lneed) / O) # Average tolerable hollowness for an output
    k       ::Int64         = 1                                     # Current output 
    J       ::Vector{Int64} = ones(Int64, R)                        # First not assigned mail of each route
    B       ::Vector{Int64} = [length(r.mail) for r in s.route]     # Number of mail in each route

    s       ::Session = Session(s.Lmax, [Route(r.id, zeros(Int64, O), r.mail) for r in s.route], zeros(Int64, O)) # Empty Session

    # ====================< Re-optimisation >====================

    stack   ::Vector{Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}, Union{Model, Nothing}}} = [] # Stored information about each filled output

    total_count_reopt     ::Int64         = 0                             # Total number of re-optimisation
    total_max_reopt       ::Int64         = ceil(Int64, 2O)               # Maximum number of re-optimisation
    count_reopt     ::Vector{Int64} = zeros(Int64, O)               # Count for each output how many time the model has been re-optimized
    max_reopt       ::Int64         = max(5, ceil(Int64, R/4))
    threshold_reopt ::Int64         = max(10, ceil(Int64, O/10))    # min output index before big backtrack step
    last_k          ::Int64         = 0                             # last output that required reoptimisation
    repeated_k      ::Int64         = 0                             # how many times last_k has been reoptimisation

    while k <= O && total_count_reopt <= total_max_reopt

        print_verbose("<k=$k,reopt=$total_count_reopt> ")
        
        if length(stack) < k
            
            V::Vector{Int64}    = [(0 < j <= length(r.mail)) ? r.mail[j] : 0 for (j, r) in zip(J, s.route)]     # Volume of the first not-assigned mail of each route
            F::Vector{Int64}    = [r for r=1:R if ((0 < J[r] <= B[r]) && (B[r] - J[r] >= O - k))]               # Set of route that need to force a mail
            C::Int64            = Lmax - sum([V[r] for r in F])                                                 # Remaining capacity after forced affectations
            M::Vector{Int64}    = [r for r=1:R if ((0 < J[r] <= B[r]) && (B[r] - J[r] < O - k) && (V[r] <= C))] # Set of route the model consider to insert mail

            print_verbose("<$Lneed> ", (ANSI_bold, ANSI_cyan))

            forced_affect_success::Bool, Lneed = Forced_mail_affectation_V4!(s, k, F, V, J, Lneed)
            
            if forced_affect_success

                print_verbose("<$Lneed> ", (ANSI_bold, ANSI_cyan))

                model, model_affect_success, Lneed = Model_mail_affectation_V4!(s, k, O, M, V, J, B, Lneed, τ, δ, nothing, tl, env)

                print_verbose("<$Lneed> ", (ANSI_bold, ANSI_cyan))

                if model_affect_success
                    push!(stack, (deepcopy(F), deepcopy(M), deepcopy(V), model))
                    (Lneed == 0) ? (k=O+1) : (k += 1)
                else
                    k -= 1
                end
            else
                k -= 1
            end
        else
            isempty(stack) && (return s, false)

            k, total_count_reopt, last_k, repeated_k, rollback = backtrack(k, total_count_reopt, total_max_reopt, count_reopt, max_reopt, threshold_reopt, last_k, repeated_k, O)

            for r=1:R
                if s.route[r].assignment[k+1] != 0
                    s.route[r].assignment[k+1] = 0
                    J[r] -= 1
                end
            end
            Lneed += s.load[k+1]
            s.load[k+1] = 0
            F, M, V, model = ([], [], [], nothing)

            for i=0:rollback
                (i != 0) && (k -= 1)

                F, M, V, model = pop!(stack)

                for r=1:R
                    if s.route[r].assignment[k] != 0
                        s.route[r].assignment[k] = 0
                        J[r] -= 1
                    end
                end
                Lneed += s.load[k]
                s.load[k] = 0
            end

            if model === nothing
                print_verbose("<No model> ", ANSI_magenta)
                k -= 1
            else

                print_verbose("<$Lneed> ", (ANSI_bold, ANSI_cyan))

                _, Lneed = Forced_mail_affectation_V4!(s, k, F, V, J, Lneed)

                print_verbose("<$Lneed> ", (ANSI_bold, ANSI_cyan))
                print_verbose("<Re-optimisation> ", ANSI_magenta)

                update_01KP_V4!(M, V, model)

                model, model_affect_success, Lneed = Model_mail_affectation_V4!(s, k, O, M, V, J, B, Lneed, τ, δ, model, tl, env)

                print_verbose("<$Lneed> ", (ANSI_bold, ANSI_cyan))

                if model_affect_success
                    (Lneed == 0) ? (k=O+1) : (k += 1)
                end

                push!(stack, (deepcopy(F), deepcopy(M), deepcopy(V), model))
            end
        end
        println_verbose()
    end

    return s, k == O+1
end

## ============================================================================================================== ##
##                                               ##    ##   ######                                                ##
##                                               ##    ##  ##                                                     ##
##                                                ##  ##    ######                                                ##
##                                                ##  ##         ##                                               ##
##                                                  ##      ######                                                ##
## ============================================================================================================== ##






