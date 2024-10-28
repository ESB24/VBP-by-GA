using Gurobi
using JuMP

begin 
    include("Instance.jl")
end


"""
```Julia
    model_01LP(instance::Instance, N::Int64, time_limit::Int64 = 600)
```
0-1 LP model for the full Vector Bin-Packing problem.
`instance` is the evaluated instance and `N` is an upper bound on the number of session.
"""

function model_01LP(instance::Instance, N::Int64, time_limit::Int64 = 600)
    # ===========< Notations: >===========
    # is the set of all rounds
    R::Vector{Int64} = collect(1:instance.nbRound)

    # is the set of batches for the round r ∈ R
    Br::Vector{Vector{Int64}} = [[k for (k, v) in enumerate(r.assignment) if v != 0] for r in instance.rounds]

    # is the volume of the j-th batch in the round r ∈ R. Here, j ∈ B(r)
    vrj::Vector{Vector{Int64}} = [[v for v in r.assignment if v != 0] for r in instance.rounds]

    # is the set of available machine outputs
    O::Vector{Int64} = collect(1:instance.nbOut) 

    # is the interval of potential outputs for the j-th batch in the round r ∈ R
    Orj::Vector{Vector{Vector{Int64}}} = [[collect(j: length(O) - length(Br[r]) + j) for j in Br[r]] for r in R]

    # is the set of mail bathes of round r, which can be potentially assigned to the output k ∈ O
    Urk::Vector{Vector{Vector{Int64}}} = [[[k for (k, v) in enumerate(Orj[r]) if o in v] for o in O] for r in R]

    # is the maximal load per output for any sorting session
    Lmax = instance.C

    # is the set of sorting sessions, where N is an upper bound on its number
    S::Vector{Int64} = collect(1:N)

    # ===========< Model: >===========

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "TimeLimit", time_limit)

    # ===========< Variables: >===========

    # 1 if the j-th batch of round r is assigned to the output k ∈ O(r)j, 0 otherwise
    @variable(model, x[r in R, j in Br[r], k in Orj[r][j]], Bin)

    # 1 if round r is located to session s, 0 otherwise
    @variable(model, y[r in R, s in S], Bin)

    # 1 if at least one round is located to session s, 0 otherwise
    @variable(model, z[s in S], Bin)

    # 1 if round r is allocated to session s and the j-th batch of round r is assigned to the output k ∈ O(r) j , 0 otherwise
    @variable(model, θ[r in R , j in Br[r] , k in Orj[r][j] , s in S], Bin)

    # ===========< Objective: >===========

    # (1) -> inimizing the number of used session
    @objective(model, Min, sum(z))

    # ===========< Constraint: >===========

    # (2) -> force each round to be allocated to exactly one session
    for r in R
        @constraint(model, sum([y[r, s] for s in S]) == 1)
    end

    # (3) -> number of allocated rounds per session can not be greater than the number of machine outputs is
    for s in S
        @constraint(model, sum([y[r, s] for r in R]) <= length(O) * z[s])
    end

    # (4) -> require that each batch j of round r be assigned to exactly oneoutput
    for r in R
        for j in Br[r]
            @constraint(model, sum([x[r, j, k] for k in Orj[r][j]]) == 1)
        end
    end

    # (5) -> for a given round, there can be at most one batch per output.
    for r in R
        for k in O
            @constraint(model, sum([x[r, j, k] for j in Urk[r][k]]) <= 1)
        end
    end

    # (6) -> link between the variables x, y and θ
    for r in R
        for s in S
            for j in Br[r]
                for k in Orj[r][j]
                    @constraint(model, x[r, j, k] + y[r, s] <= θ[r, j, k, s] + 1)
                end
            end
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
    for s in S
        for k in O
            @constraint(model, sum([sum([(vrj[r][j] * θ[r, j, k, s]) for j in Urk[r][k]]) for r in R]) <= Lmax * z[s])
        end
    end

    # (9) -> used to avoid symmetry in this problem by preventing the next session from opening when the previous one is still empty
    for s in S
        if s != N
            @constraint(model, z[s+1] <= z[s])
        end
    end

    optimize!(model)

    if termination_status(model) == OPTIMAL || MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0
        # ==========< Re-building Solution instance >==========
        perm::Vector{Int64} = zeros(Int64, instance.nbRound)
        permId::Int64 = 1
        sol::Solution = Solution(perm, [Session(instance.C, instance.nbOut) for s in S if s != 0])
        for s in S
            if value(z[s]) != 0

                roundId::Int64 = 0
                for r in R
                    if value(y[r, s]) == 1
                        perm[permId] = r
                        permId += 1
                        push!(sol.sessions[s].rounds, Round(r, zeros(Int64, length(O)), instance.rounds[r].batches))
                        roundId += 1
                        for k in O
                            val::Int64 = 0

                            for j in Urk[r][k]

                                if value(θ[r, j, k, s]) == 1
                                    val = vrj[r][j]
                                end
                            end

                            sol.sessions[s].rounds[roundId].assignment[k] = val
                            sol.sessions[s].loads[k] += val
                        end

                    end
                end

            end
        end
        return model, sol
    else
        return model, nothing
    end
end


"""
```Julia
    model_01LP_warmstart(instance::Instance, sol::Solution, time_limit::Int64 = 600)::Tuple{Model, Union{Nothing, Solution}}
```
0-1 LP model for the full Vector Bin-Packing problem. 
`instance` is the evaluated instance and `sol` is a solution used to do a warmup.
"""
function model_01LP_warmstart(instance::Instance, sol::Solution, time_limit::Int64 = 600)::Tuple{Model, Union{Nothing, Solution}}
# ===========< Parameters: >===========
    N = length(sol.sessions)

    # is the set of all rounds
    R::Vector{Int64} = collect(1:instance.nbRound)

    # is the set of batches for the round r ∈ R
    Br::Vector{Vector{Int64}} = [[k for (k, v) in enumerate(r.assignment) if v != 0] for r in instance.rounds]

    # is the volume of the j-th batch in the round r ∈ R. Here, j ∈ B(r)
    vrj::Vector{Vector{Int64}} = [[v for v in r.assignment if v != 0] for r in instance.rounds]

    # is the set of available machine outputs
    O::Vector{Int64} = collect(1:instance.nbOut) 

    # is the interval of potential outputs for the j-th batch in the round r ∈ R
    Orj::Vector{Vector{Vector{Int64}}} = [[collect(j: length(O) - length(Br[r]) + j) for j in Br[r]] for r in R]
    
    # is the set of mail bathes of round r, which can be potentially assigned to the output k ∈ O
    Urk::Vector{Vector{Vector{Int64}}} = [[[k for (k, v) in enumerate(Orj[r]) if o in v] for o in O] for r in R]

    # is the maximal load per output for any sorting session
    Lmax = instance.C

    # is the set of sorting sessions, where N is an upper bound on its number
    S::Vector{Int64} = collect(1:N)

# ===========< Model: >===========

    model = Model(Gurobi.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "TimeLimit", time_limit)

# ===========< Variables: >===========

    # 1 if the j-th batch of round r is assigned to the output k ∈ O(r)j, 0 otherwise
    @variable(model, x[r in R, j in Br[r], k in Orj[r][j]], Bin)

    # 1 if round r is located to session s, 0 otherwise
    @variable(model, y[r in R, s in S], Bin)

    # 1 if at least one round is located to session s, 0 otherwise
    @variable(model, z[s in S], Bin)

    # 1 if round r is allocated to session s and the j-th batch of round r is assigned to the output k ∈ O(r) j , 0 otherwise
    @variable(model, θ[r in R , j in Br[r] , k in Orj[r][j] , s in S], Bin)

# ===========< Warmup: >===========

    for (sId::Int64, s::Session) in enumerate(sol.sessions)
        set_start_value(z[sId], 1)

        for (rId::Int64, r::Round) in enumerate(s.rounds)
            set_start_value(y[r.id, sId], 1)

            batch::Int64 = 1
            batchMax::Int64 = length(r.batches)
            out::Int64 = 1

            while batch <= batchMax
                if r.assignment[out] != 0
                    set_start_value(x[r.id, batch, out], 1)
                    set_start_value(θ[r.id, batch, out, sId], 1)
                    batch += 1
                end
                out += 1
            end
        end
    end

# ===========< Objective: >===========

    # (1) -> inimizing the number of used session
    @objective(model, Min, sum(z))

# ===========< Constraint: >===========

    # (2) -> force each round to be allocated to exactly one session
    for r in R
        @constraint(model, sum([y[r, s] for s in S]) == 1)
    end

    # (3) -> number of allocated rounds per session can not be greater than the number of machine outputs is
    for s in S
        @constraint(model, sum([y[r, s] for r in R]) <= length(O) * z[s])
    end

    # (4) -> require that each batch j of round r be assigned to exactly oneoutput
    for r in R
        for j in Br[r]
            @constraint(model, sum([x[r, j, k] for k in Orj[r][j]]) == 1)
        end
    end

    # (5) -> for a given round, there can be at most one batch per output.
    # for r in R
    #     for k in O
    #         @constraint(model, sum([x[r, j, k] for j in Urk[r][k]]) <= 1)
    #     end
    # end

    # for r in R
    #     @constraint(model, sum([x[r, 1, k] for k in Orj[r][1]]) == 1)
    # end

    # (6) -> link between the variables x, y and θ
    for r in R
        for s in S
            for j in Br[r]
                for k in Orj[r][j]
                    @constraint(model, x[r, j, k] + y[r, s] <= θ[r, j, k, s] + 1)
                end
            end
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
    for s in S
        for k in O
            @constraint(model, sum([sum([(vrj[r][j] * θ[r, j, k, s]) for j in Urk[r][k]]) for r in R]) <= Lmax * z[s])
        end
    end

    # (9) -> used to avoid symmetry in this problem by preventing the next session from opening when the previous one is still empty
    for s in S
        if s != N
            @constraint(model, z[s+1] <= z[s])
        end
    end

    optimize!(model)

# ==========< Results >==========

    if termination_status(model) == OPTIMAL || MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0
        perm::Vector{Int64} = zeros(Int64, instance.nbRound)
        permId::Int64 = 1
        sol::Solution = Solution(perm, [Session(instance.C, instance.nbOut) for s in S if s != 0])
        for s in S
            if value(z[s]) != 0

                roundId::Int64 = 0
                for r in R
                    if value(y[r, s]) == 1
                        perm[permId] = r
                        permId += 1
                        push!(sol.sessions[s].rounds, Round(r, zeros(Int64, length(O)), instance.rounds[r].batches))
                        roundId += 1
                        for k in O
                            val::Int64 = 0

                            for j in Urk[r][k]

                                if value(θ[r, j, k, s]) == 1
                                    val = vrj[r][j]
                                end
                            end

                            sol.sessions[s].rounds[roundId].assignment[k] = val
                            sol.sessions[s].loads[k] += val
                        end

                    end
                end

            end
        end
        return model, sol
    else
        return model, nothing
    end
end

