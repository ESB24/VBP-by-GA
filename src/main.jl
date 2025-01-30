begin # Library
    # Test
	import BenchmarkTools
	import ProfileCanvas

    # Miscelaneous
    using Random
end

begin # Files
    include("parsing.jl")

    include("Round.jl")
    include("Session.jl")
    include("Solution.jl")

    include("Instance.jl")
    include("MILP.jl")
    include("Miscelaneous.jl")
    
    # include("OptiMove.jl")
    
    include("GeneticAlgorithm.jl")

    include("SessionRebuild.jl")
end

call = 0
repairedBuild = 0
improvedOverAll = 0
locked = 0
delta1 = 0.2
delta2 = 0.9

# ================================================================================= #
#  #######  ##    #   ######  #######   #####   ##    #   ######  #######   ######  #
#     #     # #   #  #           #     #     #  # #   #  #        #        #        #
#     #     #  #  #   #####      #     #######  #  #  #  #        ####      #####   #
#     #     #   # #        #     #     #     #  #   # #  #        #              #  #
#  #######  #    ##  ######      #     #     #  #    ##   ######  #######  ######   #
# ================================================================================= #

begin
    # ==========< Gathering Datas >==========
    mat = parseMailXLSX("../data/Trafic/trafic_20_200_35_1.xlsx")
    _, O = size(mat) 
    # ==========< Initialisation >==========
    # =====< Full instance >=====
    λ = 0.25
    instance::Instance = Instance(mat, λ)
    # instance.C = max(round(Int64, sum(mat) / (5 * O)), maximum(mat))
    # =====< Small instance >=====
    instanceS::Instance = Instance(mat[1:40, :], λ)
    # instanceS.C = max(round(Int64, sum(mat) / (5 * O)), maximum(mat))
end

begin
    # Oscar
    numberOfOutputs, numberOfRounds, numberOfSessions, sessions, Lmax = parseOptiInstance("TER/data/hard/test/instance_1_200_644_opt_10_50_50.txt")
    mat = vect2D_to_matrix(normalizeOptiInstance(sessions))
    λ = 0.25
    instance::Instance = Instance(mat, λ)
    instance.C = Lmax
end
begin
    # My instances
    Lmax, mat, nbSession = parseMyInstance("TER/data/Gaussian/instanceGaussian_1_O200_R1000_C300_opt_50.txt")
    λ = 0.25
    instance::Instance = Instance(mat, λ)
    instance.C = Lmax
end

# O = 20
writeInstanceBatterie_Gaussian("TER/data/Gaussian", 10, 10, 20, 20, 100, 2, 20, 2, 5, 16:20)
writeInstanceBatterie_Distrib("TER/data/Distrib", 10, 10, 60, 20, 100, 2, 20, 2, 5, 16:20)
writeInstanceBatterie_BigBatche("TER/data/BigBatch", 10, 10, 20, 20, 100, 2, 20, 2, 5, 16:20)


writeInstanceBatterie_Distrib("TER/data/Distrib", 10, 10, 60, 1, 100, 2, 20, 300, 20, distrib_1onN)
writeInstanceBatterie_BigBatche("TER/data/BigBatch", 10, 10, 60, 1:10, 80:100, 1:5, 2, 20, 300, 20)
# O = 200
writeInstanceBatterie_Gaussian("TER/data/Gaussian", 1, 1, 20, 200, 150, 2, 100, 10, 8, 16:20)
writeInstanceBatterie_Distrib("TER/data/Distrib", 100, 1, 30, 200, 40, 2, 100, 1, 40, distrib_1onN)
writeInstanceBatterie_BigBatche("TER/data/BigBatch", 100, 1, 20, 200, 80, 2, 100, 1:10, 60:70, 1:3)
begin
    call = 0
    repairedBuild = 0
    improvedOverAll = 0
    locked = 0
    sumObj = 0

    Δ               ::Int64                             = 2
    # TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD
    TAG_AddRound    ::Vector{Type{<:SimpleAddRound}}    = Type{<:SimpleAddRound}[SAVA, NFBA, OPTIMOVE_STAGE1]
    env             ::Gurobi.Env                        = Gurobi.Env()
    tl              ::Int64                             = 5

    iter::Int64 = 10
    for i=1:iter
        # Gaussian
        Lmax, mat, nbSession = parseMyInstance("TER/data/BigBatch/instanceBigBatche_$(i)_O200_R600_C300_opt_10.txt")
        λ = 0.25
        instance::Instance = Instance(mat[1:55, :], λ)
        instance.C = Lmax

        perm            ::Vector{Int64}                     = sortperm([(-fitness(r, Max), -fitness(r, NbBatches)) for r in instance.rounds])

        sol, flag, _, _ = GA(instance, TAG_FitGA=ObjGA_2)# sol = buildSolution_FFD_final(instance, perm, TAG_FitSes)
        sumObj += fitness(sol, ObjGA)
    end
    println("\n TAG -> $(OPTIMOVE_STAGE1)\n call -> $(call)\n repairedBuild -> $(repairedBuild)\n repaired % -> $((repairedBuild * 100) / call)\n improvedOverAll -> $(improvedOverAll)\n mean obj -> $(sumObj / iter)")
end

# ================================================================================= #
#                   ######   #######  ##    #   ######  #     #                     #
#                   #     #  #        # #   #  #        #     #                     #
#                   ######   ####     #  #  #  #        #######                     #
#                   #     #  #        #   # #  #        #     #                     #
#                   #######  #######  #    ##   ######  #     #                     #
# ================================================================================= #

# ==========< Benchmark 1: >==========
# best current result in average
perm            ::Vector{Int64}                     = sortperm([(-fitness(r, Max), -fitness(r, NbBatches), -fitness(r, StdBatches)) for r in instance.rounds])
Δ               ::Int64                             = 2
TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD
TAG_AddRound    ::Vector{Type{<:SimpleAddRound}}    = Type{<:SimpleAddRound}[SAVA, NFBA, OPTIMOVE_VALID_INF]
env             ::Gurobi.Env                        = Gurobi.Env()
tl              ::Int64                             = 5

call = 0
repairedBuild = 0
improvedOverAll = 0
locked = 0

s::Session = Session(instance.C, [], zeros(Int64, instance.nbOut))
sol::Solution = Solution(perm, [s])

buildProcedure!(instance, sol, sol.permutation[1])
buildProcedure!(instance, sol, sol.permutation[2])
buildProcedure!(instance, sol, sol.permutation[3])
buildProcedure!(instance, sol, sol.permutation[4])

BenchmarkTools.@benchmark begin; buildSolution_FF(instance, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRound=TAG_AddRound, env=env, tl=tl); end
BenchmarkTools.@benchmark buildSolution_FFD_final(instance, perm, TAG_FitSes)
BenchmarkTools.@benchmark buildSolution_FFD_SAVA_SPCF_OPTIMOVE(instance, perm, TAG_FitSes=TAG_FitSes, env=env)



begin
    sol = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRound=TAG_AddRound, env=env, tl=tl)
    println(sol)
end

tmp = fitness(sol, ObjGA)
ProfileCanvas.@profview begin
                            for _=1:100
                                sol = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRound=TAG_AddRound, env=env, tl=tl)
                            end
                        end
begin
    for TAG_OPTIMOVE in [OPTIMOVE_S1_NEGSTD, OPTIMOVE_STAGE1, OPTIMOVE_VALID_INF, OPTIMOVE_3S] # OPTIMOVE_1BATCH, OPTIMOVE_VALID
        println(">====================<")
        call = 0
        repairedBuild = 0
        improvedOverAll = 0
        sumObj = 0
        iter = 10000
        TAG_AddRound    ::Vector{Type{<:SimpleAddRound}}    = Type{<:SimpleAddRound}[SAVA, TAG_OPTIMOVE]
        for _=1:iter
            sol = buildSolution_FF(instance, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRound=TAG_AddRound, env=env, tl=tl)
            sumObj += fitness(sol, ObjGA)
        end
        println("\n TAG -> $(TAG_OPTIMOVE)\n call -> $(call)\n repairedBuild -> $(repairedBuild)\n repaired % -> $((repairedBuild * 100) / call)\n improvedOverAll -> $(improvedOverAll)\n mean obj -> $(sumObj / iter)")
    end
end

begin
    call = 0
    repairedBuild = 0
    improvedOverAll = 0
    TAG_AddRound    ::Vector{Type{<:SimpleAddRound}}    = Type{<:SimpleAddRound}[SAVA, OPTIMOVE_S1_NEGSTD]
    for _=1:100
        # println("|")
        buildSolution_FF(instance, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRound=TAG_AddRound, env=env, tl=tl)
    end
    # s = buildSolution_FF(instance, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRound=TAG_AddRound, env=env, tl=tl)
    # print(s)
    println("\n call -> $(call)\n repairedBuild -> $(repairedBuild)\n improvedOverAll -> $(improvedOverAll)")
end

# ==========< Benchmark 2: >==========
# current laposte behavior
BenchmarkTools.@benchmark   begin
    perm            ::Vector{Int64}                     = sortperm([(-fitness(r, Identity)) for r in instance.rounds])
    Δ               ::Int64                             = 2
    TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD
    TAG_AddRound    ::Vector{Type{<:SimpleAddRound}}    = Type{<:SimpleAddRound}[RAW]
    s = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRound=TAG_AddRound, env=env)
end

# ==========< Benchmark 3: >==========
# TER test (random permutation)
BenchmarkTools.@benchmark   begin
    perm            ::Vector{Int64}                     = sortedPerm(instance, nothing) # nothing as sorting criteria create a random permutation
    Δ               ::Int64                             = 2
    TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD
    TAG_AddRound    ::Vector{Type{<:SimpleAddRound}}    = Type{<:SimpleAddRound}[SAVANT, OPTIMOVE]
    s = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRound=TAG_AddRound)
end

# ==========< Benchmark 4: >==========
# SAVANT but select the solution with the best standard deviation not the first found
BenchmarkTools.@benchmark   begin
    perm            ::Vector{Int64}                     = sortedPerm(instance, nothing) # nothing as sorting criteria create a random permutation
    Δ               ::Int64                             = 2
    TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD
    TAG_AddRound    ::Vector{Type{<:SimpleAddRound}}    = Type{<:SimpleAddRound}[SAVANT_MIN_STD, OPTIMOVE]
    s = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRound=TAG_AddRound)
end


# ==========< Benchmark 5: >==========
# SAVANT but select the solution with the best standard deviation not the first found
BenchmarkTools.@benchmark   begin
    perm            ::Vector{Int64}                     = sortedPerm(instance, StdAssignment) # nothing as sorting criteria create a random permutation
    Δ               ::Int64                             = 2
    TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD
    TAG_AddRound    ::Vector{Type{<:SimpleAddRound}}    = Type{<:SimpleAddRound}[RAW, SPCF]
    s = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRound=TAG_AddRound)
end

# ================================================================================= #
#                                ######      #####                                  #
#                               #           #     #                                 #
#                               #   ####    #######                                 #
#                               #      #    #     #                                 #
#                                ######     #     #                                 #
# ================================================================================= #
begin # GA
    instance::Instance, _ = parseAnyInstance("instanceBigBatche_1_O200_R200_C80_opt_10.txt")

    max_time        ::Float64                       = float(300)
    max_pop         ::Int64                         = 100
    λ_core          ::Float64                       = 0.5
    λ_elite         ::Float64                       = 0.4
    max_gen         ::Int64                         = 1000
    TAG_FitGA                                       = ObjGA_2
    TAG_Select                                      = SelectElite
    TAG_Crossover   ::Tuple                         = (Edit, Point1, PointN)
    TAG_Accept                                      = Always
    maxStagnate     ::Int64                         = 10

    best_sol::Solution, flag_opti::Bool, gen_val, Runtime::Float64 = GA(
        instance, 
        max_time=max_time, 
        max_pop=max_pop, 
        max_gen=max_gen,
        λ_core=λ_core,
        λ_elite=λ_elite,
        TAG_FitGA=TAG_FitGA,
        TAG_Select=TAG_Select,
        TAG_Crossover=TAG_Crossover,
        TAG_Accept=TAG_Accept,
        maxStagnate=maxStagnate
    )

    println(best_sol)
end

# ================================================================================= #
#                           ##   ##  #######  #        ######                       #
#                           # # # #     #     #        #     #                      #
#                           #  #  #     #     #        ######                       #
#                           #     #     #     #        #                            #
#                           #     #  #######  #######  #                            #
# ================================================================================= #

begin # Construction of solution using MILP
    perm::Vector{Int64} = sortperm([(-fitness(r, Max), -fitness(r, NbBatches), -fitness(r, StdBatches)) for r in instance.rounds])

    s0::Solution = buildSolution_FFD_final(instance, perm, LoadSTD)
    model, sBest = model_01LP_warmstart(instance, s0, 20)
end

termination_status(model)
# computation time
cpu = MOI.get(model, Gurobi.ModelAttribute("RunTime"))
# number of admissible solution
nb_sol = MOI.get(model, Gurobi.ModelAttribute("SolCount"))
# gap between best two bound (lower and upper)
gap = MOI.get(model, Gurobi.ModelAttribute("MIPGap"))
# best found solution (only if nb_sol > 0)
sol = MOI.get(model, Gurobi.ModelAttribute("ObjVal"))

# ================================================================================= #
#  ##   ##  #######  #        ######            #     #   #####   ######   ##   ##  #
#  # # # #     #     #        #     #           #     #  #     #  #     #  # # # #  #
#  #  #  #     #     #        ######    #####   #  #  #  #######  ######   #  #  #  #
#  #     #     #     #        #                 # # # #  #     #  #    #   #     #  #
#  #     #  #######  #######  #                 ##   ##  #     #  #     #  #     #  #
# ================================================================================= #

begin # Construction of solution using MILP
    s0::Solution = buildSolution_FFD(instanceS, sortedPerm(instanceS, Type{<:FitnessRound}[Max, NbBatches]))
    display(s0)

    N = length(s0.sessions)
    
    model, sBest = model_01LP_warmstart(instanceS, s0, 900)
end

termination_status(model)
# computation time
cpu = MOI.get(model, Gurobi.ModelAttribute("RunTime"))
# number of admissible solution
nb_sol = MOI.get(model, Gurobi.ModelAttribute("SolCount"))
# gap between best two bound (lower and upper)
gap = MOI.get(model, Gurobi.ModelAttribute("MIPGap"))
# best found solution (only if nb_sol > 0)
sol = MOI.get(model, Gurobi.ModelAttribute("ObjVal"))

# ================================================================================= #
#                           ##   ##  #######   ######   ######                      #
#                           # # # #     #     #        #                            #
#                           #  #  #     #      #####   #                            #
#                           #     #     #           #  #                            #
#                           #     #  #######  ######    ######                      #
# ================================================================================= #
begin
    # all_δ1::Vector{Float64} = [0., 0.2, 0.4, 0.6, 0.8, 1.]
    # all_δ2::Vector{Float64} = [0., 0.2, 0.4, 0.6, 0.8, 1.]

    all_δ1::Vector{Float64} = [0.2]#, 0.2, 0.4, 0.6, 0.8, 1.0]
    all_δ2::Vector{Float64} = [0.0]#, 0.2, 0.4, 0.6, 0.8, 1.0]
    all_δ3::Vector{Float64} = [0.0]#, 0.2, 0.4, 0.6, 0.8, 1.0]
    all_δ4::Vector{Float64} = [0.0]#, 0.2, 0.4, 0.6, 0.8, 1.0]

    res::Array{Int64, 4} = zeros(Int64, length(all_δ1), length(all_δ2), length(all_δ3), length(all_δ4))
    sum_obj::Array{Int64, 4} = zeros(Int64, length(all_δ1), length(all_δ2), length(all_δ3), length(all_δ4))

    env::Gurobi.Env = Gurobi.Env()
    tl::Int64 = 10

    for i=1:100
        print(" - $i - ")
        instance::Instance, _ = parseAnyInstance("instanceDistribHard_$(i)_O100_R20_C100_opt_1.txt")
        s = Session(instance.C, instance.rounds, compute_output(instance.rounds))

        for id_δ1 = 1:length(all_δ1)
            for id_δ2 = 1:length(all_δ2)
                for id_δ3 = 1:length(all_δ3)
                    for id_δ4 = 1:length(all_δ4)
                        s0, valid = rebuildSession_knapSack_model_V3!(s, tl, env, [all_δ1[id_δ1], all_δ2[id_δ2], all_δ3[id_δ3], all_δ4[id_δ4]])
                        sum_obj[id_δ1, id_δ2, id_δ3, id_δ4] += sum(s0.loads)
                        valid && (res[id_δ1, id_δ2, id_δ3, id_δ4] += 1)
                        print("x")
                    end
                end
            end
        end
        println()
    end

    print("optimum:")
    println(res)
    print("mean:")
    println(sum_obj / (100*200) )
end

begin
    nb_inst::Int64 = 100
    
    all_std::Vector{Float64} = Vector{Float64}(undef, nb_inst)
    all_sum::Vector{Float64} = Vector{Float64}(undef, nb_inst)

    
    for i=1:nb_inst
        mat, _ = parseAnyInstanceMat("instanceBigBatche_$(i)_O200_R20_C80_opt_1.txt")
        all_std[i] = std(mat)
        all_sum[i] = sum(mat)/count(x -> x != 0, mat)
    end
end

begin
    all_path = [
        [(1, "myTrafic_$(i)_O100_R100.txt") for i=1:10];
        [(1, "myTrafic_$(i)_O40_R100.txt") for i=1:10];
        [(2, "instanceDistrib_$(i)_O100_R150_C40_opt_5.txt") for i=1:10];
        [(2, "instanceDistrib_$(i)_O40_R150_C40_opt_5.txt") for i=1:10];
        [(3, "instanceGaussian_$(i)_O100_R100_C100_opt_5.txt") for i=1:10];
        [(3, "instanceGaussian_$(i)_O40_R100_C100_opt_5.txt") for i=1:10];
        [(4, "instanceBigBatche_$(i)_O100_R100_C80_opt_5.txt") for i=1:10];
        [(4, "instanceBigBatche_$(i)_O40_R100_C80_opt_5.txt") for i=1:10];
        [(5, "instanceDistribHard_$(i)_O100_R100_C100_opt_5.txt") for i=1:10];
        [(5, "instanceDistribHard_$(i)_O40_R100_C100_opt_5.txt") for i=1:10];
    ]

    res = []
    env = Gurobi.Env()

    cpt = 0

    all_TAG = (MaxMin, Max, StdBatches, StdAssignment, NbBatches, MeanBatches, MeanAssignment)
    for TAG1 in all_TAG
        for TAG2 in all_TAG
            if TAG2 != TAG1
                #for TAG3 in all_TAG
                    #if TAG3 != TAG2
                        cpt += 1
                        println("\n - $cpt - ")

                        opti = zeros(Float64, 6)
                        obj = zeros(Float64, 6)
                        ses = zeros(Float64, 6)
                        for (i, path) in all_path
                            print("-")
                            instance::Instance, _ = parseAnyInstance(path)
                            perm::Vector{Int64} = sortperm([(-fitness(r, TAG1), -fitness(r, TAG2)) for r in instance.rounds]) # , -fitness(r, TAG3)
                            sol = buildSolution_FFD_final(instance, perm, 10, env)

                            (length(sol.sessions) == minSession(instance)) && (opti[i] += 1)
                            obj[i] += fitness(sol, ObjGA_2)
                            ses[i] += length(sol.sessions)

                            (length(sol.sessions) == minSession(instance)) && (opti[end] += 1)
                            obj[end] += fitness(sol, ObjGA_2)
                            ses[end] += length(sol.sessions)
                        end

                        push!(res, ((TAG1, TAG2), opti, obj, ses)) # , TAG3
                #     end
                # end
            end
        end
    end

    sort!(by = x -> (x[2][2], x[4][2], x[3][2]), res)

    for (tags, opt, obj, ses) in res
        println("- opt = $opt, obj = $obj, ses = $ses <- $(tags)")
    end
end

sort!(by = x -> (x[2][6], x[4][6], x[3][6]), res)
for (tags, opt, obj, ses) in res
    println("$(tags)")
end


# ================================================================================= #
#                           #######  ##    #   ######  #######                      #
#                              #     # #   #  #           #                         #
#                              #     #  #  #   #####      #                         #
#                              #     #   # #        #     #                         #
#                           #######  #    ##  ######      #                         #
# ================================================================================= #

I = 1
S = 1
O = 200
R = 10

# writeInstanceBatterie_Gaussian("TER/data/Gaussian", I, S, R, O, 100, 2, ceil(Int64, 2O/3), 10, 6, 16:20)
# writeInstanceBatterie_BigBatche("TER/data/BigBatch", I, S, R, O, 60, 2, ceil(Int64, 2O/3), 1:5, 50:55, 1:3)
# writeInstanceBatterie_Distrib("TER/data/DistribHard", I, S, R, O, 100, 2, ceil(Int64, 2O/3), 15, 40, distrib_1onN)

writeInstanceBatterie_Distrib("TER/data/DistribHard", I, S, R, O, 50, 2, ceil(Int64, 2O/3), 15, 40, distrib_1onN)

for O in [20, 40, 100, 200]
    for R in [20, 40, 100, 200, 300, 400]
        writeInstanceBatterie_Distrib_easy(I, R, O)
    end
end

begin
    env = Gurobi.Env()
    tl_perm = 100
    tl_rebuild = 10

    opt = 0
    for i=1:100
        instance, opti = parseAnyInstance("instanceGaussian_$(i)_O40_R200_C100_opt_10.txt")
        md, res, nb_ses, perm = full_partitioning_01_KP(instance.rounds, 100, 20, tl_perm, env)
        sol::Solution = buildSolution_BFD_final(instance, perm, tl_rebuild, env)
        length(sol.sessions) == opti ? (opt += 1, println("$(i) opti!")) : (println("$(i) end"))
    end
    print("|opti| = $(opt)")
end

instance, opti = parseAnyInstance("instanceGaussian_1_O200_R200_C100_opt_10.txt")
instance, opti = parseAnyInstance("myTrafic_1_O20_R20.txt")

tmp = filter(e -> !(e.id in [11, 27, 29, 40, 37, 9, 7, 14, 34, 19, 38, 35, 24, 28, 36, 23, 4, 18, 2]), instance.rounds)
instanceS = Instance(instance.C, instance.nbOut, length(tmp), tmp)

s = Session(instanceS.C, instanceS.rounds, compute_output(instanceS.rounds))
s0, valid = rebuildSession_knapSack_model_V3!(s)

sol = buildSolution_BF_final(instance)

C, mat, min_session = parseMyInstance_completed("TER/data/Distrib/instanceDistrib_1_O20_R60_C40_opt_2.txt")
R2, O2 = size(mat)
s = Session(C, O2)
for i=1:R2
    tmp = filter(x -> x != 0, mat[i, :])
    push!(s.rounds, Round(i, mat[i, :], ntuple(i -> tmp[i], length(tmp))))
end
compute_output!(s)
count(x -> x != 0, mat)/(20)
