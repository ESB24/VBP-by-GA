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

begin
    # ==========< Gathering Datas >==========
    mat = parseMailXLSX("../data/Trafic/trafic_20_200_35_1.xlsx")
    _, O = size(mat) 
    # ==========< Initialisation >==========
    # =====< Full instance >=====
    λ = 0.25
    instance::Instance = Instance(mat, λ)
    # instance.Lmax = max(round(Int64, sum(mat) / (5 * O)), maximum(mat))
    # =====< Small instance >=====
    instanceS::Instance = Instance(mat[1:40, :], λ)
    # instanceS.Lmax = max(round(Int64, sum(mat) / (5 * O)), maximum(mat))
end

buildSolution_BF(instance)

begin
    # Oscar
    numberOfOutputs, numberOfRoutes, numberOfSessions, sessions, Lmax = parseOptiInstance("TER/data/hard/test/instance_1_200_644_opt_10_50_50.txt")
    mat = vect2D_to_matrix(normalizeOptiInstance(sessions))
    λ = 0.25
    instance::Instance = Instance(mat, λ)
    instance.Lmax = Lmax
end
begin
    # My instances
    Lmax, mat, nbSession = parseMyInstance("TER/data/Gaussian/instanceGaussian_1_O200_R1000_C300_opt_50.txt")
    λ = 0.25
    instance::Instance = Instance(mat, λ)
    instance.Lmax = Lmax
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
    TAG_AddRoute    ::Vector{Type{<:SimpleAddRoute}}    = Type{<:SimpleAddRoute}[SAVA, NFBA, OPTIMOVE_STAGE1]
    env             ::Gurobi.Env                        = Gurobi.Env()
    tl              ::Int64                             = 5

    iter::Int64 = 10
    for i=1:iter
        # Gaussian
        Lmax, mat, nbSession = parseMyInstance("TER/data/BigBatch/instanceBigBatche_$(i)_O200_R600_C300_opt_10.txt")
        λ = 0.25
        instance::Instance = Instance(mat[1:55, :], λ)
        instance.Lmax = Lmax

        perm            ::Vector{Int64}                     = sortperm([(-fitness(r, Max), -fitness(r, NbBatches)) for r in instance.route])

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
perm            ::Vector{Int64}                     = sortperm([(-fitness(r, Max), -fitness(r, NbBatches), -fitness(r, StdBatches)) for r in instance.route])
Δ               ::Int64                             = 2
TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD
TAG_AddRoute    ::Vector{Type{<:SimpleAddRoute}}    = Type{<:SimpleAddRoute}[SAVA, NFBA, OPTIMOVE_VALID_INF]
env             ::Gurobi.Env                        = Gurobi.Env()
tl              ::Int64                             = 5

call = 0
repairedBuild = 0
improvedOverAll = 0
locked = 0

s::Session = Session(instance.Lmax, [], zeros(Int64, instance.nbOut))
sol::Solution = Solution(perm, [s])

buildProcedure!(instance, sol, sol.permutation[1])
buildProcedure!(instance, sol, sol.permutation[2])
buildProcedure!(instance, sol, sol.permutation[3])
buildProcedure!(instance, sol, sol.permutation[4])

BenchmarkTools.@benchmark begin; buildSolution_FF(instance, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRoute=TAG_AddRoute, env=env, tl=tl); end
BenchmarkTools.@benchmark buildSolution_FFD_final(instance, perm, TAG_FitSes)
BenchmarkTools.@benchmark buildSolution_FFD_SAVA_SPCF_OPTIMOVE(instance, perm, TAG_FitSes=TAG_FitSes, env=env)



begin
    sol = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRoute=TAG_AddRoute, env=env, tl=tl)
    println(sol)
end

tmp = fitness(sol, ObjGA)
ProfileCanvas.@profview begin
                            for _=1:100
                                sol = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRoute=TAG_AddRoute, env=env, tl=tl)
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
        TAG_AddRoute    ::Vector{Type{<:SimpleAddRoute}}    = Type{<:SimpleAddRoute}[SAVA, TAG_OPTIMOVE]
        for _=1:iter
            sol = buildSolution_FF(instance, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRoute=TAG_AddRoute, env=env, tl=tl)
            sumObj += fitness(sol, ObjGA)
        end
        println("\n TAG -> $(TAG_OPTIMOVE)\n call -> $(call)\n repairedBuild -> $(repairedBuild)\n repaired % -> $((repairedBuild * 100) / call)\n improvedOverAll -> $(improvedOverAll)\n mean obj -> $(sumObj / iter)")
    end
end

begin
    call = 0
    repairedBuild = 0
    improvedOverAll = 0
    TAG_AddRoute    ::Vector{Type{<:SimpleAddRoute}}    = Type{<:SimpleAddRoute}[SAVA, OPTIMOVE_S1_NEGSTD]
    for _=1:100
        # println("|")
        buildSolution_FF(instance, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRoute=TAG_AddRoute, env=env, tl=tl)
    end
    # s = buildSolution_FF(instance, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRoute=TAG_AddRoute, env=env, tl=tl)
    # print(s)
    println("\n call -> $(call)\n repairedBuild -> $(repairedBuild)\n improvedOverAll -> $(improvedOverAll)")
end

# ==========< Benchmark 2: >==========
# current laposte behavior
BenchmarkTools.@benchmark   begin
    perm            ::Vector{Int64}                     = sortperm([(-fitness(r, Identity)) for r in instance.route])
    Δ               ::Int64                             = 2
    TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD
    TAG_AddRoute    ::Vector{Type{<:SimpleAddRoute}}    = Type{<:SimpleAddRoute}[RAW]
    s = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRoute=TAG_AddRoute, env=env)
end

# ==========< Benchmark 3: >==========
# TER test (random permutation)
BenchmarkTools.@benchmark   begin
    perm            ::Vector{Int64}                     = sortedPerm(instance, nothing) # nothing as sorting criteria create a random permutation
    Δ               ::Int64                             = 2
    TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD
    TAG_AddRoute    ::Vector{Type{<:SimpleAddRoute}}    = Type{<:SimpleAddRoute}[SAVANT, OPTIMOVE]
    s = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRoute=TAG_AddRoute)
end

# ==========< Benchmark 4: >==========
# SAVANT but select the solution with the best standard deviation not the first found
BenchmarkTools.@benchmark   begin
    perm            ::Vector{Int64}                     = sortedPerm(instance, nothing) # nothing as sorting criteria create a random permutation
    Δ               ::Int64                             = 2
    TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD
    TAG_AddRoute    ::Vector{Type{<:SimpleAddRoute}}    = Type{<:SimpleAddRoute}[SAVANT_MIN_STD, OPTIMOVE]
    s = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRoute=TAG_AddRoute)
end


# ==========< Benchmark 5: >==========
# SAVANT but select the solution with the best standard deviation not the first found
BenchmarkTools.@benchmark   begin
    perm            ::Vector{Int64}                     = sortedPerm(instance, StdAssignment) # nothing as sorting criteria create a random permutation
    Δ               ::Int64                             = 2
    TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD
    TAG_AddRoute    ::Vector{Type{<:SimpleAddRoute}}    = Type{<:SimpleAddRoute}[RAW, SPCF]
    s = buildSolution_FFD(instance, perm, Δ=Δ, TAG_FitSes=TAG_FitSes, TAG_AddRoute=TAG_AddRoute)
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
    perm::Vector{Int64} = sortperm([(-fitness(r, Max), -fitness(r, NbBatches), -fitness(r, StdBatches)) for r in instance.route])

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
    s0::Solution = buildSolution_FFD(instanceS, sortedPerm(instanceS, Type{<:FitnessRoute}[Max, NbBatches]))
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
        s = Session(instance.Lmax, instance.route, compute_output(instance.route))

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
                            perm::Vector{Int64} = sortperm([(-fitness(r, TAG1), -fitness(r, TAG2)) for r in instance.route]) # , -fitness(r, TAG3)
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

## ============================================================================================================== ##
##            ########  ##    ##   #######  ########   ######   ##    ##   ######   ########   #######            ##
##               ##     ###   ##  ##           ##     ##    ##  ###   ##  ##    ##  ##        ##                  ##
##               ##     ## ## ##   ######      ##     ########  ## ## ##  ##        #####      ######             ##
##               ##     ##   ###        ##     ##     ##    ##  ##   ###  ##    ##  ##              ##            ##
##            ########  ##    ##  #######      ##     ##    ##  ##    ##   ######   ########  #######             ##
## ============================================================================================================== ##

# Changing the number of: output ("O"), session ("S") or instances ("I") can be done within each section without altering the instance 
# (except for the standard deviation of the mail to route distribution which is very uncontroled)

# the number of route per session impact: the capacity ("Lmax") and possibly some other parameter depending in the instance type
# (e.g. the number of generated gaussian for the contained instance)

## ============================================================================================================== ##
##                 #######        ##    #######  ########   #######             ######    ######                  ##
##                 ##    ##      ##    ##        ##        ##         ######         ##  ##   ###                 ##
##                 #######      ##      ######   #####      ######               ####    ## ## ##                 ##
##                 ##  ##      ##            ##  ##              ##   ######   ##        ###   ##                 ##
##                 ##   ##    ##       #######   ########  #######              ######    ######                  ##
## ============================================================================================================== ##

begin
    O = 200 
    n = 20
    nbGaussian = round(Int64, O/2):round(Int64, 3O/4)
    mat = generateMultipleGaussianSum(O, n, nbGaussian)
    x = 1:O
    p = plot()
    for i = 1:n
        plot!(x, mat[i, :] .+ (i*1.2), label="")
    end 
    display(p)
end

# USAGE: uncomment 1 instance generator and the corresponding parser:
# -> generate a new instance under: "../data/Tests"
# -> parse instance and display: std/mean volume/mail number + plots...

begin
    I           = 1
    S           = 1
    O           = 40
    R           = 20
    target      = "../data/Tests"
    # writeInstanceBatterie_Distrib(target, I, S, R, O, 65, 2, ceil(Int64, 2O/3), 1, 30, distrib_1onN)
    # writeInstanceBatterie_Contained(target, I, S, R, O, 75, 2, ceil(Int64, 2O/3), 4, 9, round(Int64, O/2):round(Int64, 3O/4))
    # writeInstanceBatterie_Skewed(target, I, S, R, O, 60, 2, ceil(Int64, 2O/3), 1:5, 40:44, 1:6)
    writeInstanceBatterie_Chunk(target, I, S, R, O, 95, 2, ceil(Int64, 2O/3), 10:10, 1:20, 1:20)

    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceIndus_1_O$(O)_R$(R)_C65_opt_$(S).txt")
    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceContained_1_O$(O)_R$(R)_C75_opt_$(S).txt")
    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceSkewed_1_O$(O)_R$(R)_C60_opt_$(S).txt")
    C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceChunk_1_O$(O)_R$(R)_C95_opt_$(S).txt")
    R2, O2 = size(mat)
    s = Session(C, O2)
    for i=1:R2
        tmp = filter(x -> x != 0, mat[i, :])
        push!(s.route, Route(i, mat[i, :], ntuple(i -> tmp[i], length(tmp))))
    end
    compute_output!(s)
    mail_nb = [count(x -> x != 0, mat[r, :]) for r=1:R2]
    mail_vol = [mat[i, j] for i=1:R2 for j=1:O2 if mat[i, j] != 0]
    println(s)
    println("avg mail number = $(round(mean(mail_nb), digits=3))")
    println("std mail number = $(round(std(mail_nb), digits=3))")
    println("min mail number = $(minimum(mail_nb)), max mail number = $(maximum(mail_nb))")
    println("avg mail volume = $(round(mean(mail_vol), digits=3))")
    println("std mail volume = $(round(std(mail_vol), digits=3))")
    println("min mail volume = $(minimum(mail_vol)), max mail volume = $(maximum(mail_vol))")
    # println("$(sort(mail_nb))")
    # println("$(sort(mail_vol))")

    x_number = minimum(mail_nb):maximum(mail_nb)
    y_number = [count(x->x==i, mail_nb) for i=x_number]
    p_number = plot(xlabel = "mail number", ylabel = "route number", title="mail per route distribution", xlims=(minimum(mail_nb)-.5, maximum(mail_nb)+.5))
    bar!(x_number, y_number, label="")

    x_vol = minimum(mail_vol):maximum(mail_vol)
    y_vol = [count(x->x==i, mail_vol) for i=x_vol]
    p_volume = plot(xlabel = "mail sizes", ylabel = "occurrence", title="distribution of mail volume", xlims=(minimum(mail_vol)-.5, maximum(mail_vol)+0.5))
    bar!(x_vol, y_vol, label="")

    display(plot(p_number, p_volume, layout=(2,1), size = (1250, 1600)))
end

begin
    I           = 100
    R           = 20
    target1      = "../data/HardIndus/"
    target2      = "../data/Contained/"
    target3      = "../data/Skewed/"
    target4      = "../data/Chunk/"
    for O in [200]
        for S in [1, 2, 5, 10, 20]
            writeInstanceBatterie_Distrib(target1, I, S, R, O, 65, 2, ceil(Int64, 2O/3), 1, 30, distrib_1onN)
            writeInstanceBatterie_Contained(target2, I, S, R, O, 75, 2, ceil(Int64, 2O/3), 4, 9, round(Int64, O/2):round(Int64, 3O/4))
            writeInstanceBatterie_Skewed(target3, I, S, R, O, 60, 2, ceil(Int64, 2O/3), 1:5, 40:44, 1:6)
            writeInstanceBatterie_Chunk(target4, I, S, R, O, 95, 2, ceil(Int64, 2O/3), 10:10, 1:20, 1:20)
        end
    end
end


## ============================================================================================================== ##
##                 #######        ##    #######  ########   #######               ####    ######                  ##
##                 ##    ##      ##    ##        ##        ##         ######     ## ##   ##   ###                 ##
##                 #######      ##      ######   #####      ######              ##  ##   ## ## ##                 ##
##                 ##  ##      ##            ##  ##              ##   ######   ########  ###   ##                 ##
##                 ##   ##    ##       #######   ########  #######                  ##    ######                  ##
## ============================================================================================================== ##

begin
    O = 200 
    n = 20
    nbGaussian = round(Int64, O/2):round(Int64, 3O/4)
    mat = generateMultipleGaussianSum(O, n, nbGaussian)
    x = 1:O
    p = plot()
    for i = 1:n
        plot!(x, mat[i, :] .+ (i*1.2), label="")
    end 
    display(p)
end

# USAGE: uncomment 1 instance generator and the corresponding parser:
# -> generate a new instance under: "../data/Tests"
# -> parse instance and display: std/mean volume/mail number + plots...

begin
    I           = 1
    S           = 1
    O           = 40
    R           = 40
    target      = "../data/Tests"
    # writeInstanceBatterie_Distrib(target, I, S, R, O, 120, 2, ceil(Int64, 2O/3), 1, 30, distrib_1onN)
    # writeInstanceBatterie_Contained(target, I, S, R, O, 150, 2, ceil(Int64, 2O/3), 4, 19, round(Int64, O/2):round(Int64, 3O/4))
    writeInstanceBatterie_Skewed(target, I, S, R, O, 120, 2, ceil(Int64, 2O/3), 1:5, 40:44, 1:7)
    # writeInstanceBatterie_Chunk(target, I, S, R, O, 205, 2, ceil(Int64, 2O/3), 10:10, 1:20, 1:25)

    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceIndus_1_O$(O)_R$(R)_C120_opt_$(S).txt")
    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceContained_1_O$(O)_R$(R)_C150_opt_$(S).txt")
    C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceSkewed_1_O$(O)_R$(R)_C120_opt_$(S).txt")
    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceChunk_1_O$(O)_R$(R)_C205_opt_$(S).txt")
    R2, O2 = size(mat)
    s = Session(C, O2)
    for i=1:R2
        tmp = filter(x -> x != 0, mat[i, :])
        push!(s.route, Route(i, mat[i, :], ntuple(i -> tmp[i], length(tmp))))
    end
    compute_output!(s)
    mail_nb = [count(x -> x != 0, mat[r, :]) for r=1:R2]
    mail_vol = [mat[i, j] for i=1:R2 for j=1:O2 if mat[i, j] != 0]
    println(s)
    println("avg mail number = $(round(mean(mail_nb), digits=3))")
    println("std mail number = $(round(std(mail_nb), digits=3))")
    println("min mail number = $(minimum(mail_nb)), max mail number = $(maximum(mail_nb))")
    println("avg mail volume = $(round(mean(mail_vol), digits=3))")
    println("std mail volume = $(round(std(mail_vol), digits=3))")
    println("min mail volume = $(minimum(mail_vol)), max mail volume = $(maximum(mail_vol))")
    # println("$(sort(mail_nb))")
    # println("$(sort(mail_vol))")

    x_number = minimum(mail_nb):maximum(mail_nb)
    y_number = [count(x->x==i, mail_nb) for i=x_number]
    p_number = plot(xlabel = "mail number", ylabel = "route number", title="mail per route distribution", xlims=(minimum(mail_nb)-.5, maximum(mail_nb)+.5))
    bar!(x_number, y_number, label="")

    x_vol = minimum(mail_vol):maximum(mail_vol)
    y_vol = [count(x->x==i, mail_vol) for i=x_vol]
    p_volume = plot(xlabel = "mail sizes", ylabel = "occurrence", title="distribution of mail volume", xlims=(minimum(mail_vol)-.5, maximum(mail_vol)+0.5))
    bar!(x_vol, y_vol, label="")

    display(plot(p_number, p_volume, layout=(2,1), size = (1250, 1600)))
end

begin
    I           = 100
    R           = 40
    target1      = "../data/HardIndus/"
    target2      = "../data/Contained/"
    target3      = "../data/Skewed/"
    target4      = "../data/Chunk/"
    for O in [80]
        for S in [1, 2, 5, 10, 20]
            writeInstanceBatterie_Distrib(target1, I, S, R, O, 120, 2, ceil(Int64, 2O/3), 1, 30, distrib_1onN)
            writeInstanceBatterie_Contained(target2, I, S, R, O, 150, 2, ceil(Int64, 2O/3), 4, 19, round(Int64, O/2):round(Int64, 3O/4))
            writeInstanceBatterie_Skewed(target3, I, S, R, O, 120, 2, ceil(Int64, 2O/3), 1:5, 40:44, 1:7)
            writeInstanceBatterie_Chunk(target4, I, S, R, O, 205, 2, ceil(Int64, 2O/3), 10:10, 1:20, 1:25)
        end
    end
end

begin
    I           = 100
    for R in [100]
        for O in [20, 40, 80, 120]
            writeInstanceBatterie_Distrib_easy(I, R, O)
        end
    end
end
## ============================================================================================================== ##
##                 #######        ##    #######  ########   #######             ######    ######                  ##
##                 ##    ##      ##    ##        ##        ##         ######   ##    ##  ##   ###                 ##
##                 #######      ##      ######   #####      ######              ######   ## ## ##                 ##
##                 ##  ##      ##            ##  ##              ##   ######   ##    ##  ###   ##                 ##
##                 ##   ##    ##       #######   ########  #######              ######    ######                  ##
## ============================================================================================================== ##

begin
    O = 200 
    n = 20
    nbGaussian = round(Int64, O/2):round(Int64, 3O/4)
    mat = generateMultipleGaussianSum(O, n, nbGaussian)
    x = 1:O
    p = plot()
    for i = 1:n
        plot!(x, mat[i, :] .+ (i*1.2), label="")
    end 
    display(p)
end

# USAGE: uncomment 1 instance generator and the corresponding parser:
# -> generate a new instance under: "../data/Tests"
# -> parse instance and display: std/mean volume/mail number + plots...

begin
    I           = 1
    S           = 1
    O           = 80
    R           = 80
    target      = "../data/Tests"
    # writeInstanceBatterie_Distrib(target, I, S, R, O, 240, 2, ceil(Int64, 2O/3), 1, 30, distrib_1onN)
    # writeInstanceBatterie_Contained(target, I, S, R, O, 300, 2, ceil(Int64, 2O/3), 4, 38, round(Int64, O/2):round(Int64, 3O/4))
    # writeInstanceBatterie_Skewed(target, I, S, R, O, 300, 2, ceil(Int64, 2O/3), 1:5, 40:44, 1:7)
    writeInstanceBatterie_Chunk(target, I, S, R, O, 405, 2, ceil(Int64, 2O/3), 10:10, 1:20, 1:25)

    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceIndus_1_O$(O)_R$(R)_C240_opt_$(S).txt")
    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceContained_1_O$(O)_R$(R)_C300_opt_$(S).txt")
    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceSkewed_1_O$(O)_R$(R)_C300_opt_$(S).txt")
    C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceChunk_1_O$(O)_R$(R)_C405_opt_$(S).txt")
    R2, O2 = size(mat)
    s = Session(C, O2)
    for i=1:R2
        tmp = filter(x -> x != 0, mat[i, :])
        push!(s.route, Route(i, mat[i, :], ntuple(i -> tmp[i], length(tmp))))
    end
    compute_output!(s)
    mail_nb = [count(x -> x != 0, mat[r, :]) for r=1:R2]
    mail_vol = [mat[i, j] for i=1:R2 for j=1:O2 if mat[i, j] != 0]
    # println(s)
    println("avg mail number = $(round(mean(mail_nb), digits=3))")
    println("std mail number = $(round(std(mail_nb), digits=3))")
    println("min mail number = $(minimum(mail_nb)), max mail number = $(maximum(mail_nb))")
    println("avg mail volume = $(round(mean(mail_vol), digits=3))")
    println("std mail volume = $(round(std(mail_vol), digits=3))")
    println("min mail volume = $(minimum(mail_vol)), max mail volume = $(maximum(mail_vol))")
    # println("$(sort(mail_nb))")
    # println("$(sort(mail_vol))")

    x_number = minimum(mail_nb):maximum(mail_nb)
    y_number = [count(x->x==i, mail_nb) for i=x_number]
    p_number = plot(xlabel = "mail number", ylabel = "route number", title="mail per route distribution", xlims=(minimum(mail_nb)-.5, maximum(mail_nb)+.5))
    bar!(x_number, y_number, label="")

    x_vol = minimum(mail_vol):maximum(mail_vol)
    y_vol = [count(x->x==i, mail_vol) for i=x_vol]
    p_volume = plot(xlabel = "mail sizes", ylabel = "occurrence", title="distribution of mail volume", xlims=(minimum(mail_vol)-.5, maximum(mail_vol)+0.5))
    bar!(x_vol, y_vol, label="")

    display(plot(p_number, p_volume, layout=(2,1), size = (1250, 1600)))
end

begin
    I           = 100
    R           = 80
    target1      = "../data/HardIndus/"
    target2      = "../data/Contained/"
    target3      = "../data/Skewed/"
    target4      = "../data/Chunk/"
    for O in [80, 120]
        for S in [1, 2, 5, 10]
            writeInstanceBatterie_Distrib(target1, I, S, R, O, 240, 2, ceil(Int64, 2O/3), 1, 30, distrib_1onN)
            writeInstanceBatterie_Contained(target2, I, S, R, O, 300, 2, ceil(Int64, 2O/3), 4, 38, round(Int64, O/2):round(Int64, 3O/4))
            writeInstanceBatterie_Skewed(target3, I, S, R, O, 300, 2, ceil(Int64, 2O/3), 1:5, 40:44, 1:7)
            writeInstanceBatterie_Chunk(target4, I, S, R, O, 405, 2, ceil(Int64, 2O/3), 10:10, 1:20, 1:25)
        end
    end
end

## ============================================================================================================== ##
##            #######        ##    #######  ########   #######             /###      ######    ######             ##
##            ##    ##      ##    ##        ##        ##         ######    # ##           ##  ##   ###            ##
##            #######      ##      ######   #####      ######                ##       ####    ## ## ##            ##
##            ##  ##      ##            ##  ##              ##   ######      ##     ##        ###   ##            ##
##            ##   ##    ##       #######   ########  #######              ######    ######    ######             ##
## ============================================================================================================== ##

begin
    O = 200 
    n = 20
    nbGaussian = round(Int64, O/2):round(Int64, 3O/4)
    mat = generateMultipleGaussianSum(O, n, nbGaussian)
    x = 1:O
    p = plot()
    for i = 1:n
        plot!(x, mat[i, :] .+ (i*1.2), label="")
    end 
    display(p)
end

# USAGE: uncomment 1 instance generator and the corresponding parser:
# -> generate a new instance under: "../data/Tests"
# -> parse instance and display: std/mean volume/mail number + plots...

begin
    I           = 1
    S           = 1
    O           = 120
    R           = 120
    target      = "../data/Tests"
    # writeInstanceBatterie_Distrib(target, I, S, R, O, 400, 2, ceil(Int64, 2O/3), 1, 30, distrib_1onN)
    # writeInstanceBatterie_Contained(target, I, S, R, O, 450, 2, ceil(Int64, 2O/3), 4, 60, round(Int64, O/2):round(Int64, 3O/4))
    # writeInstanceBatterie_Skewed(target, I, S, R, O, 450, 2, ceil(Int64, 2O/3), 1:5, 40:44, 1:7)
    writeInstanceBatterie_Chunk(target, I, S, R, O, 605, 2, ceil(Int64, 2O/3), 10:10, 1:20, 1:25)

    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceIndus_1_O$(O)_R$(R)_C400_opt_$(S).txt")
    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceContained_1_O$(O)_R$(R)_C450_opt_$(S).txt")
    # C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceSkewed_1_O$(O)_R$(R)_C450_opt_$(S).txt")
    C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceChunk_1_O$(O)_R$(R)_C605_opt_$(S).txt")
    C, mat, min_session = parseMyInstance_completed("../data/Tests/instanceChunk_1_O$(O)_R$(R)_C605_opt_$(S).txt")
    R2, O2 = size(mat)
    s = Session(C, O2)
    for i=1:R2
        tmp = filter(x -> x != 0, mat[i, :])
        push!(s.route, Route(i, mat[i, :], ntuple(i -> tmp[i], length(tmp))))
    end
    compute_output!(s)
    mail_nb = [count(x -> x != 0, mat[r, :]) for r=1:R2]
    mail_vol = [mat[i, j] for i=1:R2 for j=1:O2 if mat[i, j] != 0]
    # println(s)
    println("avg mail number = $(round(mean(mail_nb), digits=3))")
    println("std mail number = $(round(std(mail_nb), digits=3))")
    println("min mail number = $(minimum(mail_nb)), max mail number = $(maximum(mail_nb))")
    println("avg mail volume = $(round(mean(mail_vol), digits=3))")
    println("std mail volume = $(round(std(mail_vol), digits=3))")
    println("min mail volume = $(minimum(mail_vol)), max mail volume = $(maximum(mail_vol))")
    # println("$(sort(mail_nb))")
    # println("$(sort(mail_vol))")

    x_number = minimum(mail_nb):maximum(mail_nb)
    y_number = [count(x->x==i, mail_nb) for i=x_number]
    p_number = plot(xlabel = "mail number", ylabel = "route number", title="mail per route distribution", xlims=(minimum(mail_nb)-.5, maximum(mail_nb)+.5))
    bar!(x_number, y_number, label="")

    x_vol = minimum(mail_vol):maximum(mail_vol)
    y_vol = [count(x->x==i, mail_vol) for i=x_vol]
    p_volume = plot(xlabel = "mail sizes", ylabel = "occurrence", title="distribution of mail volume", xlims=(minimum(mail_vol)-.5, maximum(mail_vol)+0.5))
    bar!(x_vol, y_vol, label="")

    display(plot(p_number, p_volume, layout=(2,1), size = (1250, 1600)))
end

begin
    I           = 100
    R           = 120
    target1      = "../data/HardIndus/"
    target2      = "../data/Contained/"
    target3      = "../data/Skewed/"
    target4      = "../data/Chunk/"
    for O in [120]
        for S in [1, 2, 5]
            writeInstanceBatterie_Distrib(target1, I, S, R, O, 400, 2, ceil(Int64, 2O/3), 1, 30, distrib_1onN)
            writeInstanceBatterie_Contained(target2, I, S, R, O, 450, 2, ceil(Int64, 2O/3), 4, 60, round(Int64, O/2):round(Int64, 3O/4))
            writeInstanceBatterie_Skewed(target3, I, S, R, O, 450, 2, ceil(Int64, 2O/3), 1:5, 40:44, 1:7)
            writeInstanceBatterie_Chunk(target4, I, S, R, O, 605, 2, ceil(Int64, 2O/3), 10:10, 1:20, 1:25)
        end
    end
end

begin
    I           = 1
    S           = 1
    O           = 120
    R           = 120

    all_mail_distrib = zeros(Int64, 120)
    all_mail_sizes = zeros(Int64, 200)

    avg_mail_nb = 0.
    std_mail_nb = 0.
    min_mail_nb = 0.
    max_mail_nb = 0.

    avg_mail_vol = 0.
    std_mail_vol = 0.
    min_mail_vol = 0.
    max_mail_vol = 0.

    for i=1:100
        instance, nbSession = parseAnyInstance("myTrafic_$(i)_O$(O)_R$(R).txt")
        C = instance.Lmax
        mat = zeros(Int64, instance.nbRoute, instance.nbOut)
        for r = 1:instance.nbRoute
            mat[r, :] = instance.route[r].assignment
        end
        # C, mat, min_session = parseMyInstance_completed("../data/HardIndus/instanceIndus_$(i)_O$(O)_R$(R)_C400_opt_$(S).txt")
        # C, mat, min_session = parseMyInstance_completed("../data/Contained/instanceContained_$(i)_O$(O)_R$(R)_C450_opt_$(S).txt")
        # C, mat, min_session = parseMyInstance_completed("../data/Skewed/instanceSkewed_$(i)_O$(O)_R$(R)_C450_opt_$(S).txt")
        # C, mat, min_session = parseMyInstance_completed("../data/Chunk/instanceChunk_$(i)_O$(O)_R$(R)_C605_opt_$(S).txt")
        R2, O2 = size(mat)
        s = Session(C, O2)
        for i=1:R2
            tmp = filter(x -> x != 0, mat[i, :])
            push!(s.route, Route(i, mat[i, :], ntuple(i -> tmp[i], length(tmp))))
        end
        compute_output!(s)
        mail_nb = [count(x -> x != 0, mat[r, :]) for r=1:R2]
        mail_vol = [mat[i, j] for i=1:R2 for j=1:O2 if mat[i, j] != 0]

        for e in mail_nb
            all_mail_distrib[e] += 1
        end
        for e in mail_vol
            all_mail_sizes[e] += 1
        end

        avg_mail_nb += mean(mail_nb)
        std_mail_nb += std(mail_nb)
        min_mail_nb += minimum(mail_nb)
        max_mail_nb += maximum(mail_nb)

        avg_mail_vol += mean(mail_vol)
        std_mail_vol += std(mail_vol)
        min_mail_vol += minimum(mail_vol)
        max_mail_vol += maximum(mail_vol)
    end

    all_mail_distrib = [(k, e) for (k, e) in enumerate(all_mail_distrib) if e != 0]
    all_mail_sizes = [(k, e) for (k, e) in enumerate(all_mail_sizes) if e != 0]

    avg_mail_nb = round(avg_mail_nb/100, digits=3)
    std_mail_nb = round(std_mail_nb/100, digits=3)
    min_mail_nb = round(min_mail_nb/100, digits=3)
    max_mail_nb = round(max_mail_nb/100, digits=3)

    avg_mail_vol = round(avg_mail_vol/100, digits=3)
    std_mail_vol = round(std_mail_vol/100, digits=3)
    min_mail_vol = round(min_mail_vol/100, digits=3)
    max_mail_vol = round(max_mail_vol/100, digits=3)

    # println(s)
    println("avg mail number = $(avg_mail_nb)")
    println("std mail number = $(std_mail_nb)")
    println("min mail number = $(min_mail_nb), max mail number = $(max_mail_nb)")
    println("avg mail volume = $(avg_mail_vol)")
    println("std mail volume = $(std_mail_vol)")
    println("min mail volume = $(min_mail_vol), max mail volume = $(max_mail_vol)")
    println("$(join(all_mail_distrib, " "))")
    println("$(join(all_mail_sizes, " "))")
end
