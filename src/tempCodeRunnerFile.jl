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
    
    include("OptiMove.jl")
    
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
    J = [4, 5, 4, 4, 4, 1, 3, 4, 2, 3, 4, 4, 3, 3, 3, 6, 2, 1, 5, 2]
    # J = [rand(1:B[r]) for r=1:length(R)]
    V = [r.mail[J[k]] for (k, r) in enumerate(R)]

    m = maximum(V)
    b = maximum([B[r] - J[r] for r=1:length(R)])

    w1 = zeros(Float64, length(R))
    w2 = zeros(Float64, length(R))

    for r=1:length(R)
        obj1_1 = (.2/(2 * length(R))) * ((B[r]-J[r]) / (b+1))
        obj1_2 = (.0/(2 * length(R))) * (V[r]/(m+1))^2

        w1[r] += V[r] + obj1_1 + obj1_2
        
        obj2_1 = (.61) * ((B[r]-J[r]) / (1+sum(B-J)))
        obj2_2 = (.39) * ((V[r]^2) / (1+sum(V.^2)))

        w2[r] += V[r] + obj2_1 + obj2_2
    end

    println(sortperm(w1))
    println(sortperm(w2))

    (sortperm(w1) == sortperm(w2)) ? print("o") : print("-")
end

function rebuild1Session(inst_name::String = "instanceSkewed_1_O200_R20_C60_opt_1.txt"; alg=0)
    tmp = VERBOSE
    global VERBOSE = true

    instance, nbSession = parseAnyInstance(inst_name)

    glob_s = Session(instance.Lmax, instance.route)
    glob_tl = 10
    glob_env = Gurobi.Env()

    if alg == 0
        glob_s, tag, Acceptance_ratio = SimulatedAnnealing(glob_s)

        println("Acceptance_ratio -> $Acceptance_ratio")
    elseif alg == 3
        glob_s, tag = rebuildSession_knapSack_model_V3!(glob_s, glob_tl, glob_env)
    else
        glob_s, tag = rebuildSession_knapSack_model_V4!(glob_s, glob_tl, glob_env)
    end
    println_verbose("$glob_s")
    println_verbose("VALID = $tag", ANSI_cyan)
    global VERBOSE = tmp
    return glob_s, tag
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

function EM_GR(instance_name::String = "instanceSkewed_1_O200_R100_C60_opt_5.txt")
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

    start = time()
    sol_BFD_1D_GreadyRebuild, flag_BFD_1D = BFD_1D_into_GreadyRebuild(instance.route, Lmax, O, R)
    time_BFD_1D_GreadyRebuild  = round(time() - start, digits=3)

    start = time()
    sol_BP_1D_GreadyRebuild, flag_BP_1D = BP_1D_into_GreadyRebuild(instance.route, Lmax, O, R)
    time_BP_1D_GreadyRebuild  = round(time() - start, digits=3)

    println("$(ANSI_cyan) ==========< Solution BFD EMPTY-MOVE     ($(length(sol_BFD_EmptyMove)) Session(s) in $(time_BFD_EmptyMove)s, valid?= $(isSolutionValid(instance, sol_BFD_EmptyMove, false))) >========== $(ANSI_reset)")
    println("$sol_BFD_EmptyMove\n\n")

    println("$(ANSI_cyan) ==========< Solution BFD GREADY REBUILD ($(length(sol_BFD_GreadyRebuild)) Session(s) in $(time_BFD_GreadyRebuild)s, valid?= $(isSolutionValid(instance, sol_BFD_GreadyRebuild, false))) >========== $(ANSI_reset)")
    println("$sol_BFD_GreadyRebuild\n\n")

    println("$(ANSI_cyan) ==========< Solution BFD-1D into GREADY REBUILD ($(length(sol_BFD_1D_GreadyRebuild)) Session(s) in $(time_BFD_1D_GreadyRebuild)s, valid?= $(flag_BFD_1D ? isSolutionValid(instance, sol_BFD_1D_GreadyRebuild, false) : "x") - $(flag_BFD_1D)) >========== $(ANSI_reset)")
    println("$sol_BFD_1D_GreadyRebuild\n\n")

    println("$(ANSI_cyan) ==========< Solution BP-1D into GREADY REBUILD ($(length(sol_BP_1D_GreadyRebuild)) Session(s) in $(time_BP_1D_GreadyRebuild)s, valid?= $(flag_BP_1D ? isSolutionValid(instance, sol_BP_1D_GreadyRebuild, false) : "x") - $(flag_BP_1D)) >========== $(ANSI_reset)")
    println("$sol_BP_1D_GreadyRebuild\n\n")
end

function getStdInstance(inst_name::String = "instanceChunk_1_O200_R20_C100_opt_1.txt")
    i, _ = parseAnyInstance(inst_name)

    all_mail = []
    for r in i.route
        all_mail = [all_mail; collect(r.mail)]
    end

    return std(all_mail)
end

using Plots

function RunSA(inst_name::String = "instanceSkewed_1_O200_R20_C60_opt_1.txt")
    tmp = VERBOSE
    global VERBOSE = true

    instance, nbSession = parseAnyInstance(inst_name)

    glob_s = shuffle!(Session(instance.Lmax, instance.route))

    glob_s, tag, Acceptance_ratio, _ = SimulatedAnnealing(glob_s)

    println("Acceptance_ratio: $(Acceptance_ratio)")

    println_verbose("$glob_s")
    println_verbose("VALID = $tag", ANSI_cyan)
    global VERBOSE = tmp

    # x = 1:length(Acceptance_ratio)
    # plot(x, Acceptance_ratio)

    return glob_s, tag
end

function RunSA_V2(inst_name::String = "instanceIndus_1_O200_R30_C40_opt_1.txt")#"instanceContained_1_O200_R20_C150_opt_1.txt")
    instance, nbSession = parseAnyInstance(inst_name)

    glob_s = shuffle!(Session(instance.Lmax, instance.route))

    s, flag, res = SimulatedAnnealing_V2(glob_s)

    println_verbose("$s")
    println_verbose("VALID = $flag", ANSI_cyan)

    return s, flag, res
end

function RunSA_V2_mergeRatio(inst_name = "instanceSkewed_1_O200_R20_C80_opt_1.txt")
    #inst_name = "instanceContained_1_O200_R20_C150_opt_1.txt"

    instance, nbSession = parseAnyInstance(inst_name)
    glob_s = shuffle!(Session(instance.Lmax, instance.route))

    sum_sol = 0
    best_sol = nothing
    for i=1:10
        
        s, flag, res = SimulatedAnnealing_V2(deepcopy(glob_s), display_plot=false)

        sum_sol += fitness(s, OverloadVolume)
        (best_sol === nothing || best_sol > fitness(s, OverloadVolume)) && (best_sol = fitness(s, OverloadVolume))
        println("resolution: $i/10 -> $(fitness(s, OverloadVolume))")
    end
    println("merge ratio=$(best_sol/(sum_sol/10))")
end

