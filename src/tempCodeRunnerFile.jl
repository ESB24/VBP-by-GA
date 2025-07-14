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

function test_BFD_EM()
    instance_name::String = "instanceSkewed_1_O200_R100_C60_opt_5.txt"
    instance, nbSession = parseAnyInstance(instance_name) # trafic_200_200_35_1.xlsx
    Lmax = instance.Lmax
    O = instance.nbOut
    R = instance.nbRoute

    start = time()

    res = BFD_SAV3_EM(instance.route, Lmax, O, R)

    return res, time() - start 
end

function test_FFD_EM()
    instance_name::String = "instanceSkewed_1_O200_R100_C60_opt_5.txt"
    instance, nbSession = parseAnyInstance(instance_name) # trafic_200_200_35_1.xlsx
    Lmax = instance.Lmax
    O = instance.nbOut
    R = instance.nbRoute

    start = time()
    
    res = FFD_SAV3_EM(instance.route, Lmax, O, R)

    return res, time() - start 
end

function test_NFD_EM()
    instance_name::String = "instanceSkewed_3_O200_R100_C60_opt_5.txt"
    instance, nbSession = parseAnyInstance(instance_name) # trafic_200_200_35_1.xlsx
    Lmax = instance.Lmax
    O = instance.nbOut
    R = instance.nbRoute

    start = time()
    
    res = NFD_SAV3_EM(instance.route, Lmax, O, R)

    return res, time() - start 
end

function test_WFD_EM()
    instance_name::String = "instanceSkewed_1_O200_R20_C60_opt_1.txt"
    instance, nbSession = parseAnyInstance(instance_name) # trafic_200_200_35_1.xlsx
    Lmax = instance.Lmax
    O = instance.nbOut
    R = instance.nbRoute

    start = time()
    
    res = WFD_SAV3_EM(instance.route, Lmax, O, R)

    return res, time() - start 
end

function test_BFD_EM(instance_name::String = "instanceIndus_1_O200_R600_C40_opt_20.txt")
    instance, nbSession = parseAnyInstance(instance_name) # trafic_200_200_35_1.xlsx
    Lmax = instance.Lmax
    O = instance.nbOut
    R = instance.nbRoute

    start = time()
    
    res = BFD_EM(instance.route, Lmax, O, R)

    return res, time() - start 
end

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

function RunSA_V2(inst_name::String = "instanceSkewed_1_O200_R20_C80_opt_1.txt")#"instanceContained_1_O200_R20_C150_opt_1.txt")
    instance, nbSession = parseAnyInstance(inst_name)

    glob_s = shuffle!(Session(instance.Lmax, instance.route))

    s, flag, res = SimulatedAnnealing_V2(glob_s, display_plot=true, display_state=true)

    
    println_verbose("$s")
    println_verbose("VALID = $flag", ANSI_cyan)

    println("sqrt overload volume = $(fitness(s, OverloadVolume))")

    return s, flag, res
end

function RunSA_V3(inst_name::String = "instanceContained_1_O20_R20_C150_opt_1.txt")#"instanceContained_1_O200_R20_C150_opt_1.txt")
    instance, nbSession = parseAnyInstance(inst_name)

    glob_s = shuffle!(Session(instance.Lmax, instance.route))

    s, flag, res = SimulatedAnnealing_V3(glob_s, display_plot=true, display_state=true)

    
    println_verbose("$s")
    println_verbose("VALID = $flag", ANSI_cyan)

    println("sqrt overload volume = $(fitness(s, OverloadVolume))")

    return s, flag, res
end

function RunSA_V4(inst_name::String = "instanceSkewed_1_O40_R40_C120_opt_1.txt")#"instanceContained_1_O200_R20_C150_opt_1.txt")
    instance, nbSession = parseAnyInstance(inst_name)

    glob_s = shuffle!(Session(instance.Lmax, instance.route))

    s1, s2, flag1, flag2, res = SA_V4(glob_s, display_plot=true, display_state=true)

    return s1, s2, flag1, flag2, res
end

function RunSA_V4_all()
    instances = ["myTrafic_$(i)_O80_R80.txt" for i=1:100]

    result_path1::String = "../data/res_myTrafic_O80_R80.txt"

    fd1 = open(result_path1, "w+")

    for (k, i) in enumerate(instances)
        s1, s2, flag1, flag2, res = RunSA_V4(i)

        write(fd1, "### instance: $i\n")
        write(fd1, " - load: $(s1.load)\n")
        write(fd1, " - std: $(fitness(s1, LoadSTD))\n")
        write(fd1, " - matrix:\n")
        for r in s1.route
            write(fd1, " $(rpad(r.id, 3)) : $(join([rpad(e, 2) for e in r.assignment], ", "))\n")
        end
        write(fd1, "\n\n")
        flush(fd1)
    end
end

function RunSA_V2_mergeRatio(inst_name = "instanceChunk_1_O200_R20_C100_opt_1.txt")
    # "instanceChunk_1_O200_R20_C100_opt_1.txt"
    # "instanceIndus_1_O200_R30_C40_opt_1.txt"
    # "instanceContained_1_O200_R20_C150_opt_1.txt"
    # "instanceSkewed_1_O200_R20_C80_opt_1.txt"

    instance, nbSession = parseAnyInstance(inst_name)
    glob_s = shuffle!(Session(instance.Lmax, instance.route))

    sum_sol = 0
    best_sol = nothing
    run = 10
    for i=1:run
        
        s, flag, res = SimulatedAnnealing_V2(deepcopy(glob_s), display_plot=true, display_state=true)

        sum_sol += fitness(s, OverloadVolume)
        (best_sol === nothing || best_sol > fitness(s, OverloadVolume)) && (best_sol = fitness(s, OverloadVolume))
        println("\nresolution: $i/$run -> $(fitness(s, OverloadVolume))\n\n\n")
    end
    println("merge ratio=$(best_sol/(sum_sol/run))\navg solution obj = $((sum_sol/run))")
end

function RunSA_V3_mergeRatio(inst_name = "instanceIndus_1_O200_R30_C40_opt_1.txt")
    # "instanceIndus_1_O200_R30_C40_opt_1.txt"
    # "instanceContained_1_O200_R20_C150_opt_1.txt"
    # "instanceChunk_1_O200_R20_C100_opt_1.txt"
    # "instanceSkewed_1_O200_R20_C80_opt_1.txt"

    instance, nbSession = parseAnyInstance(inst_name)
    glob_s = shuffle!(Session(instance.Lmax, instance.route))

    sum_sol = 0
    best_sol = nothing
    run = 10
    for i=1:run
        
        s, flag, res = SimulatedAnnealing_V3(deepcopy(glob_s), display_plot=true, display_state=true, α=.96, τ=0.5, obj=LoadSTD)

        sum_sol += fitness(s, OverloadVolume)
        (best_sol === nothing || best_sol > fitness(s, OverloadVolume)) && (best_sol = fitness(s, OverloadVolume))
        println("\nresolution: $i/$run -> $(fitness(s, OverloadVolume))\n\n\n")
    end
    println("merge ratio=$(best_sol/(sum_sol/run))\navg solution obj = $((sum_sol/run))")
end

function Tuning_attr_α_SA_V3()
    range_α = 0.89:0.01:0.99

    plt = plot(
        1,
        xlim = ((range_α)[begin], (range_α)[end]),
        xlabel = "α",
        # zlim = (minimum(res), maximum(res)),
        ylabel = "average obj",
        title = "α parameter",
        marker = 2,
    )

    display(plt)

    res = zeros(Float64, length(range_α))
    best_sol = nothing
    
    run = 4 * 5
    
    instances = [
                ["instanceIndus_$(i)_O200_R30_C40_opt_1.txt" for i=1:round(Int64, run/4)];
                ["instanceContained_$(i)_O200_R20_C150_opt_1.txt" for i=1:round(Int64, run/4)];
                ["instanceSkewed_$(i)_O200_R20_C80_opt_1.txt" for i=1:round(Int64, run/4)];
                ["instanceChunk_$(i)_O200_R20_C100_opt_1.txt" for i=1:round(Int64, run/4)]
            ]

    for (i_α, α) in enumerate(range_α)
        for i=1:run
            instance, nbSession = parseAnyInstance(instances[i])
            glob_s = shuffle!(Session(instance.Lmax, instance.route))

            s, flag, _ = SimulatedAnnealing_V3(deepcopy(glob_s), display_plot=false, display_state=false, α=α, τ=.91)

            res[i_α] += fitness(s, OverloadVolume)
            (best_sol === nothing || best_sol > fitness(s, OverloadVolume)) && (best_sol = fitness(s, OverloadVolume))
            println("resolution: $i/$run -> $(fitness(s, OverloadVolume)) with α=$α")
        end
        println("merge ratio=$(best_sol/(res[i_α]/run))\navg solution obj = $((res[i_α]/run))")
            
        push!(plt, α, round(res[i_α]/run, digits=3))
        display(plt)
    end

    # build an animated gif by pushing new points to the plot, saving every 10th frame
    tmp = round.(res ./ run, digits=3)
    return tmp, join(enumerate(tmp), " -- ")
end

function Tuning_attr_SA_V3()
    range_τ = 0.2:0.1:2.
    range_α = 0.89:0.01:.99

    res1 = zeros(Float64, length(range_α), length(range_τ))
    res2 = zeros(Float64, length(range_α), length(range_τ))
    best1 = nothing
    best2 = nothing

    result_path1::String = "../data/results1.txt"
    result_path2::String = "../data/results2.txt"

    fd1 = open(result_path1, "a")
    fd2 = open(result_path2, "a")

    write(fd1, "\n\n < std >\n(→) τ = $(collect(range_τ))\n(↓) α = $(collect(range_α))\n")
    write(fd2, "\n\n < overflow >\n(→) τ = $(collect(range_τ))\n(↓) α = $(collect(range_α))\n")
    
    run = 10 * 5
    
    instances = [
                ["myTrafic_$(i)_O120_R120.txt" for i=1:round(Int64, run/5)];
                ["instanceIndus_$(i)_O120_R120_C400_opt_1.txt" for i=1:round(Int64, run/5)];
                ["instanceContained_$(i)_O120_R120_C450_opt_1.txt" for i=1:round(Int64, run/5)];
                ["instanceSkewed_$(i)_O120_R120_C450_opt_1.txt" for i=1:round(Int64, run/5)];
                ["instanceChunk_$(i)_O120_R120_C605_opt_1.txt" for i=1:round(Int64, run/5)]
            ]

    for (i_α, α) in enumerate(range_α)
        for (i_τ, τ) in enumerate(range_τ)
            println("====================< τ=$τ and α=$α >====================")
            for i=1:run
                instance, nbSession = parseAnyInstance(instances[i])
                glob_s = shuffle!(Session(instance.Lmax, instance.route))

                s, _, flag, _, _ = SA_V4(deepcopy(glob_s), display_plot=false, display_state=false, α=α, τ=τ, obj1=LoadSTD)

                res1[i_α, i_τ] += fitness(s, LoadSTD)
                res2[i_α, i_τ] += fitness(s, OverloadVolume)

                (best1 === nothing || best1 > fitness(s, LoadSTD)) && (best1 = fitness(s, LoadSTD))
                (best2 === nothing || best2 > fitness(s, OverloadVolume)) && (best2 = fitness(s, OverloadVolume))

                println(" - run $i/$run -> std = $(fitness(s, LoadSTD)), overload = $(fitness(s, OverloadVolume))")
            end
            println("(obj 1) -> merge ratio = $(best1/(res1[i_α, i_τ]/run))\n avg solution obj = $((res1[i_α, i_τ]/run))")
            println("(obj 2) -> merge ratio = $(best2/(res2[i_α, i_τ]/run))\n avg solution obj = $((res2[i_α, i_τ]/run))")
            
            write(fd1, "$(round(res1[i_α, i_τ]/run, digits=3)), ") # $(τ) $(α) 
            flush(fd1)
            write(fd2, "$(round(res2[i_α, i_τ]/run, digits=3)), ") # $(τ) $(α) 
            flush(fd2)
        end
        write(fd1, "\n") # $(τ) $(α) 
        flush(fd1)
        write(fd2, "\n") # $(τ) $(α) 
        flush(fd2)
    end

    # build an animated gif by pushing new points to the plot, saving every 10th frame

    println("res1: max = $(maximum(res1)), min = $(minimum(res1))")
    println("res2: max = $(maximum(res2)), min = $(minimum(res2))")
    
    return res1, res2
end

function Tuning_attr_move_SA_V3()
    range_move = 50:50:400

    display()

    res = zeros(Float64, length(range_move))
    best_sol = nothing
    
    run = 4 * 2
    
    instances = [
                ["instanceIndus_$(i)_O200_R30_C40_opt_1.txt" for i=1:round(Int64, run/4)];
                ["instanceContained_$(i)_O200_R20_C150_opt_1.txt" for i=1:round(Int64, run/4)];
                ["instanceSkewed_$(i)_O200_R20_C80_opt_1.txt" for i=1:round(Int64, run/4)];
                ["instanceChunk_$(i)_O200_R20_C100_opt_1.txt" for i=1:round(Int64, run/4)]
            ]

    for (i_move, move) in enumerate(range_move)
        for i=1:run
            instance, nbSession = parseAnyInstance(instances[i])
            glob_s = shuffle!(Session(instance.Lmax, instance.route))

            s, flag, _ = SimulatedAnnealing_V3(deepcopy(glob_s), display_plot=false, display_state=true, Emax=move, Mmax=move)

            res[i_move] += fitness(s, OverloadVolume)
            (best_sol === nothing || best_sol > fitness(s, OverloadVolume)) && (best_sol = fitness(s, OverloadVolume))
            println("\nresolution: $i/$run -> $(fitness(s, OverloadVolume)) with move=$move\n")
        end
        println("merge ratio=$(best_sol/(res[i_move]/run))\navg solution obj = $((res[i_move]/run))")
    end
end

function RunGR_V2_mergeRatio(inst_name = "instanceChunk_1_O200_R20_C100_opt_1.txt")
    # "instanceChunk_1_O200_R20_C100_opt_1.txt"
    # "instanceIndus_1_O200_R30_C40_opt_1.txt"
    # "instanceContained_1_O200_R20_C150_opt_1.txt"
    # "instanceSkewed_1_O200_R20_C80_opt_1.txt"

    instance, nbSession = parseAnyInstance(inst_name)
    glob_s = shuffle!(Session(instance.Lmax, instance.route))
    tl = 10
    env = Gurobi.Env()

    sum_sol = 0
    best_sol = nothing
    run = 10
    for i=1:run
        s, added = rebuildSession_knapSack_model_V3!(deepcopy(glob_s), tl, env)

        sum_sol += fitness(s, OverloadVolume)
        (best_sol === nothing || best_sol > fitness(s, OverloadVolume)) && (best_sol = fitness(s, OverloadVolume))
        println("\nresolution: $i/$run -> $(fitness(s, OverloadVolume))\n\n\n")
    end
    println("merge ratio=$(best_sol/(sum_sol/run))\navg solution obj = $((sum_sol/run))")
end

function attr_bench_GR_V3()
    range_d1 = .71:.01:.89
    range_d2 = .81:.01:.99

    run_per_inst = 10

    instances = [
                ["instanceIndus_$(i)_O200_R30_C40_opt_1.txt" for i=1:run_per_inst];
                ["instanceContained_$(i)_O200_R20_C150_opt_1.txt" for i=1:run_per_inst];
                ["instanceSkewed_$(i)_O200_R20_C80_opt_1.txt" for i=1:run_per_inst];
                ["instanceChunk_$(i)_O200_R20_C100_opt_1.txt" for i=1:run_per_inst]
            ]

    tl = 10
    env = Gurobi.Env()

    res_opti = zeros(Int64, length(range_d1), length(range_d2))

    for (i1, d1) in enumerate(range_d1)
        for (i2, d2) in enumerate(range_d2)
            start = time()
            print("δ = [$d1, $d2] : ")
            for instance_path in instances
                instance, nbSession = parseAnyInstance(instance_path)

                s, added = rebuildSession_knapSack_model_V3!(Session(instance.Lmax, instance.route), tl, env, [d1, d2])

                added && (res_opti[i1, i2] += 1)
            end
            println("res $(res_opti[i1, i2]), in $(round(time() - start, digits=3))")
        end
    end

    println("Gready rebuild V3:")
    println("δ1 ∈ $(collect(range_d1))")
    println("δ2 ∈ $(collect(range_d2))")
    for (i, d1) in enumerate(range_d1)
        println("δ1 = $d1 -> $(res_opti[i, :])")
    end
end

function attr_bench_GR_V4()
    range_d1 = .71:.01:.89
    range_d2 = .21:.01:.31

    run_per_inst = 10

    instances = [
                ["instanceIndus_$(i)_O200_R200_C40_opt_1.txt" for i=1:run_per_inst];
                ["instanceContained_$(i)_O200_R200_C150_opt_1.txt" for i=1:run_per_inst];
                ["instanceSkewed_$(i)_O200_R200_C80_opt_1.txt" for i=1:run_per_inst];
                ["instanceChunk_$(i)_O200_R200_C100_opt_1.txt" for i=1:run_per_inst]
            ]

    tl = 10
    env = Gurobi.Env()

    res_opti = zeros(Int64, length(range_d1), length(range_d2))

    for (i1, d1) in enumerate(range_d1)
        for (i2, d2) in enumerate(range_d2)
            start = time()
            print("δ = [$d1, $d2] : ")
            for instance_path in instances
                instance, nbSession = parseAnyInstance(instance_path)

                s, added = rebuildSession_knapSack_model_V4!(Session(instance.Lmax, instance.route), tl, env, [d1, d2])

                added && (res_opti[i1, i2] += 1)
            end
            println("res $(res_opti[i1, i2]), in $(round(time() - start, digits=3))")
        end
    end

    println("Gready rebuild V4:")
    println("δ1 ∈ $(collect(range_d1))")
    println("δ2 ∈ $(collect(range_d2))")
    for (i, d1) in enumerate(range_d1)
        println("δ1 = $d1 -> $(res_opti[i, :])")
    end
end

function get_plot()
    fd2 = open("../data/Plot/SA_res_iter.dat", "w+")

    pad_iter::Int64 = ceil(Int64, log10(length(resI[1, :]))+1)

    write(fd2, "$(rpad("m", pad_iter)) s*     sc     sb    \n")

    for i=1:length(resI[1, :])
        write(fd2, "$(rpad(i, pad_iter)) ")
        write(fd2, "$(rpad(round(resI[1, i], digits=3), 6)) ")
        write(fd2, "$(rpad(round(resI[2, i], digits=3), 6)) ")
        write(fd2, "$(rpad(round(resI[3, i], digits=3), 6)) \n")
    end

    flush(fd2)

    close(fd2)

end
