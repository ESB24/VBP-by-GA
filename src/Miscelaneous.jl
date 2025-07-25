begin
    include("Instance.jl")
    include("Bound.jl")
end

using Base.Filesystem

# ================================================================================= #
#                               ######  #     # ##    #                             #
#                               #     # #     # # #   #                             #
#                               ######  #     # #  #  #                             #
#                               #  #    #     # #   # #                             #
#                               #   #    #####  #    ##                             #
# ================================================================================= #

function run_01LP(instance::Instance, s::Solution, tl::Int64 = 900)
    model, sol = model_01LP_warmstart(instance, s, tl)

    (sol === nothing) || (MOI.get(model, Gurobi.ModelAttribute("SolCount")) <= 0) && (throw(MyCustomError("Error with the warmuped 01LP -> no sol found?!")))

    ObjVal      ::Int64     = MOI.get(model, Gurobi.ModelAttribute("SolCount")) > 0 ? MOI.get(model, Gurobi.ModelAttribute("ObjVal")) : -1 
    flag_opti   ::Bool      = termination_status(model) == MOI.OPTIMAL
    Gap         ::Float64   = MOI.get(model, Gurobi.ModelAttribute("MIPGap"))
    Runtime     ::Float64   = MOI.get(model, Gurobi.ModelAttribute("Runtime"))

    return ObjVal, flag_opti, Runtime, Gap, sol
end

function run_GA(instance::Instance, ObjBound::Int64, tl::Int64 = 600)

    max_time        ::Float64                       = float(tl)
    # max_pop         ::Int64                         = 100
    # λ_core          ::Float64                       = 0.5
    # λ_elite         ::Float64                       = 0.4
    # max_gen         ::Int64                         = 1000
    # TAG_FitGA       ::Type{<:FitnessSolution}       = ObjGA_2
    # TAG_Select      ::Type{<:SelectParent}          = SelectElite
    # TAG_Crossover   ::Tuple                         = (Edit, Point1, PointN)
    # TAG_Accept      ::Type{<:AcceptSon}             = Always
    # maxStagnate     ::Int64                         = 10

    best_sol::Solution, flag_opti::Bool, gen_val, Runtime::Float64 = GA(instance, max_time=max_time)
    #     instance, 
    #     max_time=max_time, 
    #     max_pop=max_pop, 
    #     max_gen=max_gen,
    #     λ_core=λ_core,
    #     λ_elite=λ_elite,
    #     TAG_FitGA=TAG_FitGA,
    #     TAG_Select=TAG_Select,
    #     TAG_Crossover=TAG_Crossover,
    #     TAG_Accept=TAG_Accept,
    #     maxStagnate=maxStagnate
    # )

    ObjVal::Int64 = length(best_sol.sessions)
    Gap::Float64 = abs(ObjBound - ObjVal) / abs(ObjVal)

    return ObjVal, flag_opti, Runtime, Gap, best_sol
end

function run_Build(instance::Instance, ObjBound::Int64)
    perm            ::Vector{Int64}                     = sortperm([(-fitness(r, MaxMin), -fitness(r, NbBatches)) for r in instance.rounds])
    TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD

    start_time  ::Float64   = time()
    s           ::Solution  = buildSolution_BFD_final(instance, perm)
    
    ObjVal      ::Int64     = length(s.sessions)
    flag_opti   ::Bool      = ObjVal == ObjBound
    RunTime     ::Float64   = time() - start_time
    Gap         ::Float64   = abs(ObjBound - ObjVal) / abs(ObjVal)

    return ObjVal, flag_opti, RunTime, Gap, s
end

# ================================================================================= #
#                               #######  #####  ######                              #
#                                  #    #     # #     #                             #
#                                  #    ####### ######                              #
#                                  #    #     # #     #                             #
#                                  #    #     # #######                             #
# ================================================================================= #

function parseAnyInstanceMat(path::String, λ::Float64 = 0.25)::Tuple{Union{Matrix{Int64}, Nothing}, Int64}
    println(" ====================================================================================================")
    println("                                       $(path)")
    println(" ====================================================================================================")

    instance::Union{Instance, Nothing} = nothing
    if cmp(path[1:6], "trafic") == 0
        path = "TER/data/Trafic/" * path
        mat = parseMailXLSX(path)

        return mat, minSession(instance)

    elseif cmp(path[1:8], "myTrafic") == 0
        path = "TER/data/MyTrafic/" * path
        mat = parseMyInstance_easy(path)

        return mat, 1

    elseif cmp(path[1:9], "instance_") == 0
        path = "TER/data/Hard/" * path
        numberOfOutputs, numberOfRounds, numberOfSessions, sessions, Lmax = parseOptiInstance(path)
        mat = vect2D_to_matrix(normalizeOptiInstance(sessions))
        
        return mat, numberOfSessions

    elseif cmp(path[1:10], "instanceDi") == 0 || cmp(path[1:10], "instanceBi") == 0 || cmp(path[1:10], "instanceGa") == 0
        if cmp(path[1:10], "instanceDi") == 0
            path = "TER/data/Distrib/" * path
        elseif cmp(path[1:10], "instanceBi") == 0
            path = "TER/data/BigBatch/" * path
        elseif cmp(path[1:10], "instanceGa") == 0
            path = "TER/data/Gaussian/" * path
        end
        Lmax, mat, nbSession = parseMyInstance(path)

        return mat, nbSession
    end

    return mat, -1
end

function parseAnyInstance(path::String, λ::Float64 = 0.25)::Tuple{Union{Instance, Nothing}, Int64}
    #println(" ====================================================================================================")
    #println("                                       $(path)")
    #println(" ====================================================================================================")

    instance::Union{Instance, Nothing} = nothing
    if cmp(path[1:6], "trafic") == 0
        path = "../data/Trafic/" * path
        mat = parseMailXLSX(path)
        instance = Instance(mat, λ)

        global delta1 = 0.2
        global delta2 = 0.9

        return instance, minSession(instance)

    elseif cmp(path[1:8], "myTrafic") == 0
        path = "../data/MyTrafic(Indus.)/" * path
        mat = parseMyInstance_easy(path)
        instance = Instance(mat, λ)

        global delta1 = 0.2
        global delta2 = 0.9

        return instance, minSession(instance)

    elseif cmp(path[1:9], "instance_") == 0
        path = "../data/Hard/" * path
        numberOfOutputs, numberOfRounds, numberOfSessions, sessions, Lmax = parseOptiInstance(path)
        mat = vect2D_to_matrix(normalizeOptiInstance(sessions))
        instance = Instance(mat, λ)
        instance.Lmax = Lmax

        global delta1 = 0.2
        global delta2 = 0.9
        
        return instance, numberOfSessions

    else
        if cmp(path[1:10], "instanceCh") == 0
            path = "../data/Chunk/" * path

            global delta1 = 0.2
            global delta2 = 0.0
        elseif cmp(path[1:10], "instanceIn") == 0
            path = "../data/HardIndus/" * path
            global delta1 = 0.2
            global delta2 = 0.8
        elseif cmp(path[1:10], "instanceSk") == 0
            path = "../data/Skewed/" * path
            global delta1 = 0.2
            global delta2 = 0.9
        elseif cmp(path[1:10], "instanceCo") == 0
            path = "../data/Contained/" * path
            global delta1 = 0.4
            global delta2 = 0.9
        end
        Lmax, mat, nbSession = parseMyInstance(path)
        instance = Instance(mat, λ)
        instance.Lmax = Lmax

        return instance, nbSession
    end

    return instance, -1
end

function batterie_GR3_GR4(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
            (("Indus. like", 30, 200), ["instanceIndus_$(i)_O200_R30_C40_opt_1.txt" for i=1:10]),
            # (("Indus. like", 60, 200), ["instanceIndus_$(i)_O200_R60_C40_opt_2.txt" for i=1:5]),
            (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R150_C40_opt_5.txt" for i=1:10]),
            (("Contained", 20, 200), ["instanceContained_$(i)_O200_R20_C150_opt_1.txt" for i=1:10]),
            # (("Contained", 40, 200), ["instanceContained_$(i)_O200_R40_C150_opt_2.txt" for i=1:5]),
            (("Contained", 100, 200), ["instanceContained_$(i)_O200_R100_C150_opt_5.txt" for i=1:10]),
            (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R20_C80_opt_1.txt" for i=1:10]),
            # (("Skewed", 40, 200), ["instanceSkewed_$(i)_O200_R40_C80_opt_2.txt" for i=1:5]),
            (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R100_C80_opt_5.txt" for i=1:10]),
            (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R20_C100_opt_1.txt" for i=1:10]),
            # (("Chunk", 40, 200), ["instanceChunk_$(i)_O200_R40_C100_opt_2.txt" for i=1:5]),
            (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R100_C100_opt_5.txt" for i=1:10]),
            (("Indus", 200, 200), ["myTrafic_$(i)_O200_R200.txt" for i=1:10]),
            (("Indus", 800, 200), ["myTrafic_$(i)_O200_R800.txt" for i=1:10]),
            ]   ,
        result_path         ::String            = "../data/results.txt",
        tl::Int64 = 10,
        env::Gurobi.Env = Gurobi.Env(),
    )
    # ==========< Head >==========
    fd = open(result_path, "a")

    write(fd, "\n\n    \\begin{table}[h]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\begin{tabular}{c c c c|c c c c|c c c c}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{4}{l|}{Instance} & \\multicolumn{4}{l|}{BFD (\\textsc{GREADY-REBUILD-V3})} & \\multicolumn{4}{l}{BFD (\\textsc{GREADY-REBUILD-V4})} \\\\\n")
    write(fd, "            type & \$|O|\$ & \$|R|\$ & \\textsc{OPTI} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")
    # ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanlb  ::Float64 = 0.

        # EMPTY-MOVE
        EM_opti ::Int64   = 0
        EM_cpu  ::Float64 = 0.
        EM_gap  ::Float64 = 0.
        EM_obj  ::Float64 = 0.

        # GREADY REBUILD
        GR_opti ::Int64   = 0
        GR_cpu  ::Float64 = 0.
        GR_gap  ::Float64 = 0.
        GR_obj  ::Float64 = 0.

        for path in group_path
            instance::Instance, _ = parseAnyInstance(path)
            
            Lmax        = instance.Lmax
            O           = instance.nbOut
            R           = instance.nbRoute
            lb::Int64   = minSession(instance)

            meanlb += lb  
    # ==========< solve EMPTY-MOVE >==========

            start = time()
            sol_BFD_EmptyMove = BFD_GreadyRebuild_V3(instance.route, Lmax, O, R, tl, env)

            EM_opti += (length(sol_BFD_EmptyMove) == lb) ? (1) : (0)
            EM_cpu += time() - start
            EM_gap += abs(lb - length(sol_BFD_EmptyMove)) / abs(lb)
            EM_obj += length(sol_BFD_EmptyMove)

    # ==========< GREADY REBUILD >==========
                
            start = time()
            sol_BFD_GreadyRebuild = BFD_GreadyRebuild_V4(instance.route, Lmax, O, R, tl, env)

            GR_opti += (length(sol_BFD_GreadyRebuild) == lb) ? (1) : (0)
            GR_cpu += time() - start
            GR_gap += abs(lb - length(sol_BFD_GreadyRebuild)) / abs(lb)
            GR_obj += length(sol_BFD_GreadyRebuild)

        end

    # ==========< Compute res >==========

        meanlb /= length(group_path)

        EM_cpu  = round(EM_cpu / length(group_path), digits=3)
        EM_gap  = round(EM_gap / length(group_path), digits=3)
        EM_obj  = round(EM_obj / length(group_path), digits=3)

        # GREADY REBUILD
        GR_cpu  = round(GR_cpu / length(group_path), digits=3)
        GR_gap  = round(GR_gap / length(group_path), digits=3)
        GR_obj  = round(GR_obj / length(group_path), digits=3)

    # ==========< Write line >==========

    write(fd, "$(name) & $(O) & $(R) & $(meanlb)") # Head
    write(fd, " & $(EM_opti) & $(EM_cpu) & $(EM_gap) & $(EM_obj)") # EMPTY-MOVE
    write(fd, " & $(GR_opti) & $(GR_cpu) & $(GR_gap) & $(GR_obj)") # GREADY REBUILD
    write(fd, "\\\\\n")
    flush(fd)

    # ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        \\caption{TODO}\n")
    write(fd, "        \\label{table:res}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterie_SA2_SA3(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
            (("Indus. like", 30, 200), ["instanceIndus_$(i)_O200_R30_C40_opt_1.txt" for i=1:10]),
            # (("Indus. like", 60, 200), ["instanceIndus_$(i)_O200_R60_C40_opt_2.txt" for i=1:5]),
            (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R150_C40_opt_5.txt" for i=1:10]),
            (("Contained", 20, 200), ["instanceContained_$(i)_O200_R20_C150_opt_1.txt" for i=1:10]),
            # (("Contained", 40, 200), ["instanceContained_$(i)_O200_R40_C150_opt_2.txt" for i=1:5]),
            (("Contained", 100, 200), ["instanceContained_$(i)_O200_R100_C150_opt_5.txt" for i=1:10]),
            (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R20_C80_opt_1.txt" for i=1:10]),
            # (("Skewed", 40, 200), ["instanceSkewed_$(i)_O200_R40_C80_opt_2.txt" for i=1:5]),
            (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R100_C80_opt_5.txt" for i=1:10]),
            (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R20_C100_opt_1.txt" for i=1:10]),
            # (("Chunk", 40, 200), ["instanceChunk_$(i)_O200_R40_C100_opt_2.txt" for i=1:5]),
            (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R100_C100_opt_5.txt" for i=1:10]),
            (("Indus", 200, 200), ["myTrafic_$(i)_O200_R200.txt" for i=1:10]),
            (("Indus", 800, 200), ["myTrafic_$(i)_O200_R800.txt" for i=1:10]),
            ]   ,
        result_path         ::String            = "../data/results.txt",
    )
    # ==========< Head >==========
    fd = open(result_path, "a")

    write(fd, "\n\n    \\begin{table}[h]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\begin{tabular}{c c c c|c c c c|c c c c}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{4}{l|}{Instance} & \\multicolumn{4}{l|}{BFD (\\textsc{Simulated-Annealing-V2})} & \\multicolumn{4}{l}{BFD (\\textsc{Simulated-Annealing-V3})} \\\\\n")
    write(fd, "            type & \$|O|\$ & \$|R|\$ & \\textsc{OPTI} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")
    # ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanlb  ::Float64 = 0.

        # EMPTY-MOVE
        EM_opti ::Int64   = 0
        EM_cpu  ::Float64 = 0.
        EM_gap  ::Float64 = 0.
        EM_obj  ::Float64 = 0.

        # GREADY REBUILD
        GR_opti ::Int64   = 0
        GR_cpu  ::Float64 = 0.
        GR_gap  ::Float64 = 0.
        GR_obj  ::Float64 = 0.

        for path in group_path
            instance::Instance, _ = parseAnyInstance(path)
            
            Lmax        = instance.Lmax
            O           = instance.nbOut
            R           = instance.nbRoute
            lb::Int64   = minSession(instance)

            meanlb += lb  
    # ==========< solve EMPTY-MOVE >==========

            start = time()
            sol_BFD_EmptyMove = FFD_EmptyMove(instance.route, Lmax, O, R)

            EM_opti += (length(sol_BFD_EmptyMove) == lb) ? (1) : (0)
            EM_cpu += time() - start
            EM_gap += abs(lb - length(sol_BFD_EmptyMove)) / abs(lb)
            EM_obj += length(sol_BFD_EmptyMove)

    # ==========< GREADY REBUILD >==========
                
            start = time()
            sol_BFD_GreadyRebuild = BFD_SAV3_EM(instance.route, Lmax, O, R)

            GR_opti += (length(sol_BFD_GreadyRebuild) == lb) ? (1) : (0)
            GR_cpu += time() - start
            GR_gap += abs(lb - length(sol_BFD_GreadyRebuild)) / abs(lb)
            GR_obj += length(sol_BFD_GreadyRebuild)

        end

    # ==========< Compute res >==========

        meanlb /= length(group_path)

        EM_cpu  = round(EM_cpu / length(group_path), digits=3)
        EM_gap  = round(EM_gap / length(group_path), digits=3)
        EM_obj  = round(EM_obj / length(group_path), digits=3)

        # GREADY REBUILD
        GR_cpu  = round(GR_cpu / length(group_path), digits=3)
        GR_gap  = round(GR_gap / length(group_path), digits=3)
        GR_obj  = round(GR_obj / length(group_path), digits=3)

    # ==========< Write line >==========

    write(fd, "$(name) & $(O) & $(R) & $(meanlb)") # Head
    write(fd, " & $(EM_opti) & $(EM_cpu) & $(EM_gap) & $(EM_obj)") # EMPTY-MOVE
    write(fd, " & $(GR_opti) & $(GR_cpu) & $(GR_gap) & $(GR_obj)") # GREADY REBUILD
    write(fd, "\\\\\n")
    flush(fd)

    # ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        \\caption{TODO}\n")
    write(fd, "        \\label{table:res}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterie_loadSmoothing_EM_SA(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
            # (("Indus. like", 30, 200), ["instanceIndus_$(i)_O20_R20_C65_opt_1.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O20_R20_C75_opt_1.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O20_R20_C60_opt_1.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O20_R20_C95_opt_1.txt" for i=1:10]),
            (("Indus", 20, 200), ["myTrafic_$(i)_O20_R20.txt" for i=1:10]),

            (("Indus. like", 30, 200), ["instanceIndus_$(i)_O40_R20_C65_opt_1.txt" for i=1:10]),
            (("Contained", 20, 200), ["instanceContained_$(i)_O40_R20_C75_opt_1.txt" for i=1:10]),
            (("Skewed", 20, 200), ["instanceSkewed_$(i)_O40_R20_C60_opt_1.txt" for i=1:10]),
            (("Chunk", 20, 200), ["instanceChunk_$(i)_O40_R20_C95_opt_1.txt" for i=1:10]),
            (("Indus", 20, 200), ["myTrafic_$(i)_O40_R20.txt" for i=1:10]),

            # (("Indus. like", 30, 200), ["instanceIndus_$(i)_O40_R40_C120_opt_1.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O40_R40_C150_opt_1.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O40_R40_C120_opt_1.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O40_R40_C205_opt_1.txt" for i=1:10]),
            (("Indus", 20, 200), ["myTrafic_$(i)_O40_R40.txt" for i=1:10]),

            (("Indus. like", 30, 200), ["instanceIndus_$(i)_O80_R40_C120_opt_1.txt" for i=1:10]),
            (("Contained", 20, 200), ["instanceContained_$(i)_O80_R40_C150_opt_1.txt" for i=1:10]),
            (("Skewed", 20, 200), ["instanceSkewed_$(i)_O80_R40_C120_opt_1.txt" for i=1:10]),
            (("Chunk", 20, 200), ["instanceChunk_$(i)_O80_R40_C205_opt_1.txt" for i=1:10]),
            (("Indus", 20, 200), ["myTrafic_$(i)_O80_R40.txt" for i=1:10]),

            # (("Indus. like", 30, 200), ["instanceIndus_$(i)_O80_R80_C240_opt_1.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O80_R80_C300_opt_1.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O80_R80_C300_opt_1.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O80_R80_C405_opt_1.txt" for i=1:10]),
            (("Indus", 20, 200), ["myTrafic_$(i)_O80_R80.txt" for i=1:10]),

            (("Indus. like", 30, 200), ["instanceIndus_$(i)_O120_R80_C240_opt_1.txt" for i=1:10]),
            (("Contained", 20, 200), ["instanceContained_$(i)_O120_R80_C300_opt_1.txt" for i=1:10]),
            (("Skewed", 20, 200), ["instanceSkewed_$(i)_O120_R80_C300_opt_1.txt" for i=1:10]),
            (("Chunk", 20, 200), ["instanceChunk_$(i)_O120_R80_C405_opt_1.txt" for i=1:10]),
            (("Indus", 20, 200), ["myTrafic_$(i)_O120_R80.txt" for i=1:10]),

            (("Indus. like", 30, 200), ["instanceIndus_$(i)_O120_R120_C400_opt_1.txt" for i=1:10]),
            (("Contained", 20, 200), ["instanceContained_$(i)_O120_R120_C450_opt_1.txt" for i=1:10]),
            (("Skewed", 20, 200), ["instanceSkewed_$(i)_O120_R120_C450_opt_1.txt" for i=1:10]),
            (("Chunk", 20, 200), ["instanceChunk_$(i)_O120_R120_C605_opt_1.txt" for i=1:10]),
            (("Indus", 20, 200), ["myTrafic_$(i)_O120_R120.txt" for i=1:10]),
            ]   ,
        result_path1        ::String                    = "../data/results1.txt",
        result_path2        ::String                    = "../data/results2.txt",
        obj                 ::Type{<:FitnessSession}    = OverloadVolume   ,
        run                 ::Int64                     = 10,
    )
# ==========< Head >==========
    fd1 = open(result_path1, "a")
    fd2 = open(result_path2, "a")

    write(fd1,"\n\n")
    write(fd2,"\n\n")

    obj1 = LoadSTD
    obj2 = MostLoadedOut

# ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        GRP_lb1 :: Float64 = 0.
        GRP_lb2 :: Float64 = 0.

        # EMPTY-MOVE
        GRP_EM_opti1    ::Int64   = 0
        GRP_EM_opti2    ::Int64   = 0
        GRP_EM_cpu      ::Float64 = 0.
        GRP_EM_obj1     ::Float64 = 0.
        GRP_EM_obj2     ::Float64 = 0.
        GRP_EM_improv1  ::Float64 = 0.
        GRP_EM_improv2  ::Float64 = 0.
        GRP_EM_gap1     ::Float64 = 0.
        GRP_EM_gap2     ::Float64 = 0.
        GRP_EM_gap3     ::Float64 = 0.
        GRP_EM_gap4     ::Float64 = 0.
        GRP_EM_merge1   ::Float64 = 0.
        GRP_EM_merge2   ::Float64 = 0.
        GRP_EM_var1     ::Float64 = 0.
        GRP_EM_var2     ::Float64 = 0.
        GRP_EM_all1     ::Vector{Float64} = zeros(Float64, run)
        GRP_EM_all2     ::Vector{Float64} = zeros(Float64, run)

        # SIMULATED ANNEALING
        GRP_SA_opti1    ::Int64   = 0
        GRP_SA_opti2    ::Int64   = 0
        GRP_SA_cpu      ::Float64 = 0.
        GRP_SA_obj1     ::Float64 = 0.
        GRP_SA_obj2     ::Float64 = 0.
        GRP_SA_improv1  ::Float64 = 0.
        GRP_SA_improv2  ::Float64 = 0.
        GRP_SA_gap1     ::Float64 = 0.
        GRP_SA_gap2     ::Float64 = 0.
        GRP_SA_gap3     ::Float64 = 0.
        GRP_SA_gap4     ::Float64 = 0.
        GRP_SA_merge1   ::Float64 = 0.
        GRP_SA_merge2   ::Float64 = 0.
        GRP_SA_var1     ::Float64 = 0.
        GRP_SA_var2     ::Float64 = 0.
        GRP_SA_all1     ::Vector{Float64} = zeros(Float64, run)
        GRP_SA_all2     ::Vector{Float64} = zeros(Float64, run)
# ==========< Init runs >==========
        println("\nGroup $k: $name")
        for (i, path) in enumerate(group_path)
            print("($i/$(length(group_path))) > ")
            instance::Instance, _ = parseAnyInstance(path)

            O = instance.nbOut
            R = instance.nbRoute

            meanlb1      ::Float64 = 0.
            meanlb2      ::Float64 = 0.

            EM_opti1    ::Int64   = 0
            EM_opti2    ::Int64   = 0
            EM_cpu      ::Float64 = 0.
            EM_obj1     ::Float64 = 0.
            EM_obj2     ::Float64 = 0.
            EM_improv1  ::Float64 = 0.
            EM_improv2  ::Float64 = 0.
            EM_gap1     ::Float64 = 0.
            EM_gap2     ::Float64 = 0.
            EM_gap3     ::Float64 = 0.
            EM_gap4     ::Float64 = 0.
            EM_merge1   ::Float64 = 0.
            EM_merge2   ::Float64 = 0.
            EM_var1     ::Float64 = 0.
            EM_var2     ::Float64 = 0.
            EM_all1     ::Vector{Float64} = zeros(Float64, run)
            EM_all2     ::Vector{Float64} = zeros(Float64, run)

            SA_opti1    ::Int64   = 0
            SA_opti2    ::Int64   = 0
            SA_cpu      ::Float64 = 0.
            SA_obj1     ::Float64 = 0.
            SA_obj2     ::Float64 = 0.
            SA_improv1  ::Float64 = 0.
            SA_improv2  ::Float64 = 0.
            SA_gap1     ::Float64 = 0.
            SA_gap2     ::Float64 = 0.
            SA_gap3     ::Float64 = 0.
            SA_gap4     ::Float64 = 0.
            SA_merge1   ::Float64 = 0.
            SA_merge2   ::Float64 = 0.
            SA_var1     ::Float64 = 0.
            SA_var2     ::Float64 = 0.
            SA_all1     ::Vector{Float64} = zeros(Float64, run)
            SA_all2     ::Vector{Float64} = zeros(Float64, run)

# ==========< Runs >==========
            
            for k=1:run
                print("<$k")
                glob_s = shuffle!(Session(instance.Lmax, instance.route))

                start_obj1 = fitness(glob_s, LoadSTD)
                start_obj2 = fitness(glob_s, MostLoadedOut)

                lb_obj1 = lb_std(glob_s)
                lb_obj2 = ceil(mean(glob_s.load))

                start = time()
                em_s1, em_s2, _ = EM_VE(deepcopy(glob_s))
                em_cpu = time() - start

                start = time()
                sa_s1, sa_s2, _, _, _ = SA_V4(deepcopy(glob_s), display_plot=false, display_state=false, α=0.96, τ=.5)
                sa_cpu = time() - start

# ==========< Res runs >==========

                meanlb1     += lb_std(glob_s)
                meanlb2     += ceil(mean(glob_s.load))

                EM_opti1    += (fitness(em_s1, obj1) == lb_obj1) ? (1) : (0) 
                EM_opti2    += (fitness(em_s2, obj2) == lb_obj2) ? (1) : (0) 
                EM_cpu      += em_cpu
                EM_obj1     += fitness(em_s1, obj1)
                EM_obj2     += fitness(em_s2, obj2)
                EM_improv1  += (abs(start_obj1) == 0) ? (0.) : abs(start_obj1 - fitness(em_s1, obj1)) / abs(start_obj1)
                EM_improv2  += (abs(start_obj2) == 0) ? (0.) : abs(start_obj2 - fitness(em_s2, obj2)) / abs(start_obj2)
                EM_gap1     += (abs(fitness(em_s1, obj1)) == 0) ? (0.) : abs(fitness(em_s1, obj1) - lb_obj1) / abs(fitness(em_s1, obj1))
                EM_gap2     += (abs(fitness(em_s2, obj2)) == 0) ? (0.) : abs(fitness(em_s2, obj2) - lb_obj2) / abs(fitness(em_s2, obj2))
                EM_gap3     += (abs(start_obj1 - lb_obj1) == 0) ? (0.) : abs(fitness(em_s1, obj1) - lb_obj1) / abs(start_obj1 - lb_obj1)
                EM_gap4     += (abs(start_obj2 - lb_obj2) == 0) ? (0.) : abs(fitness(em_s2, obj2) - lb_obj2) / abs(start_obj2 - lb_obj2)
                EM_all1[k]  = fitness(em_s1, obj1)
                EM_all2[k]  = fitness(em_s2, obj2)

                SA_opti1    += (fitness(sa_s1, obj1) == lb_obj1) ? (1) : (0) 
                SA_opti2    += (fitness(sa_s2, obj2) == lb_obj2) ? (1) : (0) 
                SA_cpu      += sa_cpu
                SA_obj1     += fitness(sa_s1, obj1)
                SA_obj2     += fitness(sa_s2, obj2)
                SA_improv1  += (abs(start_obj1) == 0) ? (0.) : abs(start_obj1 - fitness(sa_s1, obj1)) / abs(start_obj1)
                SA_improv2  += (abs(start_obj2) == 0) ? (0.) : abs(start_obj2 - fitness(sa_s2, obj2)) / abs(start_obj2)
                SA_gap1     += (abs(fitness(sa_s1, obj1)) == 0) ? (0.) : abs(fitness(sa_s1, obj1) - lb_obj1) / abs(fitness(sa_s1, obj1))
                SA_gap2     += (abs(fitness(sa_s2, obj2)) == 0) ? (0.) : abs(fitness(sa_s2, obj2) - lb_obj2) / abs(fitness(sa_s2, obj2))
                SA_gap3     += (abs(start_obj1 - lb_obj1) == 0) ? (0.) : abs(fitness(sa_s1, obj1) - lb_obj1) / abs(start_obj1 - lb_obj1)
                SA_gap4     += (abs(start_obj2 - lb_obj2) == 0) ? (0.) : abs(fitness(sa_s2, obj2) - lb_obj2) / abs(start_obj2 - lb_obj2)
                SA_all1[k]  = fitness(sa_s1, obj1)
                SA_all2[k]  = fitness(sa_s2, obj2)

                print(">")
            end

# ==========< Res group >==========

            GRP_lb1 += meanlb1 / run
            GRP_lb2 += meanlb2 / run


            GRP_EM_opti1    += EM_opti1
            GRP_EM_opti2    += EM_opti2
            GRP_EM_cpu      += EM_cpu / run
            GRP_EM_obj1     += EM_obj1 / run
            GRP_EM_obj2     += EM_obj2 / run
            GRP_EM_improv1  += EM_improv1 / run * 100
            GRP_EM_improv2  += EM_improv2 / run * 100
            GRP_EM_gap1     += EM_gap1 / run * 100
            GRP_EM_gap2     += EM_gap2 / run * 100
            GRP_EM_gap3     += EM_gap3 / run * 100
            GRP_EM_gap4     += EM_gap4 / run * 100
            GRP_EM_merge1   += (sum(EM_obj1) == 0) ? (0.) : minimum(EM_all1) / (EM_obj1 / run)
            GRP_EM_merge2   += (sum(EM_obj2) == 0) ? (0.) : minimum(EM_all2) / (EM_obj2 / run)
            GRP_EM_var1     += sum([(e - mean(EM_all1))^2 for e in EM_all1]) / (run-1)
            GRP_EM_var2     += sum([(e - mean(EM_all2))^2 for e in EM_all2]) / (run-1)

            GRP_SA_opti1    += SA_opti1
            GRP_SA_opti2    += SA_opti2
            GRP_SA_cpu      += SA_cpu / run
            GRP_SA_obj1     += SA_obj1 / run
            GRP_SA_obj2     += SA_obj2 / run
            GRP_SA_improv1  += SA_improv1 / run * 100
            GRP_SA_improv2  += SA_improv2 / run * 100
            GRP_SA_gap1     += SA_gap1 / run * 100
            GRP_SA_gap2     += SA_gap2 / run * 100
            GRP_SA_gap3     += SA_gap3 / run * 100
            GRP_SA_gap4     += SA_gap4 / run * 100
            GRP_SA_merge1   += (sum(SA_obj1) == 0) ? (0.) : (minimum(SA_all1) / (SA_obj1 / run))
            GRP_SA_merge2   += (sum(SA_obj2) == 0) ? (0.) : minimum(SA_all2) / (SA_obj2 / run)
            GRP_SA_var1     += sum([(e - mean(SA_all1))^2 for e in SA_all1]) / (run-1)
            GRP_SA_var2     += sum([(e - mean(SA_all2))^2 for e in SA_all2]) / (run-1)
            
            println()
        end

# ==========< Average group >==========

        GRP_lb1 /= length(group_path)
        GRP_lb2 /= length(group_path)

        GRP_EM_cpu      /= length(group_path)
        GRP_EM_obj1     /= length(group_path)
        GRP_EM_obj2     /= length(group_path)
        GRP_EM_improv1  /= length(group_path)
        GRP_EM_improv2  /= length(group_path)
        GRP_EM_gap1     /= length(group_path)
        GRP_EM_gap2     /= length(group_path)
        GRP_EM_gap3     /= length(group_path)
        GRP_EM_gap4     /= length(group_path)
        GRP_EM_merge1   /= length(group_path)
        GRP_EM_merge2   /= length(group_path)
        GRP_EM_var1     /= length(group_path)
        GRP_EM_var2     /= length(group_path)

        GRP_SA_cpu      /= length(group_path)
        GRP_SA_obj1     /= length(group_path)
        GRP_SA_obj2     /= length(group_path)
        GRP_SA_improv1  /= length(group_path)
        GRP_SA_improv2  /= length(group_path)
        GRP_SA_gap1     /= length(group_path)
        GRP_SA_gap2     /= length(group_path)
        GRP_SA_gap3     /= length(group_path)
        GRP_SA_gap4     /= length(group_path)
        GRP_SA_merge1   /= length(group_path)
        GRP_SA_merge2   /= length(group_path)
        GRP_SA_var1     /= length(group_path)
        GRP_SA_var2     /= length(group_path)

# ==========< Write line >==========

        # Instances
        write(fd1, "$(name) & $(O) & $(R) & $(round(GRP_lb1, digits=3)) & $(round(GRP_lb2, digits=3))")
        # Integer
        write(fd1, " & $(GRP_EM_opti1) & $(GRP_EM_opti2)")
        # Float
        write(fd1, " & $(join(round.([GRP_EM_cpu, GRP_EM_obj1, GRP_EM_obj2, GRP_EM_improv1, GRP_EM_improv2, GRP_EM_gap1, GRP_EM_gap2, GRP_EM_gap3, GRP_EM_gap4, GRP_EM_merge1, GRP_EM_merge2, GRP_EM_var1, GRP_EM_var2], digits=3), " & ")) \\\\\n")
        flush(fd1)

        # Instances
        write(fd2, "$(name) & $(O) & $(R) & $(round(GRP_lb1, digits=3)) & $(round(GRP_lb2, digits=3))")
        # Integer
        write(fd2, " & $(GRP_SA_opti1) & $(GRP_SA_opti2)")
        # Float
        write(fd2, " & $(join(round.([GRP_SA_cpu, GRP_SA_obj1, GRP_SA_obj2, GRP_SA_improv1, GRP_SA_improv2, GRP_SA_gap1, GRP_SA_gap2, GRP_SA_gap3, GRP_SA_gap4, GRP_SA_merge1, GRP_SA_merge2, GRP_SA_var1, GRP_SA_var2], digits=3), " & ")) \\\\\n")
        flush(fd2)
    end
# ==========< END >==========

    close(fd1)
    close(fd2)
end

function batterie_loadSmoothing_EM_SA_LP(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
            (("Indus", 20, 200),        ["myTrafic_$(i)_O20_R20.txt" for i=1:2]),

            (("Diverse", 30, 200),      ["instanceIndus_$(i)_O20_R20_C65_opt_1.txt" for i=1:10]),
            (("Contained", 20, 200),    ["instanceContained_$(i)_O20_R20_C75_opt_1.txt" for i=1:10]),
            (("Skewed", 20, 200),       ["instanceSkewed_$(i)_O20_R20_C60_opt_1.txt" for i=1:10]),
            (("Chunk", 20, 200),        ["instanceChunk_$(i)_O20_R20_C95_opt_1.txt" for i=1:10]),
            (("Indus", 20, 200),        ["myTrafic_$(i)_O20_R20.txt" for i=1:10]),

            (("Diverse", 30, 200),      ["instanceIndus_$(i)_O40_R20_C65_opt_1.txt" for i=1:10]),
            (("Contained", 20, 200),    ["instanceContained_$(i)_O40_R20_C75_opt_1.txt" for i=1:10]),
            (("Skewed", 20, 200),       ["instanceSkewed_$(i)_O40_R20_C60_opt_1.txt" for i=1:10]),
            (("Chunk", 20, 200),        ["instanceChunk_$(i)_O40_R20_C95_opt_1.txt" for i=1:10]),
            (("Indus", 20, 200),        ["myTrafic_$(i)_O40_R20.txt" for i=1:10]),

            (("Diverse", 30, 200),      ["instanceIndus_$(i)_O40_R40_C120_opt_1.txt" for i=1:10]),
            (("Contained", 20, 200),    ["instanceContained_$(i)_O40_R40_C150_opt_1.txt" for i=1:10]),
            (("Skewed", 20, 200),       ["instanceSkewed_$(i)_O40_R40_C120_opt_1.txt" for i=1:10]),
            (("Chunk", 20, 200),        ["instanceChunk_$(i)_O40_R40_C205_opt_1.txt" for i=1:10]),
            (("Indus", 20, 200),        ["myTrafic_$(i)_O40_R40.txt" for i=1:10]),


            (("Diverse", 30, 200),      ["instanceIndus_$(i)_O80_R40_C120_opt_1.txt" for i=1:10]),
            (("Contained", 20, 200),    ["instanceContained_$(i)_O80_R40_C150_opt_1.txt" for i=1:10]),
            (("Skewed", 20, 200),       ["instanceSkewed_$(i)_O80_R40_C120_opt_1.txt" for i=1:10]),
            (("Chunk", 20, 200),        ["instanceChunk_$(i)_O80_R40_C205_opt_1.txt" for i=1:10]),
            (("Indus", 20, 200),        ["myTrafic_$(i)_O80_R40.txt" for i=1:10]),

            # (("Indus. like", 30, 200), ["instanceIndus_$(i)_O80_R80_C240_opt_1.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O80_R80_C300_opt_1.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O80_R80_C300_opt_1.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O80_R80_C405_opt_1.txt" for i=1:10]),

            # (("Indus. like", 30, 200), ["instanceIndus_$(i)_O120_R80_C240_opt_1.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O120_R80_C300_opt_1.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O120_R80_C300_opt_1.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O120_R80_C405_opt_1.txt" for i=1:10]),

            # (("Indus", 20, 200), ["myTrafic_$(i)_O80_R80.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O120_R80.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O120_R120.txt" for i=1:10]),

            # (("Indus. like", 30, 200), ["instanceIndus_$(i)_O120_R120_C400_opt_1.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O120_R120_C450_opt_1.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O120_R120_C450_opt_1.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O120_R120_C605_opt_1.txt" for i=1:10]),
            ]   ,
        result_path1         ::String            = "../data/results1.txt",
        result_path2         ::String            = "../data/results2.txt",
        result_path3         ::String            = "../data/results3.txt",
        obj                 ::Type{<:FitnessSession} = SampleLoad_ABS_DEV   ,
    )
# ==========< Head >==========
    fd1 = open(result_path1, "a")
    fd2 = open(result_path2, "a")
    fd3 = open(result_path3, "a")

    write(fd1,"\n\n")
    write(fd2,"\n\n")
    write(fd3,"\n\n")

    obj1 = LoadSTD
    obj2 = MostLoadedOut

    env::Gurobi.Env = Gurobi.Env()

# ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanlb1  ::Float64 = 0.
        meanlb2  ::Float64 = 0.

        # EMPTY-MOVE
        EM_opti1    ::Int64   = 0
        EM_opti2    ::Int64   = 0
        EM_cpu      ::Float64 = 0.
        EM_obj1     ::Float64 = 0.
        EM_obj2     ::Float64 = 0.
        EM_improv1  ::Float64 = 0.
        EM_improv2  ::Float64 = 0.
        EM_gap1     ::Float64 = 0.
        EM_gap2     ::Float64 = 0.
        EM_gap3     ::Float64 = 0.
        EM_gap4     ::Float64 = 0.
        EM_gap5     ::Float64 = 0.
        EM_gap6     ::Float64 = 0.

        SA_opti1    ::Int64   = 0
        SA_opti2    ::Int64   = 0
        SA_cpu      ::Float64 = 0.
        SA_obj1     ::Float64 = 0.
        SA_obj2     ::Float64 = 0.
        SA_improv1  ::Float64 = 0.
        SA_improv2  ::Float64 = 0.
        SA_gap1     ::Float64 = 0.
        SA_gap2     ::Float64 = 0.
        SA_gap3     ::Float64 = 0.
        SA_gap4     ::Float64 = 0.
        SA_gap5     ::Float64 = 0.
        SA_gap6     ::Float64 = 0.

        LP_opti1    ::Int64   = 0
        LP_opti2    ::Int64   = 0
        LP_cpu      ::Float64 = 0.
        LP_obj1     ::Float64 = 0.
        LP_obj2     ::Float64 = 0.
        LP_improv1  ::Float64 = 0.
        LP_improv2  ::Float64 = 0.
        LP_gap1     ::Float64 = 0.
        LP_gap2     ::Float64 = 0.
        LP_gap3     ::Float64 = 0.
        LP_gap4     ::Float64 = 0.
        LP_gap5     ::Float64 = 0.

# ==========< Init Run >==========
        println("\nGroup $k: $name")
        for (i, path) in enumerate(group_path)
            print("($i/$(length(group_path))) > ")
            instance::Instance, _ = parseAnyInstance(path)

            O = instance.nbOut
            R = instance.nbRoute

            glob_s = shuffle!(Session(instance.Lmax, instance.route))

            start_obj1 = fitness(glob_s, LoadSTD)
            start_obj2 = fitness(glob_s, MostLoadedOut)

            lb_obj1 = lb_std(glob_s)
            lb_obj2 = ceil(mean(glob_s.load))
# ==========< Runs >==========

            start = time()
            em_s1, em_s2, _ = EM_VE(deepcopy(glob_s), obj1=SampleLoad_ABS_DEV)
            em_cpu = time() - start

            start = time()
            sa_s1, sa_s2, _, _, _ = SA_V4(deepcopy(glob_s), display_plot=false, display_state=false, α=0.96, τ=.5, obj1=SampleLoad_ABS_DEV)
            sa_cpu = time() - start

            start = time()
            lp_md, lp_s = MILP_load_ballance_V2_warmstart(deepcopy(glob_s), 600, env)
            lp_cpu = time() - start

# ==========< res run >==========
            
            meanlb1  += lb_obj1
            meanlb2  += lb_obj2
            
            EM_opti1    += (fitness(em_s1, obj1) == lb_obj1) ? (1) : (0) 
            EM_opti2    += (fitness(em_s1, obj2) == lb_obj2) ? (1) : (0) 
            EM_cpu      += em_cpu
            EM_obj1     += fitness(em_s1, obj1)
            EM_obj2     += fitness(em_s1, obj2)
            EM_improv1  += (abs(start_obj1) == 0) ? (0.) : abs(start_obj1 - fitness(em_s1, obj1)) / abs(start_obj1)
            EM_improv2  += (abs(start_obj2) == 0) ? (0.) : abs(start_obj2 - fitness(em_s1, obj2)) / abs(start_obj2)
            EM_gap1     += (abs(fitness(em_s1, obj1)) == 0) ? (0.) : abs(fitness(em_s1, obj1) - lb_obj1) / abs(fitness(em_s1, obj1))
            EM_gap2     += (abs(fitness(em_s1, obj2)) == 0) ? (0.) : abs(fitness(em_s1, obj2) - lb_obj2) / abs(fitness(em_s1, obj2))
            EM_gap3     += (abs(start_obj1 - lb_obj1) == 0) ? (0.) : abs(fitness(em_s1, obj1) - lb_obj1) / abs(start_obj1 - lb_obj1)
            EM_gap4     += (abs(start_obj2 - lb_obj2) == 0) ? (0.) : abs(fitness(em_s1, obj2) - lb_obj2) / abs(start_obj2 - lb_obj2)
            EM_gap5     += (lp_s == nothing || fitness(lp_s, obj1) == 0) ? (0.) : (fitness(em_s1, obj1) - fitness(lp_s, obj1)) / fitness(lp_s, obj1)
            EM_gap6     += (lp_s == nothing || fitness(lp_s, obj2) == 0) ? (0.) : (fitness(em_s1, obj2) - fitness(lp_s, obj2)) / fitness(lp_s, obj2)

            SA_opti1    += (fitness(sa_s1, obj1) == lb_obj1) ? (1) : (0) 
            SA_opti2    += (fitness(sa_s1, obj2) == lb_obj2) ? (1) : (0) 
            SA_cpu      += sa_cpu
            SA_obj1     += fitness(sa_s1, obj1)
            SA_obj2     += fitness(sa_s1, obj2)
            SA_improv1  += (abs(start_obj1) == 0) ? (0.) : abs(start_obj1 - fitness(sa_s1, obj1)) / abs(start_obj1)
            SA_improv2  += (abs(start_obj2) == 0) ? (0.) : abs(start_obj2 - fitness(sa_s1, obj2)) / abs(start_obj2)
            SA_gap1     += (abs(fitness(sa_s1, obj1)) == 0) ? (0.) : abs(fitness(sa_s1, obj1) - lb_obj1) / abs(fitness(sa_s1, obj1))
            SA_gap2     += (abs(fitness(sa_s1, obj2)) == 0) ? (0.) : abs(fitness(sa_s1, obj2) - lb_obj2) / abs(fitness(sa_s1, obj2))
            SA_gap3     += (abs(start_obj1 - lb_obj1) == 0) ? (0.) : abs(fitness(sa_s1, obj1) - lb_obj1) / abs(start_obj1 - lb_obj1)
            SA_gap4     += (abs(start_obj2 - lb_obj2) == 0) ? (0.) : abs(fitness(sa_s1, obj2) - lb_obj2) / abs(start_obj2 - lb_obj2)
            SA_gap5     += (lp_s == nothing || fitness(lp_s, obj1) == 0) ? (0.) : (fitness(sa_s1, obj1) - fitness(lp_s, obj1)) / fitness(lp_s, obj1)
            SA_gap6     += (lp_s == nothing || fitness(lp_s, obj2) == 0) ? (0.) : (fitness(sa_s1, obj2) - fitness(lp_s, obj2)) / fitness(lp_s, obj2)

            LP_opti1    += (lp_s == nothing) ? (0) : ((fitness(lp_s, obj1) == lb_obj1) ? (1) : (0)) 
            LP_opti2    += (lp_s == nothing) ? (0) : ((termination_status(lp_md) == OPTIMAL) ? (1) : (0)) 
            LP_cpu      += ((lp_s == nothing) || (termination_status(lp_md) != OPTIMAL)) ? (0.) : (lp_cpu)
            LP_obj1     += (lp_s == nothing) ? (0.) : (fitness(lp_s, obj1))
            LP_obj2     += (lp_s == nothing) ? (0.) : (fitness(lp_s, obj2))
            LP_improv1  += (lp_s == nothing || abs(start_obj1) == 0) ? (0.) : (abs(start_obj1 - fitness(lp_s, obj1)) / abs(start_obj1))
            LP_improv2  += (lp_s == nothing || abs(start_obj2) == 0) ? (0.) : (abs(start_obj2 - fitness(lp_s, obj2)) / abs(start_obj2))
            LP_gap1     += (lp_s == nothing || abs(fitness(lp_s, obj1)) == 0) ? (0.) : (abs(fitness(lp_s, obj1) - lb_obj1) / abs(fitness(lp_s, obj1)))
            LP_gap2     += (lp_s == nothing || abs(fitness(lp_s, obj2)) == 0) ? (0.) : (abs(fitness(lp_s, obj2) - lb_obj2) / abs(fitness(lp_s, obj2)))
            LP_gap3     += (lp_s == nothing || abs(start_obj1 - lb_obj1) == 0) ? (0.) : (abs(fitness(lp_s, obj1) - lb_obj1) / abs(start_obj1 - lb_obj1))
            LP_gap4     += (lp_s == nothing || abs(start_obj2 - lb_obj2) == 0) ? (0.) : (abs(fitness(lp_s, obj2) - lb_obj2) / abs(start_obj2 - lb_obj2))
            LP_gap5     += MOI.get(lp_md, Gurobi.ModelAttribute("MIPGap"))

            println()
        end

# ==========< res group >==========

        meanlb1      /= length(group_path)
        meanlb2      /= length(group_path)

        EM_cpu      /= length(group_path)
        EM_obj1     /= length(group_path)
        EM_obj2     /= length(group_path)
        EM_improv1  /= length(group_path) * 0.01
        EM_improv2  /= length(group_path) * 0.01
        EM_gap1     /= length(group_path) * 0.01
        EM_gap2     /= length(group_path) * 0.01
        EM_gap3     /= length(group_path) * 0.01
        EM_gap4     /= length(group_path) * 0.01
        EM_gap5     /= length(group_path) * 0.01
        EM_gap6     /= length(group_path) * 0.01

        SA_cpu      /= length(group_path)
        SA_obj1     /= length(group_path)
        SA_obj2     /= length(group_path)
        SA_improv1  /= length(group_path) * 0.01
        SA_improv2  /= length(group_path) * 0.01
        SA_gap1     /= length(group_path) * 0.01
        SA_gap2     /= length(group_path) * 0.01
        SA_gap3     /= length(group_path) * 0.01
        SA_gap4     /= length(group_path) * 0.01
        SA_gap5     /= length(group_path) * 0.01
        SA_gap6     /= length(group_path) * 0.01

        LP_cpu      /= max(LP_opti2, 1)
        LP_obj1     /= length(group_path)
        LP_obj2     /= length(group_path)
        LP_improv1  /= length(group_path) * 0.01
        LP_improv2  /= length(group_path) * 0.01
        LP_gap1     /= length(group_path) * 0.01
        LP_gap2     /= length(group_path) * 0.01
        LP_gap3     /= length(group_path) * 0.01
        LP_gap4     /= length(group_path) * 0.01
        LP_gap5     /= length(group_path) * 0.01

# ==========< Write line >==========

    write(fd1, "$(name) & $(O) & $(R) & $(round(meanlb1, digits=3)) & $(EM_opti1) & $(round(EM_cpu, digits=3)) & $(round(EM_obj1, digits=3)) & $(round(EM_improv1, digits=3)) & $(round(EM_gap1, digits=3)) & $(round(EM_gap3, digits=3)) & $(round(EM_gap5, digits=3)) \\\\\n") # Head
    flush(fd1)

    write(fd2, "$(name) & $(O) & $(R) & $(round(meanlb1, digits=3)) & $(SA_opti1) & $(round(SA_cpu, digits=3)) & $(round(SA_obj1, digits=3)) & $(round(SA_improv1, digits=3)) & $(round(SA_gap1, digits=3)) & $(round(SA_gap3, digits=3)) & $(round(SA_gap5, digits=3)) \\\\\n") # Head
    flush(fd2)

    write(fd3, "$(name) & $(O) & $(R) & $(round(meanlb1, digits=3)) & $(LP_opti1) & $(round(LP_cpu, digits=3)) & $(round(LP_obj1, digits=3)) & $(round(LP_improv1, digits=3)) & $(round(LP_gap1, digits=3)) & $(round(LP_gap3, digits=3)) & $(round(LP_gap5, digits=3)) \\\\\n") # Head
    flush(fd3)

# ==========< End >==========

    end

    close(fd1)
    close(fd2)
    close(fd3)
end

function batterie_SA_OneSession(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
            # (("Indus. like", 30, 200), ["instanceIndus_$(i)_O200_R30_C40_opt_1.txt" for i=1:100]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O200_R20_C150_opt_1.txt" for i=1:100]),
            (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R20_C80_opt_1.txt" for i=1:100]),
            (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R20_C100_opt_1.txt" for i=1:100]),
            ]   ,
        result_path         ::String            = "../data/results2.txt",
        obj                 ::Type{<:FitnessSession} = OverloadVolume   ,
    )
    # ==========< Head >==========
    fd = open(result_path, "a")

    write(fd, "\n\n    \\begin{table}[h]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\begin{tabular}{c c c|c c c c c}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{3}{l|}{Instance} & \\multicolumn{5}{l|}{(\\textsc{Simulated-Annealing-V3})}\\\\\n")
    write(fd, "            type & \$|O|\$ & \$|R|\$ & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt OBJ} & {\\tt VAR} & {\\tt MERGE} \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")
    # ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanlb  ::Float64 = 0.

        # EMPTY-MOVE
        EM_opti ::Int64   = 0
        EM_cpu  ::Float64 = 0.
        EM_obj  ::Float64 = 0.
        EM_all_obj::Vector{Float64} = zeros(Float64, length(group_path))
        EM_best ::Union{Nothing, Float64} = nothing

        println("\nGroup $k: $name")
        for (i, path) in enumerate(group_path)
            print("($i/$(length(group_path))) > ")
            instance::Instance, _ = parseAnyInstance(path)

            glob_s = shuffle!(Session(instance.Lmax, instance.route))
    # ==========< solve EMPTY-MOVE >==========

            start = time()

            s, flag, _ = SimulatedAnnealing_V3(deepcopy(glob_s), display_plot=false, display_state=false, obj=obj, α=0.96, τ=.5)
            # s, _, flag = improvedOptiMove_S1_V2!(deepcopy(glob_s))
            # s, _, flag = improvedOptiMove_V1(deepcopy(glob_s))
            EM_opti += (flag) ? (1) : (0)
            EM_cpu += time() - start
            EM_obj += fitness(s, LoadSTD)
            EM_all_obj[i] = fitness(s, LoadSTD)

            ((EM_best == nothing) || (EM_best > fitness(s, LoadSTD))) && (EM_best = fitness(s, LoadSTD))
            println()
        end

    # ==========< Compute res >==========

        meanlb /= length(group_path)

        EM_cpu  = round(EM_cpu / length(group_path), digits=3)
        EM_obj  = round(EM_obj / length(group_path), digits=3)

        EM_var = round(sum([(EM_all_obj[i] - EM_obj)^2 for i=length(EM_all_obj)]) / (length(EM_all_obj) - 1), digits=3)
        println("EM_all_obj = $(EM_all_obj)")
        EM_merge = round(EM_best / EM_obj, digits=3)

    # ==========< Write line >==========

    write(fd, "$(name) & $(O) & $(R)") # Head
    write(fd, " & $(EM_opti) & $(EM_cpu) & $(EM_obj) & $(EM_var) & $(EM_merge)") # EMPTY-MOVE
    write(fd, "\\\\\n")
    flush(fd)

    # ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        \\caption{TODO}\n")
    write(fd, "        \\label{table:res}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterie_MILP(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
            # (("Indus", 20, 200), ["myTrafic_$(i)_O20_R20.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O20_R40.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O40_R20.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O40_R40.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O20_R80.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O40_R80.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O80_R20.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O80_R40.txt" for i=1:10]),

            # (("Diverse", 30, 200), ["instanceIndus_$(i)_O20_R20_C65_opt_1.txt" for i=1:10]),
            # (("Diverse", 30, 200), ["instanceIndus_$(i)_O20_R40_C65_opt_2.txt" for i=1:10]),
            # (("Diverse", 30, 200), ["instanceIndus_$(i)_O40_R20_C65_opt_1.txt" for i=1:10]),
            # (("Diverse", 30, 200), ["instanceIndus_$(i)_O40_R40_C65_opt_2.txt" for i=1:10]),
            # (("Diverse", 30, 200), ["instanceIndus_$(i)_O20_R100_C65_opt_5.txt" for i=1:10]),
            # (("Diverse", 30, 200), ["instanceIndus_$(i)_O80_R40_C65_opt_2.txt" for i=1:10]),

            # (("Contained", 20, 200), ["instanceContained_$(i)_O20_R20_C75_opt_1.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O20_R40_C75_opt_2.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O40_R20_C75_opt_1.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O40_R40_C75_opt_2.txt" for i=1:10]),
            
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O20_R20_C60_opt_1.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O20_R40_C60_opt_2.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O40_R20_C60_opt_1.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O40_R40_C60_opt_2.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O20_R100_C60_opt_5.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O40_R20_C60_opt_1.txt" for i=1:10]),
            
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O20_R20_C95_opt_1.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O20_R40_C95_opt_2.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O40_R20_C95_opt_1.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O40_R40_C95_opt_2.txt" for i=1:10]),
            ]   ,
        result_path         ::String            = "../data/results1.txt",
    )
# ==========< Head >==========

    fd = open(result_path, "a")

    write(fd, "\n\n    \\begin{table}[h]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\begin{tabular}{c c c c|c c c c c c c c}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{4}{l|}{Instance} & \\multicolumn{4}{l}{MILP}\\\\\n")
    write(fd, "            type & \$|O|\$ & \$|R|\$ & \\textsc{OPTI} & \\#{\\tt OPTI}\$^{mod}\$ & \\#{\\tt OPTI}\$^{sol}\$ & {\\tt CPU} & {\\tt OBJ} & {\\tt GAP}\$^{lb}\$ & {\\tt GAP}\$^{mod}\$ \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")

    # ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanlb  ::Float64 = 0.

        # EMPTY-MOVE
        EM_opti_md  ::Int64   = 0
        EM_opti_lb  ::Int64   = 0
        EM_cpu      ::Float64 = 0.
        EM_gap_md   ::Float64 = 0.
        EM_gap_lb   ::Float64 = 0.
        EM_obj      ::Float64 = 0.

        end_success = 0

        println("\nGroup $k: $name")
        for (i, path) in enumerate(group_path)
            print("($i/$(length(group_path))) > ")
            instance::Instance, _ = parseAnyInstance(path)
            
            Lmax        = instance.Lmax
            O           = instance.nbOut
            R           = instance.nbRoute
            lb::Int64   = minSession(instance)

            meanlb += lb  
    # ==========< solve EMPTY-MOVE >==========

            sol_FD_EmptyMove = nothing
            call = 0
            success = 0

            s_init, (call, success) = WFD_SA(instance.route, Lmax, O, R)

            isSolutionValid(instance, s_init) ? (println("\x1b[32m Valid \x1b[0m")) : (println("\x1b[31m Invalid \x1b[0m"))

            md, sol = model_01LP(instance, Solution(collect(1:R), s_init), 900)
            sol_s = (sol == nothing) ? (nothing) : (sol.sessions)

            println("\x1b[36m $(sol_s) \x1b[0m")

            EM_opti_md  += (termination_status(md) == OPTIMAL) ? (1) : (0) 
            EM_opti_lb  += (sol == nothing) ? (0) : ((length(sol_s) == lb) ? (1) : (0))
            EM_cpu      += (termination_status(md) == OPTIMAL) ? (MOI.get(md, Gurobi.ModelAttribute("RunTime"))) : (0) 
            EM_gap_md   += MOI.get(md, MOI.RelativeGap())
            EM_gap_lb   += (sol == nothing) ? (0) : abs(lb - length(sol_s)) / abs(lb)
            EM_obj      += (sol == nothing) ? (0) : (length(sol_s))

            println()
        end

    # ==========< Compute res >==========

        meanlb /= length(group_path)

        EM_cpu      = round(EM_cpu / length(group_path), digits=3)
        EM_gap_md   = round(EM_gap_md / length(group_path), digits=3)
        EM_gap_lb   = round(EM_gap_lb / length(group_path), digits=3)
        EM_obj      = round(EM_obj / length(group_path), digits=3)

    # ==========< Write line >==========

    write(fd, "$(name) & $(O) & $(R) & $(meanlb)") # Head
    write(fd, " & $(EM_opti_md) & $(EM_opti_lb) & $(EM_cpu) & $(EM_gap_md) & $(EM_gap_lb) & $(EM_obj) \\\\\n")
    flush(fd)

    # ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        \\caption{TODO}\n")
    write(fd, "        \\label{table:res}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterie_EM(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
            # (("Indus. like", 30, 200), ["instanceIndus_$(i)_O200_R30_C40_opt_1.txt" for i=1:10]),
            # (("Indus. like", 60, 200), ["instanceIndus_$(i)_O200_R60_C40_opt_2.txt" for i=1:5]),
            # (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R150_C40_opt_5.txt" for i=1:10]),
            # (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R300_C40_opt_10.txt" for i=1:10]),
            # (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R600_C40_opt_20.txt" for i=1:10]),

            # (("Contained", 20, 200), ["instanceContained_$(i)_O200_R20_C150_opt_1.txt" for i=1:10]),
            # (("Contained", 40, 200), ["instanceContained_$(i)_O200_R40_C150_opt_2.txt" for i=1:5]),
            # (("Contained", 100, 200), ["instanceContained_$(i)_O200_R100_C150_opt_5.txt" for i=1:10]),
            # (("Contained", 100, 200), ["instanceContained_$(i)_O200_R200_C150_opt_10.txt" for i=1:10]),
            # (("Contained", 100, 200), ["instanceContained_$(i)_O200_R400_C150_opt_20.txt" for i=1:10]),

            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R20_C60_opt_1.txt" for i=1:10]),
            # (("Skewed", 40, 200), ["instanceSkewed_$(i)_O200_R40_C80_opt_2.txt" for i=1:5]),
            # (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R100_C60_opt_5.txt" for i=1:10]),
            # (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R200_C60_opt_10.txt" for i=1:10]),
            # (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R400_C60_opt_20.txt" for i=1:10]),

            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R20_C100_opt_1.txt" for i=1:10]),
            # (("Chunk", 40, 200), ["instanceChunk_$(i)_O200_R40_C100_opt_2.txt" for i=1:5]),
            # (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R100_C100_opt_5.txt" for i=1:10]),
            # (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R200_C100_opt_10.txt" for i=1:10]),
            # (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R400_C100_opt_20.txt" for i=1:10]),

            (("Indus", 200, 200), ["myTrafic_$(i)_O40_R200.txt" for i=1:10]),
            (("Indus", 200, 200), ["myTrafic_$(i)_O100_R200.txt" for i=1:10]),
            (("Indus", 200, 200), ["myTrafic_$(i)_O200_R200.txt" for i=1:10]),
            # (("Indus", 800, 200), ["myTrafic_$(i)_O200_R800.txt" for i=1:10]),
            ]   ,
        result_path         ::String            = "../data/results.txt",
        fit_algo :: Int64 = 1,
    )
    # ==========< Head >==========

    fd = open(result_path, "a")

    algo_name_tmp = (fit_algo == 1) ? ("BFD") : ((fit_algo == 2) ? "FFD" : ((fit_algo == 3) ? "NFD" : "WFD"))

    write(fd, "\n\n    \\begin{table}[h]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\begin{tabular}{c c c c|c c c c c c c c}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{4}{l|}{Instance} & \\multicolumn{6}{l}{$(algo_name_tmp) (\\textsc{OptiMove} all 3 stages)}\\\\\n")
    write(fd, "            type & \$|O|\$ & \$|R|\$ & \\textsc{OPTI} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} & {\\tt VAR} & least loaded at (\$\\%\$) & call & repair\\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")

    # ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanlb  ::Float64 = 0.

        # EMPTY-MOVE
        EM_opti ::Int64   = 0
        EM_cpu  ::Float64 = 0.
        EM_gap  ::Float64 = 0.
        EM_obj  ::Float64 = 0.
        EM_all_obj::Vector{Float64} = zeros(Float64, length(group_path))
        EM_least_loaded::Float64 = 0.

        EM_call ::Int64 = 0
        EM_success ::Int64 = 0

        println("\nGroup $k: $name")
        for (i, path) in enumerate(group_path)
            print("($i/$(length(group_path))) > ")
            instance::Instance, _ = parseAnyInstance(path)
            
            Lmax        = instance.Lmax
            O           = instance.nbOut
            R           = instance.nbRoute
            lb::Int64   = minSession(instance)

            meanlb += lb  
    # ==========< solve EMPTY-MOVE >==========

            sol_FD_EmptyMove = nothing
            call = 0
            success = 0
            start = time()

            if fit_algo == 1
                sol_FD_EmptyMove, (call, success) = BFD_EM(instance.route, Lmax, O, R)
            elseif fit_algo == 2
                sol_FD_EmptyMove, (call, success) = FFD_EM(instance.route, Lmax, O, R)
            elseif fit_algo == 3
                sol_FD_EmptyMove, (call, success) = NFD_EM(instance.route, Lmax, O, R)
            elseif fit_algo == 4
                sol_FD_EmptyMove, (call, success) = WFD_EM(instance.route, Lmax, O, R)
            end

            EM_opti += (length(sol_FD_EmptyMove) == lb) ? (1) : (0)
            EM_cpu += time() - start
            EM_gap += abs(lb - length(sol_FD_EmptyMove)) / abs(lb)
            EM_obj += length(sol_FD_EmptyMove)
            EM_all_obj[i] = length(sol_FD_EmptyMove)
            EM_call += call
            EM_success += success

            EM_least_loaded += (minimum([sum(s.load) for s in sol_FD_EmptyMove]) * 100) / (Lmax*O)

            println()
        end

    # ==========< Compute res >==========

        meanlb /= length(group_path)

        EM_cpu  = round(EM_cpu / length(group_path), digits=3)
        EM_gap  = round(EM_gap / length(group_path), digits=3)
        EM_obj  = round(EM_obj / length(group_path), digits=3)

        EM_least_loaded = round(EM_least_loaded / length(group_path), digits=3)

        EM_var = round(sum([(EM_all_obj[i] - EM_obj)^2 for i=length(EM_all_obj)]) / (length(EM_all_obj) - 1), digits=3)

    # ==========< Write line >==========

    write(fd, "$(name) & $(O) & $(R) & $(meanlb)") # Head
    write(fd, " & $(EM_opti) & $(EM_cpu) & $(EM_gap) & $(EM_obj) & $(EM_var) & $(EM_least_loaded) & $(EM_call) & $(EM_success)") # EMPTY-MOVE
    write(fd, "\\\\\n")
    flush(fd)

    # ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        \\caption{TODO}\n")
    write(fd, "        \\label{table:res}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterie_SA(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
            (("Diverse", 30, 200), ["instanceIndus_$(i)_O20_R20_C65_opt_1.txt" for i=1:10]),
            (("Diverse", 30, 200), ["instanceIndus_$(i)_O20_R40_C65_opt_2.txt" for i=1:10]),
            (("Diverse", 30, 200), ["instanceIndus_$(i)_O20_R100_C65_opt_5.txt" for i=1:10]),
            (("Indus", 20, 200), ["myTrafic_$(i)_O200_R200.txt" for i=1:10]),

            # (("Indus", 20, 200), ["myTrafic_$(i)_O200_R20.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O200_R40.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O200_R100.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O200_R200.txt" for i=1:10]),
            # (("Indus", 20, 200), ["myTrafic_$(i)_O200_R400.txt" for i=1:10]),

            # (("Diverse", 30, 200), ["instanceIndus_$(i)_O200_R20_C65_opt_1.txt" for i=1:10]),
            # (("Diverse", 30, 200), ["instanceIndus_$(i)_O200_R40_C65_opt_2.txt" for i=1:10]),
            # (("Diverse", 30, 200), ["instanceIndus_$(i)_O200_R100_C65_opt_5.txt" for i=1:10]),
            # (("Diverse", 30, 200), ["instanceIndus_$(i)_O200_R200_C65_opt_10.txt" for i=1:10]),
            # (("Diverse", 30, 200), ["instanceIndus_$(i)_O200_R400_C65_opt_20.txt" for i=1:10]),

            # (("Contained", 20, 200), ["instanceContained_$(i)_O200_R20_C75_opt_1.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O200_R40_C75_opt_2.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O200_R100_C75_opt_5.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O200_R200_C75_opt_10.txt" for i=1:10]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O200_R400_C75_opt_20.txt" for i=1:10]),

            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R20_C60_opt_1.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R40_C60_opt_2.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R100_C60_opt_5.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R200_C60_opt_10.txt" for i=1:10]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R400_C60_opt_20.txt" for i=1:10]),

            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R20_C95_opt_1.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R40_C95_opt_2.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R100_C95_opt_5.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R200_C95_opt_10.txt" for i=1:10]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R400_C95_opt_20.txt" for i=1:10]),
        ],
        result_path         ::String            = "../data/results2.txt",
        fit_algo :: Int64 = 4,
    )
    # ==========< Head >==========

    fd = open(result_path, "a")

    algo_name_tmp = (fit_algo == 1) ? ("BFD") : ((fit_algo == 2) ? "FFD" : ((fit_algo == 3) ? "NFD" : "WFD"))

    write(fd, "\n\n    \\begin{table}[h]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\begin{tabular}{c c c c|c c c c c c c c}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{4}{l|}{Instance} & \\multicolumn{6}{l}{$(algo_name_tmp) (\\textsc{SA})}\\\\\n")
    write(fd, "            type & \$|O|\$ & \$|R|\$ & \\textsc{OPTI} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} & {\\tt VAR} & least loaded at (\$\\%\$) & call & repair\\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")

    # ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanlb  ::Float64 = 0.

        # EMPTY-MOVE
        EM_opti ::Int64   = 0
        EM_cpu  ::Float64 = 0.
        EM_gap  ::Float64 = 0.
        EM_obj  ::Float64 = 0.
        EM_all_obj::Vector{Float64} = zeros(Float64, length(group_path))
        EM_least_loaded::Float64 = 0.

        EM_call ::Int64 = 0
        EM_success ::Int64 = 0

        println("\nGroup $k: $name")
        for (i, path) in enumerate(group_path)
            print("($i/$(length(group_path))) > ")
            instance::Instance, _ = parseAnyInstance(path)
            
            Lmax        = instance.Lmax
            O           = instance.nbOut
            R           = instance.nbRoute
            lb::Int64   = minSession(instance)

            meanlb += lb  
    # ==========< solve EMPTY-MOVE >==========

            sol_FD_EmptyMove = nothing
            call = 0
            success = 0
            start = time()

            if fit_algo == 1
                sol_FD_EmptyMove, (call, success) = BFD_SA(instance.route, Lmax, O, R)
            elseif fit_algo == 2
                sol_FD_EmptyMove, (call, success) = FFD_SA(instance.route, Lmax, O, R)
            elseif fit_algo == 3
                sol_FD_EmptyMove, (call, success) = NFD_SA(instance.route, Lmax, O, R)
            elseif fit_algo == 4
                sol_FD_EmptyMove, (call, success) = WFD_SA(instance.route, Lmax, O, R)
            end

            EM_opti += (length(sol_FD_EmptyMove) == lb) ? (1) : (0)
            EM_cpu += time() - start
            EM_gap += abs(lb - length(sol_FD_EmptyMove)) / abs(lb)
            EM_obj += length(sol_FD_EmptyMove)
            EM_all_obj[i] = length(sol_FD_EmptyMove)
            EM_call += call
            EM_success += success

            EM_least_loaded += (minimum([sum(s.load) for s in sol_FD_EmptyMove]) * 100) / (Lmax*O)

            println()
        end

    # ==========< Compute res >==========

        meanlb /= length(group_path)

        EM_cpu  = round(EM_cpu / length(group_path), digits=3)
        EM_gap  = round(EM_gap / length(group_path), digits=3)
        EM_obj  = round(EM_obj / length(group_path), digits=3)

        EM_least_loaded = round(EM_least_loaded / length(group_path), digits=3)

        EM_var = round(sum([(EM_all_obj[i] - EM_obj)^2 for i=length(EM_all_obj)]) / (length(EM_all_obj) - 1), digits=3)

    # ==========< Write line >==========

    write(fd, "$(name) & $(O) & $(R) & $(meanlb)") # Head
    write(fd, " & $(EM_opti) & $(EM_cpu) & $(EM_gap) & $(EM_obj) & $(EM_var) & $(EM_least_loaded) & $(EM_call) & $(EM_success)") # EMPTY-MOVE
    write(fd, "\\\\\n")
    flush(fd)

    # ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        \\caption{TODO}\n")
    write(fd, "        \\label{table:res}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterie_GR(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
            (("Indus. like", 30, 200), ["instanceIndus_$(i)_O200_R30_C40_opt_1.txt" for i=1:10]),
            (("Indus. like", 60, 200), ["instanceIndus_$(i)_O200_R60_C40_opt_2.txt" for i=1:5]),
            (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R150_C40_opt_5.txt" for i=1:10]),
            (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R300_C40_opt_10.txt" for i=1:10]),
            (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R600_C40_opt_20.txt" for i=1:10]),

            (("Contained", 20, 200), ["instanceContained_$(i)_O200_R20_C150_opt_1.txt" for i=1:10]),
            (("Contained", 40, 200), ["instanceContained_$(i)_O200_R40_C150_opt_2.txt" for i=1:5]),
            (("Contained", 100, 200), ["instanceContained_$(i)_O200_R100_C150_opt_5.txt" for i=1:10]),
            (("Contained", 100, 200), ["instanceContained_$(i)_O200_R200_C150_opt_10.txt" for i=1:10]),
            (("Contained", 100, 200), ["instanceContained_$(i)_O200_R400_C150_opt_20.txt" for i=1:10]),

            (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R20_C60_opt_1.txt" for i=1:10]),
            (("Skewed", 40, 200), ["instanceSkewed_$(i)_O200_R40_C80_opt_2.txt" for i=1:5]),
            (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R100_C60_opt_5.txt" for i=1:10]),
            (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R200_C60_opt_10.txt" for i=1:10]),
            (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R400_C60_opt_20.txt" for i=1:10]),

            (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R20_C100_opt_1.txt" for i=1:10]),
            (("Chunk", 40, 200), ["instanceChunk_$(i)_O200_R40_C100_opt_2.txt" for i=1:5]),
            (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R100_C100_opt_5.txt" for i=1:10]),
            (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R200_C100_opt_10.txt" for i=1:10]),
            (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R400_C100_opt_20.txt" for i=1:10]),

            (("Indus", 200, 200), ["myTrafic_$(i)_O40_R200.txt" for i=1:10]),
            (("Indus", 200, 200), ["myTrafic_$(i)_O100_R200.txt" for i=1:10]),
            (("Indus", 200, 200), ["myTrafic_$(i)_O200_R200.txt" for i=1:10]),
            # (("Indus", 800, 200), ["myTrafic_$(i)_O200_R800.txt" for i=1:10]),
            ]   ,
        result_path         ::String            = "../data/results.txt",
        fit_algo :: Int64 = 1,
    )
    # ==========< Head >==========

    fd = open(result_path, "a")

    algo_name_tmp = (fit_algo == 1) ? ("BFD") : ((fit_algo == 2) ? "FFD" : ((fit_algo == 3) ? "NFD" : "WFD"))

    write(fd, "\n\n    \\begin{table}[h]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\begin{tabular}{c c c c|c c c c c c c c}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{4}{l|}{Instance} & \\multicolumn{6}{l}{$(algo_name_tmp) (\\textsc{Gready Rebuild})}\\\\\n")
    write(fd, "            type & \$|O|\$ & \$|R|\$ & \\textsc{OPTI} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} & {\\tt VAR} & least loaded at (\$\\%\$) & call & repair\\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")

    tl::Int64 = 10
    env::Gurobi.Env = Gurobi.Env()

    # ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanlb  ::Float64 = 0.

        # EMPTY-MOVE
        EM_opti ::Int64   = 0
        EM_cpu  ::Float64 = 0.
        EM_gap  ::Float64 = 0.
        EM_obj  ::Float64 = 0.
        EM_all_obj::Vector{Float64} = zeros(Float64, length(group_path))
        EM_least_loaded::Float64 = 0.

        EM_call ::Int64 = 0
        EM_success ::Int64 = 0

        println("\nGroup $k: $name")
        for (i, path) in enumerate(group_path)
            print("($i/$(length(group_path))) > ")
            instance::Instance, _ = parseAnyInstance(path)
            
            Lmax        = instance.Lmax
            O           = instance.nbOut
            R           = instance.nbRoute
            lb::Int64   = minSession(instance)

            meanlb += lb  
    # ==========< solve EMPTY-MOVE >==========

            sol_FD_EmptyMove = nothing
            call = 0
            success = 0
            start = time()

            if fit_algo == 1
                sol_FD_EmptyMove, (call, success) = BFD_GR(instance.route, Lmax, O, R, tl=tl, env=env)
            elseif fit_algo == 2
                sol_FD_EmptyMove, (call, success) = FFD_GR(instance.route, Lmax, O, R, tl=tl, env=env)
            elseif fit_algo == 3
                sol_FD_EmptyMove, (call, success) = NFD_GR(instance.route, Lmax, O, R, tl=tl, env=env)
            elseif fit_algo == 4
                sol_FD_EmptyMove, (call, success) = WFD_GR(instance.route, Lmax, O, R, tl=tl, env=env)
            end

            EM_opti += (length(sol_FD_EmptyMove) == lb) ? (1) : (0)
            EM_cpu += time() - start
            EM_gap += abs(lb - length(sol_FD_EmptyMove)) / abs(lb)
            EM_obj += length(sol_FD_EmptyMove)
            EM_all_obj[i] = length(sol_FD_EmptyMove)
            EM_call += call
            EM_success += success

            EM_least_loaded += (minimum([sum(s.load) for s in sol_FD_EmptyMove]) * 100) / (Lmax*O)

            println()
        end

    # ==========< Compute res >==========

        meanlb /= length(group_path)

        EM_cpu  = round(EM_cpu / length(group_path), digits=3)
        EM_gap  = round(EM_gap / length(group_path), digits=3)
        EM_obj  = round(EM_obj / length(group_path), digits=3)

        EM_least_loaded = round(EM_least_loaded / length(group_path), digits=3)

        EM_var = round(sum([(EM_all_obj[i] - EM_obj)^2 for i=length(EM_all_obj)]) / (length(EM_all_obj) - 1), digits=3)

    # ==========< Write line >==========

    write(fd, "$(name) & $(O) & $(R) & $(meanlb)") # Head
    write(fd, " & $(EM_opti) & $(EM_cpu) & $(EM_gap) & $(EM_obj) & $(EM_var) & $(EM_least_loaded) & $(EM_call) & $(EM_success)") # EMPTY-MOVE
    write(fd, "\\\\\n")
    flush(fd)

    # ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        \\caption{TODO}\n")
    write(fd, "        \\label{table:res}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterie_Nothing(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
            (("Indus. like", 30, 200), ["instanceIndus_$(i)_O200_R30_C40_opt_1.txt" for i=1:10]),
            (("Indus. like", 60, 200), ["instanceIndus_$(i)_O200_R60_C40_opt_2.txt" for i=1:10]),
            (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R150_C40_opt_5.txt" for i=1:10]),
            (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R300_C40_opt_10.txt" for i=1:10]),
            (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R600_C40_opt_20.txt" for i=1:10]),

            (("Contained", 20, 200), ["instanceContained_$(i)_O200_R20_C150_opt_1.txt" for i=1:10]),
            (("Contained", 40, 200), ["instanceContained_$(i)_O200_R40_C150_opt_2.txt" for i=1:10]),
            (("Contained", 100, 200), ["instanceContained_$(i)_O200_R100_C150_opt_5.txt" for i=1:10]),
            (("Contained", 100, 200), ["instanceContained_$(i)_O200_R200_C150_opt_10.txt" for i=1:10]),
            (("Contained", 100, 200), ["instanceContained_$(i)_O200_R400_C150_opt_20.txt" for i=1:10]),

            (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R20_C60_opt_1.txt" for i=1:10]),
            (("Skewed", 40, 200), ["instanceSkewed_$(i)_O200_R40_C80_opt_2.txt" for i=1:10]),
            (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R100_C60_opt_5.txt" for i=1:10]),
            (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R200_C60_opt_10.txt" for i=1:10]),
            (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R400_C60_opt_20.txt" for i=1:10]),

            (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R20_C100_opt_1.txt" for i=1:10]),
            (("Chunk", 40, 200), ["instanceChunk_$(i)_O200_R40_C100_opt_2.txt" for i=1:10]),
            (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R100_C100_opt_5.txt" for i=1:10]),
            (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R200_C100_opt_10.txt" for i=1:10]),
            (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R400_C100_opt_20.txt" for i=1:10]),

            (("Indus", 200, 200), ["myTrafic_$(i)_O40_R200.txt" for i=1:10]),
            (("Indus", 200, 200), ["myTrafic_$(i)_O100_R200.txt" for i=1:10]),
            (("Indus", 200, 200), ["myTrafic_$(i)_O200_R200.txt" for i=1:10]),
            # (("Indus", 800, 200), ["myTrafic_$(i)_O200_R800.txt" for i=1:10]),
            ]   ,
        result_path         ::String            = "../data/results.txt",
        fit_algo::Int64=1,
    )
    # ==========< Head >==========

    algo_name_tmp = (fit_algo == 1) ? ("BFD") : ((fit_algo == 2) ? "FFD" : ((fit_algo == 3) ? "NFD" : "WFD"))

    fd = open(result_path, "a")

    write(fd, "\n\n    \\begin{table}[h]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\begin{tabular}{c c c c|c c c c c c}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{4}{l|}{Instance} & \\multicolumn{6}{l}{$(algo_name_tmp) (only left allign and smooth assign)}\\\\\n")
    write(fd, "            type & \$|O|\$ & \$|R|\$ & \\textsc{OPTI} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} & {\\tt VAR} & least loaded at (\$\\%\$) \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")

    # ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanlb  ::Float64 = 0.

        # EMPTY-MOVE
        EM_opti ::Int64   = 0
        EM_cpu  ::Float64 = 0.
        EM_gap  ::Float64 = 0.
        EM_obj  ::Float64 = 0.
        EM_all_obj::Vector{Float64} = zeros(Float64, length(group_path))
        EM_least_loaded::Float64 = 0.

        println("\nGroup $k: $name")
        for (i, path) in enumerate(group_path)
            print("($i/$(length(group_path))) > ")
            instance::Instance, _ = parseAnyInstance(path)
            
            Lmax        = instance.Lmax
            O           = instance.nbOut
            R           = instance.nbRoute
            lb::Int64   = minSession(instance)

            meanlb += lb  
    # ==========< solve EMPTY-MOVE >==========
            start = time()

            if fit_algo == 1
                sol_FD_EmptyMove = BFD_Nothing(instance.route, Lmax, O, R)
            elseif fit_algo == 2
                sol_FD_EmptyMove = FFD_Nothing(instance.route, Lmax, O, R)
            elseif fit_algo == 3
                sol_FD_EmptyMove = NFD_Nothing(instance.route, Lmax, O, R)
            elseif fit_algo == 4
                sol_FD_EmptyMove = WFD_Nothing(instance.route, Lmax, O, R)
            end

            EM_opti += (length(sol_FD_EmptyMove) == lb) ? (1) : (0)
            EM_cpu += time() - start
            EM_gap += abs(lb - length(sol_FD_EmptyMove)) / abs(lb)
            EM_obj += length(sol_FD_EmptyMove)
            EM_all_obj[i] = length(sol_FD_EmptyMove)

            EM_least_loaded += (minimum([sum(s.load) for s in sol_FD_EmptyMove]) * 100) / (Lmax*O)

            println()
        end

    # ==========< Compute res >==========

        meanlb /= length(group_path)

        EM_cpu  = round(EM_cpu / length(group_path), digits=3)
        EM_gap  = round(EM_gap / length(group_path), digits=3)
        EM_obj  = round(EM_obj / length(group_path), digits=3)

        EM_least_loaded = round(EM_least_loaded / length(group_path), digits=3)

        EM_var = round(sum([(EM_all_obj[i] - EM_obj)^2 for i=length(EM_all_obj)]) / (length(EM_all_obj) - 1), digits=3)

    # ==========< Write line >==========

    write(fd, "$(name) & $(O) & $(R) & $(meanlb)") # Head
    write(fd, " & $(EM_opti) & $(EM_cpu) & $(EM_gap) & $(EM_obj) & $(EM_var) & $(EM_least_loaded)") # EMPTY-MOVE
    write(fd, "\\\\\n")
    flush(fd)

    # ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        \\caption{TODO}\n")
    write(fd, "        \\label{table:res}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterie_EM_GR(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
            # (("Indus. like", 30, 200), ["instanceIndus_$(i)_O200_R30_C40_opt_1.txt" for i=1:2]),
            # (("Indus. like", 60, 200), ["instanceIndus_$(i)_O200_R60_C40_opt_2.txt" for i=1:5]),
            (("Indus. like", 150, 200), ["instanceIndus_$(i)_O200_R150_C40_opt_5.txt" for i=1:5]),
            # (("Contained", 20, 200), ["instanceContained_$(i)_O200_R20_C150_opt_1.txt" for i=1:30]),
            # (("Contained", 40, 200), ["instanceContained_$(i)_O200_R40_C150_opt_2.txt" for i=1:5]),
            # (("Contained", 100, 200), ["instanceContained_$(i)_O200_R100_C150_opt_5.txt" for i=1:5]),
            # (("Skewed", 20, 200), ["instanceSkewed_$(i)_O200_R20_C80_opt_1.txt" for i=1:30]),
            # (("Skewed", 40, 200), ["instanceSkewed_$(i)_O200_R40_C80_opt_2.txt" for i=1:5]),
            (("Skewed", 100, 200), ["instanceSkewed_$(i)_O200_R100_C80_opt_5.txt" for i=1:5]),
            # (("Chunk", 20, 200), ["instanceChunk_$(i)_O200_R20_C100_opt_1.txt" for i=1:30]),
            # (("Chunk", 40, 200), ["instanceChunk_$(i)_O200_R40_C100_opt_2.txt" for i=1:5]),
            # (("Chunk", 100, 200), ["instanceChunk_$(i)_O200_R100_C100_opt_5.txt" for i=1:5]),
            # (("Indus", 200, 200), ["myTrafic_$(i)_O200_R200.txt" for i=1:5]),
            # (("Indus", 800, 200), ["myTrafic_$(i)_O200_R800.txt" for i=1:5]),
            ]   ,
        result_path         ::String            = "../data/results.txt",
        tl::Int64 = 10,
        env::Gurobi.Env = Gurobi.Env(),
    )
    # ==========< Head >==========
    fd = open(result_path, "a")

    write(fd, "\n\n    \\begin{table}[h]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\begin{tabular}{c c c c|c c c c|c c c c}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{4}{l|}{Instance} & \\multicolumn{4}{l|}{BFD (\\textsc{EMPTY-MOVE})} & \\multicolumn{4}{l}{BFD (\\textsc{GREADY-REBUILD})} \\\\\n")
    write(fd, "            type & \$|O|\$ & \$|R|\$ & \\textsc{OPTI} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")
    # ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanlb  ::Float64 = 0.

        # EMPTY-MOVE
        EM_opti ::Int64   = 0
        EM_cpu  ::Float64 = 0.
        EM_gap  ::Float64 = 0.
        EM_obj  ::Float64 = 0.

        # GREADY REBUILD
        GR_opti ::Int64   = 0
        GR_cpu  ::Float64 = 0.
        GR_gap  ::Float64 = 0.
        GR_obj  ::Float64 = 0.

        for path in group_path
            instance::Instance, _ = parseAnyInstance(path)
            
            Lmax        = instance.Lmax
            O           = instance.nbOut
            R           = instance.nbRoute
            lb::Int64   = minSession(instance)

            meanlb += lb  
    # ==========< solve EMPTY-MOVE >==========

            start = time()
            sol_BFD_EmptyMove = BFD_EmptyMove(instance.route, Lmax, O, R)

            EM_opti += (length(sol_BFD_EmptyMove) == lb) ? (1) : (0)
            EM_cpu += time() - start
            EM_gap += abs(lb - length(sol_BFD_EmptyMove)) / abs(lb)
            EM_obj += length(sol_BFD_EmptyMove)

    # ==========< GREADY REBUILD >==========
                
            start = time()
            sol_BFD_GreadyRebuild = BFD_GreadyRebuild(instance.route, Lmax, O, R, tl, env)

            GR_opti += (length(sol_BFD_GreadyRebuild) == lb) ? (1) : (0)
            GR_cpu += time() - start
            GR_gap += abs(lb - length(sol_BFD_GreadyRebuild)) / abs(lb)
            GR_obj += length(sol_BFD_GreadyRebuild)

        end

    # ==========< Compute res >==========

        meanlb /= length(group_path)

        EM_cpu  = round(EM_cpu / length(group_path), digits=3)
        EM_gap  = round(EM_gap / length(group_path), digits=3)
        EM_obj  = round(EM_obj / length(group_path), digits=3)

        # GREADY REBUILD
        GR_cpu  = round(GR_cpu / length(group_path), digits=3)
        GR_gap  = round(GR_gap / length(group_path), digits=3)
        GR_obj  = round(GR_obj / length(group_path), digits=3)

    # ==========< Write line >==========

    write(fd, "$(name) & $(O) & $(R) & $(meanlb)") # Head
    write(fd, " & $(EM_opti) & $(EM_cpu) & $(EM_gap) & $(EM_obj)") # EMPTY-MOVE
    write(fd, " & $(GR_opti) & $(GR_cpu) & $(GR_gap) & $(GR_obj)") # GREADY REBUILD
    write(fd, "\\\\\n")
    flush(fd)

    # ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        \\caption{TODO}\n")
    write(fd, "        \\label{table:res}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterieHard_01LP_GA_Build(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
        # (("Indus. like", 30, 200), ["instanceDistrib_$(i)_O200_R30_C40_opt_1.txt" for i=1:30]),
        # (("Indus. like", 60, 200), ["instanceDistrib_$(i)_O200_R60_C40_opt_2.txt" for i=1:30]),
        # (("Indus. like", 150, 200), ["instanceDistrib_$(i)_O200_R150_C40_opt_5.txt" for i=1:30]),
        # (("Contained", 20, 200), ["instanceGaussian_$(i)_O200_R20_C150_opt_1.txt" for i=1:30]),
        # (("Contained", 40, 200), ["instanceGaussian_$(i)_O200_R40_C150_opt_2.txt" for i=1:30]),
        # (("Contained", 100, 200), ["instanceGaussian_$(i)_O200_R100_C150_opt_5.txt" for i=1:30]),
        # (("Skewed", 20, 200), ["instanceBigBatche_$(i)_O200_R20_C80_opt_1.txt" for i=1:30]),
        # (("Skewed", 40, 200), ["instanceBigBatche_$(i)_O200_R40_C80_opt_2.txt" for i=1:30]),
        # (("Skewed", 100, 200), ["instanceBigBatche_$(i)_O200_R100_C80_opt_5.txt" for i=1:30]),
        (("Chunk", 20, 200), ["instanceDistribHard_$(i)_O200_R20_C100_opt_1.txt" for i=1:1]),],
        # (("Chunk", 40, 200), ["instanceDistribHard_$(i)_O200_R40_C100_opt_2.txt" for i=1:30]),
        # (("Chunk", 100, 200), ["instanceDistribHard_$(i)_O200_R100_C100_opt_5.txt" for i=1:30])],
        MILP_time_limit     ::Int64             = 180,
        GA_time_limit       ::Int64             = 180,
        result_path         ::String            = "../data/results.txt"
    )
# ==========< Init >==========
    nbLine::Int64 = 0

    fd = open(result_path, "a")
# ==========< Head >==========
    write(fd, "\n\n    \\begin{table}[H]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\resizebox{\\textwidth}{!}{\n")
    write(fd, "        \\begin{tabular}{|| c|c|c|c || c|c|c|c || c|c|c|c || c|c|c|c ||}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{4}{||l||}{Instance} & \\multicolumn{4}{l||}{0-1LP} & \\multicolumn{4}{l||}{GA} & \\multicolumn{4}{l||}{Build solution} \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            type & \$|R|\$ & \$|O|\$ & \\textsc{OPTI} & \\#{\\tt OPTI} &  {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")
# ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanbound       ::Float64           = 0.

        # MILP
        MILP_opti       ::Int64             = 0 
        MILP_cput       ::Vector{Float64}   = [] 
        MILP_gap        ::Vector{Float64}   = []
        MILP_obj        ::Vector{Float64}   = []

        # GA
        GA_opti         ::Int64             = 0
        GA_cput         ::Vector{Float64}   = [] 
        GA_gap          ::Vector{Float64}   = []
        GA_obj          ::Vector{Float64}   = []

        # GA
        Build_opti      ::Int64             = 0
        Build_cput      ::Vector{Float64}   = [] 
        Build_gap       ::Vector{Float64}   = []
        Build_obj       ::Vector{Float64}   = []

        for path in group_path
            instance::Instance, _ = parseAnyInstance(path)

            bound::Int64 = minSession(instance)
            meanbound += bound
# ==========< Build >==========

            ObjVal, flag_opti, RunTime, Gap, warmup_sol = run_Build(instance, bound)

            if flag_opti
                Build_opti += 1
            end
            push!(Build_gap, Gap)
            push!(Build_cput, RunTime)
            push!(Build_obj, ObjVal)

# ==========< 01LP >==========

            ObjVal, flag_opti, RunTime, Gap, _ = run_01LP(instance, warmup_sol, MILP_time_limit)
            
            if ObjVal != -1 # at leat one valid solution
                if flag_opti 
                    MILP_opti += 1
                end
            end
            push!(MILP_cput, RunTime)
            push!(MILP_gap, Gap)
            push!(MILP_obj, ObjVal)

# ==========< GA >==========
                
            ObjVal, flag_opti, RunTime, Gap = run_GA(instance, bound, GA_time_limit)
            if flag_opti
                GA_opti += 1
            end
            push!(GA_cput, RunTime)
            push!(GA_gap, Gap)
            push!(GA_obj, ObjVal)
        end

# ==========< Compute res >==========

        meanbound /= length(group_path)

        # CPUT
        MILP_mean_cput  = isempty(MILP_cput)    ? "-" : round(mean(MILP_cput)   , digits=3)
        GA_mean_cput    = isempty(GA_cput)      ? "-" : round(mean(GA_cput)     , digits=3)
        Build_mean_cput = isempty(Build_cput)   ? "-" : round(mean(Build_cput)  , digits=3)

        # GAP
        MILP_mean_gap   = isempty(MILP_gap)     ? "-" : round(mean(MILP_gap)    , digits=3)
        GA_mean_gap     = isempty(GA_gap)       ? "-" : round(mean(GA_gap)      , digits=3)
        Build_mean_gap  = isempty(Build_gap)    ? "-" : round(mean(Build_gap)   , digits=3)

        # GAP
        MILP_mean_obj   = isempty(MILP_obj)     ? "-" : round(mean(MILP_obj)    , digits=3)
        GA_mean_obj     = isempty(GA_obj)       ? "-" : round(mean(GA_obj)      , digits=3)
        Build_mean_obj  = isempty(Build_obj)    ? "-" : round(mean(Build_obj)   , digits=3)

# ==========< Write line >==========

    write(fd, "            $(name) & $(R) & $(O) & $(meanbound)") # Head
    write(fd, " & $(MILP_opti) & $(MILP_mean_cput) & $(MILP_mean_gap) & $(MILP_mean_obj)") # 01LP
    write(fd, " & $(GA_opti) & $(GA_mean_cput) & $(GA_mean_gap) & $(GA_mean_obj)") # GA
    write(fd, " & $(Build_opti) & $(Build_mean_cput) & $(Build_mean_gap) & $(Build_mean_obj)") # Build
    write(fd, "\\\\\\hline\n") # Tail

    flush(fd)

# ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        }\n")
    write(fd, "        \\caption{Results Hard instances}\n")
    write(fd, "        \\label{table:result-table-hard-instances}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterieHard_GA_Build(;
        all_group            ::Vector{Tuple{Tuple{String, Int64, Int64}, Vector{String}}}     = [ # |O|
                                                                            (("Indus. like", 30, 200), ["instanceDistrib_$(i)_O200_R30_C40_opt_1.txt" for i=1:30]),
                                                                            (("Indus. like", 60, 200), ["instanceDistrib_$(i)_O200_R60_C40_opt_2.txt" for i=1:30]),
                                                                            (("Indus. like", 150, 200), ["instanceDistrib_$(i)_O200_R150_C40_opt_5.txt" for i=1:30]),
                                                                            (("Contained", 20, 200), ["instanceGaussian_$(i)_O200_R20_C150_opt_1.txt" for i=1:30]),
                                                                            (("Contained", 40, 200), ["instanceGaussian_$(i)_O200_R40_C150_opt_2.txt" for i=1:30]),
                                                                            (("Contained", 100, 200), ["instanceGaussian_$(i)_O200_R100_C150_opt_5.txt" for i=1:30]),
                                                                            (("Skewed", 20, 200), ["instanceBigBatche_$(i)_O200_R20_C80_opt_1.txt" for i=1:30]),
                                                                            (("Skewed", 40, 200), ["instanceBigBatche_$(i)_O200_R40_C80_opt_2.txt" for i=1:30]),
                                                                            (("Skewed", 100, 200), ["instanceBigBatche_$(i)_O200_R100_C80_opt_5.txt" for i=1:30]),
                                                                            (("Chunk", 20, 200), ["instanceDistribHard_$(i)_O200_R20_C100_opt_1.txt" for i=1:30]),
                                                                            (("Chunk", 40, 200), ["instanceDistribHard_$(i)_O200_R40_C100_opt_2.txt" for i=1:30]),
                                                                            (("Chunk", 100, 200), ["instanceDistribHard_$(i)_O200_R100_C100_opt_5.txt" for i=1:30]),
                                                                            ]   ,
        GA_time_limit       ::Int64             =20    ,
        result_path         ::String            = "TER/results/results.txt"
    )
# ==========< Init >==========
    nbLine::Int64 = 0

    fd = open(result_path, "a")
# ==========< Head >==========
    write(fd, "\n\n    \\begin{table}[H]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\resizebox{\\textwidth}{!}{\n")
    write(fd, "        \\begin{tabular}{|| c|c|c|c || c|c|c|c || c|c|c|c ||}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{4}{||l||}{Instance} & \\multicolumn{4}{l||}{GA} & \\multicolumn{4}{l||}{BFD} \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \$|O|\$ & type & \$|R|\$ & \\textsc{OPTI} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & {\\tt OBJ} \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")
# ==========< Init group >==========
    for (k::Int64, ((name::String, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanbound       ::Float64           = 0.

        # GA
        GA_opti         ::Int64             = 0
        GA_cput         ::Vector{Float64}   = [] 
        GA_gap          ::Vector{Float64}   = []
        GA_obj          ::Vector{Int64}     = []

        # GA
        Build_opti      ::Int64             = 0
        Build_cput      ::Vector{Float64}   = [] 
        Build_gap       ::Vector{Float64}   = []
        Build_obj       ::Vector{Int64}     = []

        for path in group_path
            instance::Instance, _ = parseAnyInstance(path)

            bound::Int64 = minSession(instance)
            meanbound += bound  
# ==========< Build >==========

            ObjVal, flag_opti, RunTime, Gap, sol = run_Build(instance, bound)

            flag_opti && (Build_opti += 1)
            push!(Build_gap, Gap)
            push!(Build_cput, RunTime)
            push!(Build_obj, ObjVal)

# ==========< GA >==========
                
            ObjVal, flag_opti, RunTime, Gap, sol = run_GA(instance, bound, GA_time_limit)
            
            flag_opti && (GA_opti += 1)
            push!(GA_cput, RunTime)
            push!(GA_gap, Gap)
            push!(GA_obj, ObjVal)

        end

# ==========< Compute res >==========

        meanbound /= length(group_path)

        # CPUT
        GA_mean_cput        ::String    = isempty(GA_cput)      ? "-" : string(round(mean(GA_cput)     , digits=3))
        Build_mean_cput     ::String    = isempty(Build_cput)   ? "-" : string(round(mean(Build_cput)  , digits=3))

        # GAP
        GA_mean_gap         ::String    = isempty(GA_gap)       ? "-" : string(round(mean(GA_gap)      , digits=3))
        Build_mean_gap      ::String    = isempty(Build_gap)    ? "-" : string(round(mean(Build_gap)   , digits=3))

        # OBJ
        GA_mean_obj         ::String    = isempty(GA_obj)       ? "-" : string(round(mean(GA_obj)      , digits=3))
        Build_mean_obj      ::String    = isempty(Build_obj)    ? "-" : string(round(mean(Build_obj)   , digits=3))

# ==========< Write line >==========

    write(fd, "            $(O) & $(name) & $(R) & $(meanbound)") # Head
    write(fd, " & $(GA_opti) & $(GA_mean_cput) & $(GA_mean_gap) & $(GA_mean_obj)") # GA
    write(fd, " & $(Build_opti) & $(Build_mean_cput) & $(Build_mean_gap) & $(Build_mean_obj)") # Build
    write(fd, "\\\\\\hline\n") # Tail
    flush(fd)

# ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        }\n")
    write(fd, "        \\caption{Results Hard instances}\n")
    write(fd, "        \\label{table:result-table-hard-instances}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterieHard_Build_SortingCrit(;
        all_group            ::Vector{Tuple{Tuple{Float64, Int64, Int64}, Vector{String}}}     = [ # group -> ((type, |R|, |O|), [all instances of the group])
                                                                            # ((0.25, 200, 20), ["instanceGaussian_$(i)_O20_R200_C300_opt_10.txt" for i=1:10]),
                                                                            # ((0.25, 200, 200), ["instanceGaussian_$(i)_O200_R200_C300_opt_10.txt" for i=1:10]),
                                                                            # ((0.25, 600, 20), ["instanceDistrib_$(i)_O20_R600_C300_opt_10.txt" for i=1:10]),
                                                                            # ((0.25, 600, 200), ["instanceDistrib_$(i)_O200_R600_C300_opt_10.txt" for i=1:10]),
                                                                            # ((0.25, 600, 20), ["instanceBigBatche_$(i)_O20_R600_C300_opt_10.txt" for i=1:10]),
                                                                            # ((0.25, 600, 200), ["instanceBigBatche_$(i)_O200_R600_C300_opt_10.txt" for i=1:10]),
                                                                            ((0.25, 200, 20), ["myTrafic_$(i)_O20_R200.txt" for i=1:100]),
                                                                            ((0.5, 200, 20), ["myTrafic_$(i)_O20_R200.txt" for i=1:100]),
                                                                            ((1.0, 200, 20), ["myTrafic_$(i)_O20_R200.txt" for i=1:100]),
                                                                            ((0.25, 200, 40), ["myTrafic_$(i)_O40_R200.txt" for i=1:100]),
                                                                            ((0.5, 200, 40), ["myTrafic_$(i)_O40_R200.txt" for i=1:100]),
                                                                            ((1.0, 200, 40), ["myTrafic_$(i)_O40_R200.txt" for i=1:100]),
                                                                            ((0.25, 200, 80), ["myTrafic_$(i)_O80_R200.txt" for i=1:100]),
                                                                            ((0.5, 200, 80), ["myTrafic_$(i)_O80_R200.txt" for i=1:100]),
                                                                            ((1.0, 200, 80), ["myTrafic_$(i)_O80_R200.txt" for i=1:100]),
                                                                            # ((0.25, 200, 100), ["myTrafic_$(i)_O100_R200.txt" for i=1:100]),
                                                                            # ((0.5, 200, 100), ["myTrafic_$(i)_O100_R200.txt" for i=1:100]),
                                                                            # ((1.0, 200, 100), ["myTrafic_$(i)_O100_R200.txt" for i=1:100]),
                                                                            ],
        allSortingCrit      ::Vector{Tuple{DataType, Vararg{DataType}}} = [ # sorting criteria -> NTuple of Round fitness TAG
                                                                            (NbBatches, ), 
                                                                            (StdBatches, ), 
                                                                            (MaxMin, ),
                                                                            (Max, ),
                                                                            (Max, NbBatches), 
                                                                            (NbBatches, Max, StdBatches),
                                                                            ],
        result_path         ::String            = "TER/results/results.txt"
    )

# ==========< Init >==========

    nbInstances::Int64 = length(all_group)

    fd = open(result_path, "a")

# ==========< Head >==========
    
    write(fd, "\n\n    \\begin{table}[H]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\resizebox{\\textwidth}{!}{\n")
    write(fd, "        \\begin{tabular}{|| m{25mm}|c ||$(repeat(" c |", nbInstances))}\n")
    write(fd, "            \\cline{2-$(2 + nbInstances)}\n")
    write(fd, "            \\multicolumn{1}{l|}{ } & \$|R|\$ & $(join([e[2] for (e, _) in all_group], " & ")) \\\\\\cline{2-$(2 + nbInstances)}\n")
    write(fd, "            \\multicolumn{1}{l|}{ } & \$|O|\$ & $(join([e[3] for (e, _) in all_group], " & ")) \\\\\\cline{2-$(2 + nbInstances)}\n")
    write(fd, "            \\multicolumn{1}{l|}{ } & \$\\lambda\$ & $(join([e[1] for (e, _) in all_group], " & ")) \\\\\\cline{2-$(2 + nbInstances)}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")

# ==========< Init group >==========

    for sortingCrit in allSortingCrit

        Build_opti      ::Vector{Int64}   = zeros(Int64, nbInstances)
        Build_cput      ::Vector{Float64}   = zeros(Int64, nbInstances)
        Build_gap       ::Vector{Float64}   = zeros(Int64, nbInstances)

        for (k::Int64, ((λ::Float64, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)
            for path in group_path

                instance::Instance, _ = parseAnyInstance(path, λ)
                bound::Int64 = minSession(instance)

# ==========< Build >==========

                perm            ::Vector{Int64}                     = sortperm([ntuple(i -> -fitness(r, sortingCrit[i]), length(sortingCrit)) for r in instance.rounds])
                TAG_FitSes      ::Type{<:FitnessSession}            = LoadSTD

                start_time  ::Float64   = time()
                s           ::Solution  = buildSolution_FFD_final(instance, perm, TAG_FitSes)

                ObjVal      ::Int64     = length(s.sessions)
                flag_opti   ::Bool      = ObjVal == bound
                RunTime     ::Float64   = time() - start_time
                Gap         ::Float64   = abs(bound - ObjVal) / abs(ObjVal)

                if flag_opti
                    Build_opti[k] += 1
                else
                    Build_gap[k] += 1
                end
                Build_cput[k] += RunTime

# ==========< Compute res >==========
            end

            Build_gap[k] /= length(group_path)
            Build_gap[k] = round(Build_gap[k], digits=3)
            Build_cput[k] /= length(group_path)
            Build_cput[k] = round(Build_cput[k], digits=3)


# ==========< Write line >==========
        end
    write(fd, "            \\multirow{3}{25mm}{$(join(["{\\tt $(tag)}" for tag in sortingCrit], ", "))} & \\#{\\tt OPTI} & $(join(Build_opti, " & ")) \\\\\\cline{2-$(2 + nbInstances)}\n") # OPTI
    write(fd, "            & {\\tt CPU} & $(join(Build_cput, " & ")) \\\\\\cline{2-$(2 + nbInstances)}\n") # CPU
    write(fd, "            & {\\tt GAP} & $(join(Build_gap, " & ")) \\\\\\hline\\hline\n") # GAP

# ==========< Tail >==========

    end

    write(fd, "        \\end{tabular}\n")
    write(fd, "        }\n")
    write(fd, "        \\caption{Results sorting criteria comparaison}\n")
    write(fd, "        \\label{table:result-sort-crit}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end

function batterieEasy_01LP_GA_Build(;
        all_group           ::Vector{Tuple{Tuple{Float64, Int64, Int64}, Vector{String}}}   = [ # |O|
                                                                                                ((0.25, 20, 200), ["myTrafic_$(i)_O200_R20.txt" for i=1:10]),
                                                                                                ((0.5, 20, 200), ["myTrafic_$(i)_O200_R20.txt" for i=1:10]),
                                                                                                ((1.0, 20, 200), ["myTrafic_$(i)_O200_R20.txt" for i=1:10]),
                                                                                                ((0.25, 80, 200), ["myTrafic_$(i)_O200_R80.txt" for i=1:10]),
                                                                                                ((0.5, 80, 200), ["myTrafic_$(i)_O200_R80.txt" for i=1:10]),
                                                                                                ((1.0, 80, 200), ["myTrafic_$(i)_O200_R80.txt" for i=1:10]),
                                                                                                ((0.25, 150, 200), ["myTrafic_$(i)_O200_R150.txt" for i=1:10]),
                                                                                                ((0.5, 150, 200), ["myTrafic_$(i)_O200_R150.txt" for i=1:10]),
                                                                                                ((1.0, 150, 200), ["myTrafic_$(i)_O200_R150.txt" for i=1:10]),
                                                                                                ],
        MILP_time_limit     ::Int64             = 120    ,
        GA_time_limit       ::Int64             = 120    ,
        result_path         ::String            = "TER/results/results.txt"
    )
# ==========< Init >==========
    nbLine::Int64 = 0

    fd = open(result_path, "a")
# ==========< Head >==========
    write(fd, "\n\n    \\begin{table}[H]\n")
    write(fd, "        \\centering\n")
    write(fd, "        \\resizebox{\\textwidth}{!}{\n")
    write(fd, "        \\begin{tabular}{|| c|c|c || c|c|c || c|c|c || c|c|c ||}\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\multicolumn{3}{||l||}{Instance} & \\multicolumn{3}{l||}{0-1LP} & \\multicolumn{3}{l||}{GA} & \\multicolumn{3}{l||}{Build solution} \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \$|O|\$ & \$|R|\$ & \$\\lambda\$ & \\#{\\tt OPTI} &  {\\tt CPU} & {\\tt GAP} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} & \\#{\\tt OPTI} & {\\tt CPU} & {\\tt GAP} \\\\\n")
    write(fd, "            \\hline\n")
    write(fd, "            \\hline\n")
# ==========< Init group >==========
    for (k::Int64, ((λ::Float64, R::Int64, O::Int64), group_path::Vector{String})) in enumerate(all_group)

        meanbound       ::Float64           = 0.

        # MILP
        MILP_opti       ::Int64             = 0 
        MILP_cput       ::Vector{Float64}   = [] 
        MILP_gap        ::Vector{Float64}   = []

        # GA
        GA_opti         ::Int64             = 0
        GA_cput         ::Vector{Float64}   = [] 
        GA_gap          ::Vector{Float64}   = []

        # GA
        Build_opti      ::Int64             = 0
        Build_cput      ::Vector{Float64}   = [] 
        Build_gap       ::Vector{Float64}   = []

        for path in group_path
            instance::Instance, _ = parseAnyInstance(path, λ)

            bound::Int64 = minSession(instance)
            meanbound += bound
# ==========< Build >==========

            ObjVal, flag_opti, RunTime, Gap, warmup_sol = run_Build(instance, bound)

            if flag_opti
                Build_opti += 1
            else
                push!(Build_gap, Gap)
            end
            push!(Build_cput, RunTime)

# ==========< 01LP >==========

            ObjVal, flag_opti, RunTime, Gap, _ = run_01LP(instance, warmup_sol, MILP_time_limit)
            
            if ObjVal != -1 # at leat one valid solution
                if flag_opti 
                    MILP_opti += 1
                    push!(MILP_cput, RunTime)
                else
                    push!(MILP_gap, Gap)
                end
            end

# ==========< GA >==========
                
            ObjVal, flag_opti, RunTime, Gap = run_GA(instance, bound, GA_time_limit)
            if flag_opti
                GA_opti += 1
                push!(GA_cput, RunTime)
            else
                push!(GA_gap, Gap)
            end

        end

# ==========< Compute res >==========

        meanbound /= length(group_path)

        # CPUT
        MILP_mean_cput  = isempty(MILP_cput)    ? 0.0 : round(mean(MILP_cput)   , digits=3)
        GA_mean_cput    = isempty(GA_cput)      ? 0.0 : round(mean(GA_cput)     , digits=3)
        Build_mean_cput = isempty(Build_cput)   ? 0.0 : round(mean(Build_cput)  , digits=3)

        # GAP
        MILP_mean_gap   = isempty(MILP_gap)     ? 0.0 : round(mean(MILP_gap)    , digits=3)
        GA_mean_gap     = isempty(GA_gap)       ? 0.0 : round(mean(GA_gap)      , digits=3)
        Build_mean_gap  = isempty(Build_gap)    ? 0.0 : round(mean(Build_gap)   , digits=3)

# ==========< Write line >==========

    write(fd, "               $(O) & $(R) & $(λ)") # Head
    write(fd, " & $(MILP_opti) & $(MILP_mean_cput) & $(MILP_mean_gap)") # 01LP
    write(fd, " & $(GA_opti) & $(GA_mean_cput) & $(GA_mean_gap)") # GA
    write(fd, " & $(Build_opti) & $(Build_mean_cput) & $(Build_mean_gap)") # Build
    write(fd, "\\\\\\hline\n") # Tail

# ==========< Tail >==========

    end

    write(fd, "            \\hline\n")
    write(fd, "        \\end{tabular}\n")
    write(fd, "        }\n")
    write(fd, "        \\caption{Results Hard instances}\n")
    write(fd, "        \\label{table:result-table-hard-instances}\n")
    write(fd, "    \\end{table}\n")

    close(fd)
end


# ================================================================================= #
#                           ######  #        #####  #######                         #
#                           #     # #       #     #    #                            #
#                           ######  #       #     #    #                            #
#                           #       #       #     #    #                            #
#                           #       #######  #####     #                            #
# ================================================================================= #

function genVal_to_Plot(instance::Instance, gen_val::Vector{GenVal}, result_path::String = "TER/results/results.txt")
# > INIT
    
    fd = open(result_path, "a")

# > SCALE COMPUTATION
    bound = minSession(instance)

    all_gen     ::Vector{Vector{Float64}}   = [[e.elite; e.ffr; e.cross_over] for e in gen_val]
    bestVal     ::Float64                   = minimum([gen_val[end].elite; gen_val[end].ffr; gen_val[end].cross_over])
    worstVal    ::Float64                   = maximum(maximum([gen_val[g].elite; gen_val[g].ffr; gen_val[g].cross_over]) for g=1:length(gen_val))

    xMin        ::Int64 = 0
    xMax        ::Int64 = length(gen_val)
    nbXLabel    ::Int64 = 10

    yMin        ::Int64 = round(Int64, bestVal)
    yMax        ::Int64 = ceil(Int64, worstVal)
    nbYLabel    ::Int64 = 5

    gen_val = [GenVal(round.(gen.elite, digits=3), round.(gen.ffr, digits=3), round.(gen.cross_over, digits=3)) for gen in gen_val]

# > HEAD
    write(fd, "\n\n\n\n\\begin{figure}[H]\n")
    write(fd, "    \\centering\n")
    write(fd, "    \\begin{tikzpicture}[scale=1, font=\\tiny,rounded corners=2, auto, draw=black, thick]\n")
    write(fd, "        \\begin{axis}[\n")
    write(fd, "            title={Genetic Algorithm convergence},\n")
    write(fd, "            xlabel={Generation},\n")
    write(fd, "            ylabel={Objective value},\n")
    write(fd, "            xmin=$(xMin), xmax=$(xMax),\n")
    write(fd, "            ymin=$(yMin), ymax=$(yMax),\n")
    write(fd, "            xtick={$(join([round(Int64, i) for i=xMin:(xMax - xMin)/nbXLabel:xMax], ","))},\n")
    write(fd, "            ytick={$(join([round(Int64, i) for i=yMin:yMax], ","))},\n")
    write(fd, "            legend pos=north east,\n")
    write(fd, "            xmajorgrids=true,\n")
    write(fd, "            grid style=dashed,\n")
    write(fd, "        ]\n")
    write(fd, "        \n")

# > PLOTS

    write(fd, "        \\addplot[\n")
    write(fd, "        color=green!25,]\n")
    write(fd, "        coordinates {\n")
    write(fd, "        ($(join(["$(xMin + i-1),$(minimum(gen))" for (i, gen) in enumerate(all_gen)], ")(")))\n")
    write(fd, "        };\n")
    write(fd, "        \n")

# > MARKS

    write(fd, "        % Generation 0\n")
    write(fd, "        \\addplot[only marks,mark size=0.5pt,color=blue!50,opacity=0.5]\n")
    write(fd, "        coordinates {($(join(["$(xMin), $(val)" for val in gen_val[1].cross_over], ")(")))};\n")
    write(fd, "        \\addplot[only marks,mark size=0.5pt,color=red!50,opacity=0.5]\n")
    write(fd, "        coordinates {($(join(["$(xMin), $(val)" for val in gen_val[1].ffr], ")(")))};\n")
    write(fd, "        \\addplot[only marks,mark size=0.5pt,color=green!50,opacity=0.5]\n")
    write(fd, "        coordinates {($(join(["$(xMin), $(val)" for val in gen_val[1].elite], ")(")))};\n")
    write(fd, "        \n")

    write(fd, "        \\legend{Best solution, Elite, FF, Cross Over}\n")

    for (i, gen::GenVal) in enumerate(gen_val[2:end])

    write(fd, "        % Generation $(i)\n")
        if (!isempty(gen.cross_over))
    write(fd, "        \\addplot[only marks,mark size=0.5pt,color=blue!50,opacity=0.5]\n")
    write(fd, "        coordinates {($(join(["$(xMin + i),$(val)" for val in gen.cross_over], ")(")))};\n")
        end
        if (!isempty(gen.ffr))
    write(fd, "        \\addplot[only marks,mark size=0.5pt,color=red!50,opacity=0.5]\n")
    write(fd, "        coordinates {($(join(["$(xMin + i),$(val)" for val in gen.ffr], ")(")))};\n")
        end
        if (!isempty(gen.elite))
    write(fd, "        \\addplot[only marks,mark size=0.5pt,color=green!50,opacity=0.5]\n")
    write(fd, "        coordinates {($(join(["$(xMin + i),$(val)" for val in gen.elite], ")(")))};\n")
        end
    write(fd, "        \n")
    
    end

    write(fd, "        \\end{axis}\n")
    write(fd, "    \\end{tikzpicture}\n")
    write(fd, "    \\caption{Genetic Algorithm convergence}\n")
    write(fd, "    \\label{fig:GA-converge}\n")
    write(fd, "\\end{figure}\n")

    close(fd)
end

function genVal_to_Plot3D(instance::Instance, gen_val::Vector{GenVal}, result_path::String = "TER/results/results.txt")
# > INIT
    
    fd = open(result_path, "a")

# > SCALE COMPUTATION
    bound = minSession(instance)

    all_gen     ::Vector{Vector{Float64}}   = [[e.elite; e.ffr; e.cross_over] for e in gen_val]
    bestVal     ::Float64                   = minimum([gen_val[end].elite; gen_val[end].ffr; gen_val[end].cross_over])
    worstVal    ::Float64                   = maximum(maximum([gen_val[g].elite; gen_val[g].ffr; gen_val[g].cross_over]) for g=1:length(gen_val))

    yMin        ::Int64 = 0
    yMax        ::Int64 = length(gen_val)
    nbYLabel    ::Int64 = 10

    zMin        ::Int64 = round(Int64, bestVal)
    zMax        ::Int64 = ceil(Int64, worstVal)
    nbZLabel    ::Int64 = 5

    gen_val = [GenVal(round.(gen.elite, digits=3), round.(gen.ffr, digits=3), round.(gen.cross_over, digits=3)) for gen in gen_val]

# > HEAD
    write(fd, "\n\n\n\n\\begin{figure}[H]\n")
    write(fd, "    \\centering\n")
    write(fd, "    \\begin{tikzpicture}[scale=1, font=\\tiny,rounded corners=2, auto, draw=black, thick]\n")
    write(fd, "        \\begin{axis}[\n")
    write(fd, "            title={Genetic Algorithm convergence},\n")
    write(fd, "            view={60}{45},\n")
    write(fd, "            ylabel={Generation},\n")
    write(fd, "            zlabel={Objective value},\n")
    write(fd, "            xmin=1, xmax=3,\n")
    write(fd, "            ymin=$(yMin), ymax=$(yMax),\n")
    write(fd, "            zmin=$(zMin), zmax=$(zMax),\n")
    write(fd, "            ytick={$(join([round(Int64, i) for i=yMin:(yMax - yMin)/nbYLabel:yMax], ","))},\n")
    write(fd, "            ztick={$(join([round(Int64, i) for i=zMin:zMax], ","))},\n")
    #write(fd, "            legend pos=north east,\n")
    write(fd, "            ymajorgrids=true,\n")
    write(fd, "            grid style=dashed,\n")
    write(fd, "        ]\n")
    write(fd, "        \n")

# > PLOTS

    write(fd, "        \\addplot3[\n")
    write(fd, "        color=green!25,]\n")
    write(fd, "        coordinates {\n")
    write(fd, "        ($(join(["1,$(yMin + i-1),$(minimum(gen))" for (i, gen) in enumerate(all_gen)], ")(")))\n")
    write(fd, "        };\n")
    write(fd, "        \n")

# > MARKS

    write(fd, "        % Generation 0\n")
    write(fd, "        \\addplot3[only marks,mark size=0.5pt,color=blue!50]\n")
    write(fd, "        coordinates {($(join(["$(1),$(yMin), $(val)" for val in gen_val[1].cross_over], ")(")))};\n")
    write(fd, "        \\addplot3[only marks,mark size=0.5pt,color=red!50]\n")
    write(fd, "        coordinates {($(join(["$(2),$(yMin), $(val)" for val in gen_val[1].ffr], ")(")))};\n")
    write(fd, "        \\addplot3[only marks,mark size=0.5pt,color=green!50]\n")
    write(fd, "        coordinates {($(join(["$(3),$(yMin), $(val)" for val in gen_val[1].elite], ")(")))};\n")
    write(fd, "        \n")

    write(fd, "        \\legend{Best solution, Elite, FF, Cross Over}\n")

    for x=2:3
    write(fd, "        \\addplot3[\n")
    write(fd, "        color=green!25,]\n")
    write(fd, "        coordinates {\n")
    write(fd, "        ($(join(["$(x),$(yMin + i-1),$(minimum(gen))" for (i, gen) in enumerate(all_gen)], ")(")))\n")
    write(fd, "        };\n")
    write(fd, "        \n")
    end

    for (i, gen::GenVal) in enumerate(gen_val[2:end])

    write(fd, "        % Generation $(i)\n")
        if (!isempty(gen.cross_over))
    write(fd, "        \\addplot3[only marks,mark size=0.5pt,color=blue!50]\n")
    write(fd, "        coordinates {($(join(["$(1),$(yMin + i),$(val)" for val in gen.cross_over], ")(")))};\n")
        end
        if (!isempty(gen.ffr))
    write(fd, "        \\addplot3[only marks,mark size=0.5pt,color=red!50]\n")
    write(fd, "        coordinates {($(join(["$(2),$(yMin + i),$(val)" for val in gen.ffr], ")(")))};\n")
        end
        if (!isempty(gen.elite))
    write(fd, "        \\addplot3[only marks,mark size=0.5pt,color=green!50]\n")
    write(fd, "        coordinates {($(join(["$(3),$(yMin + i),$(val)" for val in gen.elite], ")(")))};\n")
        end
    write(fd, "        \n")
    
    end

    write(fd, "        \\end{axis}\n")
    write(fd, "    \\end{tikzpicture}\n")
    write(fd, "    \\caption{Genetic Algorithm convergence}\n")
    write(fd, "    \\label{fig:GA-converge}\n")
    write(fd, "\\end{figure}\n")

    close(fd)
end
    

function Plot_LaTeX(
        path::String,
        title::String,
        all_res::Vector{Tuple{String, Vector{Float64}}},
        inst_size::Vector{String}
    )

    n = length(inst_size)
    x = collect(1:n)

    bl = [0, typemax(Int64)] # Bottom Left corner
    tr = [n, 0]              # Top Right corner

    for (_, func) in all_res
        bl[2] = min(bl[2], minimum(func))
        tr[2] = round(Int64, max(tr[2], maximum(func)))
    end

    colors = ["blue", "red", "magenta", "purple", "cyan", "green", "orange", "brown", "yellow"]

    fd = open(path, "a")

# ==========< Head >==========
    write(fd, "\\begin{figure}[H]\n")
    write(fd, "    \\centering\n")
    write(fd, "        \\begin{tikzpicture}\n")

# ==========< Body >==========
    # Axes
    write(fd, "\n        % Axes\n")
    write(fd, "        \\draw[thick,->,>=latex]($(bl[1]),$(bl[2]))--($(tr[1]*1.1),$(bl[2])) node[below] {size};\n")
    write(fd, "        \\draw[thick,->,>=latex]($(bl[1]),$(bl[2]))--($(bl[1]),$(tr[2]*1.1)) node[left] {cpu time};\n")

    # function
    for (k, (name, func)) in enumerate(all_res)
        write(fd, "\n        % $(name)\n")
        write(fd, "        \\foreach \\a/\\b/\\c/\\d in {")
        for i=1:(length(func)-1)
            write(fd, "$(i)/$(func[i])/$(i+1)/$(func[i+1]), ")
        end
        write(fd, "} {\n")
        write(fd, "            \\draw[ultra thick,$(colors[((k-1)%length(colors))+1]),-](\\a,\\b)--(\\c,\\d);}\n")
    end

    # X labels
    write(fd, "\n        % X labels\n")
    write(fd, "        \\node[below] at ($(bl[1]),$(bl[2])) {$(bl[1])};\n")
    for (i, name) in enumerate(inst_size)
        write(fd, "        \\node[below] at ($(i),$(bl[2])) {$(name)};\n")
    end

    # Y labels
    write(fd, "\n        % Y labels\n")
    write(fd, "        \\node[left] at ($(bl[1]),$(bl[2])) {$(bl[2])};\n")
    write(fd, "        \\node[left] at ($(bl[1]),$(tr[2])) {$(tr[2])};\n")
    write(fd, "        \\node[left] at ($(bl[1]),$(tr[2]/4)) {$(tr[2]/4)};\n")
    write(fd, "        \\node[left] at ($(bl[1]),$(tr[2]/2)) {$(tr[2]/2)};\n")
    write(fd, "        \\node[left] at ($(bl[1]),$(tr[2]*3/4)) {$(tr[2]*3/4)};\n")

# ==========< Tail >==========
    write(fd, "    \\end{tikzpicture}\n")
    write(fd, "    \\caption{$(title)}\\label{FunH}\n")
    write(fd, "\\end{figure}\n\n\n")

    close(fd)
end

