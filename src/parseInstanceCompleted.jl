include("Instance.jl")

function parseHardInstance_completed(path::String)

    println("Pasring -> $path")

    fd  = open(path, "r")

    _ = readline(fd) # first empty line

    O       ::Int64 = parse(Int64, readline(fd)) # Number of output
    R       ::Int64 = parse(Int64, readline(fd)) # Number of route
    Lmax    ::Int64 = parse(Int64, readline(fd)) # Maximum load of each output
    opti    ::Int64 = parse(Int64, readline(fd)) # Optimal number of session

    # println("O: $O, R: $R, Lmax: $Lmax, opti: $opti")

    _ = readline(fd) # Second empty line

    sol = [Session(Lmax, [], zeros(Int64, O))]

    routes = [] # all routes

    for i=1:R
        line = readline(fd)

        if line == ""
            # Eat all empty lines
            while line == ""
                line = readline(fd)
            end
            # Open a new session
            push!(sol, Session(Lmax, [], zeros(Int64, O)))
        end

        # println("line -> $line")

        session = sol[end] # get Last session

        id, assignment, mails = split(line, ": ")

        id          = parse(Int64, id)                                  # id of the route
        assignment  = parse.(Int64, split(assignment, ", ")[1:end-1])   # assignment vector
        mails       = parse.(Int64, split(mails, ", ")[1:end-1])        # mails of the route

        push!(session.route, Route(id, assignment, ntuple(x->mails[x], length(mails)))) # add route to the session 
        compute_output!(session)  # update outpout loads

        push!(routes, Route(id, O, ntuple(x->mails[x], length(mails))))
    end

    close(fd)

    sort!(routes, by=(x->x.id))
    return O, R, Lmax, opti, sol, Instance(Lmax, O, R, routes)
end

O, R, Lmax, opti, optimal_solution, instance = parseHardInstance_completed("../data/HardIndus/instanceIndus_1_O20_R60_C40_opt_2.txt")

println("\x1b[36mSolution optimale:\x1b[0m")
println(optimal_solution)
println("\x1b[36mInstance:\x1b[0m")
println(instance)