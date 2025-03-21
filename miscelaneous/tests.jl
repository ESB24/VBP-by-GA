# ==============================< ANSI Escape Codes >==============================

begin
    for i=0:9
        println("$i: $("-"^50)")
    end
    print("Temporary message...")
    sleep(1)
    print("\u001b[2A")  # Move up a line and clear it
    println("Final message.")
end

using Random

begin
    print("$("-"^200)")
    while true
        sleep(0.5)
        if rand() ≤ 0.05
            tmp = rand(1:5)
            print("\u001b[$(tmp+1)D")
            print("\u001b[35m")
            print("|")
            print("\u001b[31m")
            print("$("<"^(tmp+1))")
            print("\u001b[$(tmp+1)D")
        else
            print("\u001b[$(2)D\u001b[32m")
            print(">|")
            print("\u001b[35m")
            print("|")
        end
    end
end

begin
    reset = "\x1b[0m"
    print("\x1b[1;4;38;2;255;0;0;48;2;0;0;255mTest$reset")
end

Threads.@spawn begin
    println("Running on thread $(Threads.threadid())")
    print("\x1b[s")
    i = 0
    while true
        sleep(0.2)
        i = (i+1)

        if (i%4 == 0)
            print("\x1b[u")
            print("◐")
        elseif (i%4 == 1)
            print("\x1b[u")
            print("◑")
        elseif (i%4 == 2)
            print("\x1b[u")
            print("◓")
        elseif (i%4 == 3)
            print("\x1b[u")
            print("◒")
        end

        i%20 == 0 && print("-"^round(Int64, i/20))
    end
end

import Base.print

struct tFrame
    length::Int64
    height::Int64
end

# struct tRGB
#     r::Int8
#     g::Int8
#     b::Int8
# end

function init_tFrame(length::Int64 = 150, height::Int64 = 10)
    println("#$("-"^length)#")
    for i=1:height
        println("|$(" "^length)|")
    end
    println("#$("-"^length)#")
    return tFrame(length, height)
end

function print(
        frame       ::tFrame                                , 
        x           ::Int64                                 , 
        y           ::Int64                                 , 
        str         ::String                                ; 
        bold        ::Bool                      = false     , 
        italic      ::Bool                      = false     ,
        underlined  ::Bool                      = false     ,
        txt_color   ::Union{String, Nothing}    = nothing   ,
        bg_color    ::Union{String, Nothing}    = nothing   ,
    )
    # Move to location
    print("\x1b[$(frame.height - y + 1)F\x1b[$(x+2)G")

    font_opt = ""

    bold        && (font_opt = "1")
    italic      && ((font_opt == "") ? (font_opt = "2") : (font_opt = font_opt * ";2"))
    underlined  && ((font_opt == "") ? (font_opt = "3") : (font_opt = font_opt * ";3"))

    # Print string
    print("\x1b[$(font_opt)m$(str)")

    # Reset
    print("\x1b[$(frame.height - y + 1)E\x1b[0m")
end

begin
    frame = init_tFrame(200)
    print(frame, 0, 0, "Hello World!")
    print(frame, 199, 9, "H")
end


