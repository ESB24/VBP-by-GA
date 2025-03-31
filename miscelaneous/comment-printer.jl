function A(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return "########") 
    (height == 4) && (return "##    ##") 
    (height == 5) && (return "##    ##") 
end
function B(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "######  ") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return "####### ") 
    (height == 4) && (return "##    ##") 
    (height == 5) && (return "####### ") 
end

function C(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return "##      ") 
    (height == 4) && (return "##    ##") 
    (height == 5) && (return " ###### ") 
end

function D(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "######  ") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return "##     #") 
    (height == 4) && (return "##    ##") 
    (height == 5) && (return "######  ") 
end

function E(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "########") 
    (height == 2) && (return "##      ") 
    (height == 3) && (return "#####   ") 
    (height == 4) && (return "##      ") 
    (height == 5) && (return "########") 
end

function F(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "########") 
    (height == 2) && (return "##      ") 
    (height == 3) && (return "#####   ") 
    (height == 4) && (return "##      ") 
    (height == 5) && (return "##      ") 
end

function G(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "##      ") 
    (height == 3) && (return "##  ### ") 
    (height == 4) && (return "##    ##") 
    (height == 5) && (return " ###### ") 
end

function H(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "##    ##") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return "########") 
    (height == 4) && (return "##    ##") 
    (height == 5) && (return "##    ##") 
end

function I(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "########") 
    (height == 2) && (return "   ##   ") 
    (height == 3) && (return "   ##   ") 
    (height == 4) && (return "   ##   ") 
    (height == 5) && (return "########") 
end

function J(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " #######") 
    (height == 2) && (return "    ##  ") 
    (height == 3) && (return "    ##  ") 
    (height == 4) && (return "##  ##  ") 
    (height == 5) && (return "  ###   ") 
end

function K(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "##    ##") 
    (height == 2) && (return "##  ##  ") 
    (height == 3) && (return "####    ") 
    (height == 4) && (return "##  ##  ") 
    (height == 5) && (return "##    ##") 
end

function L(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "##      ") 
    (height == 2) && (return "##      ") 
    (height == 3) && (return "##      ") 
    (height == 4) && (return "##      ") 
    (height == 5) && (return "########") 
end

function M(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "##    ##") 
    (height == 2) && (return "###  ###") 
    (height == 3) && (return "## ## ##") 
    (height == 4) && (return "##    ##") 
    (height == 5) && (return "##    ##") 
end

function N(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "##    ##") 
    (height == 2) && (return "###   ##") 
    (height == 3) && (return "## ## ##") 
    (height == 4) && (return "##   ###") 
    (height == 5) && (return "##    ##") 
end

function O(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return "##    ##") 
    (height == 4) && (return "##    ##") 
    (height == 5) && (return " ###### ") 
end

function P(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "####### ") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return "####### ") 
    (height == 4) && (return "##      ") 
    (height == 5) && (return "##      ") 
end

function Q(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return "## ## ##") 
    (height == 4) && (return "##  ### ") 
    (height == 5) && (return " #### ##") 
end

function R(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "####### ") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return "####### ") 
    (height == 4) && (return "##  ##  ") 
    (height == 5) && (return "##   ## ") 
end

function S(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " #######") 
    (height == 2) && (return "##      ") 
    (height == 3) && (return " ###### ") 
    (height == 4) && (return "      ##") 
    (height == 5) && (return "####### ") 
end

function T(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "########") 
    (height == 2) && (return "   ##   ") 
    (height == 3) && (return "   ##   ") 
    (height == 4) && (return "   ##   ") 
    (height == 5) && (return "   ##   ") 
end

function U(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "##    ##") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return "##    ##") 
    (height == 4) && (return "##    ##") 
    (height == 5) && (return " ###### ") 
end

function V(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "##    ##") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return " ##  ## ") 
    (height == 4) && (return " ##  ## ") 
    (height == 5) && (return "   ##   ") 
end

function W(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "##    ##") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return "## ## ##") 
    (height == 4) && (return "###  ###") 
    (height == 5) && (return "##    ##") 
end

function X(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "##    ##") 
    (height == 2) && (return " ##  ## ") 
    (height == 3) && (return "   ##   ") 
    (height == 4) && (return " ##  ## ") 
    (height == 5) && (return "##    ##") 
end

function Y(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "##    ##") 
    (height == 2) && (return " ##  ## ") 
    (height == 3) && (return "   ##   ") 
    (height == 4) && (return "   ##   ") 
    (height == 5) && (return "   ##   ") 
end

function Z(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "########") 
    (height == 2) && (return "     ## ") 
    (height == 3) && (return "   ##   ") 
    (height == 4) && (return " ##     ") 
    (height == 5) && (return "########") 
end

function zero(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "##   ###") 
    (height == 3) && (return "## ## ##") 
    (height == 4) && (return "###   ##") 
    (height == 5) && (return " ###### ") 
end

function one(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " /###   ") 
    (height == 2) && (return " # ##   ") 
    (height == 3) && (return "   ##   ") 
    (height == 4) && (return "   ##   ") 
    (height == 5) && (return " ###### ") 
end

function two(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "      ##") 
    (height == 3) && (return "  ####  ") 
    (height == 4) && (return "##      ") 
    (height == 5) && (return " ###### ") 
end

function three(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "      ##") 
    (height == 3) && (return "  ##### ") 
    (height == 4) && (return "      ##") 
    (height == 5) && (return " ###### ") 
end

function four(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "   #### ") 
    (height == 2) && (return "  ## ## ") 
    (height == 3) && (return " ##  ## ") 
    (height == 4) && (return "########") 
    (height == 5) && (return "     ## ") 
end

function five(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "##      ") 
    (height == 3) && (return " ###### ") 
    (height == 4) && (return "      ##") 
    (height == 5) && (return " ###### ") 
end

function six(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "##      ") 
    (height == 3) && (return "####### ") 
    (height == 4) && (return "##    ##") 
    (height == 5) && (return " ###### ") 
end

function seven(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "########") 
    (height == 2) && (return "     ## ") 
    (height == 3) && (return " #####= ") 
    (height == 4) && (return " ##     ") 
    (height == 5) && (return "##      ") 
end

function eight(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return " ###### ") 
    (height == 4) && (return "##    ##") 
    (height == 5) && (return " ###### ") 
end

function nine(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return " ###### ") 
    (height == 2) && (return "##    ##") 
    (height == 3) && (return " #######") 
    (height == 4) && (return "      ##") 
    (height == 5) && (return " ###### ") 
end

function space(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "        ") 
    (height == 2) && (return "        ") 
    (height == 3) && (return "        ") 
    (height == 4) && (return "        ") 
    (height == 5) && (return "        ") 
end

function dot(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "        ") 
    (height == 2) && (return "        ") 
    (height == 3) && (return "        ") 
    (height == 4) && (return " /##\\   ") 
    (height == 5) && (return " \\##/   ") 
end

function dash(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "        ") 
    (height == 2) && (return "        ") 
    (height == 3) && (return "  ####  ") 
    (height == 4) && (return "        ") 
    (height == 5) && (return "        ") 
end

function greater_than(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "        ") 
    (height == 2) && (return "  ##    ") 
    (height == 3) && (return "    ##  ") 
    (height == 4) && (return "  ##    ") 
    (height == 5) && (return "        ") 
end

function lesser_than(height::Int64)
    (height < 1 || 5 < height) && (return "")

    (height == 1) && (return "        ") 
    (height == 2) && (return "    ##  ") 
    (height == 3) && (return "  ##    ") 
    (height == 4) && (return "    ##  ") 
    (height == 5) && (return "        ") 
end

function letter(c::Char, height::Int64)
    (c == 'a') && (return A(height))
    (c == 'b') && (return B(height))
    (c == 'c') && (return C(height))
    (c == 'd') && (return D(height))
    (c == 'e') && (return E(height))
    (c == 'f') && (return F(height))
    (c == 'g') && (return G(height))
    (c == 'h') && (return H(height))
    (c == 'i') && (return I(height))
    (c == 'j') && (return J(height))
    (c == 'k') && (return K(height))
    (c == 'l') && (return L(height))
    (c == 'm') && (return M(height))
    (c == 'n') && (return N(height))
    (c == 'o') && (return O(height))
    (c == 'p') && (return P(height))
    (c == 'q') && (return Q(height))
    (c == 'r') && (return R(height))
    (c == 's') && (return S(height))
    (c == 't') && (return T(height))
    (c == 'u') && (return U(height))
    (c == 'v') && (return V(height))
    (c == 'w') && (return W(height))
    (c == 'x') && (return X(height))
    (c == 'y') && (return Y(height))
    (c == 'z') && (return Z(height))

    (c == '0') && (return zero(height))
    (c == '1') && (return one(height))
    (c == '2') && (return two(height))
    (c == '3') && (return three(height))
    (c == '4') && (return four(height))
    (c == '5') && (return five(height))
    (c == '6') && (return six(height))
    (c == '7') && (return seven(height))
    (c == '8') && (return eight(height))
    (c == '9') && (return nine(height))
    
    (c == ' ') && (return space(height))
    (c == '.') && (return dot(height))
    (c == '-') && (return dash(height))
    (c == '>') && (return greater_than(height))
    (c == '<') && (return lesser_than(height))

    return "<Error>"
end

function print_comment(str::String; path::String = "./miscelaneous/output.txt", max_char::Int64 = 11)
    (length(str) > max_char) && (println("Warning: string length > 10!"))

    comment_width = 8 * max_char + 2 * (max_char - 1)
    word_width = 8 * length(str) + 2 * (length(str) - 1)
    tab_width = floor(Int64, (comment_width - word_width) / 2)

    fd = open(path, "a+")

    write(fd, "## $("=" ^ (comment_width + 2)) ##\n")
    
    for i=1:5
        write(fd, "##  $(" " ^ tab_width)$(join(["$(letter(c, i))" for c in str], "  "))$(" " ^ tab_width)  ##\n")
    end

    write(fd, "## $("=" ^ (comment_width + 2)) ##\n\n\n")

    close(fd)
end

# print_comment("instance")
# print_comment("heuristic")
# print_comment("01-lp")
# print_comment("01-lp warm")
# print_comment("benchmark")
# print_comment("plot")
# print_comment("constructor")
# print_comment("fitness")
# print_comment("display")
# print_comment("misc")
# print_comment("lower bound")
print_comment("1d-bp>gr")
print_comment("1d-bfd>gr")
