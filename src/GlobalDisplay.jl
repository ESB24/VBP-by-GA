fd_verbose = open("../debug/verbose_out.txt", "a+")

global VERBOSE              = false
global VERBOSE_IO           = fd_verbose                    # Default = stdout
global VERBOSE_IO_IS_STDOUT = (VERBOSE_IO == stdout)

# ==============================< ANSI Escape Codes >==============================

# ANSI Reset font reset
global ANSI_reset    = "\x1b[0m"

# ANSI Basic text colors
global ANSI_black    = "\x1b[30m"
global ANSI_red      = "\x1b[31m"
global ANSI_green    = "\x1b[32m"
global ANSI_yellow   = "\x1b[33m"
global ANSI_blue     = "\x1b[34m"
global ANSI_magenta  = "\x1b[35m"
global ANSI_cyan     = "\x1b[36m"
global ANSI_white    = "\x1b[37m"

# ANSI Basic background colors
global ANSI_bg_black    = "\x1b[40m"
global ANSI_bg_red      = "\x1b[41m"
global ANSI_bg_green    = "\x1b[42m"
global ANSI_bg_yellow   = "\x1b[43m"
global ANSI_bg_blue     = "\x1b[44m"
global ANSI_bg_magenta  = "\x1b[45m"
global ANSI_bg_cyan     = "\x1b[46m"
global ANSI_bg_white    = "\x1b[47m"

# ANSI Basic font transformation
global ANSI_bold        = "\x1b[1m"
global ANSI_italic      = "\x1b[1m"
global ANSI_underlined  = "\x1b[1m"
global ANSI_strike      = "\x1b[1m"

function println_verbose(txt::String = "", args::Union{String, NTuple{N, String}} = "") where N
    if VERBOSE
        println(VERBOSE_IO, "$(VERBOSE_IO_IS_STDOUT ? join(args) : "")$txt$(VERBOSE_IO_IS_STDOUT ? ANSI_reset : "")")

        if !VERBOSE_IO_IS_STDOUT
            flush(VERBOSE_IO)
        end
    end
end

function print_verbose(txt::String = "", args::Union{String, NTuple{N, String}} = "") where N
    if VERBOSE
        print(VERBOSE_IO, "$(VERBOSE_IO_IS_STDOUT ? join(args) : "")$txt$(VERBOSE_IO_IS_STDOUT ? ANSI_reset : "")")

        if !VERBOSE_IO_IS_STDOUT
            flush(VERBOSE_IO)
        end
    end
end
;
