using Logging
using LoggingExtras
function my_file_logger(filename)
    return FormatLogger(open(filename, "a")) do io, args
        println(io, args._module, " | ", "[", args.level, "] ", args.message)
    end
end
function my_combined_logger(filename)
    loggerFile = my_file_logger(filename)
    loggerREPL = ConsoleLogger()
    return TeeLogger(loggerFile,loggerREPL)
end
function log_structure(structure)
    type = typeof(structure)
    for i in 1:nfields(structure)
        val = getfield(structure,fieldname(type,i))
        str = string(fieldname(type,i))
        @info "$type: $str = $val"
    end
end
function start_logging(savepath;filename="log.txt")
    ispath(savepath) || mkpath(savepath)
    logfile = joinpath(savepath,filename)
    global_logger(my_combined_logger(logfile))
    @info "Threads used: $(Threads.nthreads())"
    @info "Date: $(Dates.now())"
end