using Pkg; Pkg.activate("."); Pkg.instantiate()
using UnquenchingQCD
using Parameters
using JLD
using Dates
include("scripts/IterationHelper.jl")
include("scripts/logging.jl")

function main(pathtodata)
    for group in filter(isdir,readdir(pathtodata,join=true))
        for runs in readdir(group,join=true)
            for loadpath in readdir(runs,join=true)
                
                @show loadpath
                
                savepath = joinpath("output",splitpath(loadpath)[end-2:end]...)
                ispath(savepath) && continue
                start_logging(savepath)
                
                Δg2 = 0.0
                ΔNf = 0.0
                GIR = load(joinpath(loadpath,"parameters.jld"))["GIR"]
                κGl = load(joinpath(loadpath,"theory.jld"))["κGl"]
                κGh = load(joinpath(loadpath,"theory.jld"))["κGh"]
                
                @info "Threads : $(Threads.nthreads())"
                @info "Date: $(Dates.now())"
                @info "Δg2 = $Δg2" 
                @info "ΔNf = $ΔNf" 

                
                contains(loadpath,"prop")   && mainPropagators(loadpath,savepath;Δg2,ΔNf,GIR,κGl,κGh,relax=1)
                contains(loadpath,"vertex") && mainSystem(loadpath,savepath;Δg2,ΔNf,GIR,κGl,κGh,relax=1,quarkrelax=1)
            end
        end
    end
end

pathtodata = "input"
@time main(pathtodata)