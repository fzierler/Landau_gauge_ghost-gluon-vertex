using Pkg; Pkg.activate(".")
using UnquenchingQCD
using Plots
gr()
pathtodata = "./input"

pathsSUN = [joinpath("SU$(N)") for N in 2:5 ]
pathsSON = [joinpath("SO$(N)") for N in 6:9]
pathsSpN = [joinpath("Sp$(N)") for N in 4:2:8]

paths_all = [pathsSUN,pathsSON,pathsSpN]

for paths in paths_all
    for j in paths
        paths_theory = readdir(joinpath(pathtodata,j),join=true)
        for p0 in paths_theory
            subpaths = readdir(p0,join=true)
            for p in subpaths
                isdir(p) && UnquenchingQCD.redoplots(p)
            end
        end
    end
end