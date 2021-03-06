using LinearAlgebra
# some convenient definitions
dvec(r) = cross([0,0,1], r)
function epsilon(i,j,k)
    if [i,j,k] in [[1,2,3], [3,1,2], [2,3,1]]
        return +1
    elseif [i,j,k] in [[2,1,3], [3,2,1], [1,3,2]]
        return -1
    else 
        return 0
    end
end

function generate_MPO(graph::Graph, sites::Vector{Index{Int64}}, params::Dict)
    # get parameters
    mag_field = params["-B0"]
    heis_exch = params["-J"]
    absD = params["-D"]
    K = params["-K"]
    ฯ = params["-phi"]
    ฯ = params["-theta"]


    n๐ = [sin(ฯ)*cos(ฯ), sin(ฯ)*sin(ฯ), cos(ฯ)]  # the unit vector defining the direction of the magnetic field
    ๐ = ["Sx", "Sy", "Sz"]  # vector consisting of the spin matrices 0.5ฯแตข

    # automated MPO generation by a sum of operator expressions
    ampo = OpSum()
    nodes = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])
    # uniform magnetic field on all sites
    for n in nodes
        for s in 1:length(๐)
            ampo += mag_field*n๐[s], ๐[s], n[1]
        end
    end
    # loop over all bonds of the graph
    for b in graph
        dir = b.r2 .- b.r1  # the direction vector between lattice nodes
        dist = norm(dir)  # not really used here, but for completeness...
        ๐ = absD.*dvec(dir)
        if dist โ 1
            # define uniform Heisenberg interaction
            for s in 1:length(๐)
                ampo .+= +heis_exch, ๐[s], b.s1, ๐[s], b.s2
                ampo .+= +K*n๐[s], ๐[s], b.s1, ๐[s], b.s2
            end

            # the DMI interaction
            for i in 1:3, j in 1:3, k in 1:3
                if epsilon(i,j,k) != 0  # only add nonzero terms
                    ampo .+= ๐[i]*epsilon(i,j,k), ๐[j], b.s1, ๐[k], b.s2
                end
            end
        end
    end
    # @show ampo
    return MPO(ampo,sites)
end