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
    φ = params["-phi"]
    ϑ = params["-theta"]


    n𝐁 = [sin(ϑ)*cos(φ), sin(ϑ)*sin(φ), cos(ϑ)]  # the unit vector defining the direction of the magnetic field
    𝐒 = ["Sx", "Sy", "Sz"]  # vector consisting of the spin matrices 0.5σᵢ

    # automated MPO generation by a sum of operator expressions
    ampo = OpSum()
    nodes = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])
    # uniform magnetic field on all sites
    for n in nodes
        for s in 1:length(𝐒)
            ampo += mag_field*n𝐁[s], 𝐒[s], n[1]
        end
    end
    # loop over all bonds of the graph
    for b in graph
        dir = b.r2 .- b.r1  # the direction vector between lattice nodes
        dist = norm(dir)  # not really used here, but for completeness...
        𝐃 = absD.*dvec(dir)
        if dist ≈ 1
            # define uniform Heisenberg interaction
            for s in 1:length(𝐒)
                ampo .+= +heis_exch, 𝐒[s], b.s1, 𝐒[s], b.s2
                ampo .+= +K*n𝐁[s], 𝐒[s], b.s1, 𝐒[s], b.s2
            end

            # the DMI interaction
            for i in 1:3, j in 1:3, k in 1:3
                if epsilon(i,j,k) != 0  # only add nonzero terms
                    ampo .+= 𝐃[i]*epsilon(i,j,k), 𝐒[j], b.s1, 𝐒[k], b.s2
                end
            end
        end
    end
    # @show ampo
    return MPO(ampo,sites)
end