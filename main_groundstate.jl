using ITensors, LinearAlgebra

include("bonds.jl")
include("lattices.jl")
include("generate_MPO.jl")
include("dmrg1.jl")
include("projmpo1.jl")

let
    # --------------- system settings ----------------
    Ny = 8  # lattice dimension in y-direction
    Nx = Ny  # lattice dim. in x direction (redundant for some lattices)
    absD = 1.0  # DMI amplitude
    J = -1/2*absD  # Heisenberg interaction
    K = -0  # uniaxial anisotropy
    φ = 0  # polar angle of the magnetic field
    ϑ = 0  # axial angle of the magnetic field
    params = Dict("-D"=>absD, "-J"=>J, "-K"=>K, "-phi"=>φ, "-theta"=>ϑ)
    for B in (1:20).*(-absD/10)
        params["-B0"] = B
        # --------------- MPS settings ----------------
        M = 128  # set the maximum bond dimension
        Ns = 1000  # set the maximum number of sweeps
        etresh = 1e-6  # naïve stopping criterion
        outputlevel = 1  # increase output from 0,1,2

        # --------------- initialization ----------------
        sweeps = Sweeps(Ns)  # initialize sweeps object
        maxdim!(sweeps, M)
        cutoff!(sweeps, 1e-12)
        N, graph = triangular_disk(Nx, Ny)  # see available graphs in lattices.jl
        sites = siteinds("S=1", N)
        psi = randomMPS(sites, M)
        H = generate_MPO(graph, sites, params)
        obs = DMRGObserver(; energy_tol=etresh)

        # --------------- perform 1-site DMRG ----------------
        ene, psi = dmrg1(H, psi, sweeps, weight=1000, observer=obs, outputlevel=outputlevel)

        # --------------- plot the magnetization ----------------
        include("plot_magnetization.jl")
        fn = "magnetization_B$(B)_J$(J)_K$(K)_D$(absD)_M$(M).jpg"
        plot_magnetization(psi, graph, fn, 1200)
    end
end