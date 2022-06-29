using ITensors, Observers, ITensorTDVP, DataFrames, JLD, UUIDs, CSV, PyPlot

include("bonds.jl")
include("lattices.jl")
include("generate_MPO.jl")
include("dmrg1.jl")
include("projmpo1.jl")
include("observer_overload.jl")

# --------------- system settings ----------------
dir = "./output"
if !ispath(dir)
    mkpath(dir)
end
Ny = 8  # lattice dimension in y-direction
Nx = Ny  # lattice dim. in x direction (redundant for some lattices)
absD = 1.0  # DMI amplitude
J = -1/2*absD  # Heisenberg interaction
K = -0  # uniaxial anisotropy
φ = 0  # polar angle of the magnetic field
ϑ = 0  # axial angle of the magnetic field
p = Dict("-D"=>absD, "-J"=>J, "-K"=>K, "-phi"=>φ, "-theta"=>ϑ)
etresh = 1e-6  # naïve stopping criterion
# --------------- MPS settings ----------------
M = 32  # set the maximum bond dimension
Ns = 1000  # set the maximum number of sweeps
outputlevel = 2  # increase output from 0,1,2
# --------------- initialization ----------------
sweeps = Sweeps(Ns)  # initialize sweeps object
maxdim!(sweeps, M)
cutoff!(sweeps, 1e-12)
N, graph = triangular_disk(Nx, Ny)  # see available graphs in lattices.jl
nodes = sort(unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]]))
nps = getindex.(nodes,2)
sites = siteinds("S=1/2", N)
obs_dmrg = DMRGObserver(; energy_tol=etresh)
state = ["Dn" for n=1:N]
# --------------- main loop ---------------------
for B in (1:10).*(-absD/10)
    df = DataFrame()
    p["-B0"] = B
    psi = randomMPS(sites; linkdims=M)
    H = generate_MPO(graph, sites, p)
    for s in ["Sx","Sy","Sz"]
        df[!, s] = expect(psi,s)
    end
    ene = inner(psi', H, psi)
    df[!, "X"] = getindex.(nps,1)
    df[!, "Y"] = getindex.(nps,2)
    df[!, "Z"] = getindex.(nps,3)
    df[!, "sweep"] = [0 for i=1:size(df,1)]
    df[!, "half_sweep"] = [0 for i=1:size(df,1)]
    df[!, "bond"] = [0 for i=1:size(df,1)]
    df[!, "E"] = [real(ene) for i=1:size(df,1)]

    # --------------- perform 1-site DMRG ----------------
    ene, psi = dmrg1(H, psi, sweeps, weight=1000, observer=obs_dmrg, nps=nps, df=df, outputlevel=outputlevel)

    # --------------- plot the magnetization ----------------
    include("plot_magnetization.jl")
    fn = "$dir/$(uuid1())"
    plot_magnetization(psi, graph, "$(fn).png", 2000)
    jldopen("$(fn).jld", "w") do file
        write(file, "fn", fn)
        write(file, "psi", psi)
        write(file, "graph", graph)
        write(file, "parameters", p)
    end
    CSV.write("$(fn).csv", df)
end