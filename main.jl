using ITensors, PyPlot, PyCall, LinearAlgebra
# pygui(true)

include("bonds.jl")
include("lattices.jl")
include("generate_MPO.jl")

# --------------- system settings ----------------
Ny = 8  # lattice dimension in y-direction
Nx = Ny  # lattice dim. in x direction (redundant for some lattices)
ns = 1  # how many low lying energy states to compute
absD = 1.0  # DMI amplitude
J = -1/2*absD  # Heisenberg interaction
B = -1/2*absD  # magnetic field
K = -0  # uniaxial anisotropy
φ = 0  # polar angle of the magnetic field
ϑ = 0  # axial angle of the magnetic field
params = Dict("-B0"=>B, "-D"=>absD, "-J"=>J, "-K"=>K, "-phi"=>φ, "-theta"=>ϑ)

# --------------- MPS settings ----------------
M = 32  # set the maximum bond dimension
Ns = 30  # set the maximum number of sweeps
etresh = 1e-6  # naïve stopping criterion
restart = false  # restarting from states MPS
outputlevel = 1  # increase output from 0,1,2

# --------------- initialization ----------------
sweeps = Sweeps(Ns)  # initialize sweeps object
maxdim!(sweeps, M)
N, graph = triangular_disk(Nx, Ny)  # see available graphs in lattices.jl
psi0 = nothing
if restart
    psi0 = states[1]
    sites = siteinds(psi0)
else
    sites = siteinds("S=1/2", N)
    psi0 = randomMPS(sites, M)
end
H = generate_MPO(graph, sites, params)
obs = DMRGObserver(; energy_tol=etresh)

# --------------- perform 2-site DMRG ----------------
psi = psi0
states = Vector{MPS}(undef, ns)
energy = []
es = []
for n=1:ns
    global ene, psi = (n==1 ? dmrg(H, psi, sweeps, observer=obs, outputlevel=outputlevel) : dmrg(H, states[1:n-1], psi, sweeps, weight=1000, observer=obs, outputlevel=outputlevel))
    append!(es, ene)
    global states[n] = psi
end
append!(energy, [es])

# --------------- plot the magnetization ----------------
ns = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])
xs = [n[2][1] for n in ns]
ys = [n[2][2] for n in ns]
for (i, state) in enumerate(states)
    lobs = [expect(state, lob) for lob in ["Sx", "Sy", "Sz"]]
    fig = figure(figsize=2.0.*((3+3/8), (3+3/8)/1.25))
    scatter(xs, ys, cmap="RdBu_r", c=lobs[3], marker="h", s=1600, vmin=-0.5, vmax=0.5)
    quiver(xs, ys, lobs[1],lobs[2], scale=.5, units="xy", pivot="middle", color="white")
    maxx = maximum(xs)
    maxy = maximum(ys)
    ylim(-maxy-1/sqrt(3),maxy+1/sqrt(3))
    xlim(-maxx-1,maxx+1)
    axis("off")
    savefig("magnetization_B$(B)_J$(J)_K$(K)_M$(M)_state_$(i).pdf")
end
0;