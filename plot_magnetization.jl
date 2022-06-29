using PyPlot, PyCall
function plot_magnetization(psi::MPS, graph::Graph, fn::String, ps::Int64)
    ns = unique([[(b.s1, b.r1) for b in graph]; [(b.s2, b.r2) for b in graph]])
    xs = [n[2][1] for n in ns]
    ys = [n[2][2] for n in ns]
    lobs = [expect(psi, lob) for lob in ["Sx", "Sy", "Sz"]]
    fig = plt.figure(figsize=2.0.*((3+3/8), (3+3/8)/1.2))
    plt.scatter(xs, ys, cmap="RdBu_r", c=lobs[3], marker="h", s=ps, vmin=-0.5, vmax=0.5)
    plt.quiver(xs, ys, lobs[1], lobs[2], scale=1, units="xy", pivot="middle", color="white")
    maxx = maximum(xs)
    maxy = maximum(ys)
    plt.ylim(-maxy-1/sqrt(3),maxy+1/sqrt(3))
    plt.xlim(-maxx-1,maxx+1)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(fn)
    plt.close()
end