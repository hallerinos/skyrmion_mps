# rhomboid lattice (zigzag order)
function rhomboid_zigzag(Nx::Int, Ny::Int; kwargs...) :: Tuple{Int64, Graph}
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Nz > 2)
    yperiodic == true ? error("PBC not yet implemented.") : nothing

    a1 = [1.0, 0.0, 0.]
    a2 = [0.5, sqrt(3)/2, 0.]

    lattPos = []
    for nx in 0:Nx-1, ny in 0:Ny-1
      pos = a1.*nx + a2.*ny
      append!(lattPos, [pos])
    end
    lattPos = unique(lattPos)

    nns = [a1, -a2, -a1+a2]  # consider only nearest neighbors

    latt = []
    for (idr, r) in enumerate(lattPos), (idrpr, rpr) in enumerate(lattPos)
      rmrpr = rpr .- r
      for (iddir, dir) in enumerate(nns)
        rmrpr ≈ dir ? append!(latt, [Bond(idr, idrpr, r, rpr, "$iddir")]) : nothing
      end
    end

    return length(lattPos), latt
end

# rhomboid lattice (spiral order)
function rhomboid_spiral(Nx::Int, Ny::Int; kwargs...) :: Tuple{Int64, Graph}
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Nz > 2)
    yperiodic == true ? error("PBC not yet implemented.") : nothing

    a1 = [1.0, 0.0, 0.]
    a2 = [0.5, sqrt(3)/2, 0.]
    a3 = [-0.5, sqrt(3)/2, 0.]

    lattPos = []
    pos = 0.0.*a1
    append!(lattPos, [pos])
    for n=1:Nx-1
      pos += a1
      append!(lattPos, [pos])
      for n2=1:2n-1
        pos += a2
        append!(lattPos, [pos])
      end
      for uc in [-a1, -a2, a1], n2=1:2n
        pos += uc
        append!(lattPos, [pos])
      end
    end
    # @show lattPos
    lattPos = unique(lattPos)
    # plt.scatter([b[1] for b in lattPos],[b[2] for b in lattPos])


    nns = [a1, -a2, -a1+a2]  # consider only nearest neighbors

    latt = []
    for (idr, r) in enumerate(lattPos), (idrpr, rpr) in enumerate(lattPos)
      rmrpr = rpr .- r
      for (iddir, dir) in enumerate(nns)
        rmrpr ≈ dir ? append!(latt, [Bond(idr, idrpr, r, rpr, "$iddir")]) : nothing
      end
    end

    return length(lattPos), latt
end

# triangular lattice with disk boundary conditions
function triangular_disk(Nx::Int, Ny::Int; kwargs...) :: Tuple{Int64, Graph}
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Nz > 2)
    yperiodic == true ? error("PBC not yet implemented.") : nothing

    a1 = [1.0, 0.0, 0.]
    a2 = [0.5, sqrt(3)/2, 0.]

    lattPos = []
    Nxs = -4Nx:4Nx
    for nx in Nxs, ny in Nxs
      pos = a1.*nx + a2.*ny
      if norm(pos) <= minimum(([Nx, Ny])./2)
        append!(lattPos, [pos])
      end
    end
    lattPos = unique(lattPos)
    for (idlp,lp) in enumerate(lattPos)
      x, y = getindex(lp,1), getindex(lp,2)
      plt.scatter(x,y)
      plt.text(x,y,"$idlp")
    end
    plt.savefig("nodes.png")
    plt.close()

    nns = [a1, -a2, -a1+a2]  # consider only nearest neighbors

    latt = []
    for (idr, r) in enumerate(lattPos), (idrpr, rpr) in enumerate(lattPos)
      rmrpr = rpr .- r
      for (iddir, dir) in enumerate(nns)
        rmrpr ≈ dir ? append!(latt, [Bond(idr, idrpr, r, rpr, "$iddir")]) : nothing
      end
    end

    return length(lattPos), latt
end

# square lattice
function square(Nx::Int, Ny::Int; kwargs...) :: Tuple{Int64, Graph}
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Nz > 2)
    yperiodic == true ? error("PBC not yet implemented.") : nothing

    a1 = [1.0, 0.0, 0.]
    a2 = [0.0, 1.0, 0.]

    lattPos = []
    for nx in 0:Nx-1, ny in 0:Ny-1
      pos = a1.*nx + a2.*ny
      append!(lattPos, [pos])
    end
    lattPos = unique(lattPos)

    nns = [a1, a2]  # consider only nearest neighbors

    latt = []
    for (idr, r) in enumerate(lattPos), (idrpr, rpr) in enumerate(lattPos)
      rmrpr = rpr .- r
      for (iddir, dir) in enumerate(nns)
        rmrpr ≈ dir ? append!(latt, [Bond(idr, idrpr, r, rpr, "$iddir")]) : nothing
      end
    end

    return length(lattPos), latt
end

# square lattice with disk boundary conditions
function square_disk(Nx::Int, Ny::Int; kwargs...) :: Tuple{Int64, Graph}
  yperiodic = get(kwargs, :yperiodic, false)
  yperiodic = yperiodic && (Nz > 2)
  yperiodic == true ? error("PBC not yet implemented.") : nothing

  a1 = [1.0, 0.0, 0.]
  a2 = [0.0, 1.0, 0.]

  lattPos = []
  Nxs = -4Nx:4Nx
  for nx in Nxs, ny in Nxs
    pos = a1.*nx + a2.*ny
    if norm(pos) <= minimum(([Nx, Ny])./2)
      append!(lattPos, [pos])
    end
  end
  lattPos = unique(lattPos)

  nns = [a1, -a2, -a1+a2]  # consider only nearest neighbors

  latt = []
  for (idr, r) in enumerate(lattPos), (idrpr, rpr) in enumerate(lattPos)
    rmrpr = rpr .- r
    for (iddir, dir) in enumerate(nns)
      rmrpr ≈ dir ? append!(latt, [Bond(idr, idrpr, r, rpr, "$iddir")]) : nothing
    end
  end

  return length(lattPos), latt
end

# kagome lattice
function kagome(Nx::Int, Ny::Int; kwargs...) :: Tuple{Int64, Graph}
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Nz > 2)
    yperiodic == true ? error("PBC not yet implemented.") : nothing

    a1 = [1.0, 0.0, 0.]
    a1 /= norm(a1)  # ensure normalized lattice vectors
    a2 = [0.5, sqrt(3)/2, 0.]
    a2 /= norm(a1)  # ensure normalized lattice vectors
    b1 = 2.0.*a1
    b2 = 2.0.*a2

    lattPos = []
    for nx in 0:Nx-1, ny in 0:Ny-1, uc in [0.0.*a1, a1, a2]
      if (nx==Nx-1)&&(ny<Nx-1)
        pos = b1.*nx + b2.*ny + a2
        append!(lattPos, [pos])
        pos = b1.*nx + b2.*ny
        append!(lattPos, [pos])
      elseif (ny==Ny-1)&&(nx<Nx-1)
        pos = b1.*nx + b2.*ny
        append!(lattPos, [pos])
        pos = b1.*nx + b2.*ny + a1
        append!(lattPos, [pos])
      elseif (ny==Ny-1)&&(nx==Nx-1)
        pos = b1.*nx + b2.*ny
        append!(lattPos, [pos])
      else
        pos = b1.*nx + b2.*ny + uc
        append!(lattPos, [pos])
      end
    end
    lattPos = unique(lattPos)

    nns = [a1, -a2, -a1+a2]  # consider only nearest neighbors

    latt = []
    for (idr, r) in enumerate(lattPos), (idrpr, rpr) in enumerate(lattPos)
      rmrpr = rpr .- r
      for (iddir, dir) in enumerate(nns)
        rmrpr ≈ dir ? append!(latt, [Bond(idr, idrpr, r, rpr, "$iddir")]) : nothing
      end
    end

    return length(lattPos), latt
end

# kagome lattice with disk boundary conditions
function kagome_disk(Nx::Int, Ny::Int; kwargs...) :: Tuple{Int64, Graph}
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Nz > 2)
    yperiodic == true ? error("PBC not yet implemented.") : nothing

    a1 = [1.0, 0.0, 0.]
    a1 /= norm(a1)  # ensure normalized lattice vectors
    a2 = [0.5, sqrt(3)/2, 0.]
    a2 /= norm(a1)  # ensure normalized lattice vectors
    b1 = 2.0.*a1
    b2 = 2.0.*a2

    lattPos = []
    Nxs = -4Nx:4Nx
    for nx in Nxs, ny in Nxs, uc in [0.0.*a1, a1, a2]
      pos = b1.*nx + b2.*ny + uc
      if norm(pos) <= minimum(([Nx-1, Ny-1])./2)
        append!(lattPos, [pos])
      end
    end
    lattPos = unique(lattPos)

    nns = [a1, -a2, -a1+a2]  # consider only nearest neighbors

    latt = []
    for (idr, r) in enumerate(lattPos), (idrpr, rpr) in enumerate(lattPos)
      rmrpr = rpr .- r
      for (iddir, dir) in enumerate(nns)
        rmrpr ≈ dir ? append!(latt, [Bond(idr, idrpr, r, rpr, "$iddir")]) : nothing
      end
    end

    return length(lattPos), latt
end

# honeycomb lattice
function honeycomb(Nx::Int, Ny::Int; kwargs...) :: Tuple{Int64, Graph}
    yperiodic = get(kwargs, :yperiodic, false)
    yperiodic = yperiodic && (Nz > 2)
    yperiodic == true ? error("PBC not yet implemented.") : nothing

    b1 = 0.5.*[3.0, +sqrt(3), 0.]
    b2 = 0.5.*[3.0, -sqrt(3), 0.]
    a1 = [-1, +sqrt(3), 0.]
    a1 /= norm(a1)  # ensure normalized lattice vectors
    a2 = [+1, +sqrt(3), 0.]
    a2 /= norm(a2)  # ensure normalized lattice vectors

    lattPos = []
    for nx in 0:Nx-1, ny in 0:Ny-1, uc in [a1, a2]
      if (nx==0&&ny==0)
        pos = b1.*nx + b2.*ny + a2
      elseif (nx==Nx-1&&ny==Ny-1)
        pos = b1.*nx + b2.*ny + a1
      else
        pos = b1.*nx + b2.*ny + uc
      end
      append!(lattPos, [pos])
    end
    lattPos = unique(lattPos)

    nns = [a1-a2, a1, a2]  # consider only nearest neighbors

    latt = []
    for (idr, r) in enumerate(lattPos), (idrpr, rpr) in enumerate(lattPos)
      rmrpr = rpr .- r
      for (iddir, dir) in enumerate(nns)
        rmrpr ≈ dir ? append!(latt, [Bond(idr, idrpr, r, rpr, "$iddir")]) : nothing
      end
    end

    return length(lattPos), latt
end