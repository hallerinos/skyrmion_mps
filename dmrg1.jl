import ITensors.@debug_check
import ITensors.@timeit_debug
import ITensors.@printf
import ITensors.factorize
import ITensors.leftlim
import ITensors.setleftlim!
import ITensors.rightlim
import ITensors.setrightlim!
import ITensors.orthocenter
import ITensors.check_hascommoninds
import KrylovKit.eigsolve

function dmrg1(H::MPO, psi0::MPS, sweeps::Sweeps; kwargs...)::Tuple{Number, MPS, Vector{Float64}}
    check_hascommoninds(siteinds, H, psi0)
    check_hascommoninds(siteinds, H, psi0')
    # Permute the indices to have a better memory layout
    # and minimize permutations
    H = permute(H, (linkind, siteinds, linkind))
    PH = ProjMPO1(H)

    if length(psi0) == 1
    error(
        "`dmrg` currently does not support system sizes of 1. You can diagonalize the MPO tensor directly with tools like `LinearAlgebra.eigen`, `KrylovKit.eigsolve`, etc.",
    )
    end

    @debug_check begin
    # Debug level checks
    # Enable with ITensors.enable_debug_checks()
    checkflux(psi0)
    checkflux(PH)
    end

    which_decomp::Union{String,Nothing} = get(kwargs, :which_decomp, nothing)
    svd_alg::String = get(kwargs, :svd_alg, "divide_and_conquer")
    obs = get(kwargs, :observer, NoObserver())
    outputlevel::Int = get(kwargs, :outputlevel, 1)

    write_when_maxdim_exceeds::Union{Int,Nothing} = get(
    kwargs, :write_when_maxdim_exceeds, nothing
    )

    # eigsolve kwargs
    eigsolve_tol::Float64 = get(kwargs, :eigsolve_tol, 1e-14)
    eigsolve_krylovdim::Int = get(kwargs, :eigsolve_krylovdim, 3)
    eigsolve_maxiter::Int = get(kwargs, :eigsolve_maxiter, 1)
    eigsolve_verbosity::Int = get(kwargs, :eigsolve_verbosity, 0)

    # TODO: add support for non-Hermitian DMRG
    ishermitian::Bool = get(kwargs, :ishermitian, true)

    # TODO: add support for targeting other states with DMRG
    # (such as the state with the largest eigenvalue)
    # get(kwargs, :eigsolve_which_eigenvalue, :SR)
    eigsolve_which_eigenvalue::Symbol = :SR

    # TODO: use this as preferred syntax for passing arguments
    # to eigsolve
    #default_eigsolve_args = (tol = 1e-14, krylovdim = 3, maxiter = 1,
    #                         verbosity = 0, ishermitian = true,
    #                         which_eigenvalue = :SR)
    #eigsolve = get(kwargs, :eigsolve, default_eigsolve_args)

    # Keyword argument deprecations
    if haskey(kwargs, :maxiter)
    error("""maxiter keyword has been replaced by eigsolve_krylovdim.
                Note: compared to the C++ version of ITensor,
                setting eigsolve_krylovdim 3 is the same as setting
                a maxiter of 2.""")
    end

    if haskey(kwargs, :errgoal)
    error("errgoal keyword has been replaced by eigsolve_tol.")
    end

    if haskey(kwargs, :quiet)
    error("quiet keyword has been replaced by outputlevel")
    end

    psi = copy(psi0)
    N = length(psi)

    if !isortho(psi) || orthocenter(psi) != 1
    orthogonalize!(psi, 1)
    end
    @assert isortho(psi) && orthocenter(psi) == 1

    position!(PH, psi, 1)
    energy = 0.0

    sw_times = []
    for sw in 1:nsweep(sweeps)
    sw_time = @elapsed begin
        maxtruncerr = 0.0

        if !isnothing(write_when_maxdim_exceeds) &&
        maxdim(sweeps, sw) > write_when_maxdim_exceeds
        if outputlevel >= 2
            println(
            "write_when_maxdim_exceeds = $write_when_maxdim_exceeds and maxdim(sweeps, sw) = $(maxdim(sweeps, sw)), writing environment tensors to disk",
            )
        end
        PH = disk(PH)
        end

        for (b, ha) in sweepnext(N)
        @debug_check begin
            checkflux(psi)
            checkflux(PH)
        end

        @timeit_debug timer "dmrg: position!" begin
            if ha==1
                position!(PH, psi, b)
            else
                position!(PH, psi, b+1)
            end
        end

        @debug_check begin
            checkflux(psi)
            checkflux(PH)
        end

        @timeit_debug timer "dmrg: psi[b/b+1]" begin
            if ha == 1
                phi = psi[b]
            else
                phi = psi[b+1]
            end
        end

        @timeit_debug timer "dmrg: eigsolve" begin
            vals, vecs = eigsolve(
            PH,
            phi,
            1,
            eigsolve_which_eigenvalue;
            ishermitian=ishermitian,
            tol=eigsolve_tol,
            krylovdim=eigsolve_krylovdim,
            maxiter=eigsolve_maxiter,
            )
        end
        energy::Number = vals[1]
        phi::ITensor = vecs[1]

        ortho = ha == 1 ? "left" : "right"

        @debug_check begin
            checkflux(phi)
        end

        @timeit_debug timer "dmrg: factorize" begin
            if ha==1
                indsMb = inds(psi[b])
                if b==1
                    qrlinks = indsMb[1]
                else
                    qrlinks = indsMb[1:2]
                end
            else
                indsMb = inds(psi[b+1])
                if b==length(psi)-1
                    qrlinks = indsMb[2]
                else
                    qrlinks = indsMb[2:3]
                end
            end
            Q, R, spec = factorize(phi, qrlinks; which_decomp=which_decomp, positive=true, tags=tags(linkind(psi, b)))
            if ortho == "left"
                psi[b] = Q
                psi[b+1] = R*psi[b+1]
                leftlim(psi) == b - 1 && setleftlim!(psi, leftlim(psi) + 1)
                rightlim(psi) == b + 1 && setrightlim!(psi, rightlim(psi) + 1)
                (psi[b + 1] ./= norm(psi[b + 1]))
            elseif ortho == "right"
                psi[b] = psi[b]*R
                psi[b+1] = Q
                leftlim(psi) == b && setleftlim!(psi, leftlim(psi) - 1)
                rightlim(psi) == b + 2 && setrightlim!(psi, rightlim(psi) - 1)
                (psi[b] ./= norm(psi[b]))
            else
                error(
                "In replacebond!, got ortho = $ortho, only currently supports `left` and `right`."
                )
            end
        end
        maxtruncerr = max(maxtruncerr, spec.truncerr)

        @debug_check begin
            checkflux(psi)
            checkflux(PH)
        end

        if outputlevel >= 2
            @printf(
            "Sweep %d, half %d, bond (%d,%d) energy=%.12f\n", sw, ha, b, b + 1, energy
            )
            @printf(
            "  Truncated using cutoff=%.1E maxdim=%d mindim=%d\n",
            cutoff(sweeps, sw),
            maxdim(sweeps, sw),
            mindim(sweeps, sw)
            )
            flush(stdout)
        end
        sweep_is_done = (b == 1 && ha == 2)
        measure!(
          obs;
          energy=energy,
          psi=psi,
          bond=b,
          sweep=sw,
          half_sweep=ha,
          spec=spec,
          outputlevel=outputlevel,
          sweep_is_done=sweep_is_done,
        )
    end
    end
    if outputlevel >= 1
        @printf(
        "After sweep %d energy=%.12f maxlinkdim=%d maxerr=%.2E time=%.3f\n",
        sw,
        energy,
        maxlinkdim(psi),
        maxtruncerr,
        sw_time
        )
        flush(stdout)
    end
    append!(sw_times, sw_time)
    isdone = checkdone!(obs; energy=energy, psi=psi, sweep=sw, outputlevel=outputlevel)

    isdone && break
    end
    return (energy, psi, sw_times)
end