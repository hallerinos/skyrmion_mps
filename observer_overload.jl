# checkdone overload based on changes of local observables
function ITensors.checkdone!(o::DMRGObserver;kwargs...)
    lobs = o.measurements

    diff = 0.0
    if kwargs[:sweep] > 1 && lobs.count > 0
        for key in keys(lobs)
            # update maximum change of observables
            diff = maximum([diff, maximum(abs.(lobs[key][end-1] .- lobs[key][end]))])
        end
        if diff < o.etol
            println("Difference in local observables $diff < $(o.etol). Stopping DMRG.")
            return true
        end
    elseif (length(energies(o)) > o.minsweeps && abs(energies(o)[end] - energies(o)[end - 1]) < o.etol)
        println("Energy difference less than $(o.etol). Stopping DMRG.")
        return true
    end

    # exit sweeping gracefully
    try
        if isfile("stop")
            println("Stop file in root directory. Exit gracefully.")
            return true
        end
    catch err
        @error("Could not exit gracefully.", err)
        return
    end

    return false
end

# measure overload based on changes of local observables
function ITensors.measure!(o::DMRGObserver;kwargs...)
    df = DataFrame()
    for s in ["Sx", "Sy", "Sz"]
        df[!, s] = expect(kwargs[:psi],s)
    end
    df[!, "X"] = getindex.(kwargs[:nps],1)
    df[!, "Y"] = getindex.(kwargs[:nps],2)
    df[!, "Z"] = getindex.(kwargs[:nps],3)
    df[!, "sweep"] = [kwargs[:sweep] for i=1:size(df,1)]
    df[!, "half_sweep"] = [kwargs[:half_sweep] for i=1:size(df,1)]
    df[!, "bond"] = [kwargs[:bond] for i=1:size(df,1)]
    df[!, "E"] = [kwargs[:energy] for i=1:size(df,1)]
    append!(kwargs[:df], df)
    if kwargs[:bond]==length(kwargs[:psi])÷2 && kwargs[:half_sweep]==2
        push!(o.energies, kwargs[:energy])
    end
    return
end