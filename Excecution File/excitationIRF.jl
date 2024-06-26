using ProgressMeter
using Interpolations

function excitationIRF(hydro; tEnd=100.0, nDt=1001, nDw=1001, wMin=minimum(hydro["w"]), wMax=maximum(hydro["w"]))
    # Initialize progress bar
    p = Progress(sum(hydro["dof"]) * hydro["Nh"], desc="Calculating excitation IRFs...")

    # Prepare time and frequency arrays
    t = range(-tEnd, stop=tEnd, length=nDt)
    w = range(wMin, stop=wMax, length=nDw)
    N = sum(hydro["dof"]) * hydro["Nh"]

    # Calculate the impulse response function for excitation
    n = 0
    for i in 1:sum(hydro["dof"])
        for j in 1:hydro["Nh"]
            ex_re = interpolate(hydro["w"], hydro["ex_re"][i, j, :], Gridded(Linear()))
            ex_im = interpolate(hydro["w"], hydro["ex_im"][i, j, :], Gridded(Linear()))
            hydro["ex_K"][i, j, :] = (1 / π) * trapz(w, ex_re.(w) .* cos.(w .* t) .- ex_im.(w) .* sin.(w .* t))
            n += 1
            next!(p)
        end
    end

    hydro["ex_t"] = t
    hydro["ex_w"] = w

    return hydro
end

# Helper function for trapezoidal integration
function trapz(x, y)
    return sum((x[2:end] .- x[1:end-1]) .* (y[1:end-1] .+ y[2:end]) ./ 2)
end
