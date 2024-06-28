using Interpolations, Statistics, ProgressMeter

#include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/readWAMITv2.jl")


function radiationIRF(data::HydroData, tEnd = 100, nDt = 1001, nDw = 1001, wMin = minimum(data.w), wMax = maximum(data.w))
    p = ProgressMeter.Progress(length(data.dof)^2, desc = "Calculating radiation IRFs...")  # Progress bar

    t = range(Float32(0), stop = Float32(tEnd), length = nDt)
    w = range(Float32(wMin), stop = Float32(wMax), length = nDw)

    data.ra_K = zeros(Float16, 12, 12, length(t))  # Ensure this matches your data structure

    for i in 1:sum(data.dof)
        for j in 1:sum(data.dof)
            B_slice = Float32.(data.B[i, j, :])
            if ndims(B_slice) > 1
                B_slice = dropdims(B_slice, dims=3)
            end
            ra_B = LinearInterpolation(Float32.(data.w), B_slice, extrapolation_bc=Interpolations.Flat())

            integrand = [ra_B(ω) / Float32(data.rho) * cos(Float32(ω) * t_val) for t_val in t for ω in w]
            integrand = reshape(integrand, (length(t), length(w)))
            data.ra_K[i, j, :] .= (2 / π) .* sum(integrand, dims=2) .* (w[2] - w[1])  # Correcting the cumsum part
            ProgressMeter.next!(p)
        end
    end

    ra_Ainf_temp = zeros(Float32, length(data.w))
    if isempty(data.Ainf) || data.Ainf == nothing
        data.Ainf = zeros(12, 12)  # Ensure Ainf is initialized
    end
    for i in 1:sum(data.dof)
        for j in 1:sum(data.dof)
            A_slice = Float32.(data.A[i, j, :])
            if ndims(A_slice) > 1
                A_slice = dropdims(A_slice, dims=3)
            end
            ra_A = LinearInterpolation(Float32.(data.w), A_slice, extrapolation_bc=Interpolations.Flat())
            K_slice = data.ra_K[i, j, :]
            if ndims(K_slice) > 1
                K_slice = dropdims(K_slice, dims=3)
            end
            for k in 1:length(data.w)
                integrand = K_slice .* sin.(Float32(data.w[k]) .* t) ./ Float32(data.w[k])
                ra_Ainf_temp[k] = ra_A(Float32(data.w[k])) + sum(integrand) * (t[2] - t[1])  # Use broadcasting with dot syntax
            end
            data.Ainf[i, j] = mean(ra_Ainf_temp)
        end
    end

    data.ra_t = t
    data.ra_w = w

    return data
end
