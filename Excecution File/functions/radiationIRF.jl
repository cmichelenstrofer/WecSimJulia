using Interpolations, Statistics, ProgressMeter

function radiationIRF(data::HydroData, tEnd = 100, nDt = 1001, nDw = 1001, wMin = minimum(data.w), wMax = maximum(data.w))
    p = ProgressMeter.Progress(length(data.dof)^2, desc = "Calculating radiation IRFs...")

    t = range((0), stop = (tEnd), length = nDt)
    w = range((wMin), stop = (wMax), length = nDw)

    data.ra_K = zeros(Float64, sum(data.dof), sum(data.dof), length(t))

    for i in 1:sum(data.dof)
        for j in 1:sum(data.dof)
            B_slice = Float64.(data.B[i, j, :])
            ra_B = LinearInterpolation(data.w, B_slice, extrapolation_bc=Flat())

            for ti in 1:length(t)
                integrand = [ra_B(ω) * cos(ω * t[ti]) for ω in w]
                data.ra_K[i, j, ti] = (2 / π) * trapz(w, integrand)
            end
            next!(p)
        end
    end

    ra_Ainf_temp = zeros(Float32, length(data.w))
    if isempty(data.Ainf) || data.Ainf == nothing
        data.Ainf = zeros(12, 12)
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
                    ra_Ainf_temp[k] = ra_A(Float32(data.w[k])) .+ sum(integrand) .* (t[2] - t[1])
                end
                data.Ainf[i, j] = mean(ra_Ainf_temp)
            end
        end
    end

    data.ra_t = t
    data.ra_w = w

    return data
end
