using Interpolations, Statistics

function radiationIRF(data, tEnd = 100, nDt = 1001, nDw = 1001, wMin = minimum(data[:w]), wMax = maximum(data[:w]))
    p = Progress(length(data[:dof])^2, desc = "Calculating radiation IRFs...")  # Progress bar

    t = range(Float32(0), stop = Float16(tEnd), length = nDt)
    w = range(Float32(wMin), stop = Float16(wMax), length = nDw)

    data[:ra_K] = zeros(Float32, 12, 12, length(t))  # Ensure this matches your data structure

    for i in 1:sum(data[:dof])
        for j in 1:sum(data[:dof])
            B_slice = Float16.(data[:B][i, j, :])
            if ndims(B_slice) > 1
                B_slice = dropdims(B_slice, dims=3)
            end
            ra_B = LinearInterpolation(Float16.(data[:w]), B_slice, extrapolation_bc=Interpolations.Flat())

            integrand = [ra_B(ω) / Float16(data[:rho]) * cos(Float16(ω) * t_val) for t_val in t for ω in w]
            integrand = reshape(integrand, (length(t), length(w)))
            data[:ra_K][i, j, :] .= (2 / π) * cumsum(integrand, dims=1)[:, end] * (t[2] - t[1])
            next!(p)
        end
    end

    ra_Ainf_temp = zeros(Float16, length(data[:w]))
    if isempty(data[:Ainf]) || !haskey(data, :Ainf) || data[end][:code] != "WAMIT"
        for i in 1:sum(data[:dof])
            for j in 1:sum(data[:dof])
                A_slice = Float16.(data[:A][i, j, :])
                if ndims(A_slice) > 1
                    A_slice = dropdims(A_slice, dims=3)
                end
                ra_A = LinearInterpolation(Float16.(data[:w]), A_slice, extrapolation_bc=Interpolations.Flat())
                K_slice = data[:ra_K][i, j, :]
                if ndims(K_slice) > 1
                    K_slice = dropdims(K_slice, dims=3)
                end
                for k in 1:length(data[:w])
                    integrand(t) = K_slice .* sin(Float16(data[:w][k]) * t) ./ Float16(data[:w][k])
                    ra_Ainf_temp[k] = ra_A(Float16(data[:w][k])) + cumsum(integrand.(t)) * (t[2] - t[1])
                end
                data[:Ainf][i, j] = mean(ra_Ainf_temp)
            end
        end
    end

    data[:ra_t] = t
    data[:ra_w] = w

    return data
end
