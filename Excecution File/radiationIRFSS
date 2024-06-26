using LinearAlgebra
using ProgressMeter
using Interpolations
using LinearAlgebra: norm

function radiationIRFSS(hydro, Omax::Int = 10, R2t::Float64 = 0.95, tEnd::Float64 = 100.0, nDt::Int = 1001, nDw::Int = 1001, wMin::Float64 = Float64(minimum(hydro.w)), wMax::Float64 = Float64(maximum(hydro.w)))
    # Calculates the radiation IRFs in SS format

    N = sum(hydro.dof)^2
    p = ProgressMeter.Progress(N, 1, "Calculating radiation IRFs...")  # Progress bar

    # Interpolate to the given t and w
    t = range(-tEnd, stop=tEnd, length=nDt)
    w = range(wMin, stop=wMax, length=nDw)

    # Initialize SS matrices and related data structures
    dof_sum = sum(hydro.dof)
    hydro.ss_A = zeros(Float64, dof_sum, dof_sum, Omax, Omax)
    hydro.ss_B = zeros(Float64, dof_sum, dof_sum, Omax, 1)
    hydro.ss_C = zeros(Float64, dof_sum, dof_sum, 1, Omax)
    hydro.ss_D = zeros(Float64, dof_sum, dof_sum)
    hydro.ss_K = zeros(Float64, dof_sum, dof_sum, nDt)
    hydro.ss_conv = zeros(Int, dof_sum, dof_sum)
    hydro.ss_R2 = zeros(Float64, dof_sum, dof_sum)
    hydro.ss_O = zeros(Int, dof_sum, dof_sum)

    n = 0
    for i in 1:dof_sum
        for j in 1:dof_sum
            ra_B = hydro.B[i, j, :]
            if ndims(ra_B) > 1
                ra_B = reshape(ra_B, length(ra_B))  # Ensure ra_B is a 1D array
            end
            ra_B_interp = interpolate((hydro.w,), ra_B, Gridded(Linear()))
            ra_B_vals = ra_B_interp.(w)
            
            K = zeros(Float64, nDt)
            for k in 1:nDt
                K[k] = (2 / π) * trapz(w, ra_B_vals .* cos.(w .* t[k]))
                hydro.ra_K[i, j, k] = K[k]
            end

            R2i = norm(K - mean(K))  # Initial R2

            # Hankel Singular Value Decomposition
            O = 2  # Initial state space order
            y = t * K
            h = hankel(y[2:end])
            u, svh, v = svd(h)
            svh = diagm(svh)

            while R2i != 0.0
                u1 = u[1:length(K) - 2, 1:O]
                v1 = v[1:length(K) - 2, 1:O]
                u2 = u[2:length(K) - 1, 1:O]
                sqs = sqrt(svh[1:O])
                ubar = u1' * u2

                a = ubar .* ((1 ./ sqs) * sqs')
                b = v1[1, :]' .* sqs
                c = u1[1, :] .* sqs'
                d = y[1]

                iidd = inv(t[2] / 2 * (I + a))  # (T/2*I+T/2*A)^{-1} = 2/T(I+A)^{-1)
                ac = (a - I) * iidd  # (A-I)2/T(I+A)^{-1) = 2/T(A-I)(I+A)^{-1)
                bc = t[2] * (iidd * b)  # (T/2+T/2)*2/T(I+A)^{-1}B = 2(I+A)^{-1}B
                cc = c * iidd  # C*2/T(I+A)^{-1) = 2/T(I+A)^{-1)
                dc = d - t[2] / 2 * ((c * iidd) * b)  # D-T/2C (2/T(I+A)^{-1})B = D-C(I+A)^{-1})B

                ss_K = zeros(Float64, length(t))
                for k in 1:length(t)
                    ss_K[k] = (cc * exp(ac * t[k - 1]) * bc)  # Calc SS IRF approx
                end

                R2 = 1 - (norm(K - ss_K') / R2i)^2  # Calc R2 for SS IRF approx
                if R2 >= R2t
                    status = 1
                    break  # R2 threshold
                elseif O == Omax
                    status = 2
                    break  # Max SS order threshold
                else
                    O += 1  # Increase state space order
                end
            end

            if R2i != 0.0
                hydro.ss_A[i, j, 1:O, 1:O] = ac
                hydro.ss_B[i, j, 1:O, 1] = bc
                hydro.ss_C[i, j, 1, 1:O] = cc
                hydro.ss_D[i, j] = dc
                hydro.ss_K[i, j, :] = ss_K
                hydro.ss_conv[i, j] = status
                hydro.ss_R2[i, j] = R2
                hydro.ss_O[i, j] = O
            end

            n += 1
            ProgressMeter.update!(p, n)
        end
    end

    hydro.ra_t = t
    hydro.ra_w = w

    return hydro
end
