using LinearAlgebra
using ProgressMeter

# Define the Hankel function
function hankel(c::Vector{T}, r::Vector{T}=reverse(c)) where T
    m, n = length(c), length(r)
    h = zeros(T, m, n)
    for i in 1:m
        for j in 1:n
            if i + j - 1 <= m
                h[i, j] = c[i + j - 1]
            else
                h[i, j] = r[i + j - m - 1]
            end
        end
    end
    return h
end

function radiationIRFSS(hydro::HydroData, Omax::Int=10, R2t::Float64=0.95)
    # Initialize state-space matrices if not already initialized
    num_dofs = sum(hydro.dof)
    hydro.ss_A = zeros(Float64, num_dofs, num_dofs, Omax, Omax)
    hydro.ss_B = zeros(Float64, num_dofs, num_dofs, Omax, 1)
    hydro.ss_C = zeros(Float64, num_dofs, num_dofs, 1, Omax)
    hydro.ss_D = zeros(Float64, num_dofs, num_dofs, 1)
    hydro.ss_K = zeros(Float64, num_dofs, num_dofs, length(hydro.ra_t))
    hydro.ss_conv = zeros(Float64, num_dofs, num_dofs)
    hydro.ss_R2 = zeros(Float64, num_dofs, num_dofs)
    hydro.ss_O = zeros(Float64, num_dofs, num_dofs)

    # Progress bar setup
    @showprogress "Calculating state space radiation IRFs..." for i in 1:num_dofs
        t = hydro.ra_t
        dt = t[2] - t[1]

        for i in 1:num_dofs
            for j in 1:num_dofs
                K = hydro.ra_K[i, j, :]
                if ndims(K) == 3
                    K = dropdims(K, dims=(3,))
                end
                R2i = norm(K .- mean(K))  # Initial R2

                # Hankel Singular Value Decomposition
                O = 2  # Initial state space order
                y = dt .* K
                h = hankel(y[2:end])
                u, s, v = svd(h)
                svh = diagm(s)

                while R2i != 0.0
                    u1 = u[1:length(K)-2, 1:O]
                    v1 = v[1:length(K)-2, 1:O]
                    u2 = u[2:length(K)-1, 1:O]
                    sqs = sqrt.(svh[1:O])

                    # Ensure the dimensions match for multiplication
                    if size(u1, 2) != size(u2, 2) || size(u1, 2) != length(sqs)
                        error("Dimension mismatch: u1, u2, and sqs must have compatible dimensions.")
                    end

                    ubar = u1' * u2
                    println("Dimensions of ubar: ", size(ubar))

                    a = ubar .* ((1 ./ sqs) * sqs')
                    b = v1[1, :] .* sqs
                    c = u1[1, :] .* sqs'
                    d = y[1]

                    iidd = inv(dt / 2 * (I + a))
                    ac = (a - I) * iidd
                    bc = dt * (iidd * b)
                    cc = c * iidd
                    dc = d - dt / 2 * ((c * iidd) * b)
                    ss_K = zeros(Float64, length(t))
                    for k in 1:length(t)
                        ss_K[k] = (cc * exp(ac * dt * (k - 1))) * bc
                    end

                    R2 = 1 - (norm(K - ss_K') / R2i)^2
                    if R2 >= R2t
                        status = 1
                        break
                    elseif O == Omax
                        status = 2
                        break
                    else
                        O += 1
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
            end
        end
    end
    return hydro
end
