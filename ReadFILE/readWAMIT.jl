using DelimitedFiles
using ProgressMeter
using Statistics

mutable struct HydroData
    code::String
    file::String
    h::Float64
    g::Float64
    rho::Float64
    body::Vector{String}
    Nb::Int
    Nf::Int
    Nh::Int
    Ainf::Array{Float64, 2}
    A::Array{Float16, 3}
    B::Array{Float16, 3}
    T::Vector{Float32}
    w::Vector{Float16}
    Vo::Vector{Float64}
    cg::Array{Float16, 2}
    cb::Array{Float16, 2}
    Khs::Array{Float16, 3}
    theta::Vector{Float16}
    dof::Vector{Int}
    ex_ma::Array{Float32, 3}
    ex_ph::Array{Float32, 3}
    ex_re::Array{Float32, 3}
    ex_im::Array{Float32, 3}
    sc_ma::Array{Float16, 3}
    sc_ph::Array{Float16, 3}
    sc_re::Array{Float16, 3}
    sc_im::Array{Float16, 3}
    fk_ma::Array{Float16, 3}
    fk_ph::Array{Float16, 3}
    fk_re::Array{Float16, 3}
    fk_im::Array{Float16, 3}
    md_mc::Array{Float64, 3}
    md_cs::Array{Float64, 3}
    md_pi::Array{Float64, 3}
    gbm::Array{Float64, 3}
    ra_K::Array{Float32, 3}
    ra_t::Vector{Float32}
    ra_w::Vector{Float32}
    ss_A::Array{Float64, 4}
    ss_B::Array{Float64, 3}
    ss_C::Array{Float64, 4}
    ss_D::Array{Float64, 2}
    ss_K::Array{Float64, 3}
    ss_conv::Array{Float64, 2}
    ss_R2::Array{Float64, 2}
    ss_O::Array{Float64, 2}
end

# Ensure the normalizeBEM and radiationIRFSS functions are included
include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/normalizeBEMv2.jl")

function readWAMIT(file_path::String)
    num_wave_periods = 260  # Assuming 260 wave periods
    data = HydroData(
        "WAMIT",
        splitext(basename(file_path))[1],
        Inf,           # Water depth
        9.81,          # Gravity
        1000.0,        # Water density
        String[],      # Body names
        0,             # Number of bodies
        0,             # Number of wave frequencies
        0,             # Number of wave headings
        zeros(Float64, 12, 12),  # Added Mass Coefficients at WaveNumber = infinity
        zeros(Float16, 12, 12, num_wave_periods),  # Added Mass Coefficients
        zeros(Float16, 12, 12, num_wave_periods),  # Radiation Damping Coefficients
        zeros(Float32, num_wave_periods),  # Wave period
        zeros(Float16, num_wave_periods),  # Wave frequencies
        Float64[],    # Displacement volume
        zeros(Float16, 3, 2),  # Center of Gravity
        zeros(Float16, 3, 2),  # Center of Buoyancy
        zeros(Float16, 6, 6, 2),  # Hydrostatic and gravitational restoring coefficients
        Float16[],    # Wave headings
        Int[],        # Degrees of freedom for each body
        NaN * ones(Float32, 12, 1, 260),  # Exciting force magnitudes
        NaN * ones(Float32, 12, 1, 260),  # Exciting force phases
        NaN * ones(Float32, 12, 1, 260),  # Exciting force real parts
        NaN * ones(Float32, 12, 1, 260),  # Exciting force imaginary parts
        NaN * ones(Float16, 12, 1, 260),  # Scattering force magnitudes
        NaN * ones(Float16, 12, 1, 260),  # Scattering force phases
        NaN * ones(Float16, 12, 1, 260),  # Scattering force real parts
        NaN * ones(Float16, 12, 1, 260),  # Scattering force imaginary parts
        NaN * ones(Float16, 12, 1, 260),  # Froude-Krylov force magnitudes
        NaN * ones(Float16, 12, 1, 260),  # Froude-Krylov force phases
        NaN * ones(Float16, 12, 1, 260),  # Froude-Krylov force real parts
        NaN * ones(Float16, 12, 1, 260),  # Froude-Krylov force imaginary parts
        zeros(Float64, 12, 12, 260),       # Momentum Conservation Drift Forces
        zeros(Float64, 12, 12, 260),       # Control Surface Drift Forces
        zeros(Float64, 12, 12, 260),       # Pressure Integration Drift Forces
        zeros(Float64, 0, 0, 0),  # Placeholder for gbm
        zeros(Float32, 12, 12, 1001),  # Initialize ra_K with correct dimensions
        Float32[],     # Placeholder for ra_t
        Float32[],     # Placeholder for ra_w
        zeros(Float64, 12, 12, 10, 10),   # ss_A
        zeros(Float64, 12, 12, 10),       # ss_B
        zeros(Float64, 12, 12, 1, 10),    # ss_C
        zeros(Float64, 12, 12),           # ss_D
        zeros(Float64, 12, 12, 1001),     # ss_K
        zeros(Float64, 12, 12),           # ss_conv
        zeros(Float64, 12, 12),           # ss_R2
        zeros(Float64, 12, 12)            # ss_O
    )

    raw = readdlm(file_path, '\n', String)[:, 1]
    N = length(raw)
    progress = ProgressMeter.Progress(N, 1, "Reading WAMIT output file...")  # Progress bar

    reading_added_mass_inf = false
    reading_added_mass_damping = false
    reading_exciting_forces = false
    current_wave_period = 0

    k = 0
    for i in 1:N
        ProgressMeter.update!(progress, i)
        line = raw[i]

        if occursin("Body number:", line)
            k += 1
            tmp = split(line, [' ', '.'])
            tmp = filter(x -> !isempty(x), tmp)
            push!(data.body, "Body $k")
            data.Nb += 1
            push!(data.dof, 6)
            push!(data.Vo, 0.0)
        end

        if occursin("XBODY =", line)
            values = split(line)
            xg = parse(Float64, values[3])
            yg = parse(Float64, values[6])
            zg = parse(Float64, values[9])
            data.cg[:, k] = [xg, yg, zg]
        end

        if occursin("Center of Buoyancy (Xb,Yb,Zb):", line)
            values = split(line[findfirst(':', line)+1:end])
            xb = parse(Float64, values[1])
            yb = parse(Float64, values[2])
            zb = parse(Float64, values[3])
            data.cb[:, k] = [xb, yb, zb]
        end

        if occursin("Hydrostatic and gravitational", line)
            data.Khs[:, :, k] .= 0.0
            tmp1 = parse.(Float64, split(raw[i+1][findfirst(':', raw[i+1])+1:end]))
            data.Khs[3, 3:5, k] = tmp1
            data.Khs[3:5, 3, k] = tmp1
            tmp2 = parse.(Float64, split(raw[i+2][findfirst(':', raw[i+2])+1:end]))
            data.Khs[4, 4:6, k] = tmp2
            data.Khs[4:5, 4, k] = tmp2[1:2]
            tmp3 = parse.(Float64, split(raw[i+3][findfirst(':', raw[i+3])+1:end]))
            data.Khs[5, 5:6, k] = tmp3
        end

        if occursin("Water depth:", line)
            if occursin("infinite", line)
                data.h = Inf
            else
                tmp = split(line, r"\s+")
                depth_str = tmp[end]
                depth_val = tryparse(Float64, depth_str)
                if depth_val !== nothing
                    data.h = depth_val
                else
                    println("Failed to parse water depth value: ", line)
                end
            end
        end

        if occursin("Volumes (VOLX,VOLY,VOLZ):", line)
            numbers_part = split(line, ":")[2]
            numbers_str = split(strip(numbers_part))
            numbers_str = filter(s -> !isempty(s), numbers_str)
            println("Parsing Volumes: ", numbers_str)
            numbers = parse.(Float64, numbers_str)
            data.Vo[k] = mean(numbers)
        end

        if occursin("Wave period (sec) =", line)
            split_line = split(line, "=")
            if length(split_line) >= 2
                wave_period = parse(Float64, split(split_line[2])[1])
                current_wave_period += 1
                data.T[current_wave_period] = wave_period
                data.w[current_wave_period] = Float16(2 * π / wave_period)  # Ensure wave frequency is Float16
                data.Nf += 1
                println("Parsed wave period for period $current_wave_period: $wave_period")
                println("Computed wave frequency for period $current_wave_period: $(data.w[current_wave_period])")
            end
        end

        if occursin("Wave Heading", line)
            tmp = parse(Float64, split(line[findfirst(':', line)+1:end])[1])
            push!(data.theta, tmp)
            data.Nh += 1
        end

        if occursin("ADDED-MASS COEFFICIENTS", line) && occursin("Wavenumber = infinite", raw[i-2])
            println("Found ADDED-MASS COEFFICIENTS with Wavenumber = infinite at line $i")
            reading_added_mass_inf = true
            continue
        end

        if reading_added_mass_inf
            if occursin("************************************************************************", line)
                reading_added_mass_inf = false
                println("Ending ADDED-MASS COEFFICIENTS with Wavenumber = infinite at line $i")
                continue
            end

            split_line = split(line)
            if length(split_line) == 3 && !isnothing(tryparse(Int, split_line[1])) && !isnothing(tryparse(Int, split_line[2])) && !isnothing(tryparse(Float64, split_line[3]))
                I = parse(Int, split_line[1])
                J = parse(Int, split_line[2])
                A = parse(Float64, split_line[3])
                data.Ainf[I, J] = A
                println("Ainf[$I, $J] = $A")
            end
        end

        if occursin("ADDED-MASS AND DAMPING COEFFICIENTS", line)
            reading_added_mass_damping = true
            println("Found ADDED-MASS AND DAMPING COEFFICIENTS at line $i for wave period $current_wave_period")
            continue
        end

        if reading_added_mass_damping
            if occursin("************************************************************************", line)
                reading_added_mass_damping = false
                println("Ending ADDED-MASS AND DAMPING COEFFICIENTS at line $i for wave period $current_wave_period")
                continue
            end

            split_line = split(line)
            if length(split_line) == 4 && !isnothing(tryparse(Int, split_line[1])) && !isnothing(tryparse(Int, split_line[2])) && !isnothing(tryparse(Float64, split_line[3])) && !isnothing(tryparse(Float64, split_line[4]))
                I = parse(Int, split_line[1])
                J = parse(Int, split_line[2])
                A = parse(Float64, split_line[3])
                B = parse(Float64, split_line[4])
                data.A[I, J, current_wave_period] = A
                data.B[I, J, current_wave_period] = B
                println("A[$I, $J, $current_wave_period] = $A, B[$I, $J, $current_wave_period] = $B")
            end
        end

        if occursin("HASKIND EXCITING FORCES AND MOMENTS", line) || occursin("DIFFRACTION EXCITING FORCES AND MOMENTS", line) || occursin("RESPONSE AMPLITUDE OPERATORS", line)
            reading_exciting_forces = true
            println("Found EXCITING FORCES at line $i for wave period $current_wave_period")
            continue
        end

        if reading_exciting_forces
            if occursin("************************************************************************", line)
                reading_exciting_forces = false
                println("Ending EXCITING FORCES at line $i for wave period $current_wave_period")
                continue
            end

            split_line = split(line)
            if length(split_line) == 3 && !isnothing(tryparse(Int, split_line[1])) && !isnothing(tryparse(Float64, split_line[2])) && !isnothing(tryparse(Int, split_line[3]))
                I = abs(parse(Int, split_line[1]))
                wave_period_idx = current_wave_period
                if I <= 12
                    magnitude = parse(Float64, split_line[2])
                    phase_deg = parse(Int, split_line[3])
                    phase_rad = deg2rad(phase_deg)
                    real_part = magnitude * cos(phase_rad)
                    imaginary_part = magnitude * sin(phase_rad)
                    data.ex_ma[I, 1, wave_period_idx] = magnitude
                    data.ex_ph[I, 1, wave_period_idx] = phase_rad
                    data.ex_re[I, 1, wave_period_idx] = real_part
                    data.ex_im[I, 1, wave_period_idx] = imaginary_part
                    println("Exciting Force: magnitude=$magnitude, phase=$phase_rad, real_part=$real_part, imaginary_part=$imaginary_part for mode $I, wave period $wave_period_idx")
                else
                    println("Warning: Out-of-bounds index for exciting force at line $i")
                end
            end
        end

        # Handling Drift Forces (Example for momentum conservation method)
        if occursin("SURGE, SWAY & YAW DRIFT FORCES (Momentum Conservation)", line)
            data.Nh = 0  # Number of wave headings
            i = n + 1
            data.md_mc[:, :, data.Nf] .= 0.0
            while occursin("Wave Heading", raw[i])
                data.Nh += 1
                tmp = parse(Float64, split(raw[i][findfirst(':', raw[i])+1:end])[1])
                push!(data.theta, tmp)
                i += 2
                while !occursin("*******************************", raw[i]) && 
                      !occursin("EXCITING FORCES AND MOMENTS", raw[i]) &&
                      !occursin("RESPONSE AMPLITUDE OPERATORS", raw[i]) &&
                      !occursin("SURGE, SWAY & YAW DRIFT FORCES (Momentum Conservation)", raw[i]) &&
                      !occursin("SURGE, SWAY, HEAVE, ROLL, PITCH & YAW DRIFT FORCES (Control Surface)", raw[i]) &&
                      !occursin("SURGE, SWAY, HEAVE, ROLL, PITCH & YAW DRIFT FORCES (Pressure Integration)", raw[i]) &&
                      !occursin("Wave Heading", raw[i])
                    tmp = parse.(Float64, split(raw[i]))
                    I = abs(tmp[1])
                    if I <= 12
                        data.md_mc[I, data.Nh, data.Nf] = tmp[2] * cos(deg2rad(tmp[3]))
                    end
                    i += 1
                end
            end
        end
    end

    sc_file = replace(file_path, ".out" => ".3sc")
    if isfile(sc_file)
        raw_sc = readdlm(sc_file, '\n', String)[:, 1]
        n = 1
        for i in 1:data.Nf
            for j in 1:data.Nh
                for k in 1:round(Int, size(findall(!isnan, data.ex_ma), 1) / data.Nf / data.Nh)
                    n += 1
                    values = parse.(Float64, split(raw_sc[n]))
                    I = values[3]
                    if I <= 12
                        data.sc_ma[I, 1, i] = values[4]
                        data.sc_ph[I, 1, i] = deg2rad(values[5])
                        data.sc_re[I, 1, i] = values[6]
                        data.sc_im[I, 1, i] = values[7]
                    else
                        println("Warning: Out-of-bounds index for scattering force at line $n")
                    end
                end
            end
        end
    end

    fk_file = replace(file_path, ".out" => ".3fk")
    if isfile(fk_file)
        raw_fk = readdlm(fk_file, '\n', String)[:, 1]
        n = 1
        for i in 1:data.Nf
            for j in 1:data.Nh
                for k in 1:round(Int, size(findall(!isnan, data.ex_ma), 1) / data.Nf / data.Nh)
                    n += 1
                    values = parse.(Float64, split(raw_fk[n]))
                    I = values[3]
                    if I <= 12
                        data.fk_ma[I, 1, i] = values[4]
                        data.fk_ph[I, 1, i] = deg2rad(values[5])
                        data.fk_re[I, 1, i] = values[6]
                        data.fk_im[I, 1, i] = values[7]
                    else
                        println("Warning: Out-of-bounds index for Froude-Krylov force at line $n")
                    end
                end
            end
        end
    end

    cfg_file = replace(file_path, ".out" => ".cfg")
    if isfile(cfg_file)
        raw_cfg = readdlm(cfg_file, '\n', String)[:, 1]
        for line in raw_cfg
            if occursin("NEWMDS", line) || occursin("IMODESFSP", line)
                tmp = split(line, ['(', ')', '=', ' '])
                tmp = filter(x -> !isempty(x), tmp)
                if line[7] == '('
                    body_index = parse(Int, tmp[2])
                    additional_dof = parse(Int, tmp[3])
                    data.dof[body_index] += additional_dof
                else
                    data.dof[1] += parse(Int, tmp[2])
                end
            end
        end
        if sum(data.dof) > data.Nb * 6
            mmx_file = replace(file_path, ".out" => ".mmx")
            if isfile(mmx_file)
                raw_mmx = readdlm(mmx_file, '\n', String)[:, 1]
                for line in raw_mmx
                    if occursin("External force matrices", line)
                        for k in 1:sum(data.dof) * sum(data.dof)
                            tmp = parse.(Float64, split(raw_mmx[line+k+1]))
                            data.gbm[tmp[1], tmp[2], 1] = tmp[3]  # Mass
                            data.gbm[tmp[1], tmp[2], 2] = tmp[4]  # Damping
                            data.gbm[tmp[1], tmp[2], 3] = tmp[5]  # Stiffness
                        end
                    end
                end
            end
            hst_file = replace(file_path, ".out" => ".hst")
            if isfile(hst_file)
                raw_hst = readdlm(hst_file, '\n', String)[:, 1]
                for line in raw_hst[2:end]
                    tmp = parse.(Float64, split(line))
                    data.gbm[tmp[1], tmp[2], 4] = tmp[3]  # Hydrostatic Stiffness
                end
            end
        end
    end

    # Call the normalizeBEM function to normalize the data
    data = normalizeBEM(data)
    return data
end

function print_hydro_data(data::HydroData)
    println("HydroData Structure:")
    println("code: ", data.code)
    println("file: ", data.file)
    println("h: ", data.h)
    println("g: ", data.g)
    println("rho: ", data.rho)
    println("Nb: ", data.Nb)
    println("Nf: ", data.Nf)
    println("Nh: ", data.Nh)
    println("body: ", data.body)
    
    println("Ainf (first two columns): ")
    println(data.Ainf)
    
    println("A (first two columns of each layer): ")
    println(data.A[:, :, 1:3])
    
    println("B (first two columns of each layer): ")
    println(data.B[:, :, 1:3])
    
    println("T: ", data.T)
    println("w: ", data.w)
    println("Vo: ", data.Vo)
    println("cg: ", data.cg)
    println("cb: ", data.cb)
    println("Khs: ", data.Khs)
    println("theta: ", data.theta)
    println("dof: ", data.dof)
    
    println("ex_ma (first three columns of each layer): ")
    println(data.ex_ma[:, :, 1:3])
    
    println("ex_ph (first three columns of each layer): ")
    println(data.ex_ph[:, :, 1:3])
    
    println("ex_re (first three columns of each layer): ")
    println(data.ex_re[:, :, 1:3])
    
    println("ex_im (first three columns of each layer): ")
    println(data.ex_im[:, :, 1:3])
    
    println("sc_ma (first three columns of each layer): ")
    println(data.sc_ma[:, :, 1:3])
    
    println("sc_ph (first three columns of each layer): ")
    println(data.sc_ph[:, :, 1:3])
    
    println("sc_re (first three columns of each layer): ")
    println(data.sc_re[:, :, 1:3])
    
    println("sc_im (first three columns of each layer): ")
    println(data.sc_im[:, :, 1:3])
    
    println("fk_ma (first three columns of each layer): ")
    println(data.fk_ma[:, :, 1:3])
    
    println("fk_ph (first three columns of each layer): ")
    println(data.fk_ph[:, :, 1:3])
    
    println("fk_re (first three columns of each layer): ")
    println(data.fk_re[:, :, 1:3])
    
    println("fk_im (first three columns of each layer): ")
    println(data.fk_im[:, :, 1:3])
    
    #println("md_mc: ", data.md_mc)
    #println("md_cs: ", data.md_cs)
    #println("md_pi: ", data.md_pi)
    #println("gbm: ", data.gbm)
    
    println("ra_K (first three columns of each layer): ")
    println(data.ra_K[:, :, 1:3])
    
    println("ra_t: ", data.ra_t)
    println("ra_w: ", data.ra_w)
end
