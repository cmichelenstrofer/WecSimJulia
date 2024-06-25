using DelimitedFiles
using ProgressMeter
using Statistics

include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/normalizeBEM.jl")

mutable struct WAMITData
    code::String
    file::String
    h::Float64
    g::Float64
    rho::Float64
    body::Vector{String}
    Nb::Int
    Nf::Int
    Nh::Int
    Ainf::Matrix{Float64}
    A::Array{Float16, 3}
    B::Array{Float16, 3}
    T::Vector{Float32}
    w::Vector{Float16}
    Vo::Vector{Float64}
    cg::Array{Float16, 3}
    cb::Array{Float16, 3}
    Khs::Array{Float16, 3}
    theta::Vector{Float16}
    dof::Vector{Int}
    ex_ma::Array{Float16, 3}
    ex_ph::Array{Float16, 3}
    ex_re::Array{Float16, 3}
    ex_im::Array{Float16, 3}
    sc_ma::Array{Float16, 3}
    sc_ph::Array{Float16, 3}
    sc_re::Array{Float16, 3}
    sc_im::Array{Float16, 3}
    fk_ma::Array{Float16, 3}
    fk_ph::Array{Float16, 3}
    fk_re::Array{Float16, 3}
    fk_im::Array{Float16, 3}
    ra_K::Array{Float64, 3}
    ra_t::Array{Float64, 2}
    ra_w::Array{Float64, 2}
end

function readWAMIT(; filepath::String)
    data = WAMITData(
        "WAMIT",
        splitext(basename(filepath))[1],
        Inf,
        9.81,
        1000.0,
        String[],
        0,
        0,
        0,
        zeros(Float64, 12, 12),
        zeros(Float16, 12, 12, 260),
        zeros(Float16, 12, 12, 260),
        zeros(Float32, 260),
        zeros(Float16, 260),
        Float64[],
        zeros(Float16, 3, 2),
        zeros(Float16, 3, 2),
        Array{Float16}(undef, 6, 6, 2),
        Float16[],
        Int[],
        NaN * ones(Float16, 12, 1, 260),
        NaN * ones(Float16, 12, 1, 260),
        NaN * ones(Float16, 12, 1, 260),
        NaN * ones(Float16, 12, 1, 260),
        NaN * ones(Float16, 12, 1, 260),
        NaN * ones(Float16, 12, 1, 260),
        NaN * ones(Float16, 12, 1, 260),
        NaN * ones(Float16, 12, 1, 260),
        NaN * ones(Float16, 12, 1, 260),
        NaN * ones(Float16, 12, 1, 260),
        NaN * ones(Float16, 12, 1, 260),
        zeros(Float64, 12, 12, 1001),
        zeros(Float64, 1, 1001),
        zeros(Float64, 1, 1001)
    )

    raw = readdlm(filepath, '\n', String)[:, 1]
    N = length(raw)
    progress = Progress(N, 1, "Reading WAMIT output file...")  # Progress bar
    
    ## Three flags are used to parse data on whether the specific section is being read ##
    reading_added_mass_inf = false
    reading_added_mass_damping = false
    reading_exciting_forces = false

    ## Counter for wave periods and to correctly allocate the data for each particular wave period ##
    current_wave_period = 0

    k = 0
    for i in 1:N
        update!(progress, i)
        line = raw[i]

        ## For identifying the Body Names and the Number of Bodies
        if occursin("Body number:", line)
            k += 1
            tmp = split(line, [' ', '.'])
            tmp = filter(x -> !isempty(x), tmp)
            push!(data.body, "Body $k")               # Body names (assuming body names are sequential)
            data.Nb += 1                              # Number of bodies
            push!(data.dof, 6)                        # Default degrees of freedom for each body is 6
            push!(data.Vo, 0.0)                       # Initialize displacement volume for the body
        end

        ## For identifying the Center of Gravity (Xg, Yg, Zg)
        if occursin("XBODY =", line)
            values = split(line)
            xg = parse(Float64, values[3])
            yg = parse(Float64, values[6])
            zg = parse(Float64, values[9])
            data.cg[:, k] = [xg, yg, zg]
        end

        ## For identifying the Center of Buoyancy (Xb, Yb, Zb)
        if occursin("Center of Buoyancy (Xb,Yb,Zb):", line)
            values = split(line[findfirst(':', line)+1:end])
            xb = parse(Float64, values[1])
            yb = parse(Float64, values[2])
            zb = parse(Float64, values[3])
            data.cb[:, k] = [xb, yb, zb]
        end

        ## For identifying the Hydrostatic and gravitational restoring coefficients
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

        ## For identifying the Water Depth
        if occursin("Water depth:", line)
            if occursin("infinite", line)
                data.h = Inf  # Directly assign Inf for "infinite" water depth
            else
                # Extract and parse the numerical value for water depth
                tmp = split(line, r"\s+")
                depth_str = tmp[end]  # Assuming the depth value is the last element
                depth_val = tryparse(Float64, depth_str)
                if depth_val !== nothing
                    data.h = depth_val
                else
                    println("Failed to parse water depth value: ", line)
                end
            end
        end

        ## For identifying the Water Volumes (x, y and z)
        if occursin("Volumes (VOLX,VOLY,VOLZ):", line)
            # Split the line at the colon to get the numbers part
            numbers_part = split(line, ":")[2]
            
            # Split the numbers part by whitespace and filter out empty strings
            numbers_str = split(strip(numbers_part))
            numbers_str = filter(s -> !isempty(s), numbers_str)  # Filter out empty strings
            
            println("Parsing Volumes: ", numbers_str)  # Diagnostic print
            
            # Safely parse the numbers
            numbers = parse.(Float64, numbers_str)
            
            # Storing Average values of the volumes of the current body
            data.Vo[k] = mean(numbers) 
        end

        ## For identifying the Wave Period
        if occursin("Wave period (sec) =", line)
            split_line = split(line, "=")
            if length(split_line) >= 2  # Ensure the split line has enough elements
                wave_period = parse(Float64, split(split_line[2])[1])
                current_wave_period += 1
                data.T[current_wave_period] = wave_period
                data.w[current_wave_period] = 2 * π / wave_period  # Compute wave frequency and store it in the [:w] matrix
                data.Nf += 1
                println("Parsed wave period for period $current_wave_period: $wave_period")
                println("Computed wave frequency for period $current_wave_period: $(data.w[current_wave_period])")
            end
        end

        ## For identifying the Wave Headings
        if occursin("Wave Heading", line)
            tmp = parse(Float64, split(line[findfirst(':', line)+1:end])[1])
            push!(data.theta, tmp)
            data.Nh += 1
        end
        
        ## For identifying the Added Mass Coefficients at Infinite Wave Period
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
        
        ## For identifying the Added Mass Coefficients and Damping Coefficients at a Given Wave Period
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

        ## For identifying the Exciting Forces
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
    end

    # Scattering Force
    sc_file = replace(filepath, ".out" => ".3sc")
    if isfile(sc_file)
        raw_sc = readdlm(sc_file, '\n', String)[:, 1]
        n = 1
        for i in 1:data.Nf
            for j in 1:data.Nh
                for k in 1:round(Int, size(findall(!isnan, data.ex_ma), 1) / data.Nf / data.Nh)  # Number of non-zero dof
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

    # Froude-Krylov force
    fk_file = replace(filepath, ".out" => ".3fk")
    if isfile(fk_file)
        raw_fk = readdlm(fk_file, '\n', String)[:, 1]
        n = 1
        for i in 1:data.Nf
            for j in 1:data.Nh
                for k in 1:round(Int, size(findall(!isnan, data.ex_ma), 1) / data.Nf / data.Nh)  # Number of non-zero dof
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

    # Adjust DOF if generalized body modes are used
    cfg_file = replace(filepath, ".out" => ".cfg")
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
    end
    normalizeBEM(data)
    println("Normalizing data....")
    return data
end
