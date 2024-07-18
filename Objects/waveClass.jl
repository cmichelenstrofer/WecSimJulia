module WaveClassModule

# Define the WaveClass struct
mutable struct WaveClass
    bem::Dict{Symbol, Any}
    current::Dict{Symbol, Any}
    direction::Vector{Float64}
    elevationFile::String
    gamma::Union{Nothing, Float64}
    height::Union{Nothing, Float64}
    marker::Dict{Symbol, Any}
    period::Union{Nothing, Float64}
    phaseSeed::Int
    spectrumFile::String
    spectrumType::String
    viz::Dict{Symbol, Any}
    waterDepth::Union{Nothing, Float64}
    spread::Vector{Float64}

    # Internal properties
    amplitude::Union{Nothing, Float64}
    deepWater::Union{Nothing, Bool}
    dOmega::Float64
    omega::Union{Nothing, Float64}
    phase::Float64
    power::Union{Nothing, Float64}
    spectrum::Union{Nothing, Vector{Float64}}
    type::String
    typeNum::Union{Nothing, Int}
    waveAmpTime::Union{Nothing, Matrix{Float64}}
    waveAmpTimeViz::Union{Nothing, Matrix{Float64}}
    wavenumber::Union{Nothing, Float64}
end

# Constructor for WaveClass
function WaveClass(type::String)
    bem = Dict(:option => "EqualEnergy", :count => nothing, :frequency => nothing, :range => nothing)
    current = Dict(:option => 3, :depth => 0.0, :direction => 0.0, :speed => 0.0)
    direction = [0.0]
    elevationFile = "NOT DEFINED"
    gamma = nothing
    height = nothing
    marker = Dict(:location => nothing, :size => 10, :style => 1, :graphicColor => [0.26, 0.96, 0.89])
    period = nothing
    phaseSeed = 0
    spectrumFile = "NOT DEFINED"
    spectrumType = "NOT DEFINED"
    viz = Dict(:numPointsX => 50, :numPointsY => 50)
    waterDepth = nothing
    spread = [1.0]

    # Internal properties
    amplitude = nothing
    deepWater = nothing
    dOmega = 0.0
    omega = nothing
    phase = 0.0
    power = nothing
    spectrum = nothing
    typeNum = nothing
    waveAmpTime = nothing
    waveAmpTimeViz = nothing
    wavenumber = nothing

    obj = WaveClass(bem, current, direction, elevationFile, gamma, height, marker, period, phaseSeed, spectrumFile, spectrumType, viz, waterDepth, spread, amplitude, deepWater, dOmega, omega, phase, power, spectrum, type, typeNum, waveAmpTime, waveAmpTimeViz, wavenumber)
    
    set_type!(obj, type)
    return obj
end

# Define methods for WaveClass
function set_type!(obj::WaveClass, type::String)
    obj.type = type
    obj.typeNum = get_type_num(type)
end

function get_type_num(type::String)
    if type == "noWave"
        return 0
    elseif type == "noWaveCIC"
        return 1
    elseif type == "regular"
        return 10
    elseif type == "regularCIC"
        return 11
    elseif type == "irregular"
        return 20
    elseif type == "spectrumImport"
        return 21
    elseif type == "elevationImport"
        return 30
    else
        return nothing
    end
end

# Function to set up the regular wave properties
function setup(obj::WaveClass, rampTime::Float64, time::Vector{Float64}, g::Float64, rho::Float64)
    if obj.type == "regular" || obj.type == "regularCIC"
        if isnothing(obj.omega) && isnothing(obj.period)
            obj.omega = minimum(obj.bem[:range])
            obj.period = 2 * π / obj.omega
        elseif isnothing(obj.omega)
            obj.omega = 2 * π / obj.period
        else
            obj.period = 2 * π / obj.omega
        end
        obj.amplitude = obj.height / 2
        obj.wavenumber = calcWaveNumber(obj.omega, obj.waterDepth, g, obj.deepWater)
        waveElevReg!(obj, rampTime, time)
        wavePowerReg!(obj, g, rho)
    end
end

function calcWaveNumber(omega::Float64, waterDepth::Float64, g::Float64, deepWater::Bool)
    if deepWater
        return omega^2 / g
    else
        # Simplified shallow water approximation
        return omega^2 / (g * tanh(omega^2 * waterDepth / g))
    end
end

function waveElevReg!(obj::WaveClass, rampTime::Float64, time::Vector{Float64})
    maxIt = length(time)
    rampFunction = (1 .+ cos.(π .+ π .* time / rampTime)) / 2
    rampFunction[time .>= rampTime] .= 1

    obj.waveAmpTime = zeros(Float64, maxIt, 2)
    obj.waveAmpTime[:, 1] = time
    obj.waveAmpTime[:, 2] = rampFunction .* (obj.amplitude * cos.(obj.omega .* time))

    # Wave Marker
    if !isnothing(obj.marker[:location])
        SZwaveAmpTimeViz = size(obj.marker[:location])
        obj.waveAmpTimeViz = zeros(Float64, maxIt, SZwaveAmpTimeViz[1] + 1)
        for j in 1:SZwaveAmpTimeViz[1]
            obj.waveAmpTimeViz[:, 1] = time
            obj.waveAmpTimeViz[:, j + 1] = rampFunction .* obj.amplitude .* cos.(obj.omega .* time .- obj.wavenumber .* (obj.marker[:location][j, 1] .* cosd(obj.direction) + obj.marker[:location][j, 2] .* sind(obj.direction)))
        end
    end
end

function wavePowerReg!(obj::WaveClass, g::Float64, rho::Float64)
    if obj.deepWater
        # Deepwater Approximation
        obj.power = 1 / (8 * π) * rho * g^2 * (obj.amplitude)^2 * obj.period
    else
        # Full Wave Power Equation
        obj.power = rho * g * (obj.amplitude)^2 / 4 * sqrt(g / obj.wavenumber * tanh(obj.wavenumber * obj.waterDepth)) * (1 + 2 * obj.wavenumber * obj.waterDepth / sinh(2 * obj.wavenumber * obj.waterDepth))
    end
end

end # module WaveClassModule
