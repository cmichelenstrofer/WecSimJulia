using Plots
using LinearAlgebra
using ControlSystems
using QuadGK

include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/readWAMITv2.jl")
include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/objects/waveClass.jl")
using .WaveClassModule: WaveClass

hydro = readWAMIT("C:/Users/jelope/Desktop/Git/WEC-Sim/examples/RM3/hydroData/rm3.out")
# Function for cumulative trapezoidal integration
function cumtrapz(x, y)
    n = length(x)
    z = zeros(n)
    for i in 2:n
        z[i] = z[i-1] + (x[i] - x[i-1]) * (y[i] + y[i-1]) / 2
    end
    return z
end

# Define wave conditions
waves = WaveClass("regular")
waves.height = 2.5
waves.period = 8

# Define wave properties
H = waves.height                # Wave height
A = H / 2                       # Amplitude
T = waves.period                # Wave period
omega = (1 / T) * (2 * π)       # Angular frequency

# Find wave period index
waveInd = argmin(abs.(T .- hydro.T))

println("Dimensions of A: ", size(hydro.A))
println("Dimensions of B: ", size(hydro.B))
println("Dimensions of ex_re: ", size(hydro.ex_re))
println("Dimensions of Khs: ", size(hydro.Khs))

# Define excitation force based on wave conditions
ampSpect = zeros(length(hydro.w))
closestIndOmega = argmin(abs.(omega .- hydro.w))
ampSpect[closestIndOmega] = A

# Correct access pattern for excitation force coefficients
FeRao = hydro.ex_re[3, 1, :] * hydro.rho * hydro.g .+
        hydro.ex_im[3, 1, :] * hydro.rho * hydro.g * im
Fexc = ampSpect .* FeRao

# Define the intrinsic mechanical impedance for the device
mass = hydro.rho * hydro.Vo[1]  # Assuming Vo[1] is volume of the body
addedMass = hydro.A[3, 3, :]    # Added mass for surge (3rd DOF)
radiationDamping = hydro.B[3, 3, :] # Radiation damping for surge (3rd DOF)
hydrostaticStiffness = hydro.Khs[3, 3, 1] # Hydrostatic stiffness for surge (3rd DOF), assuming single value

println("Contents of addedMass: ", addedMass)
println("Contents of radiationDamping: ", radiationDamping)
println("Contents of hydrostaticStiffness: ", hydrostaticStiffness)

Zi = -(im .* hydro.w .* (mass .+ addedMass)) .+ radiationDamping .+ hydrostaticStiffness ./ (im .* hydro.w)

# Generate excitation force time series
t = range(0, stop=400, length=4001)
F_wave = A * cos.(omega .* t) .* hydro.ex_re[3, 1, closestIndOmega] * hydro.rho * hydro.g .- 
         A * sin.(omega .* t) .* hydro.ex_im[3, 1, closestIndOmega] * hydro.rho * hydro.g

# Reshape F_wave to be a 1-row matrix
F_wave_matrix = reshape(F_wave, 1, :)
println("Contents of F_wave_matrix: ", F_wave_matrix)
println("First 100 data points of F_wave_matrix: ", F_wave_matrix[1, 1:100])

# Plot excitation time series to make sure it matches WEC-Sim
plot(t, F_wave, label="Impedance model")
# Uncomment and update the next line with actual WEC-Sim output data to compare
# plot!(output.bodies.time, output.bodies.forceExcitation[:,3], label="WEC-Sim", linestyle=:dash)
xlabel!("Time (s)")
ylabel!("Force (N)")

# Calculate the response based on timeseries excitation and impedance
# Impedance in frequency domain can be converted to transfer function
# ((m + A)s^2 + Bs + C_hs)/s - flip to calculate velocity from force!
Zi_den = [mass + addedMass[closestIndOmega], radiationDamping[closestIndOmega], hydrostaticStiffness] # At closest (wave) frequency
Zi_num = [1, 0]
Zi_tf = tf(Zi_num, Zi_den)

# Calculate velocity from impedance, then position
response = lsim(Zi_tf, F_wave_matrix, t)
vel = response.y[1, :]  # Extract the response for the first (and only) output
pos = cumtrapz(t, vel) # Integrate velocity to get position

# Plot velocity - it matches WEC-Sim output!
plot(t, vel, label="Impedance model")
# Uncomment and update the next line with actual WEC-Sim output data to compare
# plot!(output.bodies.time, output.bodies.velocity[:,3], label="WEC-Sim", linestyle=:dash)
xlabel!("Time (s)")
ylabel!("Velocity (m/s)")

# Plot position
plot(t, pos, label="Impedance model")
# Uncomment and update the next line with actual WEC-Sim output data to compare
# plot!(output.bodies.time, output.bodies.position[:,3] .- mean(output.bodies.position[:,3]), label="WEC-Sim", linestyle=:dash)
xlabel!("Time (s)")
ylabel!("Position (m)")

# Render the plots in the Julia plots pane
display(plot(t, F_wave, label="Impedance model", xlabel="Time (s)", ylabel="Force (N)"))
display(plot(t, vel, label="Impedance model", xlabel="Time (s)", ylabel="Velocity (m/s)"))
display(plot(t, pos, label="Impedance model", xlabel="Time (s)", ylabel="Position (m)"))
