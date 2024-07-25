using ModelingToolkit, DifferentialEquations, LinearAlgebra, ModelingToolkitStandardLibrary
using ModelingToolkitStandardLibrary.Mechanical.Translational
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t
using Plots
using Symbolics: scalarize

# Include the readWAMIT.jl file
include("C:/Users/jelope/Desktop/WecSimJulia-main/ReadWriteFILE/readWAMIT.jl")

hydro = readWAMIT("C:/Users/jelope/Desktop/WecSimJulia-main/ReadWriteFILE/InputFiles/RM3/rm3.out")

# User-specified value to search within the :T vector
user_input_T = 8

# Find the index of the closest value to the user-specified value within the :T vector
wave_index = argmin(abs.(hydro.T .- user_input_T))

# Print dimensions for debugging
println("Dimensions of A: ", size(hydro.A))
println("Dimensions of B: ", size(hydro.B))
println("Dimensions of ex_re: ", size(hydro.ex_re))
println("Dimensions of Khs: ", size(hydro.Khs))

# Define function to get parameters based on wave index
function get_wave_parameters(hydroData, wave_index)
    if wave_index <= size(hydro.A, 3) && wave_index <= size(hydro.B, 3) && wave_index <= size(hydro.ex_re, 3)
        A = hydro.A[3, 3, wave_index]*hydro.rho  # Added mass for heave (3,3) at the given wave index
        B = hydro.B[3, 3, wave_index]*(hydro.w[wave_index])*hydro.rho  # Damping coefficient for heave (3,3) at the given wave index
        F_ext = hydro.ex_re[3, 1, wave_index]  # External forcing function for heave (3) at the given wave index
        return A, B, F_ext
    else
        error("Wave index out of bounds for the given hydrodynamic data.")
    end
end

# Extract the 3,3 element for the heave from the first configuration of Khs
spring_constant = hydro.Khs[3, 3, 1]*hydro.rho*hydro.g

# Get wave parameters for the given wave index
A, B, F_ext = get_wave_parameters(hydro, wave_index)

# Print the wave index and parameters to verify
println("Wave Index: ", wave_index)
println("A (Added Mass): ", hydro.A[3,3, wave_index])
println("Normalized Added Mass: ", A)
println("B (Damping Coefficient): ", B)
println("F_ext (External Forcing Function): ", F_ext)
println("K_hs (Hydrostatic Stiffness):", spring_constant)

# Define all variables
m = (hydro.rho * hydro.Vo[1]) + A
s = 0.0 # initial mass position
k = spring_constant
delta_s = 0 # initial spring stretch
c = B # Damping coefficient
f =  (hydro.w[wave_index]/(2π)) # sine forcing frequency
a = 1.9049e+06 # sine forcing amplitude
J = 1
gear_ratio = 1
k_motor = .5

println("Float Volume (Vo[1]): ", hydro.Vo[1])
println("Rho (rho): ", hydro.rho)

@named mass = Translational.Mass(m=m,s=s) # using full function name avoids conflicts
@named spring = Translational.Spring(k=k,delta_s=delta_s)
@named damper = Translational.Damper(d=c)
@named frame = Translational.Fixed()
@named sine_source = Blocks.Sine(frequency = f, amplitude = a)
@named force = Translational.Force()

# Setup mass spring damper
connections = [
    ModelingToolkitStandardLibrary.Mechanical.Translational.connect(frame.flange, spring.flange_a),
    ModelingToolkitStandardLibrary.Mechanical.Translational.connect(frame.flange, damper.flange_a),
    ModelingToolkitStandardLibrary.Mechanical.Translational.connect(damper.flange_b, mass.flange),
    ModelingToolkitStandardLibrary.Mechanical.Translational.connect(spring.flange_b, mass.flange),
    ModelingToolkitStandardLibrary.Mechanical.Translational.connect(mass.flange, force.flange),
    ModelingToolkitStandardLibrary.Blocks.connect(sine_source.output, force.f)
]

@named msd_model = ODESystem(connections, t, systems = [mass, spring, damper, frame, sine_source, force])

println("Set Up Successfully")

sys = structural_simplify(msd_model)
prob = ODEProblem(sys, Pair[], (0, 400.0))
sol = solve(prob)
# Plotting
display(plot(sol, idxs=[mass.s], title="Position", labels=["Position"]))
display(plot(sol, idxs=[mass.v], title="Velocity", labels=["Velocity"]))
display(plot(sol, idxs=[mass.f], title="Force", labels=["Force"]))
