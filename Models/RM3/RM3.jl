using ModelingToolkit
using DifferentialEquations
using LinearAlgebra
using ModelingToolkitStandardLibrary.Mechanical.Translational
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using ModelingToolkit: t_nounits as t
using Plots
using Symbolics: scalarize

# Include the readWAMITv2.jl file
include("C:/Users/jaime/OneDrive/Documents/Sandia National Labs/Internship/WEC-Sim-main/WEC-Sim-main/source/functions/BEMIO/readWAMITv2.jl")

# Load the hydrodynamic data
hydroData = readWAMIT("C:/Users/jaime/OneDrive/Documents/Sandia National Labs/Internship/WEC-Sim-main/WEC-Sim-main/examples/RM3/hydroData/rm3.out")

# User-specified value to search within the :T vector
user_input_T = 8  # Example value, replace with actual user input

# Find the index of the closest value to the user-specified value within the :T vector
wave_index = argmin(abs.(hydroData.T .- user_input_T))

# Print dimensions for debugging
println("Dimensions of A: ", size(hydroData.A))
println("Dimensions of B: ", size(hydroData.B))
println("Dimensions of ex_re: ", size(hydroData.ex_re))
println("Dimensions of Khs: ", size(hydroData.Khs))

# Define function to get parameters based on wave index
function get_wave_parameters(hydroData, wave_index)
    if wave_index <= size(hydroData.A, 3) && wave_index <= size(hydroData.B, 3) && wave_index <= size(hydroData.ex_re, 3)
        A = hydroData.A[3, 3, wave_index]  # Added mass for heave (3,3) at the given wave index
        B = hydroData.B[3, 3, wave_index]  # Damping coefficient for heave (3,3) at the given wave index
        F_ext = hydroData.ex_re[3, 1, wave_index]  # External forcing function for heave (3) at the given wave index
        return A, B, F_ext
    else
        error("Wave index out of bounds for the given hydrodynamic data.")
    end
end

# Extract the 3,3 element for the heave from the first configuration of Khs
spring_constant = hydroData.Khs[3, 3, 1]

# Get wave parameters for the given wave index
A, B, F_ext = get_wave_parameters(hydroData, wave_index)

# Print the wave index and parameters to verify
println("Wave Index: ", wave_index)
println("A (Added Mass): ", A)
println("B (Damping Coefficient): ", B)
println("F_ext (External Forcing Function): ", F_ext)
println("K_hs (Hydrostatic Stiffness):", spring_constant)

# Define structure for SystemBody with moments of inertia and center of gravity
mutable struct SystemBody
    name::String
    mass::Float64
    initial_position::Float64
    spring_constant::Float64
    initial_spring_stretch::Float64
    damping_coefficient::Float64
    forcing_frequency::Float64
    forcing_amplitude::Float64
    center_of_gravity::Vector{Float64}
    inertia_tensor::Matrix{Float64}
end

# Initialize SystemBody with extracted values and given moments of inertia and center of gravity
body1 = SystemBody(
    "Float",
    725833.0 + A,    # mass in kg plus added mass A
    0.0,             # initial position in meters
    spring_constant, # spring constant in N/m
    0.0,             # initial spring stretch in meters
    B,               # damping coefficient in Ns/m
    hydroData.w[wave_index], # forcing frequency in rad/s (corresponding to the wave index)
    F_ext,           # forcing amplitude in N
    [0.0, 0.0, -0.72],         # center of gravity in meters
    [20907301.0 0.0 0.0; 0.0 21306091.0 0.0; 0.0 0.0 37085481.0]  # inertia tensor in kg*m^2
)

# Variables from the body1 structure
m = body1.mass
s = body1.initial_position
k = body1.spring_constant
delta_s = body1.initial_spring_stretch
c = body1.damping_coefficient
f = body1.forcing_frequency
a = body1.forcing_amplitude

@named mass = Translational.Mass(m=m, s=s) # using full function name avoids conflicts
@named spring = Translational.Spring(k=k, delta_s=delta_s)
@named damper = Translational.Damper(d=c)
@named frame = Translational.Fixed()
@named sine_source = Blocks.Sine(frequency=f, amplitude=a)
@named force = Translational.Force()

# setup mass spring damper
msd_eqs = [connect(frame.flange, spring.flange_a),
           connect(frame.flange, damper.flange_a),
           connect(damper.flange_b, mass.flange),
           connect(spring.flange_b, mass.flange),
           connect(mass.flange, force.flange),
           connect(sine_source.output, force.f)]

@named msd_model = ODESystem(msd_eqs, t, systems=[mass, spring, damper, frame, 
                             sine_source, force])

println("Set Up Successfully")

sys = structural_simplify(msd_model)
prob = ODEProblem(sys, Pair[], (0, 400.0))
sol = solve(prob)

# Ensure the plot is displayed
display(plot(sol, idxs=[mass.s, mass.v],
     title="Mass Spring Damper",
     labels=["Position" "Velocity"]))

