###################################################################################################
include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/readWAMIT.jl")
include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/radiationIRFv2.jl")
#include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/radiationIRFSS.jl")
#include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/excitationIRF.jl")
#include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/writeBEMIOH5.jl")
#include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/plotBEMIO.jl")
###################################################################################################

# using DelimitedFiles
hydro = Dict{Any, Any}()
# Hydro data
hydro = readWAMIT(filepath = "C:/Users/jelope/Desktop/Git/WEC-Sim/examples/RM3/hydroData/rm3.out")  # Replace 'readWAMIT' with appropriate code to read .out file

# Radiation IRF
# Call the radiationIRF function with the data dictionary
data = radiationIRF(hydro)


# Radiation IRFSS
hydro = radiationIRFSS( nothing, nothing)  # Replace 'radiationIRFSS' with appropriate code

# Excitation IRF
hydro = excitationIRF(hydro, 157, nothing, nothing, nothing, nothing)  # Replace 'excitationIRF' with appropriate code

# Write BEMIOH5
writeBEMIOH5(hydro)  # Replace 'writeBEMIOH5' with appropriate code

# Plot hydro data
plotBEMIO(hydro)

