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
option = 0  # Set option: 0 for WAMIT, 1 for OSWEC

if option == 0
    # WAMIT
    hydro = readWAMIT(filepath = "C:/Users/jelope/Desktop/Git/WEC-Sim/examples/RM3/hydroData/rm3.out")
    data = radiationIRF(hydro)
    hydro = radiationIRFSS(nothing, nothing)
    hydro = excitationIRF(hydro, 157, nothing, nothing, nothing, nothing)
    writeBEMIOH5(hydro)

    # Plot hydro data
    #plotBEMIO(hydro)
elseif option == 1
    # OSWEC
    hydro = readWAMIT(filepath = "C:/Users/jelope/Desktop/Git/WEC-Sim/examples/OSWEC/hydroData/oswec.out")
    hydro = radiationIRF(hydro, 30, nothing, nothing, nothing, nothing)
    hydro = radiationIRFSS(hydro, nothing, nothing)
    hydro = excitationIRF(hydro, 30, nothing, nothing, nothing, nothing)
end


