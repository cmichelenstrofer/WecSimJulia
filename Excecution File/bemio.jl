include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/readWAMITv2.jl")
include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/radiationIRFv2.jl")
include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/radiationIRFSSv2.jl")
include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/excitationIRF.jl")
include("C:/Users/jelope/Desktop/Git/WEC-Sim/source/functions/BEMIO/normalizeBEM.jl")

# User input
option = 0

if option == 0
    # WAMIT
    hydro = readWAMIT("C:/Users/jelope/Desktop/Git/WEC-Sim/examples/RM3/hydroData/rm3.out")
    hydro = radiationIRF(hydro, 60)
    hydro = radiationIRFSS(hydro)
    hydro = excitationIRF(hydro, 157)

    # Print field values for comparison
    print_hydro_data(hydro)

    #writeBEMIOH5(hydro)

    # Plot hydro data
    # plotBEMIO(hydro)
elseif option == 1
    # OSWEC
    hydro = readWAMIT("C:/Users/jelope/Desktop/Git/WEC-Sim/examples/OSWEC/hydroData/oswec.out")

    #hydro = radiationIRF(hydro, 30.0f0, 1001, 1001, minimum(hydro.w), maximum(hydro.w))
    #hydro = radiationIRFSS(hydro)
    #hydro = excitationIRF(hydro, 30, nothing, nothing, nothing, nothing)
    #writeBEMIOH5(hydro)

    # Print field values for comparison
    print_hydro_data(hydro)

    # Plot hydro data
    # plotBEMIO(hydro)
else
    println("Invalid option selected. Please run the program again and select either 0 or 1.")
end
